#!/usr/bin/env python
# coding: utf-8

from mmtbx.ncs.ncs import ncs
from phenix.programs import map_symmetry as map_symmetry_program
from cctbx.maptbx.segment_and_split_map import run_get_ncs_from_map
from iotbx import phil
from iotbx.data_manager import DataManager


# initialize data manager
dm = DataManager()


# file IO

# this map has origin at (0,0,0)
# dm.process_model_file("../6ui6_fit_in_corner_map.pdb")
# dm.process_real_map_file("../emd_20669_corner_zero.map")

# this map has origin at data.all()/2
dm.process_model_file("../6ui6.pdb")
dm.process_real_map_file("../emd_20669.map")


# Run Tom's map_symmetry tool
mm = dm.get_real_map()
mm.shift_origin()
mm.find_map_symmetry()
print(type(mm._ncs_object) # this was not set...

# debug find_map_symmetry()
self = mm
include_helical_symmetry = False
symmetry_center = None
min_ncs_cc = None
symmetry = None
ncs_object = None
check_crystal_symmetry = True
only_proceed_if_crystal_symmetry = False
assert self.origin_is_zero()

self._warning_message = ""
self._ncs_cc = None

from cctbx.maptbx.segment_and_split_map import \
   run_get_ncs_from_map, get_params

if symmetry is None:
  symmetry = 'ALL'


if symmetry_center is None:
  # Most likely map center is (1/2,1/2,1/2) in full grid
  full_unit_cell=self.unit_cell_crystal_symmetry(
        ).unit_cell().parameters()[:3]
  symmetry_center=[]
  for x, sc in zip(full_unit_cell, self.shift_cart()):
    
    # version 1 # original version gives incorrect symmetry
    #symmetry_center.append(0.5*x + sc)

    # version 2 # this finds I(b), the correct symmetry
    symmetry_center.append(0.5*x)

  symmetry_center = tuple(symmetry_center)

params = get_params(args=[],
  symmetry = symmetry,
  include_helical_symmetry = include_helical_symmetry,
  symmetry_center = symmetry_center,
  min_ncs_cc = min_ncs_cc,
  return_params_only = True,
  )

space_group_number = None
if check_crystal_symmetry and symmetry == 'ALL' and (not ncs_object):
  # See if we can narrow it down looking at intensities at low-res
  d_min = 0.05*self.crystal_symmetry().unit_cell().volume()**0.333
  map_coeffs = self.map_as_fourier_coefficients(d_min=d_min)
  from iotbx.map_model_manager import get_map_coeffs_as_fp_phi
  f_array_info = get_map_coeffs_as_fp_phi(map_coeffs, d_min = d_min,
    n_bins = 15)
  ampl = f_array_info.f_array
  data = ampl.customized_copy(
    data = ampl.data(),sigmas = flex.double(ampl.size(),1.))
  from mmtbx.scaling.twin_analyses import symmetry_issues
  si = symmetry_issues(data)
  cs_possibility = si.xs_with_pg_choice_in_standard_setting
  space_group_number = cs_possibility.space_group_number()
  
# # neccessary to remove or for me it will return None
#   if space_group_number < 2:
#     space_group_number = None
#   if space_group_number is None and only_proceed_if_crystal_symmetry:
#     return # skip looking further

params.reconstruction_symmetry.\
      must_be_consistent_with_space_group_number = space_group_number
new_ncs_obj, ncs_cc, ncs_score = run_get_ncs_from_map(params = params,
    map_data = self.map_data(),
    crystal_symmetry = self.crystal_symmetry(),
    out = sys.stdout,
    ncs_obj = ncs_object)

# Build cif model from ncs_obj
from iotbx import cif

model = dm.get_model()
h = model.get_hierarchy()
chains = [c.id for c in h.chains()]

n_oper = ncs_group.n_ncs_oper()

# start cif building
builder = cif.builders.cif_model_builder()
builder.add_data_block("assembly_information")

# add pdbx_struct_assembly loop
headers = ['_pdbx_struct_assembly.id',
 '_pdbx_struct_assembly.details',
 '_pdbx_struct_assembly.method_details',
 '_pdbx_struct_assembly.oligomeric_details',
 '_pdbx_struct_assembly.oligomeric_count']
columns = [["1"],["Symmetry assembly "+ncs_obj.get_ncs_name()],["?"],["?"],["?"]]
builder.add_loop(headers,columns)

# add pdbx_struct_assembly_gen loop
headers = ['_pdbx_struct_assembly_gen.assembly_id',
 '_pdbx_struct_assembly_gen.oper_expression',
 '_pdbx_struct_assembly_gen.asym_id_list']
columns = [["1"],["(1-"+str(n_oper)+")"],[','.join(chains)]]
builder.add_loop(headers,columns)

# add pdbx_struct_oper_list loop
headers = ['_pdbx_struct_oper_list.id',
 '_pdbx_struct_oper_list.type',
 '_pdbx_struct_oper_list.name',
 '_pdbx_struct_oper_list.symmetry_operation',
 '_pdbx_struct_oper_list.matrix[1][1]',
 '_pdbx_struct_oper_list.matrix[1][2]',
 '_pdbx_struct_oper_list.matrix[1][3]',
 '_pdbx_struct_oper_list.vector[1]',
 '_pdbx_struct_oper_list.matrix[2][1]',
 '_pdbx_struct_oper_list.matrix[2][2]',
 '_pdbx_struct_oper_list.matrix[2][3]',
 '_pdbx_struct_oper_list.vector[2]',
 '_pdbx_struct_oper_list.matrix[3][1]',
 '_pdbx_struct_oper_list.matrix[3][2]',
 '_pdbx_struct_oper_list.matrix[3][3]',
 '_pdbx_struct_oper_list.vector[3]']


_id = list(range(1,n_oper+1))
_type = ['point symmetry operation' for i in range(n_oper)]
_name = ['?' for i in range(n_oper)]
_symmetry_operation = ['?' for i in range(n_oper)]

info_columns = [_id,_type,_name,_symmetry_operation]

rotations  = [[r[i] for r in ncs_group.rota_matrices_inv()] for i in range(9)]
#translations = [[t[i] for t in ncs_group.translations_orth_inv()] for i in range(3)]
translations = [[0.0 for t in ncs_group.translations_orth_inv()] for i in range(3)] # debug, translations are not meaningful
numeric_columns = [rotations[0],rotations[1],rotations[2],translations[0],
                   rotations[3],rotations[4],rotations[5],translations[1],
                   rotations[6],rotations[7],rotations[8],translations[2]]

columns = info_columns+numeric_columns
builder.add_loop(headers,columns)




# Combine cif string with mmcif string. 

from StringIO import StringIO
output = StringIO()
cif_model = builder.model()
cif_model = cif_model[cif_model.keys()[0]]
cif_model.show(indent="",out=output)
filestring = "data_default\n"+output.getvalue().replace("data_information\n","")+dm.get_model().model_as_mmcif().replace("data_default\n","")


dm.write_model_file(filestring,filename="../6ui6_processed.cif",overwrite=True)

