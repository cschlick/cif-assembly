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
dm.process_model_file("../6ui6_fit_in_corner_map.pdb")
dm.process_real_map_file("../emd_20669_corner_zero.map")

# this map has origin at data.all()/2
# dm.process_model_file("../6ui6.pdb")
# dm.process_real_map_file("../emd_20669.map")


# Run Tom's map_symmetry tool
mm = dm.get_real_map()
params = phil.parse(map_symmetry_program.Program.master_phil_str).extract()
ncs_obj, cc_avg, score = run_get_ncs_from_map(params=params,
        map_data=mm.map_data(),
        crystal_symmetry=mm.crystal_symmetry(),
        ncs_obj=None)


# get ncs_group object from the map_symmetry output
ncs_groups = ncs_obj.ncs_groups()
assert(len(ncs_groups)==1)
ncs_group = ncs_groups[0]

# or if reading from file...
# ncs_object = ncs()
# ncs_object.read_ncs("../MapSymmetry_4/symmetry_from_map.ncs_spec")


# this already exists
print(ncs_group.format_for_biomt())


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
translations = [[t[i] for t in ncs_group.translations_orth_inv()] for i in range(3)]
#translations = [[0.0 for t in ncs_group.translations_orth_inv()] for i in range(3)] # debug
numeric_columns = [rotations[0],rotations[1],rotations[2],translations[0],
                   rotations[3],rotations[4],rotations[5],translations[1],
                   rotations[6],rotations[7],rotations[8],translations[2]]

columns = info_columns+numeric_columns
builder.add_loop(headers,columns)


# show cif model model
builder.model().show()


model.model_as_mmcif()








# Combine cif string with mmcif string. 

from StringIO import StringIO
output = StringIO()
cif_model = builder.model()
cif_model = cif_model[cif_model.keys()[0]]
cif_model.show(indent="",out=output)
filestring = "data_default\n"+output.getvalue().replace("data_information\n","")+dm.get_model().model_as_mmcif().replace("data_default\n","")


dm.write_model_file(filestring,filename="../6ui6_processed.cif",overwrite=True)

