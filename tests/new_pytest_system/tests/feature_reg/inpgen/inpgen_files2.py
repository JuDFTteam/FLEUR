"""
The following code was used to create some files for the inpgen tests
"""

from aiida import load_profile
load_profile()

from aiida import orm
import ase
from aiida.tools.dbimporters.plugins.cod import CodDbImporter
import spglib
from aiida.common.constants import elements
from aiida_fleur.tools.StructureData_util import supercell_ncf
#from masci_tools.io.io_fleur_inpgen import write_inpgen_file
from aiida_fleur.calculation.fleurinputgen import write_inpgen_file_aiida_struct

# One fictious sc structure with all elements fleur supports
largest_charge = 110 # coresolver cutoff at 
all_elm_struc = orm.StructureData()
alat = float(largest_charge*5.0)
cell = [[alat, 0.0, 0.0],[0.0, alat, 0.0],[0.0, 0.0, alat]]
all_elm_struc.cell = cell
for i in range(1, largest_charge + 1):# enumerate(Periodictable.keys()):
    sym = elements[i]['symbol']
    print(sym)
    position = [alat/largest_charge*i, 0.0, 0.0]
    all_elm_struc.append_atom(position=position, symbols=[sym])

print(all_elm_struc)
report = write_inpgen_file_aiida_struct(all_elm_struc, path='./inpgen_input_files/file1.in')

# One large arbitrary CuO structure with 1024 K atoms
alat = 3.2
cu_structure = orm.StructureData(cell=[[alat, 0.0, 0.0],[0.0, alat, 0.0],[0.0, 0.0, alat]])
cu_structure.append_atom(position=[0.0, 0.0, 0.0], symbols=['Cu'])
cu_structure.append_atom(position=[alat/2.0, alat/2.0, alat/2.0], symbols=['O'])

# 10 K 
structure = supercell_ncf(cu_structure, 17, 17, 17)

# 1 Mio
#structure = supercell_ncf(cu_structure, 80, 80, 80)
print(structure)
print(len(structure.sites))
report = write_inpgen_file_aiida_struct(structure, path='./inpgen_input_files/file3.in')
