"""
The following code was used to create some files for the inpgen tests
"""

from aiida import load_profile
load_profile()

from aiida import orm
import ase
from aiida.tools.dbimporters.plugins.cod import CodDbImporter
import spglib


# For test
# there are 530 hall numbers in spglib
nelement = 2
count_in  = 0
count_sec = 0
count_not = 0
for i in range(1,531):
    hall_number = i
    print(hall_number)
    spacegroup_type = spglib.get_spacegroup_type(hall_number)
    sym = spglib.get_symmetry_from_database(hall_number)
    space_group_string = spacegroup_type['international_short']
    space_group_number = spacegroup_type['number']
    print(space_group_number, spacegroup_type['arithmetic_crystal_class_number'])
    res = CodDbImporter().query(spacegroup=space_group_string)#, number_of_elements=nelement)
    # use first result only
    if len(list(res)) > 0:
        print(res[0])
        count_in = count_in + 1
    else:
        res = CodDbImporter().query(spacegroup=spacegroup_type['international_full'])
        if len(list(res)) > 0:
            print('second ', res[0])
            count_sec = count_sec + 1
        else:
            count_not = count_not + 1
            print('########## Not in cod: {} {}'.format(hall_number, space_group_string))

print(count_in, count_sec, count_not)

