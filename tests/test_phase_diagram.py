import sys
import os
sys.path.append("../")

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")

from ase.io import read

from stability.convex_hull.do_vasp_calculation import do_vasp_calculation
from stability.convex_hull.get_energy_above_hull import get_energy_above_hull

bulk = read("Ba8Zr8O24.POSCAR")
# bulk = read("Ba2Sr6Zr8O24.POSCAR")
# bulk = read("Ba4Sr4Zn7Ce1O24.POSCAR")
# bulk = read("Ba4Sr4Zr2Al6O24.POSCAR")

# bulk = bulk*[2, 2, 2]  # needed when loading BaZrO3.cif (to make Ba8Zr8O24)

do_vasp = False

if do_vasp:
    do_vasp_calculation(atoms=bulk)

# -- energy values from author's CSV file (raw energy)
energy = -333.820839   # Ba8 Zr8 O24,
# energy = -331.97628    # Ba2 Sr6 Zr8 O24
# energy = -212.11205    # Ba4 Sr4 Zn7 Ce1 O24
# energy = -281.23813    # Ba4 Sr4 Zr2 Al6 O24

e_above_hull_A, e_above_hull_C, e_above_hull_X = get_energy_above_hull(atoms=bulk, energy=energy)

print("\n")
print(f"energy above hull at anode condition (per atom, eV): {e_above_hull_A:5.3f}")
print(f"energy above hull at cathode condition (per atom, eV): {e_above_hull_C:5.3f}")
print(f"energy above hull at CO2 condition (per atom, eV): {e_above_hull_X:5.3f}")
