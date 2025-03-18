import sys
import os
sys.path.append("../")

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")

from ase.io import read

from stability.convex_hull.do_vasp_calculation import do_vasp_calculation
from stability.convex_hull.get_energy_above_hull import get_energy_above_hull

bulk = read("Ba1Sr7Nb8O24.vasp")
bulk = read("Ba1Sr7V8O24.vasp")
bulk = read("Ca1La7Fe8O24.vasp")

# bulk = bulk*[2, 2, 2]  # needed when loading BaZrO3.cif (to make Ba8Zr8O24)

do_vasp = False

if do_vasp:
    do_vasp_calculation(atoms=bulk)

# -- energy values from author's CSV file (corrected)
# energy = -306.463721  # Sm4_Y4_Ni8_O24
# energy = -339.945046  # Ba1_Sr7_Nb8_O24
# energy = -320.169744  # Ba1_Sr7_V8_O24
energy = -331.534147  # Ca1_La7_Fe8_O24

e_above_hull_A, e_above_hull_C, e_above_hull_X = get_energy_above_hull(atoms=bulk, energy=energy)

print("\n")
print(f"energy above hull at anode condition (per atom, eV): {e_above_hull_A:6.4f}")
print(f"energy above hull at cathode condition (per atom, eV): {e_above_hull_C:6.4f}")
print(f"energy above hull at CO2 condition (per atom, eV): {e_above_hull_X:6.4f}")

