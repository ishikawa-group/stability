import sys
import os
sys.path.append("../")

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")

from ase.io import read

from stability.convex_hull.do_vasp_calculation import do_vasp_calculation
from stability.convex_hull.get_energy_above_hull import get_energy_above_hull

# bulk = read("BaZrO3.cif")
bulk = read("Ba2Sr6Zr8O24.POSCAR")
# bulk = read("Ba4Sr4Zn7Ce1O24.POSCAR")
# bulk = read("Ba4Sr4Zr2Al6O24.POSCAR")

# bulk = bulk*[2, 2, 2]

do_vasp = False

if do_vasp:
    do_vasp_calculation(atoms=bulk)

# energy = -333.78748451  # Ba8 Zr8 O24
energy = -331.976283      # Ba2 Sr6 Zr8 O24
# energy = -212.112054    # Ba4 Sr4 Zn7 Ce1 O24
# energy = -281.23813     # Ba4 Sr4 Zr2 Al6 O24

e_above_hull_A, e_above_hull_C = get_energy_above_hull(atoms=bulk, energy=energy)

print(f"energy above hull at anode/cathode conditions (per atom, eV): {e_above_hull_A:5.3f}, {e_above_hull_C:5.3f}")
