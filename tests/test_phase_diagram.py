import sys
import os
sys.path.append("../")

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")

from ase.io import read

from stability.convex_hull.do_vasp_calculation import do_vasp_calculation
from stability.convex_hull.get_energy_above_hull import get_energy_above_hull

bulk = read("BaZrO3.cif")
bulk = bulk*[2, 2, 2]

do_vasp = False

if do_vasp:
    do_vasp_calculation(atoms=bulk)

energy = -333.78748451  # -333.820839  # vasp-calculated energy
e_above_hull_A, e_above_hull_C = get_energy_above_hull(atoms=bulk, energy=energy)

print(f"energy above hull at anode/cathode conditions (per atom, eV): {e_above_hull_A:5.3f},{e_above_hull_C:5.3f}")
