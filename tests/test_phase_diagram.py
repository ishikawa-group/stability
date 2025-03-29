from ase.io import read
from stability import convex_hull
import warnings

warnings.filterwarnings("ignore")

# bulk = read("Ba1Sr7Nb8O24.vasp")  # anode: 0.382792854, cathode: 0.842792854, CO2: 0.846661652
# bulk = read("Ba1Sr7V8O24.vasp")   # anode: 0.212304164, cathode: 0.528835205, CO2: 0.528835205
# bulk = read("Ca1La7Fe8O24.vasp")  # anode: 0.260603574, cathode: 0.190561059, CO2: 0.20452976
# bulk = read("Cs1Ba7Nb8O24.vasp")  # anode: 0.261439645, cathode: 0.663939645, CO2: 0.680963907
# bulk = read("Cs1Sr7Nb8O24.vasp")  # anode: 0.363453727, cathode: 0.765953727, CO2: 0.765953727
# bulk = read("Cs1Sr7V8O24.vasp")   # anode: 0.299563552, cathode: 0.561660693, CO2: 0.561660693
bulk = read("Rb1Ba7Nb8O24.vasp")  # anode: 0.216330437, cathode: 0.618830437, CO2: 0.642664404

# bulk = bulk*[2, 2, 2]  # needed when loading BaZrO3.cif (to make Ba8Zr8O24)

do_vasp = False

if do_vasp:
    convex_hull.do_vasp_calculation(atoms=bulk)

# -- energy values from author's CSV file (corrected)
# energy = -339.945046   # Ba1_Sr7_Nb8_O24
# energy = -320.169744   # Ba1_Sr7_V8_O24
# energy = -331.534147   # Ca1_La7_Fe8_O24
# energy = -338.659999   # Cs1_Ba7_Nb8_O24
# energy = -337.02237    # Cs1_Sr7_Nb8_O24
# energy = -315.246376   # Cs1_Sr7_V8_O24
energy = -339.378272   # Rb1_Ba7_Nb8_O24

e_above_hull_A, e_above_hull_C, e_above_hull_X = convex_hull.get_energy_above_hull(atoms=bulk, energy=energy)

print("\n")
print(f"energy above hull at anode condition (per atom, eV): {e_above_hull_A:8.6f}")
print(f"energy above hull at cathode condition (per atom, eV): {e_above_hull_C:8.6f}")
print(f"energy above hull at CO2 condition (per atom, eV): {e_above_hull_X:8.6f}")
