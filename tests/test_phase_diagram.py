import sys
import os
import warnings

# Ignore warnings
warnings.filterwarnings("ignore")

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from stability.convex_hull.phase_diagram import (
    initialize_global_variables,
    prepare_material_entries,
    calculate_phase_diagram_condition_A,
    calculate_phase_diagram_condition_C,
    calculate_phase_diagram_condition_X,
)

# Set test material and energy
api = os.getenv("MAPI")
TestMat_Comp = "Ba8Zr8O24"
TestMat_Ener = -331.28931146  # VASP-calculated value in eV

# Initialize global variables
initialize_global_variables()

from stability.convex_hull.phase_diagram import (
    entriesGases_A, locked_Chem_Potential_A,
    entriesGases_C, locked_Chem_Potential_C,
    entriesGases_X, locked_Chem_Potential_X,
)

# Prepare material entries
(   all_entries_A,
    all_entries_C,
    entriesTotal_X,
    TestMat_entry_A,
    TestMat_entry_C,
    TestMat_entry_X
) = prepare_material_entries(api, TestMat_Comp, TestMat_Ener)

# Calculate phase diagrams and energies for each condition
pd_A, energy_per_atom_A, hull_energy_A, energy_above_hull_A = (
    calculate_phase_diagram_condition_A(
        all_entries_A, entriesGases_A, locked_Chem_Potential_A, TestMat_entry_A
    )
)

pd_C, energy_per_atom_C, hull_energy_C, energy_above_hull_C = (
    calculate_phase_diagram_condition_C(
        all_entries_C, entriesGases_C, locked_Chem_Potential_C, TestMat_entry_C
    )
)

"""
pd_X, energy_per_atom_X, hull_energy_X, energy_above_hull_X = (
    calculate_phase_diagram_condition_X(
        entriesTotal_X, entriesGases_X, locked_Chem_Potential_X, TestMat_entry_X
    )
)
"""

# Print results for condition A (Hydrogen-rich)
print("Condition A (Hydrogen-rich):")
print(f"Energy per atom: {energy_per_atom_A:5.3f} eV/atom")
print(f"Hull energy per_atom: {hull_energy_A:5.3f} eV/atom")
print(f"Energy above hull per atom: {energy_above_hull_A:5.3f} eV/atom")

# Print results for condition C (Oxygen-rich)
print("\nCondition C (Oxygen-rich):")
print(f"Energy per atom: {energy_per_atom_C:5.3f} eV/atom")
print(f"Hull energy per_atom: {hull_energy_C:5.3f} eV/atom")
print(f"Energy above hull per atom: {energy_above_hull_C:5.3f} eV/atom")

"""
# Print results for condition X (CO2-rich)
print("\nCondition X (CO2-rich):")
print(f"Energy per atom: {energy_per_atom_X:5.3f} eV/atom")
print(f"Hull energy per_atom: {hull_energy_X:5.3f} eV/atom")
print(f"Energy above hull per atom: {energy_above_hull_X:5.3f} eV/atom")
"""
