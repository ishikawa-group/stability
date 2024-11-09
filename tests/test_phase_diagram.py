import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from stability.convex_hull.phase_diagram import (
    initialize_global_variables,
    prepare_material_entries,
    calculate_phase_diagram_condition_A,
    calculate_phase_diagram_condition_C,
    calculate_phase_diagram_condition_X,
)

# Set test material and energy
api = "kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK"
TestMat_Comp = "Ba8Zr8O24"
TestMat_Ener = -333.71584216

# Initialize global variables
initialize_global_variables()

from stability.convex_hull.phase_diagram import (
    entriesGases_A, locked_Chem_Potential_A,
    entriesGases_C, locked_Chem_Potential_C,
    entriesGases_X, locked_Chem_Potential_X
)

# Prepare material entries
(
    all_entries_A,
    all_entries_C,
    all_entries_X,
    TestMat_entry_A,
    TestMat_entry_C,
    TestMat_entry_X
) = prepare_material_entries(api, TestMat_Comp, TestMat_Ener)

# Calculate phase diagrams and energies for each condition
pd_A, energy_per_atom_A, formation_energy_A, energy_above_hull_A = (
    calculate_phase_diagram_condition_A(
        all_entries_A, entriesGases_A, locked_Chem_Potential_A, TestMat_entry_A
    )
)
pd_C, energy_per_atom_C, formation_energy_C, energy_above_hull_C = (
    calculate_phase_diagram_condition_C(
        all_entries_C, entriesGases_C, locked_Chem_Potential_C, TestMat_entry_C
    )
)
pd_X, energy_per_atom_X, formation_energy_X, energy_above_hull_X = (
    calculate_phase_diagram_condition_X(
        all_entries_X, entriesGases_X, locked_Chem_Potential_X, TestMat_entry_X
    )
)

# Print results for condition A (Hydrogen-rich)
print("Condition A (Hydrogen-rich):")
print("Phase diagram:", pd_A)
print("Energy per atom:", energy_per_atom_A)
print("Formation energy:", formation_energy_A)
print("Energy above hull:", energy_above_hull_A)

# Print results for condition C (Oxygen-rich)
print("\nCondition C (Oxygen-rich):")
print("Phase diagram:", pd_C)
print("Energy per atom:", energy_per_atom_C)
print("Formation energy:", formation_energy_C)
print("Energy above hull:", energy_above_hull_C)

# Print results for condition X (CO2-rich)
print("\nCondition X (CO2-rich):")
print("Phase diagram:", pd_X)
print("Energy per atom:", energy_per_atom_X)
print("Formation energy:", formation_energy_X)
print("Energy above hull:", energy_above_hull_X)
