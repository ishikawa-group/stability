#!/usr/bin/env python

# Headers and Imports
from pymatgen.core import Composition, Element
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram
# from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from mp_api.client import MPRester


def initialize_global_variables():
    """
    Initializes global gas entries and chemical potentials for conditions A, C, and X.
    This function sets up the required parameters for each environment (A, C, and X) by
    defining energies, compositions, and entries for relevant gases.
    """
    global H_energy_A, O_energy_A, H_energy_C, O_energy_C, O_energy_X, CO2_Ener_X, CO_energy_X
    global H2_Comp, O2_Comp, H2O_Comp, CO_Comp, CO2_Comp
    global entriesGases_A, entriesGases_C, entriesGases_X
    global chempot_A, chempot_C, chempot_X

    # Define energies for different environments
    H_energy_A = -4.024
    O_energy_A = -8.006

    H_energy_C = -4.997
    O_energy_C = -6.166

    O_energy_X = -6.166
    # O2_Entry_X = ComputedEntry(O2_Comp, O_energy_X * O2_Comp.num_atoms)
    CO2_Ener_X = -25.556
    CO_energy_X = -20.232

    # Define compositions for gases
    H2_Comp  = Composition("H2")
    O2_Comp  = Composition("O2")
    H2O_Comp = Composition("H2O")
    CO_Comp  = Composition("CO")
    CO2_Comp = Composition("X")

    # Define entries for each condition
    entriesGases_A = [
        ComputedEntry(H2_Comp, H_energy_A * H2_Comp.num_atoms),
        ComputedEntry(O2_Comp, O_energy_A * O2_Comp.num_atoms),
        ComputedEntry(H2O_Comp, H_energy_A * H2_Comp.num_atoms + O_energy_A * 0.5 * O2_Comp.num_atoms)
    ]
    entriesGases_C = [
        ComputedEntry(H2_Comp, H_energy_C * H2_Comp.num_atoms),
        ComputedEntry(O2_Comp, O_energy_C * O2_Comp.num_atoms),
        ComputedEntry(H2O_Comp, H_energy_C * H2_Comp.num_atoms + O_energy_C * 0.5 * O2_Comp.num_atoms)
    ]
    entriesGases_X = [
        ComputedEntry(O2_Comp, O_energy_X * O2_Comp.num_atoms),
        ComputedEntry(CO_Comp, CO_energy_X),
        ComputedEntry(CO2_Comp, CO2_Ener_X)
    ]

    chempot_A = {"H": H_energy_A, "O": O_energy_A}
    chempot_C = {"H": H_energy_C, "O": O_energy_C}
    chempot_X = {"O": O_energy_X, "X": CO2_Ener_X}


def prepare_material_entries(api, input_comp, input_energy):
    """
    Prepare material entries and fetch entries from the Materials Project for a given material.

    Args:
        api: API key for the Materials Project.
        input_comp (str): Composition of the material.
        input_energy (float): Energy of the material.

    Returns:
        tuple: List of entries for Condition A, list of entries for Condition C, list of entries for Condition X,
               ComputedEntry object for the material under Condition A,
               ComputedEntry object for the material under Condition C,
               ComputedEntry object for the material under Condition X.
    """
    # Define the material's composition
    input_comp = Composition(input_comp)

    # empirical correction for O atom; see SI S2.2 in the paper
    input_energy += -0.687 * input_comp["O"]

    # Define computed entries for the material under conditions A, C and X
    energy_A = input_energy - O_energy_A * input_comp["O"]
    energy_C = input_energy - O_energy_C * input_comp["O"]
    energy_X = input_energy - O_energy_X * input_comp["O"]

    input_A = ComputedEntry(input_comp, energy_A)
    input_C = ComputedEntry(input_comp, energy_C)
    input_X = ComputedEntry(input_comp, energy_X)

    # Initialize compatibility module
    # compat = MaterialsProject2020Compatibility()
    compat = MaterialsProjectCompatibility()

    # Automatically extract elements in test_mat_comp
    elements = [str(element) for element in input_comp.elements]

    # Initialize MPRester to get material entries from the Materials Project
    with MPRester(api) as api:
        entries_MP_Org_AC = api.get_entries_in_chemsys(elements + ["O", "H"])
        entries_MP_Org_X  = api.get_entries_in_chemsys(elements + ["O", "C"])

    # Process entries with the compatibility module
    entries_MP_Org_AC = compat.process_entries(entries_MP_Org_AC)
    entries_MP_Org_X  = compat.process_entries(entries_MP_Org_X)

    # Combine the fetched entries with the VASP computed entries for Conditions A and C
    entries_VASP_A = [input_A]
    entries_VASP_C = [input_C]

    all_entries_A = entries_MP_Org_AC + entries_VASP_A
    all_entries_C = entries_MP_Org_AC + entries_VASP_C
    entriesTotal_X = entries_MP_Org_X

    return all_entries_A, all_entries_C, entriesTotal_X, input_A, input_C, input_X


def calculate_phase_diagram_condition(all_entries, entriesGases, chempot, input_entry):
    """
    Calculate the phase diagram and energy above the hull for Condition A (Hydrogen-rich environment).

    Args:
        all_entries (list): List of ComputedEntry objects for Condition.
        entriesGases (list): List of ComputedEntry objects for gases under Condition.
        chempot (dict): Locked chemical potentials for Condition.
        input_entry (ComputedEntry): ComputedEntry object for the material under Condition.

    Returns:
        tuple: Phase diagram for Condition A, energy per atom, formation energy, and energy above hull.
    """
    # Define species to eliminate
    eliminate_AC = ["H2", "O2", "H2O"]

    # Filter out these species from all_entries_A
    all_entries = list(filter(lambda e: e.composition.reduced_formula not in eliminate_AC, all_entries))
    all_entries += entriesGases

    pd = GrandPotentialPhaseDiagram(all_entries, chempot)

    num_O = input_entry.composition["O"]
    num_wo_O = input_entry.composition.num_atoms - num_O

    comp = Composition({el: input_entry.composition[el]
                        for el in input_entry.composition.elements if el != Element("O")})

    input_wo_O = ComputedEntry(comp, input_entry.energy + chempot["O"] * num_O)

    e_per_atom = input_entry.energy / input_entry.composition.num_atoms
    e_hull = pd.get_hull_energy(input_wo_O.composition)
    e_hull_per_atom = pd.get_hull_energy_per_atom(input_wo_O.composition)
    e_above_hull_per_atom = input_entry.energy / num_wo_O - e_hull_per_atom

    return pd, e_per_atom, e_hull, e_hull_per_atom, e_above_hull_per_atom


"""
# Condition C
def calculate_phase_diagram_condition_C(all_entries_C, entriesGases_C, chempot_C, test_mat_entry_C):
    # Calculate the phase diagram and energy above the hull for Condition C (Oxygen-rich environment).

    Args:
        all_entries_C (list): List of ComputedEntry objects for Condition C.
        entriesGases_C (list): List of ComputedEntry objects for gases under Condition C.
        chempot_C (dict): Locked chemical potentials for Condition C.
        test_mat_entry_C (ComputedEntry): ComputedEntry object for the material under Condition C.

    Returns:
        tuple: Phase diagram for Condition C, energy per atom, formation energy, and energy above hull.
        
    # Define species to eliminate
    eliminate_AC = ['H2', 'O2', 'H2O']

    # Filter out these species from all_entries_C
    all_entries_C = list(filter(lambda e: e.composition.reduced_formula not in eliminate_AC, all_entries_C))
    all_entries_C += entriesGases_C

    # Create phase diagram for Condition C
    pd_C = GrandPotentialPhaseDiagram(all_entries_C, chempot_C)

    # One should use number of atoms without O, because O atom energy is subtracted from entry
    num_atoms_without_O = test_mat_entry_C.composition.num_atoms - test_mat_entry_C.composition["O"]

    energy_per_atom_C = test_mat_entry_C.energy / num_atoms_without_O
    hull_energy_C = pd_C.get_hull_energy(test_mat_entry_C.composition) / num_atoms_without_O
    energy_above_hull_C = (test_mat_entry_C.energy / num_atoms_without_O
                           - pd_C.get_hull_energy_per_atom(test_mat_entry_C.composition))

    return pd_C, energy_per_atom_C, hull_energy_C, energy_above_hull_C
    
"""


# Condition X
def calculate_phase_diagram_condition_X(entriesTotal_X, entriesGases_X, chempot_X, test_mat_entry_X):
    """
    Calculate the phase diagram and energy above the hull for Condition X (CO2-rich environment).

    Args:
        entriesTotal_X (list): List of ComputedEntry objects for Condition X.
        entriesGases_X (list): List of ComputedEntry objects for gases under Condition X.
        chempot_X (dict): Locked chemical potentials for Condition X.
        test_mat_entry_X (ComputedEntry): ComputedEntry object for the material under Condition X.

    Returns:
        tuple: Phase diagram for Condition X, energy per atom, formation energy, and energy above hull.
    """
    # TODO: We need to split materials by using "split_composition.sh".

    # Filter out CO, X, and O2 from entriesTotal_X
    eliminate_X = ['CO', 'X', 'O2']
    all_entries_X = list(filter(lambda e: e.composition.reduced_formula not in eliminate_X, entriesTotal_X))

    # Add back entries for Condition X
    entries_VASP_X = [test_mat_entry_X]
    all_entries_X += entries_VASP_X + entriesGases_X

    # Create phase diagram for Condition X
    pd_X = GrandPotentialPhaseDiagram(all_entries_X, chempot_X)

    # One should use number of atoms without O, because O atom energy is subtracted from entry
    num_atoms_without_O = test_mat_entry_X.composition.num_atoms - test_mat_entry_X.composition["O"]

    energy_per_atom_X = test_mat_entry_X.energy / num_atoms_without_O
    hull_energy_X = pd_X.get_hull_energy(test_mat_entry_X.composition) / num_atoms_without_O
    energy_above_hull_X = (test_mat_entry_X.energy / num_atoms_without_O
                           - pd_X.get_hull_energy_per_atom(test_mat_entry_C.composition))

    return pd_X, energy_per_atom_X, hull_energy_X, energy_above_hull_X
