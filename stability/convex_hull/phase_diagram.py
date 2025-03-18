#!/usr/bin/env python

# Headers and Imports
from pymatgen.core import Composition, Element, periodic_table
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from mp_api.client import MPRester


def initialize_global_variables():
    """
    Initializes global gas entries and chemical potentials for conditions A, C, and X.
    This function sets up the required parameters for each environment (A, C, and X) by
    defining energies, compositions, and entries for relevant gases.
    """
    global H_energy_A, O_energy_A, H_energy_C, O_energy_C, O_energy_X, CO2_energy_X, CO_energy_X
    global H2_Comp, O2_Comp, H2O_Comp, CO_Comp, CO2_Comp
    global gas_entries_A, gas_entries_C, gas_entries_X
    global chempot_A, chempot_C, chempot_X

    # Define energies for different environments
    H_energy_A = -4.024
    O_energy_A = -8.006

    H_energy_C = -4.997
    O_energy_C = -6.166

    O_energy_X = -8.006
    CO2_energy_X = -24.388
    CO_energy_X = -20.232

    # Define dummy species for X
    X = periodic_table.DummySpecie('X')

    # Define compositions for gases
    H2_Comp = Composition("H2")
    O2_Comp = Composition("O2")
    H2O_Comp = Composition("H2O")
    CO_Comp = Composition("CO")
    CO2_Comp = Composition("X")

    # Define entries for each condition
    gas_entries_A = [
        ComputedEntry(H2_Comp, H_energy_A * H2_Comp.num_atoms),
        ComputedEntry(O2_Comp, O_energy_A * O2_Comp.num_atoms),
        ComputedEntry(H2O_Comp, H_energy_A * H2_Comp.num_atoms + O_energy_A * 0.5 * O2_Comp.num_atoms)
    ]
    gas_entries_C = [
        ComputedEntry(H2_Comp, H_energy_C * H2_Comp.num_atoms),
        ComputedEntry(O2_Comp, O_energy_C * O2_Comp.num_atoms),
        ComputedEntry(H2O_Comp, H_energy_C * H2_Comp.num_atoms + O_energy_C * 0.5 * O2_Comp.num_atoms)
    ]
    gas_entries_X = [
        ComputedEntry(O2_Comp, O_energy_X * O2_Comp.num_atoms),
        ComputedEntry(CO_Comp, CO_energy_X),
        ComputedEntry(CO2_Comp, CO2_energy_X)
    ]

    chempot_A = {"H": H_energy_A, "O": O_energy_A}
    chempot_C = {"H": H_energy_C, "O": O_energy_C}
    # chempot_X = {"X": CO2_energy_X, "O2": O_energy_X*2}
    chempot_X = {"X": CO2_energy_X, "O": O_energy_X}


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
    # input_energy += -0.687 * input_comp["O"]

    # Define computed entries for the material under conditions A, C and X
    energy_A = input_energy - O_energy_A * input_comp["O"]
    energy_C = input_energy - O_energy_C * input_comp["O"]
    energy_X = input_energy - O_energy_X * input_comp["O"]

    input_A = ComputedEntry(input_comp, energy_A)
    input_C = ComputedEntry(input_comp, energy_C)
    input_X = ComputedEntry(input_comp, energy_X)

    # Initialize compatibility module
    compat = MaterialsProject2020Compatibility()

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
    all_entries_X = entries_MP_Org_X

    return all_entries_A, all_entries_C, all_entries_X, input_A, input_C, input_X


def calculate_phase_diagram_condition(all_entries, gas_entries, chempot, input_entry):
    """
    Calculate the phase diagram and energy above the hull for Condition A (Hydrogen-rich environment).

    Args:
        all_entries (list): List of ComputedEntry objects for Condition.
        gas_entries (list): List of ComputedEntry objects for gases under Condition.
        chempot (dict): Locked chemical potentials for Condition.
        input_entry (ComputedEntry): ComputedEntry object for the material under Condition.

    Returns:
        tuple: Phase diagram for Condition A, energy per atom, formation energy, and energy above hull.
    """
    # Define species to eliminate
    eliminate_AC = ["H2", "O2", "H2O"]

    # Filter out these species from all_entries_A
    all_entries = list(filter(lambda e: e.composition.reduced_formula not in eliminate_AC, all_entries))
    all_entries += gas_entries

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


def calculate_phase_diagram_condition_X(all_entries=None, gas_entries=None, chempot=None, input_entry=None):
    import re
    import os
    import subprocess
    import collections
    import warnings
    import pathlib

    # Ignore warnings
    warnings.filterwarnings("ignore")

    from pymatgen.core import Composition, periodic_table
    from pymatgen.entries.computed_entries import ComputedEntry
    from pymatgen.analysis.phase_diagram import GrandPotPDEntry, PDEntry, PDPlotter, PhaseDiagram
    from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
    from pathlib import Path

    # Modify composition to account for CO2 (as X)
    current_dir = Path(__file__).parent.absolute()
    split_composition = current_dir / "split_composition.sh"

    for j in range(0, len(all_entries)):
        CurrentEntry = all_entries[j]
        get_Composition = CurrentEntry.composition

        # Modify composition text
        with open('composition.txt', 'w') as out_file:
            originalLine = str(get_Composition)
            fixedLine1 = originalLine.replace('\n', 'Ba8Zr8O24')
            fixedLine2 = re.sub("[A-Za-z]+", lambda ele: " " + ele[0] + " ", fixedLine1)
            print(fixedLine2, file=out_file)

        # Execute external script to split composition
        subprocess.run(split_composition)

        with open('./num_Lines.txt', 'r') as out_file:
            numLinesO = out_file.readlines()
            numLines = numLinesO[0]
            numLines = int(numLines)

        numO = 0
        numC = 0
        holder = [0 for i in range(numLines)]

        with open('composition.txt', 'r') as out_file:
            holderOut = out_file.readlines()
            for k in range(0, numLines):
                holder[k] = holderOut[k]
                holder[k] = holder[k].replace('\n', 'Ba8Zr8O24')

                if holder[k] == 'C':
                    numCs = holderOut[k + 1]
                    numCs = str(numCs.replace('\n', 'Ba8Zr8O24'))
                    numCs = int(numCs)
                elif holder[k] == 'O':
                    numOs = holderOut[k + 1]
                    numOs = str(numOs.replace('\n', 'Ba8Zr8O24'))
                    numOs = int(numOs)

        C_exist_check = 'numC' in locals()
        O_exist_check = 'numO' in locals()

        if not C_exist_check:
            numC = 0

        if not O_exist_check:
            numO = 0

        if numC == 0 and numO == 0:
            all_entries[j] = CurrentEntry
        elif numC > 0 and numO == 0:
            all_entries[j] = CurrentEntry
        elif numC == 0 and numO > 0:
            all_entries[j] = CurrentEntry
        elif numC > 0 and numO < 2:
            all_entries[j] = CurrentEntry
        elif numC > 0 and numO > 0:
            numX = 0
            Mod_Comp = 'Ba'
            if 2 * numC == numO:
                numX = numC
                g = 1
                while g < numLines + 1:
                    if holder[g-1] == 'C':
                        Mod_Comp = Mod_Comp + 'X' + str(numX)
                        g += 2
                    elif holder[g-1] == 'O':
                        g += 2
                    else:
                        Mod_Comp = holder[g-1] + holder[g] + Mod_Comp
                        g += 2
            elif 2 * numC > numO:
                while True:
                    numO -= 2
                    if numO < 0:
                        numO += 2
                        break
                    numX += 1
                    numC -= 1
                g = 1
                while g < numLines + 1:
                    if holder[g-1] == 'C':
                        Mod_Comp = Mod_Comp + 'X' + str(numX)
                        Mod_Comp = Mod_Comp + 'C' + str(numC)
                        g += 2
                    elif holder[g-1] == 'O':
                        if numO == 0:
                            g += 2
                        elif numO > 0:
                            Mod_Comp = Mod_Comp + 'O' + str(numO)
                            g += 2
                    else:
                        Mod_Comps = holder[g-1] + holder[g] + Mod_Comp
                        g += 2
            elif 2 * numC < numO:
                numX = numC
                numO = numO - 2 * numC
                g = 1
                while g < numLines + 1:
                    if holder[g-1] == 'C':
                        Mod_Comp = Mod_Comp + 'X'
                        Mod_Comp = Mod_Comp + str(numC)
                        g += 2
                    elif holder[g-1] == 'O':
                        Mod_Comp = Mod_Comp + 'O'
                        Mod_Comp = Mod_Comp + str(numO)
                        g += 2
                    else:
                        Mod_Comp = holder[g-1] + holder[g] + Mod_Comp
                        g += 2

                energyEntry = CurrentEntry.energy
                make_Mod_Entry = ComputedEntry(Mod_Comp, energyEntry)
                all_entries[j] = make_Mod_Entry

    ################################### CO(Element Z) Checks ###################################################

    for j in range(0, len(all_entries)):

        CurrentEntry = all_entries[j]
        get_Composition = CurrentEntry.composition

        with open('composition.txt', 'w') as out_file:
            originalLine = str(get_Composition)
            fixedLine1 = originalLine.replace('\n', 'Ba8Zr8O24')
            fixedLine2 = re.sub("[A-Za-z]+", lambda ele: " " + ele[0] + " ", fixedLine1)
            print(fixedLine2, file=out_file)

        subprocess.run(split_composition)

        with open('./num_Lines.txt', 'r') as out_file:
            numLinesO = out_file.readlines()
            numLines = numLinesO[0]
            numLines = int(numLines)

        numO = 0
        numC = 0
        holder = [0 for chi in range(numLines)]

        with open('composition.txt', 'r') as out_file:
            holderOut = out_file.readlines()

            for k in range(0, numLines):
                holder[k] = holderOut[k]
                holder[k] = holder[k].replace('\n', 'Ba8Zr8O24')

                if holder[k] == 'C':
                    numCs = holderOut[k + 1]
                    numCs = str(numCs.replace('\n', 'Ba8Zr8O24'))
                    numCs = int(numCs)
                elif holder[k] == 'O':
                    numOs = holderOut[k + 1]
                    numOs = str(numOs.replace('\n', 'Ba8Zr8O24'))
                    numOs = int(numOs)

            C_exist_check = 'numC' in locals()
            O_exist_check = 'numO' in locals()

            if not C_exist_check:
                numC = 0

            if not O_exist_check:
                numO = 0

            if numC > 0 and numO == 0:
                numZ = 0
                Mod_Comp = 'Ba8Zr8O24'
                if numC == numO:
                    numZ = numC
                    g = 1
                    while g < numLines+1:
                        if holder[g-1] == 'C':
                            Mod_Comp = Mod_Comp + 'Z' + str(numZ)
                            g += 2
                        elif holder[g-1] == 'O':
                            g += 2
                        else:
                            Mod_Comp = holder[g-1]+holder[g]+Mod_Comp
                            g += 2
                elif numC > numO:
                    while True:
                        numC -= 1
                        if numC < 0:
                            numC += 1
                            break
                        numZ += 1
                        numO -= 1
                    g = 1
                    while g < numLines+1:
                        if holder[g-1] == 'C':
                            Mod_Comp = Mod_Comp + 'Z' + str(numZ)
                            Mod_Comp = Mod_Comp + 'C' + str(numC)
                            g += 2
                        elif holder[g-1] == 'O':
                            if numO == 0:
                                g += 2
                            elif numO > 0:
                                Mod_Comp = Mod_Comp + 'O' + str(numO)
                                g += 2
                        else:
                            Mod_Comps = holder[g-1]+holder[g]+Mod_Comp
                            g += 2

                    energyEntry = CurrentEntry.energy
                    make_Mod_Entry = ComputedEntry(Mod_Comp, energyEntry)
                    all_entries[j] = make_Mod_Entry
                    print(all_entries[j])

    # --- Gas Phase Correction CO2 ---

    CO_entries = [e for e in all_entries if e.composition.reduced_formula == 'CO']
    CO2_entries = [e for e in all_entries if e.composition.reduced_formula == 'X']
    O2_entries = [e for e in all_entries if e.composition.reduced_formula == 'O2']

    eliminate = ["CO", "X", "O2"]
    all_entries = list(filter(lambda e: e.composition.reduced_formula not in eliminate, all_entries))

    input_entry = ComputedEntry(input_entry.composition, input_entry.energy)

    entries_VASP = [input_entry]
    all_entries = all_entries + entries_VASP + gas_entries

    # --- Calculate Convex Hull under Environmental Conditions ---
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
