def get_energy_above_hull(atoms=None, energy=None):
    import sys
    import os
    import warnings

    # Ignore warnings
    warnings.filterwarnings("ignore")

    sys.path.append("../")

    from stability.convex_hull.phase_diagram import (
        initialize_global_variables,
        prepare_material_entries,
        calculate_phase_diagram_condition,
        calculate_phase_diagram_condition_X,
    )

    # Set test material and energy
    api = os.getenv("MAPI")
    input_comp = atoms.get_chemical_formula()  # "Ba8Zr8O24"
    input_energy = energy

    print(f"input_comp: {input_comp}")
    print(f"input_energy: {input_energy:8.6f}")

    # Initialize global variables
    initialize_global_variables()

    from stability.convex_hull.phase_diagram import (
        entriesGases_A, chempot_A,
        entriesGases_C, chempot_C,
        entriesGases_X, chempot_X,
    )

    # Prepare material entries
    (all_entries_A, all_entries_C, entriesTotal_X, input_A, input_C, input_X) \
        = prepare_material_entries(api, input_comp, input_energy)

    # Calculate phase diagrams and energies for each condition
    pd_A, e_per_atom_A, e_hull_A, e_hull_per_atom_A, e_above_hull_per_atom_A \
        = calculate_phase_diagram_condition(all_entries_A, entriesGases_A, chempot_A, input_A)

    pd_C, e_per_atom_C, e_hull_C, e_hull_per_atom_C, e_above_hull_per_atom_C \
        = calculate_phase_diagram_condition(all_entries_C, entriesGases_C, chempot_C, input_C)

    """
    pd_X, energy_per_atom_X, hull_energy_X, energy_above_hull_X = (
        calculate_phase_diagram_condition_X(
            entriesTotal_X, entriesGases_X, chempot_X, TestMat_entry_X))
    """
    # --- Print results for condition A (Hydrogen-rich)
    # print("Condition A (Hydrogen-rich):")
    # print(f"Energy per atom: {e_per_atom_A:6.4f} eV/atom")
    # print(f"Hull energy per_atom: {e_hull_per_atom_A:6.4f} eV/atom")
    # print(f"Energy above hull per atom: {e_above_hull_per_atom_A:6.4f} eV/atom")

    # --- Print results for condition C (Oxygen-rich)
    # print("\nCondition C (Oxygen-rich):")
    # print(f"Energy per atom: {e_per_atom_C:6.4f} eV/atom")
    # print(f"Hull energy per_atom: {e_hull_per_atom_C:6.4f} eV/atom")
    # print(f"Energy above hull per atom: {e_above_hull_per_atom_C:6.4f} eV/atom")

    """
    # Print results for condition X (CO2-rich)
    print("\nCondition X (CO2-rich):")
    print(f"Energy per atom: {energy_per_atom_X:6.4f} eV/atom")
    print(f"Hull energy per_atom: {hull_energy_X:6.4f} eV/atom")
    print(f"Energy above hull per atom: {energy_above_hull_X:6.4f} eV/atom")
    """

    return e_above_hull_per_atom_A, e_above_hull_per_atom_C
