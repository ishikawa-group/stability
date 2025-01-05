def get_energy_above_hull(atoms=None, energy=None):
    import sys
    import os
    import warnings
    from stability.convex_hull.phase_diagram import initialize_global_variables
    from stability.convex_hull.phase_diagram import prepare_material_entries
    from stability.convex_hull.phase_diagram import calculate_phase_diagram_condition
    from stability.convex_hull.phase_diagram import calculate_phase_diagram_condition_X

    # Ignore warnings
    warnings.filterwarnings("ignore")

    sys.path.append("../")

    # Set test material and energy
    api = os.getenv("MAPI")
    input_comp = atoms.get_chemical_formula()
    input_energy = energy

    print(f"input_composition: {input_comp}")
    print(f"input_energy: {input_energy:8.6f}")

    # Initialize global variables
    initialize_global_variables()

    from stability.convex_hull.phase_diagram import gas_entries_A, chempot_A
    from stability.convex_hull.phase_diagram import gas_entries_C, chempot_C
    from stability.convex_hull.phase_diagram import gas_entries_X, chempot_X

    # Prepare material entries
    (all_entries_A, all_entries_C, all_entries_X, input_A, input_C, input_X) \
        = prepare_material_entries(api, input_comp, input_energy)

    # --- Calculate phase diagrams and energies for each condition
    # anode
    pd_A, e_per_atom_A, e_hull_A, e_hull_per_atom_A, e_above_hull_per_atom_A \
        = calculate_phase_diagram_condition(all_entries_A, gas_entries_A, chempot_A, input_A)

    # cathode
    pd_C, e_per_atom_C, e_hull_C, e_hull_per_atom_C, e_above_hull_per_atom_C \
        = calculate_phase_diagram_condition(all_entries_C, gas_entries_C, chempot_C, input_C)

    # CO2
    pd_X, e_per_atom_X, e_hull_X, e_hull_per_atom_X, e_above_hull_per_atom_X \
        = calculate_phase_diagram_condition_X(all_entries_X, gas_entries_X, chempot_X, input_X)

    # --- Print results for condition A (Hydrogen-rich)
    print("Condition A (Hydrogen-rich):")
    print(f"Energy per atom: {e_per_atom_A:5.3f} eV/atom")
    print(f"Hull energy per_atom: {e_hull_per_atom_A:5.3f} eV/atom")
    print(f"Energy above hull per atom: {e_above_hull_per_atom_A:5.3f} eV/atom")

    # --- Print results for condition C (Oxygen-rich)
    print("\nCondition C (Oxygen-rich):")
    print(f"Energy per atom: {e_per_atom_C:5.3f} eV/atom")
    print(f"Hull energy per_atom: {e_hull_per_atom_C:5.3f} eV/atom")
    print(f"Energy above hull per atom: {e_above_hull_per_atom_C:5.3f} eV/atom")

    # Print results for condition X (CO2-rich)
    print("\nCondition X (CO2-rich):")
    print(f"Energy per atom: {e_per_atom_X:5.3f} eV/atom")
    print(f"Hull energy per_atom: {e_hull_per_atom_X:5.3f} eV/atom")
    print(f"Energy above hull per atom: {e_above_hull_per_atom_X:5.3f} eV/atom")

    return e_above_hull_per_atom_A, e_above_hull_per_atom_C, e_above_hull_per_atom_X
