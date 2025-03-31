def get_energy_above_hull(atoms=None, calculator=None, energy=None):
    import sys
    import os
    import warnings
    from stability.convex_hull.phase_diagram import initialize_global_variables
    from stability.convex_hull.phase_diagram import prepare_material_entries
    from stability.convex_hull.phase_diagram import calculate_phase_diagram
    from stability.convex_hull.phase_diagram import calculate_phase_diagram_CO2

    # Ignore warnings
    warnings.filterwarnings("ignore")

    # Set test material and energy
    api = os.getenv("MAPI")
    input_comp = atoms.get_chemical_formula()

    # calculator energy here
    if calculator is None:
        if energy is not None:
            input_energy = energy
        else:
            raise ValueError("Error: either calculator or energy should be given.")
    else:
        if "vasp" in calculator:
            raise ValueError("vasp not implemented")
        elif "m3gnet" in calculator:
            import matgl
            from matgl.ext.ase import PESCalculator

            potential = matgl.load_model("M3GNet-MP-2021.2.8-PES")
            atoms.calc = PESCalculator(potential=potential)
            input_energy = atoms.get_potential_energy()

        input_energy = atoms.get_potential_energy()

    print(f"input_composition: {input_comp}")
    print(f"input_energy: {input_energy:8.6f}")

    # Initialize global variables
    initialize_global_variables()

    from stability.convex_hull.phase_diagram import gas_entries_A, chempot_A
    from stability.convex_hull.phase_diagram import gas_entries_C, chempot_C
    from stability.convex_hull.phase_diagram import gas_entries_X, chempot_X

    # Prepare material entries
    (all_entries_A, all_entries_C, all_entries_X, input_A, input_C, input_X) \
        = prepare_material_entries(api, input_comp, input_energy, correct_O_energy=True)

    # --- Calculate phase diagrams and energies for each condition
    # anode
    pd_A, e_per_atom_A, e_hull_A, e_hull_per_atom_A, e_above_hull_per_atom_A \
        = calculate_phase_diagram(all_entries_A, gas_entries_A, chempot_A, input_A)

    # cathode
    pd_C, e_per_atom_C, e_hull_C, e_hull_per_atom_C, e_above_hull_per_atom_C \
        = calculate_phase_diagram(all_entries_C, gas_entries_C, chempot_C, input_C)

    # CO2
    pd_X, e_per_atom_X, e_hull_X, e_hull_per_atom_X, e_above_hull_per_atom_X \
        = calculate_phase_diagram_CO2(all_entries_X, gas_entries_X, chempot_X, input_X)

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
