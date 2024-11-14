import sys
import os
import warnings

# Ignore warnings
warnings.filterwarnings("ignore")

sys.path.append("../")
from ase.io import read
from ase.calculators.vasp import Vasp

cif = "BaZrO3.cif"
atoms = read(cif)
atoms *= [2, 2, 2]

# Set up VASP calculator with standard settings
tmpdir = "tmpdir_convex_hull"
atoms.calc = Vasp(prec="normal", xc="pbe", ispin=2, lorbit=10,
                  ibrion=2, nsw=10, isif=8,
                  encut=520, ediff=1e-6, algo="Normal", nelm=50, nelmin=5,
                  kpts=[4, 4, 4], kgamma=True,
                  ismear=0, sigma=0.05,
                  lwave=False, lcharg=False,
                  npar=4, nsim=4,
                  directory=tmpdir,
                  lreal=False,
                  )

# Calculate total energy (needed for band structure calculation)
energy = atoms.get_potential_energy()

print(f"Energy: {energy} eV")

quit()

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
# TestMat_Ener = -331.28931146  # VASP-calculated value in eV
TestMat_Ener = energy

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
