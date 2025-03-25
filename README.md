# Stability analysis tools
* A toolset for materials stability calculations, focusing on energy analysis under different environmental conditions.
* This package provides advanced functionality for analyzing material stability through formation energy calculations, convex hull analysis, and phase diagram construction.

## Features
### 1. Formation Energy Calculations
- Computes formation energies from first-principles calculations
- Supports VASP integration for energy calculations
- Handles multiple atomic species and complex chemical compositions

### 2. Energy Above Hull Calculations
* Analyzes stability under three distinct environmental conditions:
- **Hydrogen-rich condition (A)**: Suitable for analyzing materials in reducing environments
- **Oxygen-rich condition (C)**: For oxidizing conditions analysis

### 3. Phase Diagram Analysis
- Constructs grand potential phase diagrams
- Handles chemical potential variations
- Supports multiple component systems
- Integrates with Materials Project database for comprehensive phase analysis

## Installation

```bash
pip install git+https://github.com/ishikawa-group/stability.git
```

## Dependencies
### Core Dependencies
- **pymatgen**: For materials analysis and manipulation
- **ase**: Atomic Simulation Environment for structure handling
- **mp-api**: Materials Project API interface
- **numpy**: Numerical computations

### Additional Requirements
- Materials Project API key (set as MAPI environment variable)
- VASP (optional, for DFT calculations)

## Detailed Usage
### 1. Energy Above Hull Calculation
* Calculate the energy above hull for a material under different environmental conditions:

```python
from stability.convex_hull import get_energy_above_hull
from ase.io import read

# Load structure from VASP format
atoms = read("structure.vasp")

# Specify total energy from DFT calculations
energy = -100.0  # eV

# Calculate energy above hull for all conditions
e_above_hull_A, e_above_hull_C, e_above_hull_X = get_energy_above_hull(
    atoms=atoms,
    energy=energy
)

# Print results
print(f"Energy above hull (H-rich): {e_above_hull_A:5.3f} eV/atom")
print(f"Energy above hull (O-rich): {e_above_hull_C:5.3f} eV/atom") 
print(f"Energy above hull (CO2-rich): {e_above_hull_X:5.3f} eV/atom")
```

The energy above hull values indicate:
- Positive values: Material is metastable
- Zero or negative values: Material is thermodynamically stable
- Larger positive values: Higher degree of instability

### 2. Phase Diagram Construction
* Create and analyze phase diagrams under specific conditions:

```python
from stability.convex_hull import (
    initialize_global_variables,
    prepare_material_entries,
    calculate_phase_diagram_condition
)
import os

# Materials Project API key must be set in environment
api = os.getenv("MAPI")

# Initialize global variables for chemical potentials and gas entries
initialize_global_variables()

# Example composition and energy
composition = "Ba8Zr8O24"
energy = -212.11205  # eV

# Prepare entries for all conditions
entries_A, entries_C, input_A, input_C = prepare_material_entries(
    api=api,
    input_comp=composition,
    input_energy=energy
)

# Calculate phase diagram for hydrogen-rich condition
pd_A, e_per_atom_A, e_hull_A, e_hull_per_atom_A, e_above_hull_per_atom_A = calculate_phase_diagram_condition(
    all_entries=entries_A,
    gas_entries=gas_entries_A,
    chempot=chempot_A,
    input_entry=input_A
)
```

* The phase diagram calculation:
  - Fetches relevant entries from Materials Project
  - Applies appropriate chemical potential corrections
  - Constructs grand potential phase diagram
  - Calculates stability metrics

## Technical Details
### Chemical Potential Definitions
- **Condition A (H-rich)**: μ(H) = -4.024 eV, μ(O) = -8.006 eV
- **Condition C (O-rich)**: μ(H) = -4.997 eV, μ(O) = -6.166 eV
