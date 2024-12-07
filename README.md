# stability
* Assemble of codes related to the stability analysis.

## convex_hull
* Given the material composition and the total (internal) energy, 
  following properties are calculated.
  1. formation energy
  2. energy above the convex hull

### Usage
```python
from ase.io import read

from stability.convex_hull.do_vasp_calculation import do_vasp_calculation
rom stability.convex_hull.get_energy_above_hull import get_energy_above_hull

bulk = read("BaZrO3.cif")
bulk = bulk*[2, 2, 2]

energy = -333.78748451  # VASP-calculateion should be done beforehand
e_above_hull_A, e_above_hull_C = get_energy_above_hull(atoms=bulk, energy=energy)
print(f"energy above hull at anode/cathode conditions (per atom, eV): {e_above_hull_A:5.3f},{e_above_hull_C:5.3f}")
```

## pourbaix
* Drawing the Pourbaix diagram (pH vs. V) using ab-initio etc. calculations.
