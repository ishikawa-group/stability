# stability
* Assemble of codes related to the stability analysis.

## convex_hull
* Given the material composition and the total (internal) energy, following properties are calculated.
  1. formation energy
  2. energy above the convex hull

### Usage
```python
from ase.io import read

from stability.convex_hull.do_vasp_calculation import do_vasp_calculation
from stability.convex_hull.get_convex_hull import get_convex_hull

bulk = read("BaZrO3.cif")
bulk = bulk*[2, 2, 2]

energy = -333.78748451  # VASP-calculateion should be done beforehand
get_convex_hull(atoms=bulk, energy=energy)
```

## pourbaix
* Drawing the Pourbaix diagram (pH vs. V) using ab-initio etc. calculations.
