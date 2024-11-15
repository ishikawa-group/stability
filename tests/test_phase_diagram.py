import sys
import os
import warnings

from stability.convex_hull import do_vasp_calculation
from stability.convex_hull import get_convex_hull

from ase.io import read

# Ignore warnings
warnings.filterwarnings("ignore")

sys.path.append("../")

bulk = read("BaZrO3.cif")
bulk = bulk*[2, 2, 2]

do_vasp_calculation(atoms=bulk)
