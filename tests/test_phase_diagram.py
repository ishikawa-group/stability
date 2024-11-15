import sys
import os
sys.path.append("../")

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")

from ase.io import read

from stability.convex_hull.do_vasp_calculation import do_vasp_calculation
from stability.convex_hull.get_convex_hull import get_convex_hull

bulk = read("BaZrO3.cif")
bulk = bulk*[2, 2, 2]

do_vasp = False

if do_vasp:
    do_vasp_calculation(atoms=bulk)

get_convex_hull(atoms=bulk, energy=-333.78748451)
