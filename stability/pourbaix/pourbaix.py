import os
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter

API = os.getenv["MY_API"]

mpr = MPRester(API)
entries = mpr.get_pourbaix_entries(["Cu"])
pbx = PourbaixDiagram(entries)
plotter = PourbaixPlotter(pbx)
plotter.get_pourbaix_plot().show()
