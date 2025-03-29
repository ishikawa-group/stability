"""
Stabilityは、材料の安定性を解析するためのPythonパッケージです。

主な機能：
* Convex hull解析
* Pourbaix図の生成と解析

本パッケージは以下のモジュールを提供します：
* convex_hull - 構造安定性の解析
* pourbaix - 電気化学的安定性の解析
"""

__version__ = "0.1.0"

from . import convex_hull
from .convex_hull import (
    do_vasp_calculation,
    get_energy_above_hull,
    initialize_global_variables,
    prepare_material_entries,
    calculate_phase_diagram,
    calculate_phase_diagram_CO2
)

try:
    from . import pourbaix
    from .pourbaix.pourbaix import PourBaixDiagram
    __all__ = [
        "convex_hull",
        "pourbaix",
        "PourBaixDiagram",
        "do_vasp_calculation",
        "get_energy_above_hull",
        "initialize_global_variables",
        "prepare_material_entries",
        "calculate_phase_diagram",
        "calculate_phase_diagram_CO2"
    ]
except (ImportError, ValueError):
    __all__ = [
        "convex_hull",
        "do_vasp_calculation",
        "get_energy_above_hull",
        "initialize_global_variables",
        "prepare_material_entries",
        "calculate_phase_diagram",
        "calculate_phase_diagram_CO2"
    ]
