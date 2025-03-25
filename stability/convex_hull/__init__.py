"""
Convex hull analysis tools for stability calculations.

This module provides functionality for:
1. Formation energy calculations
2. Energy above hull calculations
"""

from .do_vasp_calculation import do_vasp_calculation
from .get_energy_above_hull import get_energy_above_hull
from .phase_diagram import (initialize_global_variables, prepare_material_entries,
                            calculate_phase_diagram, calculate_phase_diagram_CO2)

__all__ = [
    "do_vasp_calculation",
    "get_energy_above_hull",
    "initialize_global_variables", "prepare_material_entries", "calculate_phase_diagram",
    "calculate_phase_diagram_CO2"
]
