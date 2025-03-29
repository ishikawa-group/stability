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
from . import pourbaix
from .convex_hull import get_energy_above_hull

__all__ = [
    "convex_hull",
    "pourbaix",
    "get_energy_above_hull"
]
