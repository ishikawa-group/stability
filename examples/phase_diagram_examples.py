"""
Phase diagram examples for the stability package.
This script calculates energy above hull for various materials
under different conditions (anode, cathode, and CO2).
"""

from ase.io import read
from stability.convex_hull.get_energy_above_hull import get_energy_above_hull
import warnings
import os
from pathlib import Path
import yaml

# 警告を無視
warnings.filterwarnings("ignore")

def analyze_material(material):
    """指定された材料のhull energyを解析

    Args:
        material (dict): 材料情報（名前、エネルギー、ファイル名）

    Returns:
        tuple: (anode, cathode, CO2)のhull energy
    """
    structure_file = Path(material["file"])
    if not structure_file.exists():
        print(f"Warning: {structure_file} not found")
        return None

    bulk = read(str(structure_file))
    
    e_above_hull_A, e_above_hull_C, e_above_hull_X = get_energy_above_hull(
        atoms=bulk, calculator=None, energy=material["energy"]
    )
    return e_above_hull_A, e_above_hull_C, e_above_hull_X


if __name__ == "__main__":
    """
    複数の材料に対してhull energyを計算し、結果を表示
    """
    # ディレクトリのパスを取得
    current_dir = Path(__file__).parent
    materials_file = current_dir / "materials.yaml"

    # YAMLファイルから材料データを読み込む
    with open(materials_file, "r") as f:
        data = yaml.safe_load(f)
        materials = data["materials"]

    print("\nEnergy above hull analysis:\n")
    print(f"{'Material':15} {'Anode':>10} {'Cathode':>10} {'CO2':>10}")
    print("-" * 47)

    for material in materials:
        result = analyze_material(material)
        if result:
            e_above_hull_A, e_above_hull_C, e_above_hull_X = result
            print(f"{material['name']:15} {e_above_hull_A:10.6f} {e_above_hull_C:10.6f} {e_above_hull_X:10.6f}")