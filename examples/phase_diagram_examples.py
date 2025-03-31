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

# 警告を無視
warnings.filterwarnings("ignore")

# 材料データの定義
MATERIALS = [
    {
        "name": "Ba1Sr7Nb8O24",
        "energy": -339.945046,
        "file": "Ba1Sr7Nb8O24.vasp",
    },
    {
        "name": "Ba1Sr7V8O24",
        "energy": -320.169744,
        "file": "Ba1Sr7V8O24.vasp",
    },
    {
        "name": "Ca1La7Fe8O24",
        "energy": -331.534147,
        "file": "Ca1La7Fe8O24.vasp",
    },
    {
        "name": "Cs1Ba7Nb8O24",
        "energy": -338.659999,
        "file": "Cs1Ba7Nb8O24.vasp",
    },
    {
        "name": "Cs1Sr7Nb8O24",
        "energy": -337.02237,
        "file": "Cs1Sr7Nb8O24.vasp",
    },
    {
        "name": "Cs1Sr7V8O24",
        "energy": -315.246376,
        "file": "Cs1Sr7V8O24.vasp",
    },
    {
        "name": "Rb1Ba7Nb8O24",
        "energy": -339.378272,
        "file": "Rb1Ba7Nb8O24.vasp",
    }
]


def analyze_material(material, tests_dir):
    """指定された材料のhull energyを解析

    Args:
        material (dict): 材料情報（名前、エネルギー、ファイル名）
        tests_dir (Path): 構造ファイルのあるディレクトリ

    Returns:
        tuple: (anode, cathode, CO2)のhull energy
    """
    structure_file = tests_dir / material["file"]
    if not structure_file.exists():
        print(f"Warning: {structure_file} not found")
        return None

    bulk = read(str(structure_file))
    e_above_hull_A, e_above_hull_C, e_above_hull_X = get_energy_above_hull(
        atoms=bulk,
        energy=material["energy"]
    )
    return e_above_hull_A, e_above_hull_C, e_above_hull_X


def main():
    """
    複数の材料に対してhull energyを計算し、結果を表示
    """
    # testsディレクトリのパスを取得
    current_dir = Path(__file__).parent
    tests_dir = current_dir.parent / "tests"

    print("\nEnergy above hull analysis:\n")
    print(f"{'Material':15} {'Anode':>10} {'Cathode':>10} {'CO2':>10}")
    print("-" * 47)

    for material in MATERIALS:
        result = analyze_material(material, tests_dir)
        if result:
            e_above_hull_A, e_above_hull_C, e_above_hull_X = result
            print(f"{material['name']:15} {e_above_hull_A:10.6f} {e_above_hull_C:10.6f} {e_above_hull_X:10.6f}")


if __name__ == "__main__":
    main()