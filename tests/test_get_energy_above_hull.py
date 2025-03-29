import pytest
from ase.io import read
from stability import get_energy_above_hull
import numpy as np

# テストデータ
TEST_STRUCTURE = "tests/Rb1Ba7Nb8O24.vasp"
TEST_ENERGY = -339.378272

@pytest.fixture
def test_atoms():
    """テスト用の構造データを読み込む"""
    return read(TEST_STRUCTURE)

def test_input_validation():
    """入力パラメータのバリデーションテスト"""
    with pytest.raises(AttributeError):
        # atomsがNoneの場合
        get_energy_above_hull(atoms=None, energy=TEST_ENERGY)

    atoms = read(TEST_STRUCTURE)
    with pytest.raises(TypeError):
        # energyがNoneの場合
        get_energy_above_hull(atoms=atoms, energy=None)

def test_get_energy_above_hull_basic(test_atoms):
    """基本的な機能のテスト"""
    # 既知のRb1Ba7Nb8O24での計算結果と比較
    e_above_hull_A, e_above_hull_C, e_above_hull_X = get_energy_above_hull(
        atoms=test_atoms,
        energy=TEST_ENERGY
    )

    # 期待される値
    expected_A = 0.216330437
    expected_C = 0.618830437
    expected_X = 0.618830437  # CO2条件はcathode条件と同じ値になる（O_energy_X = O_energy_C）

    # 結果の検証（小数点以下6桁まで）
    np.testing.assert_almost_equal(e_above_hull_A, expected_A, decimal=6)
    np.testing.assert_almost_equal(e_above_hull_C, expected_C, decimal=6)
    np.testing.assert_almost_equal(e_above_hull_X, expected_X, decimal=6)

def test_return_type(test_atoms):
    """戻り値の型とサイズのテスト"""
    result = get_energy_above_hull(atoms=test_atoms, energy=TEST_ENERGY)
    
    # 戻り値はタプルで、長さは3
    assert isinstance(result, tuple)
    assert len(result) == 3
    
    # 各値はfloat型
    for value in result:
        assert isinstance(value, float)

def test_energy_above_hull_positive(test_atoms):
    """hull energyが非負であることの確認"""
    e_above_hull_A, e_above_hull_C, e_above_hull_X = get_energy_above_hull(
        atoms=test_atoms,
        energy=TEST_ENERGY
    )
    
    assert e_above_hull_A >= 0
    assert e_above_hull_C >= 0
    assert e_above_hull_X >= 0

def test_chemical_formula(test_atoms):
    """化学式の取得が正しく動作することの確認"""
    expected_formula = "Ba7Nb8O24Rb"
    assert test_atoms.get_chemical_formula() == expected_formula