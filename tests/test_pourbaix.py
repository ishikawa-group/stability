import os
import pytest
from stability.pourbaix import PourBaixDiagram

# APIキーが設定されているかどうかをチェックする
def check_api_key():
    return os.getenv("MAPI") is not None

# APIキーが設定されていない場合はスキップする
@pytest.mark.skipif(not check_api_key(), reason="Please supply an API key. See https://materialsproject.org/api for details.")
def test_pourbaix_diagram():
    """Pourbaix図の生成と操作のテスト"""
    # 銅のPourbaix図エントリーの取得
    entries = PourBaixDiagram.get_pourbaix_entries(["Cu"])
    assert entries is not None, "Pourbaix図エントリーの取得に失敗"
    assert len(entries) > 0, "Pourbaix図エントリーが空"

    # Pourbaix図の生成
    pbx = PourBaixDiagram(entries)
    assert pbx is not None, "PourBaixDiagramの初期化に失敗"

    # プロッターの生成と図の取得
    plotter = pbx.get_plotter()
    assert plotter is not None, "プロッターの取得に失敗"

    # 安定相の取得テスト（pH=7, V=0での安定相）
    stable_entry = pbx.get_stable_entry(pH=7, voltage=0)
    assert stable_entry is not None, "安定相の取得に失敗"

def test_api_key_error():
    """APIキーが設定されていない場合のエラーテスト"""
    os.environ["MAPI"] = ""  # 一時的にAPIキーを空に設定
    with pytest.raises(ValueError, match="Please supply an API key. See https://materialsproject.org/api for details."):
        PourBaixDiagram.get_pourbaix_entries(["Cu"])
    # テスト後に環境変数を元に戻す（もし設定されていた場合）
    original_api = os.getenv("MAPI")
    if original_api:
        os.environ["MAPI"] = original_api