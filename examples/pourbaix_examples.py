import os
from stability import PourBaixDiagram
import matplotlib.pyplot as plt

def plot_pourbaix_diagram(elements, title):
    """指定された元素のPourbaix図を生成して表示

    Args:
        elements (list): 元素記号のリスト（例：["Cu"]）
        title (str): プロットのタイトル
    """
    # Pourbaix図エントリーの取得
    entries = PourBaixDiagram.get_pourbaix_entries(elements)
    
    # Pourbaix図の生成
    pbx = PourBaixDiagram(entries)
    
    # プロッターの取得とプロット
    plotter = pbx.get_plotter()
    plt.figure(figsize=(10, 8))
    plot = plotter.get_pourbaix_plot()
    plt.title(title)
    
    # 安定相の取得と表示（pH=7, V=0での例）
    stable_entry = pbx.get_stable_entry(pH=7, voltage=0)
    if stable_entry:
        print(f"\n{title} at pH=7, V=0:")
        print(f"Stable phase: {stable_entry.name}")
    
    return plot

def main():
    """主な金属のPourbaix図を生成する例"""
    # 環境変数のチェック
    if not os.getenv("MAPI"):
        print("Please set the MAPI environment variable with your Materials Project API key")
        return

    # 1. 銅のPourbaix図
    plot_pourbaix_diagram(["Cu"], "Copper Pourbaix Diagram")

    # 2. 鉄のPourbaix図
    plot_pourbaix_diagram(["Fe"], "Iron Pourbaix Diagram")

    # 3. ニッケルのPourbaix図
    plot_pourbaix_diagram(["Ni"], "Nickel Pourbaix Diagram")

    # 4. 複合系の例（Fe-Cr系）
    plot_pourbaix_diagram(["Fe", "Cr"], "Fe-Cr Pourbaix Diagram")

    plt.show()

if __name__ == "__main__":
    main()