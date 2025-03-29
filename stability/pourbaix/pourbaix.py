import os
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter

class PourBaixDiagram:
    def __init__(self, entries):
        """PourBaixDiagramクラスの初期化

        Args:
            entries: Pourbaix図のエントリーリスト
        """
        self.diagram = PourbaixDiagram(entries)
        self._plotter = None

    @staticmethod
    def get_pourbaix_entries(elements):
        """指定された元素のPourbaix図エントリーを取得

        Args:
            elements (List[str]): 元素記号のリスト（例：["Cu"]）

        Returns:
            List: Pourbaix図のエントリーリスト
        """
        API = os.getenv("MAPI")
        if API is None:
            raise ValueError("MAPI environment variable is not set")

        mpr = MPRester(API)
        return mpr.get_pourbaix_entries(elements)

    def get_plotter(self):
        """Pourbaix図のプロッターを取得

        Returns:
            PourbaixPlotter: プロッターオブジェクト
        """
        if self._plotter is None:
            self._plotter = PourbaixPlotter(self.diagram)
        return self._plotter

    def show_plot(self):
        """Pourbaix図を表示"""
        plotter = self.get_plotter()
        plotter.get_pourbaix_plot().show()

    def get_stable_entry(self, pH, voltage):
        """指定されたpHと電位での安定相を取得

        Args:
            pH (float): pH値
            voltage (float): 電位（V vs. SHE）

        Returns:
            Entry: 安定相のエントリー
        """
        return self.diagram.get_stable_entry(pH, voltage)
