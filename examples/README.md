# Examples

このディレクトリには、`stability`パッケージの使用例が含まれています。

## Pourbaix diagram examples

`pourbaix_examples.py`は、Materials Project APIを使用してPourbaix図を生成する例を示しています。

### セットアップ

1. Materials Project APIキーの取得
   * [Materials Project](https://materialsproject.org/)にアカウントを作成
   * [Dashboard](https://materialsproject.org/dashboard)からAPIキーを取得

2. 環境変数の設定
   ```bash
   export MAPI="your-api-key-here"
   ```

3. 必要なパッケージのインストール
   ```bash
   pip install matplotlib mpcontribs-client
   ```

### 実行方法

```bash
python examples/pourbaix_examples.py
```

### 機能

* 単一元素のPourbaix図生成（Cu, Fe, Ni）
* 複合系のPourbaix図生成（Fe-Cr系）
* pH=7, V=0での安定相の表示
* プロット結果の自動表示

### 出力例

* 各元素のPourbaix図が順番に表示されます
* コンソールには各条件での安定相が出力されます

## Phase diagram examples

`phase_diagram_examples.py`は、異なる材料のconvex hull解析を実行する例を示しています。

### 実行方法

```bash
python examples/phase_diagram_examples.py
```

### 機能

* 複数の材料（Ba-Sr-Nb-O系、Cs-Ba-Nb-O系など）のhull energy解析
* 3つの条件（anode、cathode、CO2）での安定性評価
* 結果の表形式での表示

### 出力例

```
Energy above hull analysis:

Material         Anode   Cathode       CO2
-----------------------------------------------
Ba1Sr7Nb8O24  0.382793  0.842793  0.846662
Ba1Sr7V8O24   0.212304  0.528835  0.528835
Ca1La7Fe8O24  0.260604  0.190561  0.204530
...
```