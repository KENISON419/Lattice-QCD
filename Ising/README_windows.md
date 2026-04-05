# Ising 相転移プロット (Windows PowerShell) — 手順を丁寧に解説

このフォルダのコードで、次の 2 種類のグラフを作成できます。

1. **縦軸 m、横軸 h**（温度を臨界温度 `Tc` の前後 2 本で同時表示）
2. **縦軸 m、横軸 T**（`h=0`）

---

## まず「新しく作成したファイル」が何をしているか

### 1) `run_phase_transition_windows.ps1`
Windows PowerShell から実行する**入口スクリプト**です。以下を順番に自動実行します。

1. C++コードをコンパイルして `ising.x` を作る
2. Python の必要パッケージ（`matplotlib`, `numpy`）を入れる
3. `plot_phase_transition.py` を呼び出して計算と作図を実行する

つまり、**この `.ps1` 1本を実行すれば最終画像まで生成**されます。

### 2) `plot_phase_transition.py`
計算・データ整理・作図を担当する Python スクリプトです。  
やっている内容は次の通りです。

- `ising.x` を複数条件で実行
- 標準出力から `T`, `H`, `mag` を正規表現で抽出
- `m-h`（`T=0.9Tc`, `1.1Tc`）を 1 枚に重ね描き
- `m-T`（`h=0`）を描画
- `.png` 画像と `.dat` データを `results/` に保存

### 3) `README_windows.md`（このファイル）
Windowsユーザー向けに、前提環境・実行手順・トラブル対応をまとめた説明書です。

---

## 前提環境（最初に確認）

- PowerShell 5+ または PowerShell 7+
- C++ コンパイラ（`g++` または `clang++`）
- Python 3

### PowerShell で確認するコマンド

```powershell
$PSVersionTable.PSVersion
g++ --version
python --version
```

> `g++` が無い場合は MSYS2/MinGW などで導入してください。  
> `python` が無い場合は Python 3 をインストールしてください。

---

## 実行方法（コピペで OK）

### 手順 1: フォルダへ移動

```powershell
cd .\Ising
```

### 手順 2: スクリプトを実行

```powershell
powershell -ExecutionPolicy Bypass -File .\run_phase_transition_windows.ps1
```

### 手順 3: 出力を確認

```powershell
dir .\results
```

生成されるファイル:

- `results/m_vs_h_around_Tc.png`
- `results/m_vs_T.png`
- `results/m_vs_h_around_Tc.dat`
- `results/m_vs_T.dat`

---

## 計算条件を変えたい場合

以下の引数を `.ps1` に渡せます。

- `-L`：格子サイズ（2D なら全体サイズは `L x L`）
- `-NTherm`：熱平衡化ステップ数
- `-NSweep`：測定間スイープ数
- `-NMc`：測定回数
- `-PythonExe`：Python 実行コマンド（例: `python`, `py -3`）

### 例（高精度寄り）

```powershell
powershell -ExecutionPolicy Bypass -File .\run_phase_transition_windows.ps1 `
  -L 40 -NTherm 5000 -NSweep 30 -NMc 7000
```

---

## うまく動かない時（よくある原因）

1. **`g++` が見つからない**
   - C++ コンパイラ未導入。MSYS2/MinGW を導入し、PATH を通す。

2. **`pip install` で失敗する**
   - ネットワークや社内プロキシ制限の可能性あり。
   - 先に手動で `matplotlib` と `numpy` を導入してから再実行。

3. **PowerShell の実行ポリシーで止まる**
   - 上記コマンドのように `-ExecutionPolicy Bypass` を付ける。

---

## 物理モデル（使っている式）

- 2次元 Ising モデル（周期境界条件）
- ハミルトニアン  
  \(E = -J \sum_{\langle i,j \rangle} s_i s_j - \mu_0 h \sum_i s_i\)
- 既定値: `J=1`, `mu0=1`
- 臨界温度（2次元厳密解）  
  `Tc = 2J / ln(1+sqrt(2))`
