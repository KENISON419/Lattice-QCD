# Ising相転移の実行手順（手取り足取り版）

このドキュメントは、`Ising/plot_phase_transition.py` を使って **実際に相転移グラフを作る** ための具体的な手順です。

---

## 0. どこで作業するか

```bash
cd /workspace/Lattice-QCD/Ising
pwd
```

`pwd` の表示が `/workspace/Lattice-QCD/Ising` ならOKです。

---

## 1. 実行ファイルをビルドする

> このリポジトリの makefile は Intel コンパイラ設定が初期値なので、ここでは `g++/gcc` を明示します。

```bash
make clean
make CXX=g++ CC=gcc
```

成功すると `ising.x` が生成されます。

確認:

```bash
ls -l ising.x
```

---

## 2. （推奨）まずは短い試運転

いきなり重い計算を回す前に、軽い設定で動作確認します。

```bash
python3 plot_phase_transition.py \
  --exe ./ising.x \
  --Ls 12 16 \
  --T-min 1.8 --T-max 3.0 --n-T 6 \
  --n-therm 800 --n-sweep 8 --n-mc 400 \
  --out phase_transition_smoke.png
```

### 実行後の見方

- `matplotlib` が入っていれば `phase_transition_smoke.png` が生成されます。
- `matplotlib` がない場合は、
  - `phase_scan_data/scan_L*.tsv`
  - `phase_scan_data/plot_phase_transition.gp`
  が生成されます（gnuplot用）。

---

## 3. 本番設定で相転移を描く（2D Ising, H=0）

```bash
python3 plot_phase_transition.py \
  --exe ./ising.x \
  --Ls 16 24 32 \
  --T-min 1.6 --T-max 3.2 --n-T 16 \
  --n-therm 4000 --n-sweep 20 --n-mc 4000 \
  --out ising_phase_transition.png
```

この設定で、一般に次が見えてきます。

- `|m|`：低温で大きく、高温で小さくなる
- `chi`（帯磁率）：臨界近傍でピーク
- `C`（比熱）：臨界近傍でピーク

---

## 4. matplotlib が無い場合（gnuplotで描画）

`plot_phase_transition.py` 実行後に `phase_scan_data/plot_phase_transition.gp` が作られていれば、

```bash
gnuplot phase_scan_data/plot_phase_transition.gp
```

で `--out` で指定した PNG が生成されます。

---

## 5. よくある詰まりポイント

### 5-1. `ModuleNotFoundError: No module named 'matplotlib'`

- 対処A: matplotlibを入れる
- 対処B: そのまま gnuplot フォールバックを使う（このスクリプトは自動でTSV+gpを書き出します）

### 5-2. `gnuplot: command not found`

- gnuplot 未導入です。matplotlib ルートを使うか、gnuplot を導入してください。

### 5-3. 実行が遅い

以下を下げると速くなります（精度は下がります）。

- `--n-therm`
- `--n-sweep`
- `--n-mc`
- `--Ls`（格子サイズ）

---

## 6. まずはこの3コマンドだけ覚えればOK

```bash
cd /workspace/Lattice-QCD/Ising
make clean && make CXX=g++ CC=gcc
python3 plot_phase_transition.py --exe ./ising.x --Ls 16 24 32 --T-min 1.6 --T-max 3.2 --n-T 16 --n-therm 4000 --n-sweep 20 --n-mc 4000 --out ising_phase_transition.png
```

これで「相転移を示す図」を作るところまで進めます。
