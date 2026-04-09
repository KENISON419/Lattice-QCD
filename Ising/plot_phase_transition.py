#!/usr/bin/env python3
import argparse
import math
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

LINE_RE = re.compile(r"T=\s*([\-0-9.eE+]+)\s+H=\s*([\-0-9.eE+]+)\s+mag=\s*([\-0-9.eE+]+)")


def run_sim(exe: Path, ndim: int, L: int, J: float, mu0: float, flg_init: int,
            n_therm: int, n_sweep: int, n_mc: int,
            T_i: float, T_f: float, n_T: int,
            H_i: float, H_f: float, n_H: int):
    cmd = [
        str(exe),
        "-n_dim", str(ndim),
        "-L", str(L),
        "-J", str(J),
        "-mu0", str(mu0),
        "-flg_init", str(flg_init),
        "-n_therm", str(n_therm),
        "-n_sweep", str(n_sweep),
        "-n_mc", str(n_mc),
        "-T", str(T_i), str(T_f), str(n_T),
        "-H", str(H_i), str(H_f), str(n_H),
    ]
    out = subprocess.check_output(cmd, text=True)

    rows = []
    for line in out.splitlines():
        m = LINE_RE.search(line)
        if not m:
            continue
        rows.append((float(m.group(1)), float(m.group(2)), float(m.group(3))))
    return rows


def main():
    ap = argparse.ArgumentParser(description="Generate Ising phase transition plots.")
    ap.add_argument("--exe", default="./ising.x", help="Path to executable")
    ap.add_argument("--L", type=int, default=32)
    ap.add_argument("--n-therm", type=int, default=3000)
    ap.add_argument("--n-sweep", type=int, default=20)
    ap.add_argument("--n-mc", type=int, default=5000)
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    exe = Path(args.exe)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ndim = 2
    J = 1.0
    mu0 = 1.0
    flg_init = 1

    Tc = 2.0 * J / math.log(1.0 + math.sqrt(2.0))
    Ts = [0.90 * Tc, 1.10 * Tc]

    # m-h curves for T below/above Tc
    h_rows = []
    for T in Ts:
        rows = run_sim(
            exe, ndim, args.L, J, mu0, flg_init,
            args.n_therm, args.n_sweep, args.n_mc,
            T, T, 0,
            -0.6, 0.6, 24,
        )
        h_rows.append((T, rows))

    plt.figure(figsize=(8, 6))
    for T, rows in h_rows:
        H = [r[1] for r in rows]
        m = [r[2] for r in rows]
        plt.plot(H, m, marker="o", ms=3, lw=1.3, label=f"T={T:.3f}")
    plt.axvline(0.0, color="k", ls="--", lw=0.8)
    plt.xlabel("h")
    plt.ylabel("m")
    plt.title("Ising model: m-h curves around Tc")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    mh_png = outdir / "m_vs_h_around_Tc.png"
    plt.savefig(mh_png, dpi=160)

    # m-T curve at h=0
    t_rows = run_sim(
        exe, ndim, args.L, J, mu0, flg_init,
        args.n_therm, args.n_sweep, args.n_mc,
        1.2, 3.4, 22,
        0.0, 0.0, 0,
    )
    T_arr = np.array([r[0] for r in t_rows])
    m_arr = np.abs(np.array([r[2] for r in t_rows]))

    plt.figure(figsize=(8, 6))
    plt.plot(T_arr, m_arr, marker="o", ms=3, lw=1.3)
    plt.axvline(Tc, color="r", ls="--", lw=1.0, label=f"Tc={Tc:.3f}")
    plt.xlabel("T")
    plt.ylabel("m")
    plt.title("Ising model: m-T at h=0")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    mt_png = outdir / "m_vs_T.png"
    plt.savefig(mt_png, dpi=160)

    np.savetxt(
        outdir / "m_vs_h_around_Tc.dat",
        np.array([(t, h, m) for t, rows in h_rows for (_, h, m) in rows]),
        header="T h m",
    )
    np.savetxt(outdir / "m_vs_T.dat", np.column_stack([T_arr, m_arr]), header="T |m|")

    print(f"Saved: {mh_png}")
    print(f"Saved: {mt_png}")


if __name__ == "__main__":
    main()
