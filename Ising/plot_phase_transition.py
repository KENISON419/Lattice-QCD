#!/usr/bin/env python3
import argparse
import math
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

LINE_RE = re.compile(
    r"T=\s*([\-0-9.eE+]+)\s+H=\s*([\-0-9.eE+]+)\s+"
    r"mag=\s*([\-0-9.eE+]+)\s+mag_err=\s*([\-0-9.eE+]+)\s+"
    r"abs_mag=\s*([\-0-9.eE+]+)\s+abs_mag_err=\s*([\-0-9.eE+]+)\s+"
    r"mag_chi=\s*([\-0-9.eE+]+)\s+mag_chi_err=\s*([\-0-9.eE+]+)"
)
THERM_RE = re.compile(
    r"THERM\s+T=\s*([\-0-9.eE+]+)\s+H=\s*([\-0-9.eE+]+)\s+step=\s*(\d+)\s+m=\s*([\-0-9.eE+]+)"
)


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
    therm_rows = []
    for line in out.splitlines():
        m = LINE_RE.search(line)
        if m:
            rows.append({
                "T": float(m.group(1)),
                "H": float(m.group(2)),
                "mag": float(m.group(3)),
                "mag_err": float(m.group(4)),
                "abs_mag": float(m.group(5)),
                "abs_mag_err": float(m.group(6)),
                "chi": float(m.group(7)),
                "chi_err": float(m.group(8)),
            })
            continue
        t = THERM_RE.search(line)
        if t:
            therm_rows.append({
                "T": float(t.group(1)),
                "H": float(t.group(2)),
                "step": int(t.group(3)),
                "m": float(t.group(4)),
            })
    return rows, therm_rows


def solve_mf_m(T: float, H: np.ndarray, J: float, mu0: float, z: float):
    m = np.zeros_like(H)
    m_prev = -1.0
    for i, h in enumerate(H):
        x = m_prev
        for _ in range(2000):
            x_new = np.tanh((z * J * x + mu0 * h) / T)
            if abs(x_new - x) < 1.0e-12:
                x = x_new
                break
            x = x_new
        m[i] = x
        m_prev = x
    return m


def theory_abs_m_mf(T: np.ndarray, J: float, z: float):
    m = np.zeros_like(T)
    for i, t in enumerate(T):
        if t >= z * J:
            m[i] = 0.0
            continue
        x = 1.0
        for _ in range(2000):
            x_new = np.tanh((z * J * x) / t)
            if abs(x_new - x) < 1.0e-12:
                x = x_new
                break
            x = x_new
        m[i] = abs(x)
    return m


def theory_chi_mf(T: np.ndarray, J: float, mu0: float, z: float):
    chi = np.zeros_like(T)
    for i, t in enumerate(T):
        beta = 1.0 / t
        if t >= z * J:
            chi[i] = beta * (mu0 ** 2) / max(1.0e-12, (1.0 - beta * z * J))
        else:
            # spontaneous solution at h=0
            m = 1.0
            for _ in range(2000):
                m_new = np.tanh((z * J * m) / t)
                if abs(m_new - m) < 1.0e-12:
                    m = m_new
                    break
                m = m_new
            denom = 1.0 - beta * z * J * (1.0 - m * m)
            chi[i] = beta * (mu0 ** 2) * (1.0 - m * m) / max(1.0e-12, denom)
    return chi


def theory_abs_m_2d(T: np.ndarray, J: float):
    Tc = 2.0 * J / math.log(1.0 + math.sqrt(2.0))
    m = np.zeros_like(T)
    mask = T < Tc
    x = np.sinh(2.0 * J / T[mask])
    m[mask] = np.power(1.0 - np.power(x, -4.0), 1.0 / 8.0)
    return m, Tc


def main():
    ap = argparse.ArgumentParser(description="Generate Ising phase transition plots.")
    ap.add_argument("--exe", default="./ising.x", help="Path to executable")
    ap.add_argument("--ndim", type=int, default=2, choices=[1, 2])
    ap.add_argument("--L", type=int, default=32)
    ap.add_argument("--n-therm", type=int, default=3000)
    ap.add_argument("--n-sweep", type=int, default=20)
    ap.add_argument("--n-mc", type=int, default=5000)
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    exe = Path(args.exe)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ndim = args.ndim
    J = 1.0
    mu0 = 1.0
    flg_init = 1
    z = 2.0 * ndim

    if ndim == 2:
        Tc_ref = 2.0 * J / math.log(1.0 + math.sqrt(2.0))
    else:
        Tc_ref = None
    T_low = 1.2 if ndim == 1 else 0.90 * Tc_ref
    T_high = 2.8 if ndim == 1 else 1.10 * Tc_ref
    Ts = [T_low, T_high]

    # m-h curves for T below/above Tc
    h_rows = []
    for T in Ts:
        rows, _ = run_sim(
            exe, ndim, args.L, J, mu0, flg_init,
            args.n_therm, args.n_sweep, args.n_mc,
            T, T, 0,
            -0.6, 0.6, 24,
        )
        h_rows.append((T, rows))

    plt.figure(figsize=(8, 6))
    for T, rows in h_rows:
        H = np.array([r["H"] for r in rows])
        m = np.array([r["mag"] for r in rows])
        m_err = np.array([r["mag_err"] for r in rows])
        plt.errorbar(H, m, yerr=m_err, fmt="none", capsize=2, label=f"MC T={T:.3f}")
        H_dense = np.linspace(H.min(), H.max(), 400)
        plt.plot(H_dense, solve_mf_m(T, H_dense, J, mu0, z), lw=1.0, ls="--", label=f"Mean-field T={T:.3f}")
    plt.axvline(0.0, color="k", ls="--", lw=0.8)
    plt.xlabel("h")
    plt.ylabel("m")
    plt.title("Ising model: m-h curves (jackknife error bars)")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    mh_png = outdir / "m_vs_h.png"
    plt.savefig(mh_png, dpi=160)

    # m-T and chi-T curve at h=0
    t_min = 0.4 if ndim == 1 else 1.2
    t_max = 3.4
    t_rows, _ = run_sim(
        exe, ndim, args.L, J, mu0, flg_init,
        args.n_therm, args.n_sweep, args.n_mc,
        t_min, t_max, 22,
        0.0, 0.0, 0,
    )
    T_arr = np.array([r["T"] for r in t_rows])
    m_arr = np.array([r["abs_mag"] for r in t_rows])
    m_err = np.array([r["abs_mag_err"] for r in t_rows])
    chi_arr = np.array([r["chi"] for r in t_rows])
    chi_err = np.array([r["chi_err"] for r in t_rows])

    plt.figure(figsize=(8, 6))
    plt.errorbar(T_arr, m_arr, yerr=m_err, fmt="none", capsize=2, label="MC |m|")
    T_dense = np.linspace(T_arr.min(), T_arr.max(), 500)
    if ndim == 2:
        m_theory, Tc = theory_abs_m_2d(T_dense, J)
        plt.plot(T_dense, m_theory, lw=1.2, ls="--", label="2D exact |m| (h=0)")
        plt.axvline(Tc, color="r", ls="--", lw=1.0, label=f"2D exact Tc={Tc:.3f}")
    else:
        Tc_mf = z * J
        plt.plot(T_dense, theory_abs_m_mf(T_dense, J, z), lw=1.0, ls="--", label="Mean-field |m| (h=0)")
        plt.axvline(Tc_mf, color="r", ls="--", lw=1.0, label=f"Mean-field Tc={Tc_mf:.3f}")
    plt.xlabel("T")
    plt.ylabel("|m|")
    plt.title("Ising model: |m|-T at h=0 (jackknife error bars)")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    mt_png = outdir / "m_vs_T.png"
    plt.savefig(mt_png, dpi=160)

    plt.figure(figsize=(8, 6))
    plt.errorbar(T_arr, chi_arr, yerr=chi_err, fmt="none", capsize=2, label="MC χ")
    chi_theory = theory_chi_mf(T_dense, J, mu0, z)
    Tc_mf = z * J
    plt.plot(T_dense, chi_theory, lw=1.1, ls="--", label="Mean-field χ (h=0)")
    plt.axvline(Tc_mf, color="r", ls="--", lw=1.0, label=f"Mean-field Tc={Tc_mf:.3f}")
    plt.xlabel("T")
    plt.ylabel("χ")
    plt.title("Ising model: χ-T at h=0 (jackknife error bars)")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    chi_png = outdir / "chi_vs_T.png"
    plt.savefig(chi_png, dpi=160)

    therm_target_T = Tc_ref if ndim == 2 else 1.0
    _, therm_rows = run_sim(
        exe, ndim, args.L, J, mu0, flg_init,
        args.n_therm, args.n_sweep, max(100, args.n_mc // 10),
        therm_target_T, therm_target_T, 0,
        0.0, 0.0, 0,
    )
    therm_x = np.array([r["step"] for r in therm_rows], dtype=int)
    therm_m = np.array([r["m"] for r in therm_rows])
    plt.figure(figsize=(8, 6))
    plt.plot(therm_x, therm_m, lw=1.0)
    plt.xlabel("update count during thermalization")
    plt.ylabel("m")
    plt.title("Thermalization check: update count vs m")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    therm_png = outdir / "thermalization_m_vs_update.png"
    plt.savefig(therm_png, dpi=160)

    np.savetxt(
        outdir / "m_vs_h.dat",
        np.array([
            (t, r["H"], r["mag"], r["mag_err"]) for t, rows in h_rows for r in rows
        ]),
        header="T h m m_err",
    )
    np.savetxt(outdir / "m_vs_T.dat", np.column_stack([T_arr, m_arr, m_err]), header="T |m| |m|_err")
    np.savetxt(outdir / "chi_vs_T.dat", np.column_stack([T_arr, chi_arr, chi_err]), header="T chi chi_err")
    np.savetxt(outdir / "thermalization_m_vs_update.dat", np.column_stack([therm_x, therm_m]), header="update m")

    print(f"Saved: {mh_png}")
    print(f"Saved: {mt_png}")
    print(f"Saved: {chi_png}")
    print(f"Saved: {therm_png}")


if __name__ == "__main__":
    main()
