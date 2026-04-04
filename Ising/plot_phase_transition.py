#!/usr/bin/env python3
"""Run 2D Ising scans and generate phase-transition plots.

If matplotlib is available, writes a PNG directly.
If not, writes data TSV files and a gnuplot script.
"""

from __future__ import annotations

import argparse
import math
import re
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot Ising phase-transition observables")
    p.add_argument("--exe", default="./ising.x", help="Path to Ising executable")
    p.add_argument("--Ls", nargs="+", type=int, default=[16, 24, 32], help="Lattice sizes")
    p.add_argument("--J", type=float, default=1.0)
    p.add_argument("--mu0", type=float, default=1.0)
    p.add_argument("--H", type=float, default=0.0)
    p.add_argument("--T-min", type=float, default=1.6)
    p.add_argument("--T-max", type=float, default=3.2)
    p.add_argument("--n-T", type=int, default=16)
    p.add_argument("--flg-init", type=int, default=1, help="0: cold, 1: hot")
    p.add_argument("--n-therm", type=int, default=4000)
    p.add_argument("--n-sweep", type=int, default=20)
    p.add_argument("--n-mc", type=int, default=4000)
    p.add_argument("--out", default="ising_phase_transition.png", help="Output image path")
    p.add_argument("--data-dir", default="phase_scan_data", help="Directory for exported TSV data")
    return p.parse_args()


def run_scan(exe: str, L: int, args: argparse.Namespace) -> dict[str, list[float]]:
    cmd = [
        exe,
        "-n_dim", "2",
        "-L", str(L),
        "-J", str(args.J),
        "-mu0", str(args.mu0),
        "-flg_init", str(args.flg_init),
        "-n_therm", str(args.n_therm),
        "-n_sweep", str(args.n_sweep),
        "-n_mc", str(args.n_mc),
        "-T", str(args.T_min), str(args.T_max), str(args.n_T),
        "-H", str(args.H), str(args.H), "0",
    ]
    proc = subprocess.run(cmd, check=True, text=True, capture_output=True)

    pattern = re.compile(
        r"T=\s*([-+0-9.eE]+)\s+H=\s*([-+0-9.eE]+)\s+"
        r"mag=\s*([-+0-9.eE]+)\s+mag_err=\s*([-+0-9.eE]+)\s+"
        r"mag_chi=\s*([-+0-9.eE]+)\s+mag_chi_err=\s*([-+0-9.eE]+)\s+"
        r"E=\s*([-+0-9.eE]+)\s+E_err=\s*([-+0-9.eE]+)\s+"
        r"C_heat=\s*([-+0-9.eE]+)\s+C_heat_err=\s*([-+0-9.eE]+)"
    )

    data = {k: [] for k in ["T", "mag", "mag_err", "chi", "chi_err", "E", "E_err", "C", "C_err"]}
    for line in proc.stdout.splitlines():
        m = pattern.search(line)
        if not m:
            continue
        data["T"].append(float(m.group(1)))
        data["mag"].append(float(m.group(3)))
        data["mag_err"].append(float(m.group(4)))
        data["chi"].append(float(m.group(5)))
        data["chi_err"].append(float(m.group(6)))
        data["E"].append(float(m.group(7)))
        data["E_err"].append(float(m.group(8)))
        data["C"].append(float(m.group(9)))
        data["C_err"].append(float(m.group(10)))

    if not data["T"]:
        raise RuntimeError(f"No observables parsed for L={L}.\nOutput:\n{proc.stdout}")
    return data


def write_tsv(scans: dict[int, dict[str, list[float]]], data_dir: Path) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    header = "T\tmag\tmag_err\tchi\tchi_err\tE\tE_err\tC\tC_err\n"
    for L, d in scans.items():
        p = data_dir / f"scan_L{L}.tsv"
        with p.open("w", encoding="utf-8") as f:
            f.write(header)
            for i in range(len(d["T"])):
                f.write(
                    f"{d['T'][i]:.8f}\t{d['mag'][i]:.8f}\t{d['mag_err'][i]:.8f}\t"
                    f"{d['chi'][i]:.8f}\t{d['chi_err'][i]:.8f}\t{d['E'][i]:.8f}\t"
                    f"{d['E_err'][i]:.8f}\t{d['C'][i]:.8f}\t{d['C_err'][i]:.8f}\n"
                )


def write_gnuplot_script(ls: list[int], data_dir: Path, out: Path, tc: float) -> Path:
    gp = data_dir / "plot_phase_transition.gp"
    mag_terms = ", ".join(
        [f"'{data_dir}/scan_L{L}.tsv' using 1:2:3 with yerrorlines title 'L={L}'" for L in ls]
    )
    chi_terms = ", ".join(
        [f"'{data_dir}/scan_L{L}.tsv' using 1:4:5 with yerrorlines title 'L={L}'" for L in ls]
    )
    c_terms = ", ".join(
        [f"'{data_dir}/scan_L{L}.tsv' using 1:8:9 with yerrorlines title 'L={L}'" for L in ls]
    )

    lines = [
        "set terminal pngcairo size 1000,1200",
        f"set output '{out}'",
        "set multiplot layout 3,1 title '2D Ising phase transition (H=0)'",
        "set grid",
        "set key left top",
        f"Tc={tc:.10f}",
        "set ylabel '<|m|>'",
        f"plot {mag_terms}",
        "set arrow from Tc,graph 0 to Tc,graph 1 nohead dt 2 lc rgb 'black'",
        "set ylabel 'chi'",
        f"plot {chi_terms}",
        "set ylabel 'C'",
        "set xlabel 'T'",
        f"plot {c_terms}",
        "unset multiplot",
    ]

    gp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return gp


def try_matplotlib(scans: dict[int, dict[str, list[float]]], ls: list[int], out: Path, tc: float) -> bool:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False

    fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
    for L in ls:
        d = scans[L]
        axs[0].errorbar(d["T"], d["mag"], yerr=d["mag_err"], marker="o", capsize=2, lw=1, label=f"L={L}")
        axs[1].errorbar(d["T"], d["chi"], yerr=d["chi_err"], marker="o", capsize=2, lw=1, label=f"L={L}")
        axs[2].errorbar(d["T"], d["C"], yerr=d["C_err"], marker="o", capsize=2, lw=1, label=f"L={L}")

    axs[0].set_ylabel(r"$\langle |m| \rangle$")
    axs[1].set_ylabel(r"$\chi$")
    axs[2].set_ylabel(r"$C$")
    axs[2].set_xlabel("T")

    for ax in axs:
        ax.axvline(tc, color="k", linestyle="--", alpha=0.5, label=fr"$T_c\approx{tc:.3f}$")
        ax.grid(alpha=0.3)
        ax.legend()

    fig.suptitle("2D Ising phase transition (H=0)")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    print(f"Saved plot: {out}")
    return True


def main() -> None:
    args = parse_args()
    exe = str(Path(args.exe).resolve())
    out = Path(args.out)
    data_dir = Path(args.data_dir)

    scans = {L: run_scan(exe, L, args) for L in args.Ls}
    write_tsv(scans, data_dir)

    tc_2d = 2.0 / math.log(1.0 + math.sqrt(2.0))
    if not try_matplotlib(scans, args.Ls, out, tc_2d):
        gp = write_gnuplot_script(args.Ls, data_dir, out, tc_2d)
        print("matplotlib not found; wrote data + gnuplot script instead.")
        print(f"Data dir: {data_dir}")
        print(f"Gnuplot script: {gp}")
        print(f"Run: gnuplot {gp}")


if __name__ == "__main__":
    main()
