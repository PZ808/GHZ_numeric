"""Make presentation figures from the C++ Schwarzschild Bondi diagnostics.

Run ``build/bondi_plot_data tests/data/psi0_schwarzschild_r0_10_lmax20_mostly_minus.csv
plots/bondi_schwarzschild`` first, then pass that output directory to this script.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "bondi_mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "bondi_xdg_cache"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def save(fig: plt.Figure, output: Path) -> None:
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def redshift_figure(data_dir: Path) -> None:
    data = rows(data_dir / "redshift_modes.csv")
    ells = sorted({int(row["ell"]) for row in data})
    ms = list(range(-max(ells), max(ells) + 1))
    values = np.full((len(ells), len(ms)), np.nan)
    helical = np.full_like(values, np.nan)
    detuned = np.full_like(values, np.nan)
    for row in data:
        i, j = ells.index(int(row["ell"])), ms.index(int(row["m"]))
        values[i, j] = float(row["abs_u_dot_zeta_helical"])
        helical[i, j] = float(row["abs_lie_uu_helical"])
        detuned[i, j] = float(row["abs_lie_uu_detuned"])

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.1), constrained_layout=True)
    fig.suptitle(r"Schwarzschild circular orbit: gauge contribution to $h_{uu}$",
                 fontsize=16, fontweight="semibold")
    cmap = plt.colormaps["viridis"].copy()
    cmap.set_bad("#e9e9e9")

    for ax, array, title, color_label in (
        (axes[0], values, r"Gauge-vector contraction $|u\cdot\zeta|$",
         r"$\log_{10}|u\cdot\zeta|$"),
        (axes[2], detuned, r"Control: $\omega=m\Omega+0.004$",
         r"$\log_{10}|(\mathcal{L}_\zeta g)_{uu}|$"),
    ):
        positive = np.log10(np.where(array > 0, array, np.nan))
        image = ax.imshow(positive, origin="lower", aspect="auto", cmap=cmap)
        fig.colorbar(image, ax=ax, shrink=0.78, label=color_label)
        ax.set_title(title, fontsize=12)

    ax = axes[1]
    helical_log = np.where(np.isfinite(helical),
                           np.log10(np.maximum(helical, 1e-22)), np.nan)
    image = ax.imshow(helical_log, origin="lower", aspect="auto", cmap=cmap,
                      vmin=-22, vmax=-18)
    fig.colorbar(image, ax=ax, shrink=0.78,
                 label=r"$\log_{10}\max(|(\mathcal{L}_\zeta g)_{uu}|,10^{-22})$")
    ax.set_title(r"Helical: $\omega=m\Omega$", fontsize=12)
    maximum = np.nanmax(helical)
    ax.text(0.5, 0.96, f"max = {maximum:.2e} across {len(data)} modes",
            transform=ax.transAxes, ha="center", va="top", fontsize=9,
            color="white", bbox={"facecolor": "#182338", "edgecolor": "none",
                                 "alpha": 0.9, "pad": 4})

    for ax in axes:
        ax.set_xticks(np.arange(len(ms)), ms, fontsize=7)
        ax.set_yticks(np.arange(len(ells)), ells)
        ax.set_xlabel("m")
        ax.set_ylabel(r"$\ell$")
        ax.grid(False)
    fig.supxlabel(
        r"$M=1$, $r_0=10$; $(\mathcal{L}_\zeta g)_{uu}="
        r"2i u^u(m\Omega-\omega)(u\cdot\zeta)e^{i(m\phi-\omega u)}$",
        fontsize=10, y=-0.04,
    )
    save(fig, data_dir / "lie_zeta_redshift_modes")


def m_mode_figure(data_dir: Path) -> None:
    grid = rows(data_dir / "m_mode_grid.csv")
    summary = {int(row["m"]): row for row in rows(data_dir / "m_mode_summary.csv")}
    modes = sorted(summary)
    fig, axes = plt.subplots(2, len(modes), figsize=(15, 7.1), sharex="col",
                             constrained_layout=True, height_ratios=[1.15, 1])
    fig.suptitle(r"Held collocation in $m$ modes vs. analytic Schwarzschild $\ell$ sum",
                 fontsize=16, fontweight="semibold")
    colors = ["#155e75", "#a34a2a", "#6950a1"]
    for column, m in enumerate(modes):
        subset = sorted((row for row in grid if int(row["m"]) == m),
                        key=lambda row: float(row["z"]))
        z = np.array([float(row["z"]) for row in subset])
        numeric = np.array([complex(float(row["seed_numeric_real"]),
                                    float(row["seed_numeric_imag"])) for row in subset])
        exact = np.array([complex(float(row["seed_exact_real"]),
                                  float(row["seed_exact_imag"])) for row in subset])
        scale = np.max(np.abs(exact))
        error = np.abs(numeric - exact) / scale
        color = colors[column % len(colors)]
        top, bottom = axes[:, column]
        top.plot(z, np.abs(exact), color=color, lw=2, label="analytic ℓ sum")
        top.scatter(z, np.abs(numeric), facecolors="white", edgecolors=color,
                    s=28, zorder=3, label="m-mode collocation")
        top.set_title(rf"$m={m}$, $N_z={summary[m]['nz']}$", fontsize=12)
        top.set_ylabel(r"$|\widetilde\Phi_0^H|$")
        top.grid(alpha=0.25)
        if column == 0:
            top.legend(loc="best", frameon=False, fontsize=9)
        bottom.semilogy(z, np.maximum(error, 1e-17), "o-", color=color,
                        markersize=3.5, lw=1.5)
        bottom.set_ylim(1e-17, max(1e-6, np.max(error) * 20))
        bottom.set_xlabel(r"$z=\cos\theta$")
        bottom.set_ylabel("complex seed error / max |exact|")
        bottom.grid(alpha=0.25, which="both")
        bottom.text(0.04, 0.95,
                    "max error = {:.2e}\noperator residual = {:.2e}".format(
                        float(summary[m]["relative_seed_error"]),
                        float(summary[m]["relative_exact_equation_residual"])),
                    transform=bottom.transAxes, va="top", fontsize=9,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85})
    fig.supxlabel(r"Input: fitted, $\ell$-summed $\psi_0$ coefficients; mostly-minus signature",
                 fontsize=10, y=-0.03)
    save(fig, data_dir / "m_mode_held_collocation")


def metric_falloff_figure(data_dir: Path) -> None:
    data = rows(data_dir / "metric_falloff.csv")
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True,
                             constrained_layout=True)
    fig.suptitle(r"Schwarzschild adjustment: $2\,\mathrm{Re}\,S^\dagger"
                 r"(\Phi^{\mathrm{Adj}})+\mathcal{L}_{\zeta}g$",
                 fontsize=16, fontweight="semibold")
    for i, m in enumerate((2, 3)):
        for j, component in enumerate(("nn", "nm", "mm")):
            subset = sorted((row for row in data if int(row["m"]) == m
                             and row["component"] == component),
                            key=lambda row: float(row["r"]))
            r = np.array([float(row["r"]) for row in subset])
            reconstructed = np.array([float(row["abs_reconstructed"]) for row in subset])
            lie = np.array([float(row["abs_lie_zeta"]) for row in subset])
            adjusted = np.array([float(row["abs_direct_sum"]) for row in subset])
            leading = float(subset[0]["abs_adjusted_coefficient_1_over_r"])
            ax = axes[i, j]
            ax.loglog(r, reconstructed, color="#a34a2a", lw=1.7,
                      label=r"$|2\,\mathrm{Re}\,S^\dagger|$")
            ax.loglog(r, lie, color="#386c91", lw=1.7, linestyle="--",
                      label=r"$|\mathcal{L}_{\zeta}g|$")
            ax.loglog(r, adjusted, color="#142e40", lw=2.5,
                      label="|sum| (direct evaluation)")
            ax.loglog(r, leading / r, color="#788a24", lw=1.2,
                      linestyle=":", label=r"$|A_{-1}|/r$")
            slope = np.polyfit(np.log(r[-20:]), np.log(adjusted[-20:]), 1)[0]
            ax.text(0.04, 0.07, f"large-r slope {slope:.3f}",
                    transform=ax.transAxes, fontsize=9,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85})
            ax.set_title(rf"$m={m}$, $h_{{{component}}}$", fontsize=12)
            ax.set_ylabel("component magnitude")
            ax.set_xlabel(r"$r/M$")
            ax.grid(alpha=0.22, which="both")
            if i == 0 and j == 0:
                ax.legend(loc="upper right", frameon=False, fontsize=8)
    fig.supxlabel(r"Fixed $z=0.37$, $u=\phi=0$; fitted $\psi_0$ data summed to $\ell=20$.",
                 y=-0.025, fontsize=10)
    save(fig, data_dir / "adjusted_metric_radial_falloff")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", type=Path)
    args = parser.parse_args()
    redshift_figure(args.data_dir)
    m_mode_figure(args.data_dir)
    if (args.data_dir / "metric_falloff.csv").exists():
        metric_falloff_figure(args.data_dir)


if __name__ == "__main__":
    main()
