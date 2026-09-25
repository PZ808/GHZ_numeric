"""Plot the SN-seeded Kerr m-mode radial-decay diagnostic; requires matplotlib."""
import csv
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ghz_mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader((ROOT / "Data/gsn_kerr_m2_decay_n13.csv").open()))
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
for col, component in enumerate(("nn", "nm", "mm")):
    data = [row for row in rows if row["component"] == component]
    r = np.array([float(row["r"]) for row in data])
    adjusted = np.array([float(row["adjusted"]) for row in data])
    for key, label, color, style in (
        ("reconstructed", r"$\mathcal{S}^{\dagger}\Phi+\mathrm{c.c.}$", "#6d87b6", "-"),
        ("lie", r"$\mathcal{L}_{\zeta}g$", "#de9b57", "--"),
        ("adjusted", "Sum", "#166b52", "-"),
    ):
        axes[0, col].loglog(r, [float(row[key]) for row in data], style,
                           color=color, label=label, lw=2)
    axes[0, col].set_title(r"$h_{" + component + r"}$")
    axes[0, col].grid(alpha=.18, which="both")
    axes[1, col].semilogx(r, r * adjusted, color="#166b52", lw=2)
    axes[1, col].grid(alpha=.18, which="both")
    axes[1, col].set_xlabel(r"$r/M$")
    i, j = 20, 40
    slope = np.log(adjusted[j]/adjusted[i]) / np.log(r[j]/r[i])
    axes[0, col].text(.05, .06, f"slope (200–2000 M): {slope:.4f}",
                      transform=axes[0, col].transAxes, fontsize=10)
axes[0, 0].set_ylabel("Max magnitude on angular grid")
axes[1, 0].set_ylabel(r"$r\,\max_z |h^{\mathrm{adjusted}}|$")
axes[0, 0].legend(fontsize=9)
fig.suptitle(r"Circular Kerr orbit: $a/M=0.5$, $r_0/M=10$, $m=2$" + "\n"
             + r"SN source $\ell=2\ldots8$; angular m-mode solve, $N_z=13$; $M=\mu=1$", fontsize=14)
folder = ROOT / "plots/gsn_kerr"
folder.mkdir(parents=True, exist_ok=True)
fig.savefig(folder / "radial_decay.png", dpi=180)
fig.savefig(folder / "radial_decay.pdf")
print(folder / "radial_decay.png")
