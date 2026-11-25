#!/usr/bin/env python3
"""Generate report figures from the CSV outputs; figures never contain fabricated data."""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "hasil"
FIG = ROOT / "gambar"
FIG.mkdir(exist_ok=True)


def read_csv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def geometry_plot():
    rows = read_csv(OUT / "koefisien_variasi_geometri.csv")
    fig, ax = plt.subplots(figsize=(9, 5.5))
    specs = [
        ("porosity", "Porositas"),
        ("specific_surface_corrected_m_inv", "Luas permukaan spesifik terkoreksi"),
    ]
    for key, label in specs:
        data = [r for r in rows if r["property"] == key]
        x = np.array([float(r["L_voxel"]) for r in data])
        y = np.array([100 * float(r["Cv"]) for r in data])
        line, = ax.plot(x, y, marker="o", linewidth=2, label=label)
        fit = np.polyfit(np.log(x), np.log(y), 1)
        xf = np.geomspace(x.min(), x.max(), 100)
        yf = np.exp(fit[1]) * xf ** fit[0]
        ax.plot(xf, yf, linestyle="--", color=line.get_color(), alpha=.75,
                label=f"Fit {label}: b={fit[0]:.2f}")
        # The −3/2 reference is anchored to the first observed point.
        yref = y[0] * (xf / x[0]) ** (-1.5)
        ax.plot(xf, yref, linestyle=":", color=line.get_color(), alpha=.65,
                label=f"Acuan −3/2 ({label})")
    ax.axhline(5, linestyle="-.", color="black", linewidth=1.2, label="Ambang REV 5%")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Ukuran sisi subvolume, L (voxel)")
    ax.set_ylabel("Koefisien variasi, Cv (%)")
    ax.set_title("Variasi sifat geometri terhadap ukuran subvolume — Kirby")
    ax.grid(True, which="both", alpha=.25)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "cv_geometri.png", dpi=200)
    plt.close(fig)


def permeability_plot():
    p = OUT / "koefisien_variasi_permeabilitas.csv"
    if not p.exists():
        return
    rows = read_csv(p)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for axis in "xyz":
        data = [r for r in rows if r["axis"] == axis and r.get("Cv", "") not in ("", "nan", "None")]
        if not data: continue
        x = np.array([float(r["L_voxel"]) for r in data])
        y = np.array([100 * float(r["Cv"]) for r in data])
        ax.plot(x, y, marker="o", linewidth=2, label=f"Arah {axis.upper()}")
    ax.axhline(5, linestyle="-.", color="black", linewidth=1.2, label="Ambang REV 5%")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Ukuran sisi subvolume, L (voxel)")
    ax.set_ylabel("Koefisien variasi permeabilitas, Cv (%)")
    ax.set_title("Variasi permeabilitas terhadap ukuran subvolume — Kirby")
    ax.grid(True, which="both", alpha=.25); ax.legend()
    fig.tight_layout(); fig.savefig(FIG / "cv_permeabilitas.png", dpi=200); plt.close(fig)


if __name__ == "__main__":
    geometry_plot()
    permeability_plot()
    print("Gambar laporan diperbarui di", FIG)
