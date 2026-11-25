#!/usr/bin/env python3
"""Pemeriksaan ulang domain acuan dan ketergantungan spasial subvolume."""
from pathlib import Path
import numpy as np
import pandas as pd
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "program"))
import rev_anisotropi as rev

def main():
    pore, meta, phi_total = rev.load_image()
    print(f"Ukuran citra: {pore.shape}; porositas total: {phi_total:.8f}")
    print("Porositas P3 (tanpa periodisitas):")
    for axis, name in enumerate("xyz"):
        part = rev.porosity_partition_p3(pore, axis)
        print(f"  {name}: phi_eff={part['effective'].mean():.8f}; "
              f"phi_flow={part['flow'].mean():.8f}; "
              f"phi_buntu={part['dead_end'].mean():.8f}; "
              f"phi_isolasi={part['isolated'].mean():.8f}")
    print("Porositas P11 (syarat periodik untuk domain LBM):")
    for axis, name in enumerate("xyz"):
        active = rev.through_pore_cluster(pore, axis)
        print(f"  {name}: phi_periodic_through={active.mean():.8f}; n={active.sum()}")

    coords = pd.read_csv(ROOT / "hasil/subvolume_coordinates.csv")
    rows = []
    for L, group in coords.groupby("L_voxel"):
        starts = group[["x0", "y0", "z0"]].to_numpy()
        overlaps = []
        for i in range(len(starts)):
            for j in range(i):
                intersection = np.maximum(0, int(L) - np.abs(starts[i] - starts[j])).prod()
                overlaps.append(intersection / float(L ** 3))
        rows.append({"L_voxel": int(L), "n_pairs": len(overlaps),
                     "mean_pairwise_overlap_fraction": float(np.mean(overlaps)),
                     "median_pairwise_overlap_fraction": float(np.median(overlaps))})
    result = pd.DataFrame(rows)
    result.to_csv(ROOT / "hasil/diagnostik_tumpang_tindih.csv", index=False)
    print("Diagnostik tumpang tindih:")
    print(result.to_string(index=False))

if __name__ == "__main__":
    main()
