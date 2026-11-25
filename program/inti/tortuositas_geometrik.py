#!/usr/bin/env python3
"""Compute 6-neighbour geometric tortuosity on one representative REV-scale Kirby subvolume.
Uses the same through-pore connectivity definition as the course P11 notebook.
"""
from pathlib import Path
import csv, importlib.util, json, sys, time
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("rev_anisotropi", ROOT / "program/inti/rev_anisotropi.py")
rev = importlib.util.module_from_spec(spec); sys.modules[spec.name] = rev; spec.loader.exec_module(rev)


def graph_for_axis(active, axis):
    shape = active.shape
    coords = np.argwhere(active)
    n = len(coords)
    idx = -np.ones(shape, dtype=np.int32)
    idx[tuple(coords.T)] = np.arange(n, dtype=np.int32)
    rows, cols = [], []
    for direction in range(3):
        q = coords.copy()
        q[:, direction] += 1
        valid = (q[:, axis] >= 0) & (q[:, axis] < shape[axis])
        for a in range(3):
            if a != axis:
                q[:, a] %= shape[a]
        q = q[valid]; p = coords[valid]
        connected = active[tuple(q.T)]
        p = p[connected]; q = q[connected]
        rows.append(idx[tuple(p.T)])
        cols.append(idx[tuple(q.T)])
    row = np.concatenate(rows); col = np.concatenate(cols)
    data = np.ones(row.size, dtype=np.float32)
    return sparse.coo_matrix((data, (row, col)), shape=(n, n)).tocsr(), coords, idx


def main():
    pore, meta, phi = rev.load_image()
    L = 128
    rng = np.random.default_rng(rev.SEED + 101)
    starts = rev.rand_subvolumes(pore.shape, L, rev.N_FLOW, rng)
    start = starts[0]
    sub = rev.extract(pore, start, L)
    results = []
    for axis, name in enumerate("xyz"):
        t0 = time.time()
        active = rev.through_pore_cluster(sub, axis)
        if not active.any():
            results.append(dict(L_voxel=L, L_um=L*2.25, x0=int(start[0]), y0=int(start[1]), z0=int(start[2]),
                                axis=name, through_nodes=0, source_nodes=0, target_nodes=0,
                                tau_geometric_6n=None, status="no through path", elapsed_s=time.time()-t0))
            continue
        G, coords, idx = graph_for_axis(active, axis)
        inlet = np.zeros(active.shape, dtype=bool); outlet = np.zeros(active.shape, dtype=bool)
        s0 = [slice(None)]*3; s1 = [slice(None)]*3
        s0[axis] = 0; s1[axis] = -1
        inlet[tuple(s0)] = True; outlet[tuple(s1)] = True
        sources = idx[inlet & active]; targets = idx[outlet & active]
        sources = sources[sources >= 0]; targets = targets[targets >= 0]
        dist = csgraph.dijkstra(G, directed=False, indices=sources, min_only=True)
        d = dist[targets]
        d = d[np.isfinite(d)]
        tau = float(np.mean(d)/(L-1)) if d.size else float("nan")
        results.append(dict(L_voxel=L, L_um=L*2.25, x0=int(start[0]), y0=int(start[1]), z0=int(start[2]),
                            axis=name, through_nodes=int(active.sum()), source_nodes=int(len(sources)),
                            target_nodes=int(len(targets)), tau_geometric_6n=tau,
                            status="ok" if np.isfinite(tau) else "no finite geodesic", elapsed_s=time.time()-t0))
        print(name, results[-1], flush=True)
    out = ROOT / "hasil/tortuositas_geometrik.csv"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0])); w.writeheader(); w.writerows(results)
    (ROOT / "hasil/tortuositas_metadata.json").write_text(json.dumps({
        "method":"Dijkstra shortest path, 6 face-neighbours; periodic transverse axes; no periodicity along flow axis",
        "sample":"Kirby", "sample_size_voxel":L, "seed":rev.SEED+101,
        "subvolume_start_xyz":start.tolist(), "voxel_size_um":2.25,
        "note":"Geometric tortuosity for one largest-tested subvolume; not flow-weighted tortuosity."}, indent=2))
    print("Saved", out)

if __name__ == "__main__":
    main()
