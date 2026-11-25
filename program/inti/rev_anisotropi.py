#!/usr/bin/env python3
"""REV and anisotropy analysis for Kirby sandstone (course project).

The script is deliberately staged: run the cheap geometric statistics first, then
run the optional flow solver. All random coordinates and all numeric results are
persisted. Raw image convention follows the course bundle: uint8, C order, 0=pore.
"""
from __future__ import annotations
import argparse, csv, json, math, os, sys, time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = ROOT / "hasil"
FIG = ROOT / "gambar"
SIZES = [32, 48, 64, 96, 128]
VOXEL_M = 2.25e-6
SEED = 20261125
N_GEOM = 30
N_FLOW = 6

def load_image():
    meta = json.loads((DATA/"metadata.json").read_text())
    raw = DATA/meta["raw_file"]
    a = np.fromfile(raw, dtype=np.uint8)
    expected = int(np.prod(meta["size_voxels_xyz"]))
    if a.size != expected: raise ValueError(f"Raw size mismatch: {a.size} != {expected}")
    pore = (a.reshape(tuple(meta["size_voxels_xyz"]), order=meta["order"]) == 0)
    phi = float(pore.mean())
    ref = json.loads((DATA/"metadata.json").read_text())
    return pore, meta, phi

def periodic_components(mask, axes):
    """Label connected components while joining opposite faces on selected axes."""
    lab, n = ndi.label(mask, structure=ndi.generate_binary_structure(3, 1))
    if n == 0:
        return lab
    parent = np.arange(n + 1, dtype=np.int64)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for ax in axes:
        lo = np.take(lab, 0, axis=ax)
        hi = np.take(lab, -1, axis=ax)
        valid = (lo > 0) & (hi > 0)
        pairs = np.unique(np.column_stack((lo[valid], hi[valid])), axis=0)
        for x, y in pairs:
            rx, ry = find(int(x)), find(int(y))
            if rx != ry:
                parent[rx] = ry
    roots = np.arange(n + 1)
    for i in range(1, n + 1):
        roots[i] = find(i)
    return roots[lab]


def through_pore_cluster(pore, axis):
    """Return all pore components periodic transversely and spanning the flow axis.

    This follows the course P11 notebook definition. All spanning clusters are kept;
    choosing only the largest would incorrectly discard parallel flow paths.
    """
    transverse = [a for a in range(3) if a != axis]
    lab_trans = periodic_components(pore, transverse)
    first = set(np.unique(np.take(lab_trans, 0, axis=axis))) - {0}
    last = set(np.unique(np.take(lab_trans, -1, axis=axis))) - {0}
    spanning_labels = first & last
    if not spanning_labels:
        return np.zeros_like(pore, dtype=bool)
    selected = np.isin(lab_trans, list(spanning_labels))
    # Exact P11 convention: use components periodic on all three axes to
    # reunite transverse-spanning paths that may wrap around the periodic box.
    lab_periodic = periodic_components(pore, [0, 1, 2])
    periodic_ids = set(np.unique(lab_periodic[selected])) - {0}
    if not periodic_ids:
        return np.zeros_like(pore, dtype=bool)
    return np.isin(lab_periodic, list(periodic_ids))



def porosity_partition_p3(pore, axis=0):
    """P3 definition: 6-connected components without periodic face joining.

    phi_tot  : all pore voxels / bulk voxels
    phi_eff  : pore voxels that touch either inlet or outlet face
    phi_flow : pore voxels in components touching both faces
    """
    lab, n = ndi.label(pore)  # default 3-D connectivity: 6 face neighbours
    if n == 0:
        zero = np.zeros_like(pore, dtype=bool)
        return {"effective": zero, "flow": zero, "dead_end": zero,
                "isolated": zero, "n_components": 0}
    front = set(np.unique(np.take(lab, 0, axis=axis))) - {0}
    back = set(np.unique(np.take(lab, -1, axis=axis))) - {0}
    through_ids = front & back
    touching_ids = front | back
    all_ids = set(range(1, n + 1))
    isolated_ids = all_ids - touching_ids
    dead_end_ids = touching_ids - through_ids
    through = np.isin(lab, list(through_ids)) if through_ids else np.zeros_like(pore, bool)
    effective = pore & ~np.isin(lab, list(isolated_ids))
    dead_end = np.isin(lab, list(dead_end_ids)) if dead_end_ids else np.zeros_like(pore, bool)
    isolated = np.isin(lab, list(isolated_ids)) if isolated_ids else np.zeros_like(pore, bool)
    return {"effective": effective, "flow": through, "dead_end": dead_end,
            "isolated": isolated, "n_components": int(n)}


def domain_reference():
    """Write the three P3 porosity definitions and compare them with course CSV."""
    pore, meta, phi_total = load_image()
    rows = []
    for axis, name in enumerate("xyz"):
        part = porosity_partition_p3(pore, axis)
        rows.append({
            "axis": name,
            "phi_total": phi_total,
            "phi_effective": float(part["effective"].mean()),
            "phi_flow": float(part["flow"].mean()),
            "phi_dead_end": float(part["dead_end"].mean()),
            "phi_isolated": float(part["isolated"].mean()),
            "n_components": part["n_components"],
        })
    save_csv(OUT / "porositas_domain_penuh.csv", rows)
    reference = list(csv.DictReader(open(DATA / "referensi_permeabilitas.csv", encoding="utf-8")))
    ref = [r for r in reference if r["batuan"] == "Kirby"]
    domain = {
        "sample": "Kirby", "shape_voxel": list(pore.shape),
        "voxel_size_um": 2.25, "domain_side_um": 256 * 2.25,
        "porosity_definitions": "P3: 6-connected components without periodic face joining",
        "porosity_by_axis": rows,
        "course_csv_reference": [{"axis": r["arah"], "phi_total": float(r["phi"]),
                                   "phi_connected_reference": float(r["phi_terhubung"]),
                                   "K_mD": float(r["K_mD"]), "converged": r["konvergen"]}
                                  for r in ref],
        "note": "Untuk domain acuan, phi_flow mengikuti definisi P3. Saat LBM subvolume, penyaringan pori tembus menggunakan kondisi periodik P11; kedua takrif dibedakan.",
    }
    (OUT / "domain_acuan.json").write_text(json.dumps(domain, indent=2))
    print(json.dumps(domain, indent=2))
    return domain

def largest_percolating_cluster(pore, axis):
    """Compatibility alias: all topologically through-flow pore clusters."""
    return through_pore_cluster(pore, axis)

def surface_area_specific(full_pore, start, L):
    """Specific pore-solid area from voxel faces, including interfaces at crop boundaries.

    Internal faces are counted once. At each crop face, the adjacent voxel from the
    parent image is used as a one-voxel halo (when available), avoiding both omission
    of true boundary-crossing interfaces and artificial periodic crop interfaces.
    """
    x0, y0, z0 = map(int, start)
    sub = extract(full_pore, start, L)
    faces = 0
    for ax, lo in enumerate((x0, y0, z0)):
        # Interfaces wholly inside the sampled cube.
        a = [slice(None)] * 3; b = [slice(None)] * 3
        a[ax] = slice(0, -1); b[ax] = slice(1, None)
        faces += int(np.count_nonzero(sub[tuple(a)] != sub[tuple(b)]))
        # Negative crop boundary, if a parent voxel exists there.
        if lo > 0:
            boundary = [slice(None)] * 3; boundary[ax] = 0
            outside = [slice(x0, x0 + L), slice(y0, y0 + L), slice(z0, z0 + L)]
            outside[ax] = lo - 1
            faces += int(np.count_nonzero(sub[tuple(boundary)] != full_pore[tuple(outside)]))
        # Positive crop boundary, if a parent voxel exists there.
        hi = lo + L
        if hi < full_pore.shape[ax]:
            boundary = [slice(None)] * 3; boundary[ax] = -1
            outside = [slice(x0, x0 + L), slice(y0, y0 + L), slice(z0, z0 + L)]
            outside[ax] = hi
            faces += int(np.count_nonzero(sub[tuple(boundary)] != full_pore[tuple(outside)]))
    voxel_area = VOXEL_M**2
    bulk_volume = sub.size * VOXEL_M**3
    raw = (faces * voxel_area) / bulk_volume
    return raw, raw * (2.0 / 3.0)

def rand_subvolumes(shape, L, n, rng):
    starts = rng.integers(0, np.array(shape)-L+1, size=(n,3))
    # reject duplicate starts
    starts=np.unique(starts,axis=0)
    while len(starts)<n:
        extra=rng.integers(0,np.array(shape)-L+1,size=(n-len(starts),3))
        starts=np.unique(np.vstack([starts,extra]),axis=0)
    return starts[:n]

def extract(a, start, L):
    x,y,z=map(int,start); return a[x:x+L,y:y+L,z:z+L]

def cv_stats(vals):
    x=np.asarray(vals,float); x=x[np.isfinite(x)]
    if x.size==0: return (float("nan"),)*4
    mu=float(np.mean(x)); sd=float(np.std(x,ddof=0)); cv=sd/abs(mu) if mu else float("nan")
    return mu,sd,cv,int(x.size)

def save_csv(path, rows, fields=None):
    if not rows: return
    fields=fields or list(rows[0])
    with open(path,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)

def physical_checks_pass(row):
    """True only when all P11-style diagnostic fields are present and pass."""
    try:
        checks=json.loads(row.get("physical_checks_json") or "{}")
    except (TypeError, json.JSONDecodeError):
        return False
    return bool(checks) and all(bool(v) for v in checks.values())

def geometry():
    pore,meta,phi_full=load_image()
    rng=np.random.default_rng(SEED)
    pos_rows=[]; geom_rows=[]
    for L in SIZES:
        starts=rand_subvolumes(pore.shape,L,N_GEOM,rng)
        for j,st in enumerate(starts):
            sub=extract(pore,st,L)
            phi=float(sub.mean())
            raw,corrected=surface_area_specific(pore,st,L)
            pos_rows.append(dict(L_voxel=L,L_um=L*2.25,replicate=j,
                                 x0=int(st[0]),y0=int(st[1]),z0=int(st[2])))
            geom_rows.append(dict(L_voxel=L,L_um=L*2.25,replicate=j,
                                  porosity=phi,specific_surface_raw_m_inv=raw,
                                  specific_surface_corrected_m_inv=corrected))
    save_csv(OUT/"subvolume_coordinates.csv",pos_rows)
    save_csv(OUT/"geometri_subvolume.csv",geom_rows)
    # Ringkasan statistik per ukuran. Nilai luas permukaan mentah dan terkoreksi
    # disimpan bersama agar faktor koreksi tangga voxel dapat diaudit.
    summary=[]
    wide=[]
    for L in SIZES:
        rr=[r for r in geom_rows if r["L_voxel"]==L]
        row={"L_voxel":L,"L_um":L*2.25,"n_subvolume":len(rr)}
        for key,label in [("porosity","porosity"),
                          ("specific_surface_raw_m_inv","surface_raw_m_inv"),
                          ("specific_surface_corrected_m_inv","surface_corrected_m_inv")]:
            mu,sd,cv,n=cv_stats([r[key] for r in rr])
            row[f"{label}_mean"]=mu; row[f"{label}_std"]=sd; row[f"{label}_Cv"]=cv
            if label in ("porosity","surface_corrected_m_inv"):
                summary.append(dict(property=key,L_voxel=L,L_um=L*2.25,mean=mu,std=sd,Cv=cv,n=n,
                                    L_REV_threshold_5pct=(cv<.05)))
        wide.append(row)
    save_csv(OUT/"ringkasan_statistik_geometri.csv",wide)
    slopes={}
    for key in ["porosity","specific_surface_corrected_m_inv"]:
        sub=[r for r in summary if r["property"]==key and r["Cv"]>0]
        slopes[key]=float(np.polyfit(np.log([r["L_voxel"] for r in sub]),
                                    np.log([r["Cv"] for r in sub]),1)[0]) if len(sub)>=2 else None
    save_csv(OUT/"koefisien_variasi_geometri.csv",summary)
    metaout={"sample":"Kirby","seed":SEED,"n_geometry_per_size":N_GEOM,
             "sizes_voxel":SIZES,"sizes_um":[L*2.25 for L in SIZES],
             "full_image_porosity":phi_full,"slope_logCv_logL":slopes,
             "REV":{"porosity":None,"specific_surface":None},
             "surface_reference":"Bulk volume; voxel-face count including crop-boundary interfaces via parent-image halo; 2/3 staircase bias correction"}
    for prop,short in [("porosity","porosity"),("specific_surface_corrected_m_inv","specific_surface")]:
        rows=[r for r in summary if r["property"]==prop]
        hits=[r for r in rows if r["Cv"]<.05]
        metaout["REV"][short]=({"L_voxel":hits[0]["L_voxel"],"L_um":hits[0]["L_um"]} if hits else
                               {"status":"not_reached","largest_L_voxel":SIZES[-1],"Cv_at_largest":rows[-1]["Cv"]})
    (OUT/"ringkasan_geometri.json").write_text(json.dumps(metaout,indent=2))
    print(json.dumps(metaout,indent=2))
    return pore,metaout

# --- LBM D3Q19 TRT solver: periodic transverse directions, pressure/body-force drive.
def lattice():
    c=[[0,0,0]]
    for a in range(3):
        for s in (1,-1):
            v=[0,0,0]; v[a]=s; c.append(v)
    for a in range(3):
        for b in range(a+1,3):
            for s1 in (1,-1):
                for s2 in (1,-1):
                    v=[0,0,0]; v[a]=s1; v[b]=s2; c.append(v)
    c=np.asarray(c,dtype=np.int8)
    w=np.array([1/3]+[1/18]*6+[1/36]*12,float)
    opp=np.array([np.where(np.all(c==-v,axis=1))[0][0] for v in c])
    return c,w,opp

from numba import njit

@njit
def lbm_step_global(f, fb, nbr, c, w, opp, force, wp, wm):
    Q,n=f.shape
    for j in range(n):
        rho=0.0
        for qi in range(Q): rho+=f[qi,j]
        u0=u1=u2=0.0
        for qi in range(Q):
            u0+=c[qi,0]*f[qi,j]; u1+=c[qi,1]*f[qi,j]; u2+=c[qi,2]*f[qi,j]
        u0=(u0+0.5*force[0])/rho; u1=(u1+0.5*force[1])/rho; u2=(u2+0.5*force[2])/rho
        uu=u0*u0+u1*u1+u2*u2
        uF=u0*force[0]+u1*force[1]+u2*force[2]
        for qi in range(Q):
            o=opp[qi]
            if o<qi: continue
            cu=c[qi,0]*u0+c[qi,1]*u1+c[qi,2]*u2
            cF=c[qi,0]*force[0]+c[qi,1]*force[1]+c[qi,2]*force[2]
            ei=w[qi]*rho*(1+3*cu+4.5*cu*cu-1.5*uu)
            eo=w[o]*rho*(1-3*cu+4.5*cu*cu-1.5*uu)
            fp=0.5*(f[qi,j]+f[o,j]); fm=0.5*(f[qi,j]-f[o,j])
            ep=0.5*(ei+eo); em=0.5*(ei-eo)
            Sg=w[qi]*(-3*uF+9*cu*cF); Sa=w[qi]*3*cF
            gp=-wp*(fp-ep)+(1-0.5*wp)*Sg
            gm=-wm*(fm-em)+(1-0.5*wm)*Sa
            fb[qi,j]=f[qi,j]+gp+gm
            if o!=qi: fb[o,j]=f[o,j]+gp-gm
    for qi in range(Q):
        for j in range(n):
            src=nbr[qi,j]
            if src>=0: f[qi,j]=fb[qi,src]
            else: f[qi,j]=fb[opp[qi],j]

def run_lbm(pore, axis, g=1e-6, tau=1.0, max_steps=100000, check_every=100, tol=1e-7, flux_tol=1e-3):
    """TRT + bounce-back. Stops only when temporal residual and slice-flux CV both meet tolerance."""
    c,w,opp=lattice()
    nu=(tau-.5)/3
    active=largest_percolating_cluster(pore,axis)
    if not active.any():
        return dict(k_lattice=0.0,converged=True,steps=0,history=[],through=False)
    sh=active.shape
    idx=-np.ones(sh,np.int64); idx[active]=np.arange(active.sum())
    n=int(active.sum()); Q=19
    nbr=np.empty((Q,n),np.int64)
    for qi in range(Q):
        shift=tuple(int(v) for v in c[qi])
        rolled=np.roll(idx,shift=shift,axis=(0,1,2))
        nbr[qi]=rolled[active]
    force=np.zeros(3); force[axis]=g
    wp=1.0/tau; wm=1.0/(0.5+3/16/(tau-.5))
    f=np.repeat(w[:,None],n,axis=1).copy(); fb=np.empty_like(f)
    history=[]; prev=None; converged=False; mean_u=0.
    coords=np.array(np.nonzero(active)).T
    transverse_area=pore.size//sh[axis]
    for it in range(1,max_steps+1):
        lbm_step_global(f,fb,nbr,c,w,opp,force,wp,wm)
        if it%check_every==0:
            rho=f.sum(axis=0)
            mom=np.zeros((3,n))
            for q in range(Q): mom += c[q,:,None]*f[q]
            ux=(mom[axis]+.5*force[axis])/rho
            mean_u=float(ux.sum()/pore.size)
            residual=float("nan") if prev is None else abs(mean_u-prev)/max(abs(mean_u),1e-30)
            slice_sums=np.bincount(coords[:,axis],weights=ux,minlength=sh[axis])
            slice_means=slice_sums/transverse_area
            flux_cv=float(np.std(slice_means)/max(abs(float(np.mean(slice_means))),1e-30))
            history.append([it,mean_u,residual,flux_cv])
            if prev is not None and residual<tol and flux_cv<flux_tol:
                converged=True; break
            prev=mean_u
    # Reconstruct final velocity and physical-regime diagnostics (as in course P11/P13).
    rho = f.sum(axis=0)
    mom = np.zeros((3, n), dtype=np.float64)
    for qi in range(Q):
        mom += c[qi, :, None] * f[qi]
    vel = (mom + 0.5 * force[:, None]) / rho[None, :]
    speed = np.sqrt(np.sum(vel * vel, axis=0))
    umax = float(speed.max()) if speed.size else 0.0
    mach = umax * math.sqrt(3.0)
    # Recompute superficial velocity from the final post-streaming field. The
    # mean is divided by total bulk voxels, matching the course P11 convention.
    mean_u = float(np.sum(vel[axis]) / pore.size)
    # Slice-flux variation is a separate P11 validity diagnostic, not the
    # temporal convergence residual. Each slice is averaged over its transverse plane.
    coords = np.array(np.nonzero(active)).T
    slice_sums = np.bincount(coords[:, axis], weights=vel[axis], minlength=sh[axis])
    transverse_area = pore.size // sh[axis]
    slice_means = slice_sums / transverse_area
    slice_flux_cv = float(np.std(slice_means) / max(abs(float(np.mean(slice_means))), 1e-30))
    # Hydraulic diameter from the voxel-face area of the active pore space.
    muka = sum(int(np.count_nonzero(active != np.roll(active, 1, axis=a))) for a in range(3))
    dh = 4.0 * float(n) / max(muka, 1)
    # P11 defines pore Reynolds number from the mean axial velocity on active pores.
    u_pore_axial = float(np.mean(vel[axis])) if n else 0.0
    reynolds = u_pore_axial * dh / nu
    tau_u = float(np.sum(speed) / max(np.sum(vel[axis]), 1e-30))
    mass_error = abs(float(rho.mean()) - 1.0)
    # Lattice intrinsic permeability k = nu <u_superficial> / g.
    k=nu*mean_u/g
    physical_checks = {
        "mach_lt_0_01": bool(mach < 0.01),
        "reynolds_lt_1": bool(reynolds < 1.0),
        "mass_error_lt_1e_6": bool(mass_error < 1e-6),
        "slice_flux_cv_lt_1e_3": bool(slice_flux_cv < 1e-3),
    }
    return dict(k_lattice=float(k),converged=converged,steps=it,history=history,through=True,
                active_pore_voxels=n,mean_velocity=mean_u,mach=mach,
                reynolds_pore=reynolds,tau_u=tau_u,mass_error=mass_error,
                slice_flux_cv=slice_flux_cv,slice_flux_profile=slice_means.tolist(),physical_checks=physical_checks)


def run_lbm_cupy(pore, axis, g=1e-6, tau=1.0, max_steps=400000, check_every=100, tol=1e-7, flux_tol=1e-3):
    """GPU implementation of the same D3Q19 TRT/Guo/half-way-bounce-back scheme.

    Requires CuPy and an available CUDA device. It intentionally keeps the full
    subvolume grid (rather than compacting active nodes) to use efficient GPU rolls.
    """
    try:
        import cupy as cp
    except Exception as exc:
        raise RuntimeError("CuPy belum terpasang. Pilih runtime GPU dan pasang cupy-cuda12x.") from exc
    if cp.cuda.runtime.getDeviceCount() < 1:
        raise RuntimeError("Tidak ada perangkat CUDA yang aktif.")

    c_np, w_np, opp_np = lattice()
    c = cp.asarray(c_np, dtype=cp.int32)
    w = cp.asarray(w_np, dtype=cp.float64)
    active_np = through_pore_cluster(pore, axis)
    if not active_np.any():
        return dict(k_lattice=0.0, converged=True, steps=0, history=[], through=False,
                    active_pore_voxels=0, mean_velocity=0.0, mach=0.0)
    active = cp.asarray(active_np)
    shape = pore.shape
    Q = len(w_np)
    nu = (tau - .5) / 3.0
    wp = 1.0 / tau
    wm = 1.0 / (0.5 + 3.0 / 16.0 / (tau - .5))
    F = np.zeros(3, dtype=np.float64); F[axis] = g
    transverse = tuple(a for a in range(3) if a != axis)
    # Precompute which source node is active for each pull-streaming direction.
    # Dua salinan fungsi distribusi sudah cukup: setelah collision selesai ditulis
    # ke fb, hasil streaming dapat ditulis kembali ke f karena seluruh sumber
    # streaming dibaca dari fb (bukan dari f). Ini menjaga asumsi memori 304 byte/node.
    f = cp.empty((Q,) + shape, dtype=cp.float64)
    fb = cp.empty_like(f)
    for q in range(Q):
        f[q] = w_np[q] * active

    history = []
    previous = None
    converged = False
    mean_u = 0.0
    for step in range(1, max_steps + 1):
        rho = cp.sum(f, axis=0)
        rho_safe = cp.where(active, rho, 1.0)
        m0 = cp.zeros(shape, dtype=cp.float64)
        m1 = cp.zeros(shape, dtype=cp.float64)
        m2 = cp.zeros(shape, dtype=cp.float64)
        for q in range(Q):
            if c_np[q, 0]: m0 += int(c_np[q, 0]) * f[q]
            if c_np[q, 1]: m1 += int(c_np[q, 1]) * f[q]
            if c_np[q, 2]: m2 += int(c_np[q, 2]) * f[q]
        ux = (m0 + .5 * F[0]) / rho_safe
        uy = (m1 + .5 * F[1]) / rho_safe
        uz = (m2 + .5 * F[2]) / rho_safe
        u2 = ux * ux + uy * uy + uz * uz
        uF = ux * F[0] + uy * F[1] + uz * F[2]
        for q in range(Q):
            o = int(opp_np[q])
            if o < q:
                continue
            cx, cy, cz = map(int, c_np[q])
            cu = cx * ux + cy * uy + cz * uz
            cF = cx * F[0] + cy * F[1] + cz * F[2]
            ei = w_np[q] * rho_safe * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * u2)
            eo = w_np[o] * rho_safe * (1 - 3 * cu + 4.5 * cu * cu - 1.5 * u2)
            fp = .5 * (f[q] + f[o])
            fm = .5 * (f[q] - f[o])
            ep = .5 * (ei + eo)
            em = .5 * (ei - eo)
            sg = w_np[q] * (-3 * uF + 9 * cu * cF)
            sa = w_np[q] * 3 * cF
            gp = -wp * (fp - ep) + (1 - .5 * wp) * sg
            gm = -wm * (fm - em) + (1 - .5 * wm) * sa
            fb[q] = cp.where(active, f[q] + gp + gm, 0.0)
            if o != q:
                fb[o] = cp.where(active, f[o] + gp - gm, 0.0)

        for q in range(Q):
            o = int(opp_np[q])
            shift = tuple(int(v) for v in c_np[q])
            source_active_q = cp.roll(active, shift=shift, axis=(0, 1, 2))
            shifted = cp.roll(fb[q], shift=shift, axis=(0, 1, 2))
            f[q] = cp.where(active, cp.where(source_active_q, shifted, fb[o]), 0.0)
            del source_active_q, shifted

        if step % check_every == 0:
            mean_u = float((cp.sum(cp.where(active, ux if axis == 0 else uy if axis == 1 else uz, 0.0))
                            / pore.size).get())
            residual = float("nan") if previous is None else abs(mean_u - previous) / max(abs(mean_u), 1e-30)
            check_slice_means = cp.mean(cp.where(active, ux if axis == 0 else uy if axis == 1 else uz, 0.0), axis=transverse)
            check_slice_cv = float((cp.std(check_slice_means) / max(abs(float(cp.mean(check_slice_means).get())), 1e-30)).get())
            history.append([step, mean_u, residual, check_slice_cv])
            if previous is not None and np.isfinite(residual) and residual < tol and check_slice_cv < flux_tol:
                converged = True
                break
            previous = mean_u

    # Reconstruct the final velocity from the post-streaming populations.
    rho = cp.sum(f, axis=0)
    rho_safe = cp.where(active, rho, 1.0)
    mom = [cp.zeros(shape, dtype=cp.float64) for _ in range(3)]
    for q in range(Q):
        for a in range(3):
            if c_np[q, a]:
                mom[a] += int(c_np[q, a]) * f[q]
    vel = [(mom[a] + .5 * F[a]) / rho_safe for a in range(3)]
    u_axis = vel[axis]
    mean_u = float((cp.sum(cp.where(active, u_axis, 0.0)) / pore.size).get())
    speed2 = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]
    umax = float(cp.max(cp.where(active, cp.sqrt(speed2), 0.0)).get())
    mach = umax * math.sqrt(3.0)
    transverse = tuple(a for a in range(3) if a != axis)
    slice_flux = cp.mean(cp.where(active, u_axis, 0.0), axis=transverse)
    slice_mean = float(cp.mean(slice_flux).get())
    slice_cv = float((cp.std(slice_flux) / max(abs(slice_mean), 1e-30)).get())
    # Hydraulic diameter and pore Reynolds number are diagnostic, not a replacement
    # for the convergence history or the anisotropy comparison.
    muka = sum(int(np.count_nonzero(active_np != np.roll(active_np, 1, axis=a))) for a in range(3))
    dh_voxel = 4.0 * float(active_np.sum()) / max(muka, 1)
    u_pore_axial = float((cp.sum(cp.where(active, u_axis, 0.0)) / max(int(active_np.sum()), 1)).get())
    reynolds = u_pore_axial * dh_voxel / nu
    sum_speed = float(cp.sum(cp.where(active, cp.sqrt(speed2), 0.0)).get())
    sum_axial = float(cp.sum(cp.where(active, u_axis, 0.0)).get())
    tau_u = sum_speed / max(sum_axial, 1e-30)
    rho_active_mean = float((cp.sum(cp.where(active, rho, 0.0)) / max(int(active_np.sum()), 1)).get())
    mass_error = abs(rho_active_mean - 1.0)
    physical_checks = {
        "mach_lt_0_01": bool(mach < 0.01),
        "reynolds_lt_1": bool(reynolds < 1.0),
        "mass_error_lt_1e_6": bool(mass_error < 1e-6),
        "slice_flux_cv_lt_1e_3": bool(slice_cv < 1e-3),
    }
    return dict(k_lattice=float(nu * mean_u / g), converged=converged, steps=step,
                history=history, through=True, active_pore_voxels=int(active_np.sum()),
                mean_velocity=mean_u, mach=mach, reynolds_pore=reynolds,tau_u=tau_u,
                slice_flux_cv=slice_cv, mass_error=mass_error, physical_checks=physical_checks)


def choose_flow_backend(requested):
    if requested == "cpu":
        return "cpu"
    if requested == "cupy":
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount() < 1:
            raise RuntimeError("Backend CuPy dipilih, tetapi CUDA tidak tersedia.")
        return "cupy"
    try:
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount() > 0:
            return "cupy"
    except Exception:
        pass
    return "cpu"

def permeability(max_steps=400000, backend="auto", selected_sizes=None):
    pore,meta,phi_full=load_image()
    rng=np.random.default_rng(SEED+101)
    backend_used=choose_flow_backend(backend)
    print(f"Backend aliran: {backend_used}",flush=True)
    posrows=[]; rows=[]
    selected_sizes=set(SIZES if selected_sizes is None else selected_sizes)
    checkpoint=OUT/"permeabilitas_subvolume_checkpoint.csv"
    if checkpoint.exists():
        with open(checkpoint,encoding="utf-8") as f:
            rows=list(csv.DictReader(f))
        row_index={(int(r["L_voxel"]),int(r["replicate"]),r["axis"]):i for i,r in enumerate(rows)}
        completed=set()
        for r in rows:
            key=(int(r["L_voxel"]),int(r["replicate"]),r["axis"])
            same_backend=(r.get("backend")==backend_used)
            is_disconnected=(str(r.get("through")).lower() in ("false","0"))
            is_converged=(str(r.get("converged")).lower() in ("true","1"))
            checks_present=bool(r.get("physical_checks_json"))
            if same_backend and (is_disconnected or (is_converged and checks_present)):
                completed.add(key)
        print(f"Ripanjut: {len(completed)} kasus lengkap pada backend {backend_used}.",flush=True)
    else:
        row_index={}; completed=set()
    for L in SIZES:
        starts=rand_subvolumes(pore.shape,L,N_FLOW,rng)
        for j,st in enumerate(starts):
            posrows.append(dict(L_voxel=L,replicate=j,seed=SEED+101,x0=int(st[0]),y0=int(st[1]),z0=int(st[2])))
        if L not in selected_sizes:
            continue
        for j,st in enumerate(starts):
            sub=extract(pore,st,L)
            for axis,name in enumerate("xyz"):
                if (L,j,name) in completed: continue
                t0=time.time()
                solver = run_lbm_cupy if backend_used == "cupy" else run_lbm
                r=solver(sub,axis,max_steps=max_steps,check_every=100,tol=1e-7)
                # 1 mD = 0.986923e-15 m^2 (not 0.986923e-12 m^2, which is 1 D)
                k_md=r["k_lattice"]*(VOXEL_M**2)/(0.986923e-15)
                new_row=dict(L_voxel=L,L_um=L*2.25,replicate=j,axis=name,
                    k_lattice=r["k_lattice"],k_mD=k_md,through=r["through"],
                    converged=r["converged"],steps=r["steps"],
                    final_residual=(r["history"][-1][2] if len(r["history"])>1 else None),
                    final_slice_flux_cv=(r["history"][-1][3] if r["history"] else None),
                    mean_velocity=r.get("mean_velocity",None),mach=r.get("mach",None),
                    reynolds_pore=r.get("reynolds_pore",None),tau_u=r.get("tau_u",None),
                    mass_error=r.get("mass_error",None),slice_flux_cv=r.get("slice_flux_cv",None),
                    physical_checks_json=json.dumps(r.get("physical_checks",{})),backend=backend_used,
                    elapsed_s=time.time()-t0,history_json=json.dumps(r["history"]))
                key=(L,j,name)
                if key in row_index:
                    rows[row_index[key]]=new_row
                else:
                    row_index[key]=len(rows); rows.append(new_row)
                print(f"L={L:3} rep={j} axis={name} k={k_md:.5g} mD "
                      f"through={r['through']} conv={r['converged']} steps={r['steps']}",flush=True)
                completed.add((L,j,name))
                save_csv(checkpoint,rows)
    save_csv(OUT/"permeabilitas_subvolume.csv",rows)
    save_csv(OUT/"koordinat_subvolume_aliran.csv",posrows)
    # Disconnected, nonconverged, and physically invalid cases remain in the raw
    # table. Means/Cv use only through-flow runs that converged and pass the
    # P11-style Mach, Reynolds, mass, and slice-flux checks.
    summary=[]
    for L in SIZES:
        per_axis_valid = {}
        for axis in "xyz":
            rr=[r for r in rows if int(r["L_voxel"])==L and r["axis"]==axis]
            through=[r for r in rr if str(r["through"]).lower() in ("true","1")]
            converged=[r for r in through if str(r["converged"]).lower() in ("true","1") and float(r["k_mD"])>0]
            good=[r for r in converged if physical_checks_pass(r)]
            mu,sd,cv,n=cv_stats([float(r["k_mD"]) for r in good])
            per_axis_valid[axis]={int(r["replicate"]):float(r["k_mD"]) for r in good}
            tau_vals=[float(r["tau_u"]) for r in good if r.get("tau_u") not in (None, "")]
            summary.append(dict(property="permeability_directional",L_voxel=L,L_um=L*2.25,axis=axis,
                n_total=len(rr),n_through=len(through),n_not_through=len(rr)-len(through),
                n_nonconverged=len(through)-len(converged),
                n_physical_checks_passed=len(good),
                n_physical_checks_incomplete_or_failed=len(converged)-len(good),n_usable=n,
                mean_k_mD=mu,std_k_mD=sd,Cv=cv,mean_tau_u=(float(np.mean(tau_vals)) if tau_vals else float("nan")),n_for_Cv=n))
        # The scalar permeability for REV is the geometric mean k_geo within each
        # replicate, requiring all three directions to be through and converged.
        kgeo=[]
        n_missing_nonthrough=0
        n_missing_nonconverged=0
        n_missing_physical=0
        for rep in range(N_FLOW):
            vals=[per_axis_valid[a].get(rep) for a in "xyz"]
            if all(v is not None and v>0 for v in vals):
                kgeo.append(float(np.exp(np.mean(np.log(vals)))))
            else:
                replicate_rows=[r for r in rows if int(r["L_voxel"])==L and int(r["replicate"])==rep]
                if any(str(r["through"]).lower() not in ("true", "1") for r in replicate_rows):
                    n_missing_nonthrough += 1
                elif any(str(r["converged"]).lower() not in ("true", "1") for r in replicate_rows):
                    n_missing_nonconverged += 1
                elif any(not physical_checks_pass(r) for r in replicate_rows):
                    n_missing_physical += 1
        mu,sd,cv,n=cv_stats(kgeo)
        summary.append(dict(property="permeability_geometric_mean",L_voxel=L,L_um=L*2.25,axis="geo",
            n_total=N_FLOW,n_valid_all_axes=n,n_incomplete_any_axis=N_FLOW-n,
            n_missing_due_to_nonthrough=n_missing_nonthrough,
            n_missing_due_to_nonconvergence=n_missing_nonconverged,
            n_missing_due_to_physical_checks=n_missing_physical,
            mean_k_mD=mu,std_k_mD=sd,Cv=cv,n_for_Cv=n))
    save_csv(OUT/"koefisien_variasi_permeabilitas.csv",summary)
    # Extrapolate at Cv=5% for the isotropic scalar k_geo and retain directional checks.
    estimates={}
    for axis in ["geo", "x", "y", "z"]:
        vals=[r for r in summary if r["axis"]==axis and np.isfinite(float(r["Cv"])) and float(r["Cv"])>0]
        hits=[r for r in vals if float(r["Cv"])<0.05]
        if hits:
            estimates[axis]={"status":"reached","L_voxel":hits[0]["L_voxel"],"L_um":hits[0]["L_um"]}
        elif vals:
            last=max(vals,key=lambda r:r["L_voxel"])
            Lest=last["L_voxel"]*(float(last["Cv"])/0.05)**(2/3)
            mem_bytes=304*(Lest**3)
            estimates[axis]={"status":"not_reached_extrapolated","Cv_at_largest":float(last["Cv"]),
                "L_est_voxel_using_L_minus_3_over_2":float(Lest),"L_est_um":float(Lest*2.25),
                "memory_GB_decimal":float(mem_bytes/1e9),"memory_formula":"304 * L^3 bytes"}
        else:
            estimates[axis]={"status":"insufficient_valid_subvolumes"}
    (OUT/"perkiraan_REV_permeabilitas.json").write_text(json.dumps(estimates,indent=2))
    return rows,summary

def anisotropy():
    """Compare directions using the same paired subvolumes at each L.

    A replicate enters the directional comparison only when all three axes are
    through-flow and converged. Counts of through/non-through/nonconverged
    simulations remain visible in the output so exclusions are never hidden.
    """
    p=OUT/"permeabilitas_subvolume.csv"
    if not p.exists():
        raise SystemExit("Jalankan tahap permeabilitas lebih dahulu.")
    rows=list(csv.DictReader(open(p,encoding="utf-8")))
    out=[]
    for L in SIZES:
        by_rep={}
        for r in rows:
            if int(r["L_voxel"]) != L:
                continue
            by_rep.setdefault(int(r["replicate"]),{})[r["axis"]]=r
        paired=[]
        for rep, directions in by_rep.items():
            ok=True
            for a in "xyz":
                r=directions.get(a)
                if r is None or str(r["through"]).lower() not in ("true","1") or \
                   str(r["converged"]).lower() not in ("true","1") or float(r["k_mD"]) <= 0 or \
                   not physical_checks_pass(r):
                    ok=False; break
            if ok:
                paired.append((rep,{a:float(directions[a]["k_mD"]) for a in "xyz"}))
        directional={a:[v[a] for _,v in paired] for a in "xyz"}
        means={a:(float(np.mean(v)) if v else float("nan")) for a,v in directional.items()}
        sds={a:(float(np.std(v,ddof=0)) if v else float("nan")) for a,v in directional.items()}
        counts={a:sum(1 for r in rows if int(r["L_voxel"])==L and r["axis"]==a and
                      str(r["through"]).lower() in ("true","1")) for a in "xyz"}
        nonconv={a:sum(1 for r in rows if int(r["L_voxel"])==L and r["axis"]==a and
                       str(r["through"]).lower() in ("true","1") and
                       str(r["converged"]).lower() not in ("true","1")) for a in "xyz"}
        physfail={a:sum(1 for r in rows if int(r["L_voxel"])==L and r["axis"]==a and
                        str(r["through"]).lower() in ("true","1") and
                        str(r["converged"]).lower() in ("true","1") and
                        not physical_checks_pass(r)) for a in "xyz"}
        if len(paired)>=2:
            amax=max("xyz",key=lambda a:means[a]); amin=min("xyz",key=lambda a:means[a])
            delta=means[amax]-means[amin]
            ratio=delta/sds[amax] if sds[amax]>0 else (float("inf") if delta>0 else 0.0)
            anis=means[amax]/means[amin] if means[amin]>0 else float("nan")
            test=bool(ratio>2)
        else:
            amax=amin=None; delta=ratio=anis=float("nan"); test=None
        out.append(dict(L_voxel=L,L_um=L*2.25,
            mean_kx_mD=means["x"],mean_ky_mD=means["y"],mean_kz_mD=means["z"],
            std_kx_mD=sds["x"],std_ky_mD=sds["y"],std_kz_mD=sds["z"],
            n_through_x=counts["x"],n_through_y=counts["y"],n_through_z=counts["z"],
            n_nonconverged_x=nonconv["x"],n_nonconverged_y=nonconv["y"],n_nonconverged_z=nonconv["z"],
            n_physical_invalid_x=physfail["x"],n_physical_invalid_y=physfail["y"],n_physical_invalid_z=physfail["z"],
            n_paired_valid=len(paired),paired_replicates=";".join(str(r) for r,_ in paired),
            kmax_axis=amax,kmin_axis=amin,anisotropy_ratio=anis,
            delta_kmax_min_mD=delta,delta_over_sd_kmax=ratio,
            distinguishable_over_2sd=test,
            comparison_note="Rerata arah dihitung pada replikasi berpasangan yang menembus, konvergen, dan lolos pemeriksaan fisis pada X,Y,Z."))
    save_csv(OUT/"penimbangan_anisotropi.csv",out)
    return out

def fit_slopes():
    geom=list(csv.DictReader(open(OUT/"koefisien_variasi_geometri.csv",encoding="utf-8")))
    per=list(csv.DictReader(open(OUT/"koefisien_variasi_permeabilitas.csv",encoding="utf-8"))) if (OUT/"koefisien_variasi_permeabilitas.csv").exists() else []
    slopes={}
    for key,rows in [("porosity",[r for r in geom if r["property"]=="porosity"]),
                     ("surface",[r for r in geom if r["property"]=="specific_surface_corrected_m_inv"]),
                     ("permeability_geometric_mean",[r for r in per if r.get("property")=="permeability_geometric_mean"])]:
        vals=[]
        for r in rows:
            cv=float(r["Cv"])
            if np.isfinite(cv) and cv>0: vals.append((float(r["L_voxel"]),cv))
        slopes[key]=float(np.polyfit(np.log([x[0] for x in vals]),np.log([x[1] for x in vals]),1)[0]) if len(vals)>=2 else None
    return slopes

def combined_cv_plot():
    """Final three-property Cv(L) plot required by the Bab 6 project."""
    geom = list(csv.DictReader(open(OUT/"koefisien_variasi_geometri.csv", encoding="utf-8")))
    per_path = OUT/"koefisien_variasi_permeabilitas.csv"
    per = list(csv.DictReader(open(per_path, encoding="utf-8"))) if per_path.exists() else []
    fig, ax = plt.subplots(figsize=(9, 5.8))
    specs = [
        ("Porositas", [r for r in geom if r["property"]=="porosity"]),
        ("Luas permukaan spesifik", [r for r in geom if r["property"]=="specific_surface_corrected_m_inv"]),
        ("Permeabilitas rerata geometrik", [r for r in per if r.get("property")=="permeability_geometric_mean"]),
    ]
    for label, rows in specs:
        vals=[]
        for r in rows:
            try:
                L=float(r["L_voxel"]); cv=float(r["Cv"])
                if np.isfinite(cv) and cv>0: vals.append((L,100*cv))
            except (ValueError,KeyError):
                continue
        if len(vals)<1: continue
        vals=sorted(vals)
        x=np.array([v[0] for v in vals]); y=np.array([v[1] for v in vals])
        line,=ax.plot(x,y,marker="o",linewidth=2,label=label)
        yref=y[0]*(x/x[0])**(-1.5)
        ax.plot(x,yref,linestyle=":",color=line.get_color(),alpha=.65,
                label=f"Acuan −3/2 ({label})")
    ax.axhline(5,linestyle="--",color="black",linewidth=1.1,label="Ambang REV 5%")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Ukuran sisi subvolume, L (voxel)")
    ax.set_ylabel("Koefisien variasi, Cv (%)")
    ax.set_title("Koefisien variasi tiga sifat terhadap ukuran subvolume — Kirby")
    ax.grid(True,which="both",alpha=.25); ax.legend(fontsize=8,ncol=2)
    fig.tight_layout(); fig.savefig(FIG/"cv_tiga_sifat.png",dpi=200); plt.close(fig)
    return fit_slopes()


def main():
    OUT.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
    p=argparse.ArgumentParser()
    p.add_argument("--stage",choices=["domain","geometry","flow","anisotropy","all"],default="geometry")
    p.add_argument("--max-steps",type=int,default=400000,help="batas langkah tiap kasus LBM; run yang tidak konvergen tetap dicatat")
    p.add_argument("--backend",choices=["auto","cpu","cupy"],default="auto",help="auto memilih CuPy bila CUDA tersedia")
    p.add_argument("--flow-sizes",nargs="+",type=int,choices=SIZES,default=None,
                   help="ukuran sisi untuk tahap flow; jalankan bertahap agar checkpoint mudah diperiksa")
    args=p.parse_args()
    if args.stage in ("domain","all"): domain_reference()
    if args.stage in ("geometry","all"): geometry()
    if args.stage in ("flow","all"):
        permeability(max_steps=args.max_steps,backend=args.backend,selected_sizes=args.flow_sizes)
    if args.stage in ("anisotropy","all"):
        anisotropy()
    if args.stage in ("anisotropy","all"):
        slopes=combined_cv_plot()
        print("Kemiringan log(Cv)-log(L):",json.dumps(slopes,indent=2))
    print("Finished:",args.stage)

if __name__=="__main__": main()
