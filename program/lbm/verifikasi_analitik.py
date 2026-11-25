"""Verification helper: flow through straight circular channel (D3Q19 TRT)."""
import csv, sys, time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "program" / "inti"))
from rev_anisotropi import run_lbm

N=24
Y,Z=np.meshgrid(np.arange(N),np.arange(N),indexing='ij')
c=(N-1)/2
rows=[]
for R in (4,5,6):
    cross=(Y-c)**2+(Z-c)**2 <= R**2
    pore=np.broadcast_to(cross[None,:,:],(N,N,N)).copy()
    t=time.time()
    r=run_lbm(pore,axis=0,g=2e-6,tau=1.0,max_steps=12000,check_every=200,tol=1e-6)
    # Continuum pipe solution, averaged over the entire lattice cross-section.
    phi=float(cross.mean())
    exact=phi*R**2/8.0
    rows.append(dict(radius_voxel=R,porosity=phi,k_lbm=r['k_lattice'],k_analytical=exact,
        relative_error=abs(r['k_lattice']/exact-1),converged=r['converged'],steps=r['steps'],
        final_residual=r['history'][-1][2] if len(r['history'])>1 else None,elapsed_s=time.time()-t))
    print(rows[-1],flush=True)
with open(ROOT/'hasil'/'verifikasi_kanal_lurus.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
