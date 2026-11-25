#!/usr/bin/env python3
"""Penilaian resolusi Kirby menggunakan granulometri pada notebook P2.

Agar hemat memori, r yang diuji dibatasi di sekitar lintasan fraksi pori 50%.
Ukuran r=4.0, 4.5, dan 5.0 voxel digunakan untuk interpolasi lokal r50.
"""
from pathlib import Path
import json, csv, time
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
DATA, OUT, FIG = ROOT / 'data', ROOT / 'hasil', ROOT / 'gambar'

def main():
    meta = json.loads((DATA / 'metadata.json').read_text())
    raw = np.fromfile(DATA / meta['raw_file'], dtype=np.uint8)
    dims = tuple(meta['size_voxels_xyz'])
    if raw.size != int(np.prod(dims)):
        raise ValueError('Jumlah voxel pada RAW tidak cocok dengan metadata')
    pore = raw.reshape(dims, order=meta['order']) == 0
    dx_um = float(meta['voxel_size_m']) * 1e6
    total_pore = int(pore.sum())
    t0 = time.time()
    distance = ndi.distance_transform_edt(pore)
    rows = []
    for r in (4.0, 4.5, 5.0):
        seeds = distance >= r
        outside = ndi.distance_transform_edt(~seeds)
        accessible = (outside <= r) & pore
        fraction = float(np.count_nonzero(accessible) / total_pore)
        rows.append({'radius_voxel': r, 'radius_um': r*dx_um,
                     'accessible_pore_fraction': fraction})
        print(f'r={r:.1f} voxel; fraksi terakses={fraction:.6f}', flush=True)
        del seeds, outside, accessible
    # Interpolasi lokal di antara titik yang mengapit fraksi 0,5.
    bracket = next(((rows[i], rows[i+1]) for i in range(len(rows)-1)
                    if rows[i]['accessible_pore_fraction'] >= .5 >= rows[i+1]['accessible_pore_fraction']), None)
    if bracket is None:
        raise RuntimeError('Fraksi 50% tidak terapit oleh radius uji; perluas radius dengan pemeriksaan memori.')
    a, b = bracket
    r50 = a['radius_voxel'] + (a['accessible_pore_fraction']-.5) * (b['radius_voxel']-a['radius_voxel']) / (a['accessible_pore_fraction']-b['accessible_pore_fraction'])
    d_vox = 2*r50
    d_um = d_vox*dx_um
    if d_vox < 2: decision='tidak_layak'
    elif d_vox < 4: decision='layak_dengan_catatan'
    elif d_vox < 8: decision='layak'
    else: decision='layak_dan_baik'
    OUT.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
    with (OUT/'granulometri_leher_pori.csv').open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    result={
        'method':'P2 granulometry: EDT(pore); seeds=(EDT>=r); accessible=(EDT(~seeds)<=r); local linear interpolation at 50% accessible pore volume',
        'sample':'Kirby','voxel_size_um':dx_um,'porosity_total':float(pore.mean()),
        'r50_voxel':float(r50),'r50_um':float(r50*dx_um),
        'estimated_controlling_throat_diameter_um':float(d_um),
        'estimated_controlling_throat_diameter_voxel':float(d_vox),
        'resolution_ratio_d_t_over_dx':float(d_vox),'decision':decision,
        'radius_points':rows,'runtime_seconds':time.time()-t0,
        'caveat':'Indikator granulometri dari citra biner, bukan pengukuran langsung tiap leher pori. Artefak grayscale tidak dapat dinilai dari citra biner.'}
    (OUT/'granulometri_ringkasan.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
    fig,ax=plt.subplots(figsize=(7,4.5))
    ax.plot([x['radius_um'] for x in rows],[x['accessible_pore_fraction'] for x in rows],'o-')
    ax.axhline(.5,color='black',linestyle='--',linewidth=1,label='Fraksi 50%')
    ax.axvline(r50*dx_um,color='gray',linestyle=':',linewidth=1,label=f'r50={r50:.2f} voxel')
    ax.set_xlabel('Jari-jari bola uji (µm)'); ax.set_ylabel('Fraksi volume pori terakses')
    ax.set_title('Granulometri pori Kirby — metode P2'); ax.grid(alpha=.25); ax.legend(); fig.tight_layout()
    fig.savefig(FIG/'granulometri_leher_pori.png',dpi=200); plt.close(fig)
    print(json.dumps(result,indent=2,ensure_ascii=False))
    return result

if __name__=='__main__': main()
