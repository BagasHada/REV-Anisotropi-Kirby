# Program dan kode analisis

Kode dipisahkan berdasarkan fungsi. `program/inti/rev_anisotropi.py` tetap menjadi mesin utama bertahap (domain, geometri, aliran LBM, dan anisotropi); skrip lain memanggil atau melengkapi tahapan tersebut. Pemisahan folder tidak berarti tiap skrip merupakan model yang independen.

## `inti/`
- `rev_anisotropi.py` — pemrosesan domain, konektivitas, sampling subvolume, statistik geometri, solver aliran LBM, checkpoint, dan ringkasan anisotropi.
- `granulometri_kelayakan.py` — indikator granulometri / resolusi pori.
- `tortuositas_geometrik.py` — tortuositas geometri dengan Dijkstra 6-tetangga.
- `visualisasi.py` — grafik dari hasil CSV.
- `audit_reproduksibilitas.py` — audit metadata, koordinat, dan diagnostik overlap.

## `lbm/`
- `verifikasi_analitik.py` — uji solver pada kanal lurus terhadap solusi Poiseuille.

Jalankan dari root repository agar path relatif `data/`, `hasil/`, dan `gambar/` ter-resolve dengan benar. Parameter tahap LBM dan opsi backend ditampilkan dengan `python program/inti/rev_anisotropi.py --help`.
