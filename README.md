# Volume Elementer Representatif dan Anisotropi — Citra Batuan Kirby

**Penulis:** Bagas Hada  
**NPM:** 251020260008  
**Program:** Magister Inovasi Regional, Universitas Padjadjaran  
**Mata kuliah:** Simulasi Berbasis Data Citra  
**Sampel:** Kirby, citra biner micro-CT 256 × 256 × 256 voxel  
**Ukuran voxel:** 2,25 µm

[![Buka workflow di Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/BagasHada/REV-Anisotropi-Kirby/blob/main/notebooks/02_workflow_revisi_REV_anisotropi_Kirby_Colab.ipynb)

Repositori ini mendokumentasikan alur komputasi proyek REV dan anisotropi pada citra batuan Kirby. Isi repositori mencakup notebook, kode analisis, metadata dan data pembanding, keluaran numerik, grafik, serta catatan metodologis. **Dokumen laporan akhir disimpan terpisah dan tidak menjadi bagian dari repositori ini.**

## Pertanyaan penelitian

1. Pada ukuran subvolume berapa porositas dan luas permukaan spesifik menunjukkan variasi statistik yang rendah?
2. Apakah ukuran REV geometri juga berlaku untuk permeabilitas?
3. Apakah perbedaan permeabilitas arah X, Y, dan Z lebih besar daripada variasi antarposisi subvolume?
4. Bagaimana tumpang tindih subvolume, konektivitas pori, dan pemeriksaan fluks memengaruhi interpretasi?

## Ringkasan hasil komputasi

- Porositas total domain penuh: **20,5284%**; porositas alir P3: **20,3326%**.
- Granulometri citra biner: `r50 = 4,821 voxel` (10,847 µm); diameter pengendali perkiraan 21,693 µm; `d_t/Δx = 9,64 voxel`. Ini merupakan indikator granulometri, bukan pengukuran langsung setiap leher pori.
- Lima ukuran subvolume geometri: 32, 48, 64, 96, dan 128 voxel; 30 subvolume per ukuran.
- CV porositas: **64,18%; 27,01%; 20,23%; 10,36%; 4,84%**.
- CV luas permukaan spesifik terkoreksi: **36,55%; 22,10%; 13,27%; 5,53%; 3,82%**.
- Kemiringan log-CV terhadap log-L: **−1,768** (porositas) dan **−1,703** (SSA terkoreksi). Secara nominal kedua sifat melewati ambang CV < 5% pada L = 128 voxel (288 µm).
- Pada L = 128, rata-rata tumpang tindih pasangan subvolume **33,55%** (median 32,35%). Karena itu, indikasi REV geometri pada ukuran maksimum harus dibaca dengan hati-hati; subvolume tidak independen secara kuat.
- Batch LBM terdahulu memiliki **90 catatan**: 34 lolos pemeriksaan versi notebook, 46 memiliki variasi fluks antaririsan tinggi, dan 10 non-percolating. Pada L = 128, 18 dari 18 kasus lolos pemeriksaan versi notebook terdahulu.
- Rerata permeabilitas batch historis pada L = 128: **482,25 mD (X)**, **475,20 mD (Y)**, dan **560,11 mD (Z)**. CV masing-masing **13,06%; 10,79%; 28,85%**, sehingga REV permeabilitas belum tercapai pada ukuran yang diuji.
- Rasio permeabilitas maksimum/minimum: **1,179**; indeks anisotropi: **0,526**. Perbedaan arah belum terpisah jelas dari variasi antarsubvolume menurut patokan analisis yang digunakan.

**Catatan validasi:** batch LBM historis belum menyimpan seluruh diagnostik tambahan program revisi (khususnya Reynolds pori dan kesalahan massa). Karena itu, batch historis dan pilot program revisi dipisahkan; statistik keduanya tidak digabung. Keluaran pilot dengan variasi fluks yang tidak memenuhi ambang disimpan sebagai diagnostik, bukan hasil yang lolos.

## Data dan keterulangan

Metadata citra tersedia di `data/metadata.json`, sedangkan nilai pembanding permeabilitas domain penuh tersedia di `data/referensi_permeabilitas.csv`. Citra RAW sumber tidak disertakan dalam repositori publik karena berasal dari bundel mata kuliah dan izin redistribusinya belum dipastikan. Untuk menjalankan analisis, peroleh citra dari bundel resmi mata kuliah dan letakkan sebagai `data/citra/Kirby_2d25um_binary_256.raw`. Spesifikasi format dan pemeriksaan ukuran dijelaskan di [`data/citra/README.md`](data/citra/README.md).

## Cara menjalankan

### Google Colab

1. Buka notebook [`02_workflow_revisi_REV_anisotropi_Kirby_Colab.ipynb`](notebooks/02_workflow_revisi_REV_anisotropi_Kirby_Colab.ipynb) di GitHub, lalu pilih **Open in Colab**.
2. Pada sel pertama, unggah ZIP repository yang diunduh melalui **Code → Download ZIP** jika folder proyek belum ada di Google Drive.
3. Pada permintaan berikutnya, unggah `Kirby_2d25um_binary_256.raw` dari bundel resmi mata kuliah. Notebook menyalinnya ke `data/citra/` pada folder kerja Google Drive.
4. Jalankan sel berurutan. Untuk LBM, jalankan ukuran satu per satu dan pertahankan checkpoint ketika runtime terputus.

### Python lokal

Python 3.10 atau lebih baru disarankan. Jalankan perintah dari root repository:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt

python program/inti/rev_anisotropi.py --stage domain
python program/inti/granulometri_kelayakan.py
python program/inti/rev_anisotropi.py --stage geometry
python program/lbm/verifikasi_analitik.py
python program/inti/tortuositas_geometrik.py
python program/inti/visualisasi.py
```

Tahap LBM subvolume dijalankan terpisah melalui `program/inti/rev_anisotropi.py --stage flow`. Periksa opsi dengan `--help`. Backend GPU memerlukan lingkungan CUDA dan CuPy yang sesuai; untuk uji kecil dapat menggunakan CPU. Jangan memulai semua ukuran sekaligus pada satu runtime Colab.

## Peta program dan keluaran

| Luaran / angka | Program utama | Berkas keluaran |
|---|---|---|
| Takrif porositas dan domain | `program/inti/rev_anisotropi.py --stage domain` | `hasil/porositas_domain_penuh.csv`, `hasil/domain_acuan.json` |
| Granulometri dan rasio resolusi | `program/inti/granulometri_kelayakan.py` | `hasil/granulometri_ringkasan.json`, `hasil/granulometri_leher_pori.csv` |
| Koordinat dan statistik geometri REV | `program/inti/rev_anisotropi.py --stage geometry` | `hasil/subvolume_coordinates.csv`, `hasil/geometri_subvolume.csv`, `hasil/koefisien_variasi_geometri.csv`, `hasil/ringkasan_statistik_geometri.csv` |
| Diagnostik tumpang tindih | `program/inti/audit_reproduksibilitas.py` | `hasil/diagnostik_tumpang_tindih.csv` |
| Verifikasi kanal lurus | `program/lbm/verifikasi_analitik.py` | `hasil/verifikasi_kanal_lurus.csv` |
| Tortuositas geometri | `program/inti/tortuositas_geometrik.py` | `hasil/tortuositas_geometrik.csv`, `hasil/tortuositas_metadata.json` |
| Permeabilitas subvolume dan anisotropi | `program/inti/rev_anisotropi.py --stage flow` dan `--stage anisotropy` | checkpoint serta ringkasan di `hasil/pilot_program_revisi/` dan `hasil/` |
| Grafik | `program/inti/visualisasi.py` | `gambar/` |

## Status keterjalanan dan reproduksi

Panduan lengkap prasyarat data, pemasangan dependensi, perintah eksekusi, serta status tahapan yang sudah diuji/yang belum tervalidasi ada di [`docs/PANDUAN_REPRODUKSI.md`](docs/PANDUAN_REPRODUKSI.md). Citra RAW resmi tidak termasuk di repository; berkas itu harus diperoleh dari bundel mata kuliah terlebih dahulu. Tahap domain, geometri, granulometri, verifikasi kanal lurus, tortuositas, dan visualisasi telah diuji ulang. Simulasi LBM 90 kasus versi revisi belum dinyatakan selesai; workbook yang disertakan memuat batch historis dan pilot parsial dipisahkan.

## Isi repositori

```text
.
├── README.md
├── CATATAN.md
├── requirements.txt
├── data/                 # metadata dan referensi; citra RAW diperoleh dari bundel kuliah
├── notebooks/            # notebook historis dan workflow revisi
├── program/
│   ├── inti/              # pipeline domain, geometri, granulometri, tortuositas, audit, grafik
│   └── lbm/               # verifikasi solver LBM
├── hasil/                # hasil numerik, batch historis, pilot, dan arsip suplemen
├── gambar/               # grafik utama
└── docs/                 # audit, metode, dan catatan batasan
```

## Interpretasi dan batasan

- Citra yang tersedia sudah berbentuk biner. Artefak berbasis intensitas grayscale (misalnya beam hardening, ring artefact, dan partial-volume) tidak dapat dievaluasi dari citra ini.
- Porositas alir P3 tidak identik dengan masker pori aktif P11 pada simulasi LBM; definisi konektivitasnya berbeda.
- Nilai tortuositas adalah tortuositas geometri Dijkstra pada satu subvolume 128³, bukan tortuositas aliran berbobot fluks.
- Hasil LBM batch terdahulu dan pilot program revisi memiliki set pemeriksaan yang berbeda. Perbedaan status harus tetap ditampilkan.
- CSV suplemen yang belum direkonsiliasi ditempatkan pada `hasil/arsip_suplemen_perlu_rekonsiliasi/`, tidak digunakan untuk mengganti sumber utama geometri. Audit perbedaannya ada pada `docs/AUDIT_DATA_TAMBAHAN.md`.
- Keluaran pilot yang variasi fluks irisan-nya melebihi 10⁻³ bukan hasil permeabilitas yang lolos pemeriksaan.

## Referensi utama

- *Proyek Akhir dan Bedah Paper: Simulasi Berbasis Data Citra*, Program Magister Inovasi Regional, Universitas Padjadjaran, Bab 2 dan Bab 6.
- Nurcahya, A. & Dharmawan, I. A. (2025). Enhanced understanding of fluid flow and effective properties in three-dimensional porous rock using lattice Boltzmann simulation and image-based correlation studies. *Physics of Fluids*, 37(5), 053607. https://doi.org/10.1063/5.0268487
- Nurcahya, A., Alexandra, A., Akmal, F. & Dharmawan, I. A. (2024). The lattice Boltzmann method and image processing techniques for effective parameter estimation of digital rock. *Applied Sciences*, 14(17), 7509. https://doi.org/10.3390/app14177509
- Dharmawan, I. A., Kartika, S. R., Alexandra, A., Melkisedek, B. L. & Emeliana, Y. (2026). Explainable machine learning for sandstone rock typing: triple-anchor convergence with closed-form permeability laws. *Chemical Engineering Journal Advances*, 101323. https://doi.org/10.1016/j.ceja.2026.101323

## Penggunaan alat bantu

ChatGPT digunakan sebagai pendamping untuk memahami materi, menata notebook dan penjelasan sel Google Colab, membantu menyusun dan meninjau sebagian kode pendukung, serta merapikan dokumentasi. Eksekusi program, pemeriksaan keluaran numerik, keputusan analisis, dan interpretasi akhir diperiksa serta dipertanggungjawabkan oleh penulis.

## Lisensi dan distribusi

Repositori ini belum menetapkan lisensi distribusi untuk citra sumber atau data milik pihak lain. Jangan menganggap isi bundel mata kuliah bebas didistribusikan ulang.
