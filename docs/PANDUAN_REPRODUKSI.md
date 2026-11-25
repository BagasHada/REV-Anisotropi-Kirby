# Panduan reproduksi dan status keterjalanan

Dokumen ini menjelaskan prasyarat agar orang lain dapat menjalankan kembali analisis dari kode di repositori. Hasil yang sudah tersimpan dan tahapan yang belum selesai dibedakan secara eksplisit.

## 1. Prasyarat data

Repositori ini **tidak menyertakan** citra RAW `Kirby_2d25um_binary_256.raw`. Citra berasal dari bundel resmi mata kuliah dan tidak disertakan karena izin redistribusi publik belum dipastikan.

Sebelum menjalankan program:

1. Peroleh `Kirby_2d25um_binary_256.raw` dari bundel resmi mata kuliah.
2. Letakkan di `data/citra/Kirby_2d25um_binary_256.raw`.
3. Pastikan ukuran file tepat 16.777.216 byte (256³ voxel × 1 byte).
4. Pertahankan nilai 0 = pori, 1 = padatan, urutan C, dan metadata di `data/metadata.json`.

Tanpa berkas RAW, program yang membaca citra akan berhenti dengan `FileNotFoundError`. Notebook operasional Colab meminta pengguna mengunggah berkas RAW tersebut.

## 2. Lingkungan Python

Disarankan Python 3.10 atau lebih baru. Dari direktori utama repositori:

```bash
python -m venv .venv
# Windows PowerShell:
.venv\\Scripts\\Activate.ps1
# Linux/macOS:
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Notebook Colab sudah mencakup pemasangan dependensi. CuPy bersifat opsional dan dipasang oleh notebook jika runtime GPU NVIDIA tersedia. Untuk analisis geometri dan uji kecil, backend CPU dapat dipakai; simulasi LBM lengkap dapat memerlukan waktu komputasi yang panjang.

## 3. Menjalankan ulang tahapan yang telah diuji

Perintah dijalankan dari root repository, setelah citra RAW diletakkan pada lokasi di atas:

```bash
python program/inti/rev_anisotropi.py --stage domain
python program/inti/rev_anisotropi.py --stage geometry
python program/inti/granulometri_kelayakan.py
python program/lbm/verifikasi_analitik.py
python program/inti/tortuositas_geometrik.py
python program/inti/audit_reproduksibilitas.py
python program/inti/visualisasi.py
```

Tahapan domain, geometri, granulometri, verifikasi kanal lurus, tortuositas, dan visualisasi telah diuji ulang pada lingkungan penyusunan paket dengan citra RAW yang diunggah pengguna. Hasil domain, CV geometri, granulometri, galat kanal, dan tortuositas cocok dengan nilai yang dipakai dalam dokumentasi proyek.

## 4. Menjalankan simulasi LBM subvolume

LBM versi revisi dijalankan terpisah per ukuran agar checkpoint mudah diperiksa. Contoh:

```bash
python program/inti/rev_anisotropi.py --stage flow --backend auto --flow-sizes 32 --max-steps 400000
```

Setelah ukuran tersebut selesai atau checkpoint tersimpan, lanjutkan ukuran 48, 64, 96, dan 128 voxel dengan perintah serupa, satu ukuran pada satu waktu. Jangan menjalankan beberapa proses bersamaan karena semuanya menulis checkpoint yang sama. Periksa status konektivitas, konvergensi, bilangan Mach, Reynolds pori, kesalahan massa, serta variasi fluks irisan. Sesudah kasus yang dimaksud tercatat, jalankan:

```bash
python program/inti/rev_anisotropi.py --stage anisotropy
```

**Penting:** seluruh 90 kasus LBM versi revisi belum dinyatakan selesai atau tervalidasi dalam paket ini. Folder `hasil/pilot_program_revisi/` berisi pilot parsial; beberapa kasus tidak lolos kriteria variasi fluks. Jangan gunakan pilot itu sebagai pengganti batch historis.

## 5. Cara membaca notebook

- `notebooks/02_workflow_revisi_REV_anisotropi_Kirby_Colab.ipynb` adalah workflow operasional yang disarankan. Buka di Google Colab, unggah ZIP repositori bila diminta, lalu unggah citra RAW resmi ketika diminta. Jalankan sel berurutan dan pertahankan checkpoint.
- `notebooks/01_analisis_awal_dan_batch_LBM_historis.ipynb` adalah rekam kerja historis. Notebook ini mengasumsikan beberapa berkas dari bundel mata kuliah diunggah terpisah, termasuk citra RAW, metadata, dan CSV acuan. Gunakan notebook operasional untuk workflow reproduksi saat ini.

## 6. Status dan asal hasil

- Hasil domain/geometri/granulometri/tortuositas/verifikasi yang tersimpan merupakan keluaran dari pipeline proyek dan berkas rinciannya ada di `hasil/`.
- `hasil/Hasil_Final_REV_Anisotropi_Kirby.xlsx` adalah workbook yang pengguna unggah. Isinya mencakup 90 catatan LBM historis beserta sheet ringkasan tambahan.
- Pada batch LBM historis, 34 kasus berstatus lolos pemeriksaan versi notebook lama, 46 memiliki variasi fluks irisan tinggi, dan 10 non-percolating. Pemeriksaan tambahan versi program revisi (Re pori dan kesalahan massa) tidak tersedia lengkap pada batch itu.
- Sheet `Ringkasan_Porositas`, `Ringkasan_SSA`, dan `Data_SSA_150` di workbook memiliki hasil yang berbeda dari dataset geometri utama. Sheet tersebut tetap dipertahankan sebagai rekaman suplemen, tetapi bukan sumber angka utama laporan. Perbandingan rinci ada di `AUDIT_DATA_TAMBAHAN.md`.
- Workbook bukan bukti bahwa 90 kasus LBM versi revisi telah dihitung ulang. Batch historis dan pilot program revisi harus tetap dipisahkan.
