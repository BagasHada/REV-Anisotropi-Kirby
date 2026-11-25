# Audit CSV tambahan terhadap dataset geometri utama

Tiga CSV berikut diunggah sebagai data tambahan. Pemeriksaan menunjukkan nilainya tidak sama dengan keluaran geometri yang menjadi dasar tabel dalam ringkasan hasil utama. Karena asal/batch dan koordinat pembentuknya belum dikonfirmasi, file dipertahankan sebagai arsip dan tidak menggantikan angka dataset geometri utama.

## 1. Ringkasan porositas tambahan

CV dari `ringkasan_porositas_kirby.csv`:

|   ukuran_voxel |   CV_persen |
|---------------:|------------:|
|             32 |    47.3292  |
|             48 |    36.444   |
|             64 |    26.3874  |
|             96 |    11.0308  |
|            128 |     5.06656 |

CV porositas yang dipakai di dataset geometri utama (Tabel 3): 64,18%; 27,01%; 20,23%; 10,36%; 4,84% untuk L=32, 48, 64, 96, 128 voxel.

Perbedaan paling jelas terdapat pada L=32 (47,33% vs 64,18%) dan L=48 (36,44% vs 27,01%). Karena itu kedua dataset tidak boleh digabung atau dipertukarkan.

## 2. Data SSA tambahan

Dari `data_ssa_subvolume_kirby.csv`, statistik sederhana kolom `SSA_koreksi_um_inv` adalah:

|   ukuran_voxel |   count |      mean |        std |   cv_pct |
|---------------:|--------:|----------:|-----------:|---------:|
|             32 |      30 | 0.0304531 | 0.00944054 | 31.0003  |
|             48 |      30 | 0.0307105 | 0.00723215 | 23.5494  |
|             64 |      30 | 0.0316181 | 0.0050877  | 16.0911  |
|             96 |      30 | 0.0311105 | 0.00242803 |  7.80454 |
|            128 |      30 | 0.0320935 | 0.00150646 |  4.69396 |

CV SSA terkoreksi yang dipakai di dataset geometri utama (Tabel 3): 36,55%; 22,10%; 13,27%; 5,53%; 3,82%. Nilai CSV tambahan juga berbeda sehingga dipertahankan terpisah.

## 3. Diagnostik overlap tambahan

Rata-rata overlap dari `diagnostik_overlap_kirby.csv`:

|   ukuran_voxel |   rata_rata_overlap_persen |
|---------------:|---------------------------:|
|             32 |                   0.128805 |
|             48 |                   0.642353 |
|             64 |                   3.37478  |
|             96 |                  13.1904   |
|            128 |                  31.5465   |

Rata-rata irisan pasangan yang dipakai di dataset geometri utama (Tabel 5): 0,36%; 1,02%; 2,89%; 12,99%; 33,55%. Karena angkanya tidak sama, CSV tambahan ini tidak menjadi sumber Gambar 3/Tabel 5 pada laporan.

## Keputusan dokumentasi

- File CSV tambahan tidak dihapus; berkas disimpan di `hasil/arsip_suplemen_perlu_rekonsiliasi/`.
- Nilai dataset geometri utama tetap mengikuti keluaran pada `hasil/` yang telah menjadi dasar laporan.
- Untuk menggabungkan atau mengganti salah satu sumber, asal batch, seed, koordinat subvolume, dan metode penghitungan perlu dicocokkan terlebih dahulu.

## 4. Sheet geometri tambahan dalam workbook hasil final

Workbook `hasil/Hasil_Final_REV_Anisotropi_Kirby.xlsx` mencakup sheet `Ringkasan_Porositas`, `Ringkasan_SSA`, dan `Data_SSA_150`. Nilai CV dari sheet ringkasan workbook juga berbeda dari dataset geometri utama. Contohnya:

- CV porositas pada L=32 adalah 47,33%, sedangkan sumber utama 64,18%; pada L=128 adalah 5,07%, sedangkan sumber utama 4,84%.
- CV SSA terkoreksi pada L=32 adalah 31,00%, sedangkan sumber utama 36,55%; pada L=128 adalah 4,69%, sedangkan sumber utama 3,82%.

Sheet tersebut dipertahankan apa adanya sebagai data suplemen dari workbook, bukan dipakai untuk mengganti keluaran geometri utama. Koordinat/seed dan asal batch perlu direkonsiliasi sebelum angka-angkanya dapat disatukan.
