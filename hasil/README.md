# Hasil dan keluaran

Folder ini memisahkan tiga kelompok hasil agar sumber angka tetap dapat ditelusuri:

1. **Geometri dan domain** — hasil domain penuh, porositas, koordinat subvolume, CV, SSA, granulometri, overlap, tortuositas, dan verifikasi kanal yang dihasilkan oleh kode proyek.
2. **Batch LBM historis** — workbook 90 kasus dari notebook awal. Statusnya adalah 34 lolos pemeriksaan notebook lama, 46 variasi fluks antaririsan tinggi, dan 10 non-percolating.
3. **Pilot program revisi** — hasil uji dari program dengan pemeriksaan tambahan. Pilot disimpan terpisah dan tidak digabung dengan batch historis.

CSV suplemen yang statistiknya berbeda dari dataset geometri utama berada di `arsip_suplemen_perlu_rekonsiliasi/`. Rincian perbandingan ada di `docs/AUDIT_DATA_TAMBAHAN.md`.

Jangan menghapus kasus non-percolating atau kasus dengan diagnostik tidak lolos. Pertahankan statusnya saat menghitung statistik dan menafsirkan hasil.

## Workbook hasil final yang diunggah

`Hasil_Final_REV_Anisotropi_Kirby.xlsx` memuat sheet `Data_LBM_90` (90 catatan), `Data_LBM_Valid` (80 kasus berangka/percolating, bukan berarti semuanya lolos ambang fluks), ringkasan permeabilitas/anisotropi, serta sheet porositas dan SSA suplemen. Pada data LBM, kolom `status` adalah penentu: 34 `Lolos pemeriksaan`, 46 `Variasi irisan tinggi`, dan 10 `Non-percolating`. Sheet porositas/SSA di workbook tidak sama dengan dataset geometri utama; lihat `docs/AUDIT_DATA_TAMBAHAN.md`.
