# Ringkasan hasil dan status proyek

## Data dan domain
- Sampel Kirby: 256³ voxel, `uint8`, nilai 0 = pori, nilai 1 = padatan, ukuran voxel 2,25 µm.
- Porositas total domain: 20,5284%; porositas alir P3: 20,3326%.
- Porositas alir P3 tanpa periodisitas berbeda dari masker pori aktif P11 pada LBM.

## Granulometri dan geometri REV
- Granulometri: r50 = 4,821 voxel (10,847 µm); diameter pengendali perkiraan = 21,693 µm; d_t/Δx = 9,64 voxel.
- Geometri: 30 subvolume per ukuran pada L = 32, 48, 64, 96, dan 128 voxel.
- CV porositas: 64,18%; 27,01%; 20,23%; 10,36%; 4,84%.
- CV SSA terkoreksi: 36,55%; 22,10%; 13,27%; 5,53%; 3,82%.
- Kemiringan log CV–log L: −1,768 (porositas) dan −1,703 (SSA terkoreksi).
- Secara nominal kedua sifat melewati ambang CV < 5% pada L = 128 voxel (288 µm), tetapi rata-rata overlap pasangan 33,55% (median 32,35%) membatasi independensi subvolume. Kesimpulan REV geometri bersifat indikatif.

## Tortuositas dan verifikasi solver
- Tortuositas geometri Dijkstra pada satu subvolume 128³: X = 1,401; Y = 1,445; Z = 1,427.
- Verifikasi kanal lurus: galat relatif 0,83%; 3,77%; dan 5,13%; ketiganya mencapai kriteria residual pada uji tersebut.

## Batch LBM historis
- Total 90 catatan: 34 lolos pemeriksaan versi notebook, 46 memiliki variasi fluks irisan tinggi, dan 10 non-percolating.
- Pada L = 128, seluruh 18 kasus lolos pemeriksaan versi notebook terdahulu.
- Rerata permeabilitas L = 128: X = 482,25 mD; Y = 475,20 mD; Z = 560,11 mD.
- CV L = 128: X = 13,06%; Y = 10,79%; Z = 28,85%; REV permeabilitas belum tercapai.
- Rasio k maksimum/minimum = 1,179; indeks anisotropi = 0,526. Perbedaan arah belum terpisah jelas dari variasi posisi menurut patokan yang digunakan.
- Batch historis belum menyimpan semua diagnostik tambahan program revisi, termasuk Reynolds pori dan kesalahan massa.

## Pilot program revisi
Pilot LBM terbaru tersimpan di `hasil/pilot_program_revisi/`. Salah satu uji L = 32 arah X memiliki residual temporal kecil, tetapi variasi fluks irisan 0,953 (ambang 0,001), sehingga tidak lolos pemeriksaan fluks. Pilot tersebut tidak digabung dengan statistik batch historis.

## CSV tambahan
Tiga CSV suplemen yang dikirim terpisah memiliki nilai statistik berbeda dari dataset geometri utama. File dipertahankan dalam folder arsip dan dijelaskan di `AUDIT_DATA_TAMBAHAN.md`; jangan digabung sebelum asal batch, seed, koordinat, dan metode diperiksa.
