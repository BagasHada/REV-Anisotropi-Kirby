# Catatan proyek

## Data dan konvensi
- Sampel: Kirby, citra biner 256³ voxel, ukuran voxel 2,25 µm.
- Nilai 0 = pori; nilai 1 = padatan.
- Domain penuh memiliki sisi 576 µm.
- Citra RAW berasal dari bundel resmi mata kuliah dan tidak disertakan dalam repository publik; lihat `data/citra/README.md`.
- `phi_flow` P3 (tanpa periodisitas) berbeda dari masker pori aktif P11 pada LBM. Jangan memperlakukannya sebagai definisi identik.

## Hasil geometri
- Porositas total domain: 20,5284%; porositas alir P3: 20,3326%.
- CV porositas per L=32,48,64,96,128: 64,18%; 27,01%; 20,23%; 10,36%; 4,84%.
- CV SSA terkoreksi: 36,55%; 22,10%; 13,27%; 5,53%; 3,82%.
- Kemiringan log-CV/log-L: −1,768 dan −1,703.
- Rata-rata overlap pada L=128: 33,55%; interpretasi REV geometri tetap nominal/indikatif.

## LBM dan anisotropi
- Batch historis: 90 catatan; 34 lolos kriteria notebook lama, 46 variasi fluks irisan tinggi, 10 non-percolating.
- Pada L=128, 18/18 kasus lolos kriteria notebook lama.
- CV permeabilitas pada L=128 (X/Y/Z): 13,06% / 10,79% / 28,85%; REV permeabilitas belum tercapai.
- Rasio k maksimum/minimum 1,179; indeks anisotropi 0,526. Belum ada bukti kuat bahwa beda arah terpisah dari variasi posisi.
- Hasil batch lama belum memiliki seluruh diagnostik program revisi. Pilot baru disimpan terpisah.

## Aturan membaca keluaran
1. Jangan menghapus kasus non-percolating. Itu merupakan informasi konektivitas, bukan otomatis kegagalan solver.
2. Jangan memasukkan kasus variasi fluks irisan tinggi ke statistik hasil lolos.
3. Jangan mencampur CSV suplemen yang belum direkonsiliasi dengan dataset utama.
4. Jangan menyebut ekstrapolasi ukuran REV sebagai hasil simulasi langsung.
5. Baca `docs/STATUS_PROYEK.md`, `docs/CATATAN_KETIDAKPASTIAN.md`, dan `docs/AUDIT_DATA_TAMBAHAN.md` sebelum menginterpretasikan hasil.
