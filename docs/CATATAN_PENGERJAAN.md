# Catatan Pengerjaan

Repositori ini mengarsipkan hasil proyek REV dan anisotropi Kirby. Alasan metodologis utama:

1. Sampel Kirby dipilih karena nisbah permeabilitas domain penuh antararah relatif dekat satu (sekitar 1,0695).
2. Citra yang tersedia sudah biner; tidak dilakukan segmentasi ulang dan tidak dapat dilakukan pemeriksaan artefak grayscale dari file ini.
3. Ukuran subvolume geometri adalah 32, 48, 64, 96, dan 128 voxel, dengan 30 sampel per ukuran.
4. Ambang operasional REV adalah CV < 5%, tetapi kesimpulan juga membaca kemiringan log-CV–log-L dan diagnostik tumpang tindih.
5. L=128 secara nominal memenuhi ambang CV untuk porositas dan SSA terkoreksi, tetapi overlap rata-rata 33,55% menandakan dependensi spasial yang kuat.
6. LBM historis berisi 90 kasus; kasus dengan variasi fluks tinggi dan non-percolating tidak disembunyikan. Batch ini tidak memiliki semua diagnostik yang ada di program revisi.
7. Pilot program revisi yang belum lolos pemeriksaan fluks tidak digabung dengan batch historis.
8. CSV suplemen berbeda dari dataset geometri utama; audit rinci berada di `AUDIT_DATA_TAMBAHAN.md`.

Eksekusi, validasi keluaran, keputusan analisis, dan interpretasi tetap menjadi tanggung jawab penulis.
