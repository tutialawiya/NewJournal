# Hasil deteksi dan analisis BAMNET 2025

Folder ini berisi hasil yang dibuat Codex dan diunggah dengan persetujuan pengguna. Data sumber tetap di folder BAMNET01–05 dan tidak diubah.

## Unduh berkas utama

Klik nama file, lalu klik **Download raw file** (ikon unduh) di kanan atas tampilan file GitHub. Login ke akun GitHub yang memiliki akses jika repository bersifat privat.

- [Excel full detection 5 menit](BAMNET_2025_Deteksi_5menit.xlsx): **197 cold pool, 2.294 wind-only, 178 temperature-only**.
- [Excel full detection 1 menit](BAMNET_2025_Deteksi_1menit.xlsx): **187 cold pool, 4.117 wind-only, 148 temperature-only**.
- [PDF bahan meeting, 10 halaman](Bahan_Meeting_BAMNET.pdf).
- [English meeting briefing, 11 pages](BAMNET_Meeting_Briefing_English.pdf): all headings and explanatory text in English, including a dedicated explanation of onset-centered normalization, X(t) − X(0).
- [Excel analisis overlap, temperature-only, arah, baseline, dan komposit](Analisis_Meeting_BAMNET.xlsx).
- [ZIP lengkap analisis meeting: kode, tabel, grafik, dan laporan](BAMNET_Analisis_Meeting_Lengkap.zip).
- [Grafik komposit 6 panel, 187 cold pool](Komposit_6panel_187.png).
- [Laporan penjelasan lengkap](LAPORAN_MEETING.md).
- [Kode deteksi fleksibel 1/5 menit](coldpool_flexible.py) dan [dependensi](requirements.txt).

## Hal penting saat membaca hasil

Excel deteksi mempertahankan mode `legacy` sesuai permintaan pengguna. Masing-masing resolusi menyertakan **57 temperature-only dari mask kosong akhir bulan**, yang ditandai di sheet `QA_MASK_KOSONG`. Dengan demikian 148 temperature-only 1 menit terdiri dari 91 event termal valid dan 57 QA; 178 temperature-only 5 menit terdiri dari 121 valid dan 57 QA. Kelompok QA tidak digunakan untuk interpretasi meteorologi utama dalam analisis meeting.

Cold pool masih kandidat berdasarkan stasiun; verifikasi satelit belum dilakukan. Overlap utama memakai toleransi onset ±20 menit, stasiun sama, dan pasangan one-to-one. Komposit memakai jendela −40..+120 menit, garis mean hitam, ±1 SD antar-event, serta normalisasi X(t)−X(0). Gust menunjukkan respons puncak lebih besar, tetapi superioritas akurasi belum terbukti karena sampel dipilih menggunakan gust.

File `SHA256.json` menyimpan checksum berkas hasil. Penjelasan metode, cakupan data, dan batas interpretasi tersedia pada laporan dan PDF.

## Jalankan detektor 5 menit

Dengan Python 3.12, dari folder ini:

```powershell
python -m pip install -r requirements.txt
python coldpool_flexible.py --data-root "C:\Users\bamnet\Downloads\TUTI\THESIS\OlahDataBamnet\Level2Data_New" --output "hasil_deteksi_5menit" --resolution-minutes 5 --boundary-policy legacy
```

Untuk data 1 menit, ganti `--resolution-minutes 1` dan gunakan folder output lain.
