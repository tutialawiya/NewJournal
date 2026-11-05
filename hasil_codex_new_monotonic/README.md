# hasil_codex_new_monotonic — BAMNET 2025

Eksperimen baru menggunakan **data 1 menit**, dengan hanya syarat monoton yang diganti. Katalog dan laporan lama tidak ditimpa.

## Unduh berkas utama

- [ZIP semua hasil baru](BAMNET_New_Monotonic_All_Results.zip).
- [Excel semua deteksi baru](detection_1min/BAMNET_semua_event.xlsx).
- [PDF meeting bahasa Indonesia](Meeting_New_Monotonic_ID.pdf).
- [English meeting PDF](Meeting_New_Monotonic_EN.pdf).
- [Excel analisis meeting](Analisis_Meeting_BAMNET.xlsx).
- [Komposit 310 cold pool, enam panel](plots/composite_6panel_310.png).
- [Kode deteksi](coldpool_flexible.py), [kode analisis](analyze.py), [kode finalisasi dan PDF](finalize_reports.py).
- [English explanation](README_EN.md).

Pada halaman file GitHub, pilih **Download raw file** untuk mengunduh Excel/PDF/ZIP.

## Definisi mean dan empat langkah monoton

Untuk setiap calon onset t pada grid **1 menit**, bentuk lima mean dari suhu asli:

| Pusat mean | Sampel suhu 1 menit yang dirata-ratakan |
|---|---|
| t−5 | t−7, t−6, t−5, t−4, t−3 |
| t | t−2, t−1, t, t+1, t+2 |
| t+5 | t+3, t+4, t+5, t+6, t+7 |
| t+10 | t+8, t+9, t+10, t+11, t+12 |
| t+15 | t+13, t+14, t+15, t+16, t+17 |

Keempat selisih mean berturut-turut harus **≤0**; nilai tetap diperbolehkan. Mean harus memiliki lima nilai valid, tanpa interpolasi. Inilah konvensi **mean terpusat 5 menit** yang dipakai dalam eksperimen ini. Mean tidak diambil dari variabel suhu 5 menit yang sudah ada; input tetap `air_temperature_1min`. Onset tetap dievaluasi setiap menit, bukan hanya tiap kelipatan 5 menit.

**Catatan penting:** empat selisih memerlukan lima mean. Lima mean dengan lima sampel yang tidak overlap memerlukan 25 sampel sumber, bukan hanya 21 sampel. Dukungan mean ini meluas dari **t−7 hingga t+17**. Ini pilihan implementasi eksplisit; bukan klaim bahwa seluruh penghalusan berada di dalam t−5..t+15. Grafik definisi dan contoh numeriknya tersedia pada `plots/new_monotonic_definition_ID.png` dan `five_mean_example_values.csv`.

## Bagian yang tidak berubah

- Penurunan total suhu **raw 1 menit**: T(t+15)−T(t−5) ≤−2°C.
- Penurunan awal suhu **raw 1 menit**: T(t)−T(t−5) ≤−0,5°C.
- Gust maksimum data 1 detik pada bin **1 menit**, Eq. (2) dan Eq. (3) bersyarat, jendela dan ambang angin, pengelompokan episode, overlap kategori, pemrosesan per bulan, dan mode `legacy` tetap sama.
- Semua durasi/shift tetap dikonversi melalui `Timing.steps()`.
- Karena hanya monoton yang diganti, **mean lima menit tidak digunakan untuk kedua ambang penurunan suhu raw**.

## Hasil deteksi

| Stasiun | Cold pool | Wind-only | Temperature-only legacy |
|---|---:|---:|---:|
| BAMNET01 | 92 | 763 | 46 |
| BAMNET02 | 66 | 479 | 23 |
| BAMNET03 | 79 | 877 | 43 |
| BAMNET04 | 42 | 934 | 48 |
| BAMNET05 | 31 | 947 | 45 |
| **Total** | **310** | **4.000** | **205** |

Total **4.515 event**. Dari 205 suhu saja, **57 adalah artefak mask kosong mode legacy** dan 148 memiliki kriteria termal valid. Semuanya tetap ditampilkan dalam Excel, dengan kolom flag dan sheet QA; interpretasi meteorologi utama hanya menggunakan kelompok valid.

Sebagai pembanding:

| Metode | Cold pool | Wind-only | Temperature-only legacy |
|---|---:|---:|---:|
| 1 menit lama, 20 selisih raw | 187 | 4.117 | 148 |
| **1 menit baru, 4 selisih mean** | **310** | **4.000** | **205** |
| 5 menit asli | 197 | 2.294 | 178 |

Jumlah lebih besar **tidak otomatis lebih akurat**; ini hasil sensitivitas terhadap definisi monoton.

## Perbandingan event dan analisis angin

Matching satu-ke-satu memerlukan stasiun sama dan toleransi onset. Pada ±20 menit:

- Baru versus **5 menit asli**: **189** pasangan cold pool, yaitu 61,0% dari 310 baru dan 95,9% dari 197 versi 5 menit.
- Baru versus **1 menit lama**: **187** pasangan cold pool. Seluruh event lama memiliki pasangan waktu dalam versi baru, tetapi ini bukan klaim semua timestamp persis sama.
- Sensitivitas toleransi dan perpindahan kategori disimpan dalam CSV/Excel; QA dikeluarkan dari matching fisik utama.

Untuk **148 temperature-only valid**, mean angin pada t..t+20 dibanding background t−40..t−20, dengan minimal 17 dari 21 sampel valid per jendela:

- 81 naik, 61 turun, 1 tetap, 5 data kurang.
- Median Δspeed **+0,110 m/s**; median dihitung pada 143 event dengan data cukup.
- 131 event memiliki arah vektor pre/post andal (magnitude masing-masing ≥0,5 m/s); median rotasi absolut **96,2°**. Pembagian arah: 45 hampir searah, 50 menyilang, 36 hampir berlawanan untuk **arah angin sesudah versus sebelum**, bukan arah propagasi cold pool.
- Spearman baseline dengan perubahan mean speed **−0,489**; untuk perubahan puncak **−0,117**. Perubahan mean mengurangkan baseline, sehingga mathematical coupling/regression to the mean perlu diperhatikan. Korelasi bukan bukti sebab-akibat.

Tabel per event mencakup mean speed/gust, arah dari rata-rata u/v, delta komponen, rotasi, serta orientasi vektor perubahan terhadap background. Arah lokal dan vektor perubahan tidak boleh langsung dianggap arah gerak cold pool.

## Komposit dan normalisasi

Komposit baru memakai **310 event**, jendela −40..+120 menit, enam panel: speed, gust, gust−speed; masing-masing absolut dan versi **X_i(t)−X_i(0)**. Normalisasi dilakukan per event sebelum menghitung mean hitam dan ±1 SD populasi antar-event; satuan tetap m/s, bukan persentase/z-score.

Nilai kosong tidak diinterpolasi. Jumlah data tersedia per menit tersimpan dalam `composite_mean_std_coverage.csv`. Versi zoom hanya membatasi tampilan, tidak membuang nilai dari statistik. Flag ekstrem dalam `extreme_sample_review_flags.csv` memakai speed>20 atau gust>30 m/s sebagai penanda tinjauan, bukan batas penolakan otomatis.

Pada 308 cold pool dengan jendela berpasangan yang layak, median kenaikan puncak speed **+2,265 m/s**, gust **+3,888 m/s**. Sampel dipilih memakai gust, sehingga amplitudo lebih besar belum membuktikan gust lebih akurat. Pilihan sementara tetap gust untuk konsistensi metode, dibantu speed/u/v untuk respons dan arah; evaluasi akurasi memerlukan label independen.

## Validasi

- Tiga uji terkendali lulus, termasuk pemeriksaan AST bahwa algoritma yang berubah hanya fungsi suhu; fungsi run hanya memperoleh metadata tambahan.
- **458 event termal valid** (310 cold pool +148 suhu saja) diperiksa ulang langsung pada data NetCDF: kedua penurunan raw tetap lolos dan keempat selisih mean tidak positif.
- Array komposit berukuran 310×161. Cakupan, normalisasi onset, jumlah QA, dan hasil matching tersedia pada `validation_report.json`.
- Excel deteksi dibuka ulang: total 4.515, cold pool 310, suhu saja 205.

## Jalankan ulang

Gunakan Python 3.12 dan `requirements.txt`:

```powershell
python -m pip install -r requirements.txt
python -m unittest -v test_new_monotonic
python coldpool_flexible.py --data-root "C:\Users\bamnet\Downloads\TUTI\THESIS\OlahDataBamnet\Level2Data_New" --output "detection_1min" --resolution-minutes 1 --boundary-policy legacy
python analyze.py --data-root "C:\Users\bamnet\Downloads\TUTI\THESIS\OlahDataBamnet\Level2Data_New" --catalog-root "catalogs" --output "."
python finalize_reports.py --data-root "C:\Users\bamnet\Downloads\TUTI\THESIS\OlahDataBamnet\Level2Data_New"
```

`catalogs/1min/all_events.csv` adalah snapshot deteksi baru yang dipakai untuk analisis. Jika mengubah/mengulang deteksi, perbarui snapshot tersebut dengan `detection_1min/all_events.csv` sebelum menjalankan analisis. Data NetCDF tidak digandakan dalam folder hasil.

## Batas interpretasi

Belum ada verifikasi satelit Hoeller; seluruh cold pool masih kandidat stasiun. Cakupan data tidak sama antarstasiun, BAMNET02 hanya memiliki sembilan file bulanan. Rolling/grouping tetap berbasis langkah sampel seperti versi awal sehingga dapat melewati celah data. Metrik `duration_min` tetap dari fungsi lama dan belum tervalidasi sebagai durasi fisik. Tidak ada perubahan data sumber, interpolasi, ataupun penghapusan outlier otomatis.
