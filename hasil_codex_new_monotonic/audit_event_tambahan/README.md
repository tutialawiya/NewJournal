# Audit 121 cold pool tambahan: 1 menit monoton baru vs 5 menit

## Unduh atau baca

- [Excel daftar lengkap dan nilai setiap kriteria](New_121_ColdPool_Event_Audit.xlsx).
- [Daftar 121 event langsung di GitHub](DAFTAR_121_EVENT.md).
- [CSV seluruh audit](new_121_event_audit.csv).
- [Flag suhu/angin setiap menit di sekitar onset](criterion_flags_near_new_onsets.csv).
- `summary.json`, `validation_report.json`, dan `audit.py`: ringkasan, pemeriksaan, dan kode audit.

## Definisi kandidat tambahan

Dari 310 cold pool 1 menit monoton baru, 189 berpasangan dengan cold pool 5 menit pada stasiun sama dalam toleransi onset ±20 menit. **121 sisanya** adalah objek audit ini. Kata tambahan tidak membuktikan sebuah kejadian fisik baru; pergeseran onset, pengelompokan, dan aturan pencocokan juga berpengaruh.

Label versi 5 menit diambil dari **pencocokan lintas kategori satu-ke-satu** yang sudah digunakan pada analisis sebelumnya. Hasilnya:

| Label 5 menit yang berpasangan | Jumlah | Bagian pembeda |
|---|---:|---|
| Temperature-only | **44** | Suhu sudah lolos pada 5 menit, tetapi konfirmasi angin/overlap tidak lolos. Pada 1 menit, episode suhu kini dikonfirmasi angin. |
| Wind-only | **5** | Angin sudah lolos pada 5 menit, tetapi tidak overlap kandidat suhu. Pada 1 menit, kandidat suhu lolos dan overlap angin. |
| Tidak ada pasangan event dalam ±20 menit | **72** | Tidak ada kandidat suhu 5 menit yang lolos dalam ±20 menit sekitar onset baru. Resolusi suhu, grid onset, ambang raw, monoton, atau data hilang dapat menjadi pembeda. |

Tidak ada pasangan bukan berarti stasiun tidak merekam perubahan cuaca, atau tidak ada episode angin yang mulai lebih awal. Kolom `nearest_5min_*_context` memberikan event lama terdekat sebagai konteks; event tersebut **bukan otomatis pasangan fisik yang sama**.

Per stasiun: BAMNET01 30, BAMNET02 21, BAMNET03 30, BAMNET04 22, BAMNET05 18.

## Di bagian mana event baru menjadi cold pool?

Setiap baris berisi:

1. `new_onset_UTC`: awal episode suhu pada katalog baru.
2. `first_joint_confirmation_UTC`: anchor pertama di episode tersebut yang benar-benar memenuhi suhu dan angin bersamaan. Onset awal dan konfirmasi tidak selalu sama.
3. `new_confirm_T20_C`: T(t+15)−T(t−5), harus ≤−2°C.
4. `new_confirm_T5_C`: T(t)−T(t−5), harus ≤−0,5°C.
5. `new_confirm_max_mean_increment_C`: nilai terbesar dari empat selisih mean, harus ≤0. Mean 5 menit berpusat pada t−5,t,t+5,t+10,t+15, dari suhu 1 menit.
6. `new_confirmation_equation`: Eq. (2) atau Eq. (3) yang mengonfirmasi angin.
7. `new_confirm_wind_max_anomaly` dan `new_confirm_wind_threshold_eq2`: untuk Eq. (2), maksimum anomali harus ≥ ambang.
8. `new_confirm_wind_max_speed` dan `new_confirm_wind_threshold_eq3`: untuk Eq. (3), maksimum gust harus ≥ ambang fallback; Eq. (3) tetap bersyarat anomali hilang.

Seluruh **121 event memiliki anchor nyata yang memenuhi kedua kriteria**, bukan hanya overlap rentang tanpa anchor bersama. Pada anchor konfirmasi pertama, 115 memakai Eq. (2), 6 memakai Eq. (3). Nilai raw suhu, mean monoton, dan perbandingan ambang angin semuanya diperiksa numerik.

## Contoh angka

**NEW-001 — BAMNET03, 4 Januari 2025 17:52 UTC**:

- Pada versi 5 menit, pasangan terdekat yang ditetapkan adalah **temperature-only 17:50**.
- Suhu versi 5 menit sudah lolos: ΔT20=−2,58°C, ΔT5=−0,70°C, monoton lolos.
- Tetapi maksimum anomali angin 5 menit **2,386 m/s < ambang 2,586 m/s**, sehingga angin gagal di anchor tersebut.
- Pada anchor konfirmasi versi 1 menit, ΔT20=−2,72°C, ΔT5=−0,62°C dan monoton mean lolos.
- Maksimum anomali angin **3,150 m/s ≥ ambang 2,761 m/s**, sehingga Eq. (2) mengonfirmasi cold pool.

Ini menunjukkan tambahan event tidak semuanya disebabkan oleh penggantian monoton: resolusi gust, statistik angin, dan alignment onset juga berbeda ketika membandingkan **5 menit versus 1 menit**.

Contoh tanpa pasangan, **NEW-003 — BAMNET02, 7 Februari 2025 06:34 UTC**: pada grid 5 menit terdekat 06:35, penurunan awal hanya −0,40°C sehingga gagal syarat 0,5°C, sedangkan versi 1 menit di anchor konfirmasi memiliki −0,61°C dan angin juga lolos.

## Cara membaca pembandingan lama

`diagnostic_5min_anchor_UTC` menggunakan onset pasangan lama jika ada; bila tidak ada pasangan, memakai grid 5 menit terdekat. Ini perbandingan diagnostik pada anchor yang **dapat berbeda waktunya**, bukan eksperimen terkontrol yang hanya mengubah satu komponen.

Pada 72 baris tanpa pasangan, di anchor 5 menit diagnostik terdapat 49 kegagalan penurunan total 2°C, 15 kegagalan penurunan awal 0,5°C, dan 21 kegagalan monoton. Kategori kegagalan dapat bertumpang tindih. Jumlah ini tidak berarti 49/15/21 penyebab unik atau pembuktian kausal; indikator paling kuat adalah tidak adanya kandidat suhu 5 menit yang lolos di jendela lokal ±20 menit.

Untuk benar-benar mengisolasi efek monoton saja, gunakan pembandingan **1 menit lama vs 1 menit monoton baru**, bukan mengatribusikan seluruh perbedaan 5 menit vs 1 menit pada monoton.

## Validasi dan batasan

Sheet Excel: `ALL_121_ADDITIONAL`, `PREV_TEMPERATURE_44`, `PREV_WIND_5`, `NO_PAIRED_EVENT_72`, `SUMMARY`, dan `METHOD`. Semua waktu UTC, dengan onset WIB juga tersedia. Tidak ada data sumber atau katalog lama yang diubah. Ini audit algoritma dan pencocokan waktu, bukan verifikasi meteorologi independen atau bukti bahwa semua kandidat tambahan merupakan cold pool fisik baru.
