# Bahan meeting: perbandingan resolusi dan variasi angin BAMNET 2025

Analisis memakai hasil detektor terbaru untuk BAMNET01–05 tahun 2025, tanpa mengubah katalog deteksi. Seluruh waktu memakai UTC. Semua hasil merupakan kandidat berbasis stasiun, belum terverifikasi satelit.

## 1. Arti angka di output

| Resolusi | Cold pool | Temperature-only, legacy | Wind-only | Total |
|---|---:|---:|---:|---:|
| 1 menit | 187 | 148 | 4.117 | 4.452 |
| 5 menit | 197 | 178 | 2.294 | 2.669 |

- **Cold pool**: episode suhu dan episode angin overlap sesuai kode pengguna; angin lolos Eq. (2) atau Eq. (3) bersyarat. Bukan jumlah stasiun atau jumlah sampel waktu.
- **Temperature-only**: episode suhu tidak overlap dengan episode angin. Ini berarti gagal konfirmasi angin algoritmik, **bukan berarti angin tidak berubah**.
- **Wind-only**: episode angin tidak overlap dengan episode suhu. Ini **bukan 4.117 cold pool tambahan**.
- Mode legacy menyertakan perilaku mask kosong akhir bulan. Dari 148 temperature-only 1 menit, **57 adalah event QA dari mask kosong dan 91 memiliki kriteria suhu valid**. Pada resolusi 5 menit, 178 terdiri dari 57 QA dan 121 valid. Hitungan asli tetap ditampilkan; inferensi meteorologi utama memakai kelompok valid.

## 2. Overlap 1 menit versus 5 menit

Pencocokan memerlukan stasiun sama dan selisih onset dalam toleransi. Algoritma assignment memilih jumlah pasangan maksimum, lalu meminimalkan total selisih waktu. Setiap event boleh dipakai sekali; event satu tidak dipasangkan ke banyak event lain. Toleransi utama **±20 menit** dipilih sebagai skala asosiasi praktis yang sejalan dengan jarak penggabungan event, bukan sebagai kepastian identitas fisik. Sensitivitas 0/5/10/20/30 menit juga dihitung.

Untuk kategori yang sama dan QA dikeluarkan:

| Kategori | N 1 menit | N 5 menit | Pasangan ±20 menit | % event 1 menit berpasangan | % event 5 menit berpasangan |
|---|---:|---:|---:|---:|---:|
| Cold pool | 187 | 197 | **122** | **65,2%** | **61,9%** |
| Temperature-only valid | 91 | 121 | **56** | **61,5%** | **46,3%** |
| Wind-only | 4.117 | 2.294 | **1.472** | **35,8%** | **64,2%** |

Cold pool: 65 event 1 menit dan 75 event 5 menit tidak memiliki pasangan **cold-pool kategori sama** pada toleransi ini; bukan otomatis bukti false positive. Jaccard cold pool = 122/(187+197−122) = **0,466**. Timestamp persis hanya cocok pada **26** cold pool; toleransi ±5 menit memberi 110, ±10 menit memberi 115, dan ±30 menit tetap 122. Karena itu membandingkan angka total 187 versus 197 saja tidak cukup. Artefak QA sangat memengaruhi overlap temperature-only: 64 pasangan timestamp persis pada katalog legacy terdiri dari 57 pasangan QA dan hanya 7 pasangan termal valid.

Pencocokan terpisah yang memperbolehkan perubahan kategori menghasilkan tabel berikut (baris = kategori 1 menit, kolom = kategori 5 menit). Ini assignment baru, bukan tabel yang harus dijumlahkan dengan pasangan kategori sama:

| 1 menit → 5 menit | Cold pool | Temperature-only valid | Wind-only |
|---|---:|---:|---:|
| Cold pool | 122 | 29 | 3 |
| Temperature-only valid | 5 | 56 | 0 |
| Wind-only | 2 | 1 | 1.472 |

Artinya sejumlah kejadian muncul pada kedua resolusi dengan label berbeda. Penghalusan suhu, definisi maksimum gust 1/5 menit, dan pergeseran onset semuanya dapat berperan. Hasil assignment hanya menunjukkan kedekatan waktu, bukan pembuktian kejadian fisik yang sama.

Gambar: `plots/overlap_sensitivity.png`. Seluruh pasangan tersedia dalam Excel dan CSV, termasuk signed offset.

## 3. Menjelaskan prompt dan kode monoton

Permintaan awal adalah menerapkan metode Hoeller. Tidak ada instruksi terpisah dari pengguna yang menentukan toleransi kenaikan suhu kecil. Implementasi menafsirkan **monotonic decrease** sebagai **setiap selisih suhu berturut-turut ≤0**, sesuai pembacaan metode dan operator yang juga ada pada kode asli yang kemudian dikirim pengguna.

Kode aktual hasil refactor:

```python
window = temperature.rolling(time=before + after + 1).construct('window')
monotonic = (window.diff('window') <= 0).all('window')
```

Maknanya **non-increasing**: nilai tetap diperbolehkan, kenaikan sekecil apa pun menggagalkan monoton. Ini berbeda dari regresi berslope negatif, suhu akhir lebih rendah dari awal, atau tren menurun secara umum. Tidak ada smoothing atau toleransi kenaikan 0,01°C yang diterapkan.

Di 5 menit, uji mencakup 5 suhu dan 4 selisih dalam 20 menit. Di 1 menit, uji mencakup 21 suhu dan 20 selisih. Ambang lain tetap: penurunan ≥2°C dari t−5 ke t+15, dan penurunan ≥0,5°C dari t−5 ke t. Kode 1 menit membandingkan selisih **5 menit**, bukan menuntut 0,5°C setiap menit.

Contoh nyata: **BAMNET01, 6 Januari 2025 07:15 UTC**, memakai anchor hasil 5 menit. Suhu 1 menit naik dari **28,91 ke 28,92°C** pada menit −2 ke −1. Ini hanya +0,01°C, tetapi gagal syarat monoton pada anchor tersebut. Seri 5 menit lolos. Seri 5 menit merupakan agregasi tersendiri, bukan sekadar mengambil setiap titik kelima dari seri 1 menit.

Tunjukkan `plots/monotonic_code_explanation.png` dan `plots/monotonic_actual_example.png`. Jangan menyebut contoh ini bukti event fisiknya tidak terdeteksi sama sekali pada resolusi 1 menit; anchor lain mungkin lolos.

Jika dosen bermaksud tren turun yang mengizinkan fluktuasi kecil, itu perubahan definisi ilmiah yang perlu disepakati dan diuji sensitivitasnya. Analisis ini belum mengganti definisi atau mengubah hasil deteksi.

## 4. Analisis temperature-only

Semua **148** event dianalisis dan tersedia di sheet `Temp_all_148`; interpretasi utama memakai **91** event dengan mask suhu valid, tersedia di `Temp_valid_91`.

Definisi pengukuran yang dipakai:

- Background: mean dari t−40 hingga t−20 menit, 21 nilai inklusif.
- Kondisi event: mean dari t hingga t+20 menit, 21 nilai inklusif.
- Δspeed = mean selama event − mean background, satuan m/s.
- Kenaikan puncak sebagai ukuran tambahan: maksimum t..t+20 − mean background.
- Setiap mean membutuhkan sedikitnya **17 dari 21 nilai valid** (≥80%); tanpa interpolasi. Data padding antarbulan digabung, timestamp overlap identik diambil satu kali.
- Naik/turun utama memakai tanda perubahan dengan toleransi numerik 1e−8. Sensitivitas ambang praktis ±0,1 m/s juga dilaporkan.

### Kecepatan naik atau turun?

Untuk 91 event valid:

| Status perubahan mean speed | Jumlah |
|---|---:|
| Naik | **45** |
| Turun | **42** |
| Tetap | **1** |
| Data angin tidak cukup | **3** |

Dari 88 yang dapat dievaluasi, 51,1% naik dan 47,7% turun. Median Δspeed **+0,081 m/s**, mean **+0,143 m/s**, kuartil 25–75% **−0,415 sampai +0,591 m/s**, rentang **−1,589 sampai +4,371 m/s**. Dengan deadband ±0,1 m/s: 44 naik, 37 turun, 7 berubah sangat kecil, 3 data kurang.

Jika semua 148 dihitung tanpa memisahkan QA, hasilnya 68 naik, 68 turun, 1 tetap, dan 11 data kurang. Angka ini diberikan untuk transparansi, tetapi mencampurkan event termal valid dengan artefak batas file.

### Perubahan arah

Mean arah tidak dihitung dengan rata-rata aritmetika derajat: 359° dan 1° seharusnya dekat 0°, bukan 180°. Analisis memakai rata-rata komponen `u_wind_1min` dan `v_wind_1min`, kemudian dikonversi menjadi arah meteorologis. Selisih sudut dibungkus ke [−180°,180°).

Untuk mengurangi arah tidak stabil saat angin lemah, interpretasi arah utama memerlukan magnitude vektor rata-rata background dan event masing-masing ≥0,5 m/s. Ini filter analisis arah, bukan kriteria deteksi baru. Ada **83 event** yang memenuhi syarat tersebut. Median perubahan arah absolut **80,4°**.

| Rotasi arah angin event relatif background | Jumlah | Median Δspeed |
|---|---:|---:|
| 0–<45° | 29 | −0,001 m/s |
| 45–<135° | 37 | +0,110 m/s |
| 135–180° | 17 | +0,337 m/s |

Ini **rotasi angin yang terukur**, bukan sudut arah propagasi cold pool. Kenaikan kecepatan tidak identik dengan sedikit atau banyak rotasi.

Untuk mendekati diskusi vektor, dihitung juga Δvektor = vektor event − vektor background dan sudutnya terhadap vektor background. Dengan background dan magnitude Δvektor ≥0,5 m/s: 4 event hampir searah, 27 menyilang, dan 45 hampir berlawanan; median Δspeed masing-masing **+0,621**, **+0,498**, dan **−0,312 m/s**. Ini diagnosis perubahan lokal, belum ukuran angin outflow murni atau arah gerak cold pool. Efek lain juga dapat menghasilkan Δvektor tersebut.

Secara geometri:

\[
\Delta |U|=\sqrt{|U_0|^2+|C|^2+2|U_0||C|\cos\theta}-|U_0|.
\]

Searah memberi penambahan kecepatan; menyilang juga dapat menaikkan speed tetapi responsnya melemah saat background lebih kuat; berlawanan bisa menurunkan speed, menghentikan, atau membalik arah. Jika perturbasi berlawanan sangat kuat, speed akhirnya bisa naik kembali. Karena itu kategori temperature-only tidak otomatis berarti cold pool lemah atau tidak ada perubahan angin. Tunjukkan `plots/vector_geometry.png`; grafik itu skematik, bukan hasil estimasi arah propagasi.

### Efek baseline

Pada 88 event yang memiliki background dan event-window cukup, korelasi Spearman background speed dengan Δmean speed adalah **ρ=−0,510**, p deskriptif **3,97×10⁻⁷**. Hubungan negatif juga muncul pada ringkasan per stasiun, walaupun ukuran beberapa kelompok kecil.

| Kuartil background speed | N | Median Δspeed | Proporsi naik |
|---|---:|---:|---:|
| Terendah, sekitar 0,409–1,224 m/s | 22 | +0,606 m/s | 86,4% |
| 1,224–1,719 m/s | 22 | +0,169 m/s | 54,5% |
| 1,719–2,267 m/s | 22 | −0,205 m/s | 36,4% |
| Tertinggi, sekitar 2,267–3,917 m/s | 22 | −0,376 m/s | 27,3% |

Namun Δspeed secara matematis sudah mengurangkan baseline. Karena itu korelasi negatif dapat dipengaruhi mathematical coupling dan regression to the mean, selain proses fisik. Untuk perubahan puncak, hubungan jauh lebih lemah: **ρ=−0,202**, p **0,059**. Jadi ada pola yang konsisten dengan efek baseline, tetapi belum pembuktian sebab-akibat atau pengaruh independen arah. P-value menganggap event independen dan bukan koreksi multiple testing; kelompok stasiun kecil perlu kehati-hatian.

Gambar utama: `plots/temperature_only_wind_analysis.png`. Data tiap event, arah, Δu/Δv, cakupan, dan statistik ada dalam Excel/CSV.

## 5. Komposit wind speed versus wind gust

Komposit memakai **187 cold pool 1 menit** dengan jendela **−40 sampai +120 menit**, sehingga setiap seri berisi 161 posisi waktu. Setiap panel memiliki 187 garis event tipis, garis mean hitam tebal, dan shading ±1 simpangan baku populasi antar-event (`ddof=0`), bukan confidence interval mean.

Layout tiga baris × dua kolom:

| Baris | Kiri: absolut | Kanan: dinolkan pada onset |
|---|---|---|
| 1 | Mean wind speed 1 menit | speed(t) − speed(0) |
| 2 | Gust maksimum 1 menit | gust(t) − gust(0) |
| 3 | Gust − speed | [gust−speed](t) − [gust−speed](0) |

Normalisasi di sini adalah **pengurangan nilai onset**, bukan pembagian atau z-score. Satuan tetap m/s. Karena angin dapat memuncak sebelum onset termal dan nilai t=0 sudah naik, kurva setelah dinolkan tidak harus meningkat lagi setelah t=0. Untuk keputusan speed versus gust, metrik background-window juga disertakan.

Data kosong dibiarkan kosong, bukan diinterpolasi. Seluruh 187 event tetap dimasukkan; jumlah observasi yang tersedia pada tiap menit bervariasi **180–187**. Semua nilai onset tersedia untuk normalisasi. File `composite_mean_std_coverage.csv` menyimpan N per menit.

Terdapat **8 sampel ekstrem pada 5 event** yang memicu flag tinjauan deskriptif speed >20 atau gust >30 m/s. Batas ini hanya flag pemeriksaan, bukan batas fisik atau alasan otomatis membuang data. Semua sampel tetap dipakai pada mean/SD. Gambar utama full-range mempertahankan seluruh rentang. Versi `_zoom` hanya membatasi tampilan sumbu dan mencatat jumlah titik di luar tampilan; mean dan SD tetap dihitung dari semua data.

Gambar wajib: `plots/composite_6panel_187.png` dan versi PDF. Versi fokus: `plots/composite_6panel_187_zoom.png`. Array 187×161 untuk reproduksi tersedia di NPZ.

## 6. Mana yang lebih baik: speed atau gust?

Pada **186 cold pool** dengan kedua jendela speed/gust yang memenuhi cakupan:

| Ukuran perubahan, m/s | Mean speed | Gust |
|---|---:|---:|
| Median perubahan mean event terhadap background | +0,387 | +0,701 |
| Median kenaikan puncak terhadap background | +1,872 | +3,135 |
| Proporsi mean event meningkat | 67,7% | 68,3% |

Gust menunjukkan respons amplitudo lebih besar pada sampel ini. Namun gust memang maksimum dalam interval, dan **187 event dipilih menggunakan kriteria gust**; perbandingan ini mengandung selection bias. Hampir samanya proporsi kenaikan mean juga tidak mendukung klaim sederhana bahwa gust selalu jauh lebih baik. Tanpa label independen dan kasus non-event, kita belum bisa menghitung sensitivity, false-alarm rate, specificity, atau menyatakan akurasi gust lebih tinggi.

**Rekomendasi sementara:** pertahankan gust sebagai indikator utama gust-front agar konsisten dengan metode Hoeller dan amplitudo respons terpilih, sambil memakai mean speed serta komponen u/v untuk analisis respons dan arah. Ini pilihan operasional sementara, bukan hasil pembuktian superioritas.

Untuk menguji mana benar-benar lebih baik, tahap berikutnya adalah mengevaluasi kedua indikator pada kandidat termal yang sama (termasuk temperature-only valid), menggunakan verifikasi independen seperti satelit/radar atau konsistensi lintas stasiun, dan memisahkan tuning dari evaluasi. Periksa pula sampel ekstrem. Mengubah ambang agar gust tampak lebih baik pada 187 event yang sudah dipilih gust tidak memberikan validasi independen.

## 7. Kalimat singkat untuk meeting

1. “187 adalah episode yang lolos konfirmasi suhu–angin; 4.117 adalah wind-only, bukan tambahan cold pool. Dari 148 temperature-only, 57 merupakan QA sehingga 91 dipakai untuk interpretasi meteorologi.”
2. “Overlap cold pool 1 dan 5 menit adalah 122 pasangan pada stasiun sama dan toleransi ±20 menit; saya juga menunjukkan sensitivitas terhadap toleransi.”
3. “Monoton pada kode berarti setiap selisih ≤0, bukan sekadar tren turun. Satu kenaikan 0,01°C dapat menggagalkan anchor 1 menit.”
4. “Pada 91 temperature-only valid, 45 mengalami kenaikan mean speed, 42 penurunan, satu tetap, tiga kurang data. Arah berubah cukup besar pada banyak kasus.”
5. “Background yang lebih kuat berkaitan dengan Δmean speed lebih kecil, tetapi hubungan ini juga dipengaruhi pengurangan baseline secara matematis.”
6. “Gust memberi respons puncak lebih besar, tetapi belum bisa dinyatakan lebih akurat karena event dipilih berdasarkan gust. Arah angin lokal juga belum sama dengan arah gerak cold pool.”

## Reproduksi dan batasan

`analyze.py` berisi analisis utama. `statistics.json` dan `validation_report.json` menyimpan statistik dan pemeriksaan. Excel memuat 148 temperature-only lengkap, subset valid, pasangan matching, metrik angin, dan komposit. Grafik dan tabel tidak mengganti katalog deteksi.

Snapshot katalog 1/5 menit ada pada `catalogs/`; data NetCDF besar tidak dimasukkan ulang ke ZIP. Dengan Python 3.12 dan dependensi `requirements.txt`, analisis utama dapat dijalankan menggunakan:

```powershell
python analyze.py --data-root "C:\Users\bamnet\Downloads\TUTI\THESIS\OlahDataBamnet\Level2Data_New" --catalog-root "catalogs" --output "hasil_analisis_meeting"
```

Perintah tersebut menghasilkan ulang analisis utama, tabel, Excel, dan grafik komposit/monoton/temperature-only. PDF briefing serta grafik tambahan geometri/zoom disertakan sebagai artefak siap pakai dalam paket.

Data diukur oleh lima stasiun dengan cakupan berbeda, sehingga hasil deskriptif ini bukan estimasi frekuensi pada lima stasiun selama 365 hari lengkap. Tidak ada interpolasi nilai hilang, penghapusan outlier secara otomatis, perubahan ambang monoton, atau klaim verifikasi satelit. Penilaian arah propagasi cold pool memerlukan informasi spasial/timing tambahan; satu stasiun hanya memberi arah angin lokal.
