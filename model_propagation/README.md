# model_propagation — model konseptual propagasi cold pool

Model dikembangkan bertahap sesuai tugas dosen: **lingkaran tanpa angin → lingkaran dengan angin latar → elips dengan angin latar**, lalu inversi untuk memperkirakan asal dari waktu kedatangan di stasiun.

**Input yang digunakan hanya katalog asli 5 menit berisi 197 kandidat cold pool**, dengan dukungan Eq. (2) pada 188 episode dan Eq. (3) pada 9 episode. Katalog 1 menit maupun monoton baru tidak dipakai. Deteksi dan data sumber sebelumnya tidak diubah.

## Berkas utama

- [ZIP semua hasil dan kode](BAMNET_Model_Propagation_All_Results.zip).
- [PDF meeting Indonesia](Propagation_Meeting_ID.pdf), [English PDF](Propagation_Meeting_EN.pdf), masing-masing 12 halaman.
- [Excel seluruh hasil](Propagation_Results.xlsx).
- [Animasi tiga model konseptual](animations/three_models.gif).
- [Animasi contoh nyata bersyarat](animations/real_conditional_example.gif).
- [Gambar tiga model](plots/conceptual_three_models.png).
- [Peta galat kandidat asal](plots/real_origin_loss_profile.png).
- [Kode model](propagation.py), [studi lengkap](run_study.py), [validasi dan PDF](finish_study.py), [pengujian](test_propagation.py).
- [English explanation](README_EN.md).

Pada GitHub, buka berkas lalu klik **Download raw file**. GIF dapat diputar di browser. PDF menjelaskan cara membedakan demonstrasi simulasi dan hasil bersyarat data BAMNET.

## 1. Model maju: bagaimana front mencapai stasiun?

### Lingkaran, tanpa angin

Dengan asal horizontal r₀=(x₀,y₀), waktu mulai t₀, dan kecepatan radial c:

`|rᵢ − r₀| = c (tᵢ − t₀)`.

Front membesar sama cepat ke semua arah. Empat parameter belum diketahui: x₀,y₀,t₀,c. Kecepatan c bukan wind speed atau gust yang diukur di stasiun; c adalah kecepatan penyebaran batas cold pool relatif terhadap pusat model.

### Lingkaran dengan angin latar

`pusat(t) = r₀ + U (t−t₀)` dan `radius(t) = c(t−t₀)`.

Stasiun tercapai ketika jaraknya dari pusat bergerak sama dengan radius. Kode menyelesaikan persamaan kuadrat untuk waktu kedatangan. Dalam studi ini U **ditentukan**, bukan dicari bersama asal. Menyamakan U dengan angin permukaan sebelum event adalah **asumsi proxy**, bukan fakta bahwa cold pool benar-benar beradveksi persis seperti angin tersebut.

### Elips dengan angin latar

Pusat tetap bergerak seperti di atas. Semi-sumbu berkembang sebagai:

- `a(t) = c√r (t−t₀)`.
- `b(t) = c/√r (t−t₀)`.

r≥1 adalah rasio sumbu mayor/minor dan theta orientasinya. c adalah mean geometrik laju pertumbuhan kedua sumbu. Setelah koordinat diputar ke sumbu elips, stasiun berada pada front ketika `(d₁/a)²+(d₂/b)²=1`. Saat r=1, model ini **tepat kembali ke lingkaran**, terlepas dari theta; hal itu diuji pada unit test.

Ini model **kinematik horizontal**, bukan simulasi dinamika fluida atmosfer. Tidak ada terrain following, perubahan c akibat waktu, shear vertikal, interaksi cold pool, atau mekanisme pembentukan ulang front.

## 2. Inversi dan kemampuan menentukan parameter

Kode mencari parameter yang meminimalkan selisih waktu kedatangan model dan pengamatan, memakai bounded nonlinear least squares dan banyak titik awal. Multi-start mengurangi risiko solusi lokal, tetapi tidak membuktikan solusi unik atau optimum global.

| Kasus | Parameter bebas | Jumlah minimum berdasarkan hitungan parameter |
|---|---:|---:|
| Lingkaran tanpa angin | x,y,t₀,c = 4 | 4 kedatangan |
| Lingkaran, U ditentukan | x,y,t₀,c = 4 | 4 kedatangan |
| Elips bebas, U ditentukan | x,y,t₀,c,r,theta = 6 | 6 kedatangan |
| Skenario nyata dengan c/U/r/theta tetap | x,y,t₀ = 3 | 3 kedatangan |

Hitungan minimum **bukan jaminan keunikan**. Empat data dengan empat parameter dapat menghasilkan residual nol tanpa menyisakan derajat bebas untuk validasi. Jika U juga harus ditentukan, tambahkan dua parameter. Parameter elips tidak teridentifikasi bebas dari hanya lima waktu kedatangan satu lintasan; jaringan aktual bahkan hanya memiliki kelompok maksimal empat stasiun pada asosiasi yang diuji.

## 3. Uji simulasi: asal diketahui

Tiga model diuji pada sembilan stasiun simulasi yang tersebar dua dimensi:

- Asal sejati (−12,−8) km, waktu mulai −60 menit, c=6 m/s.
- Kasus angin memakai U=(2,1) m/s.
- Elips sejati r=1,6 dan theta=30°; pada jaringan sembilan stasiun, r dan theta ikut dicari.

Pada data tanpa noise, ketiga inversi memulihkan asal hingga presisi numerik (galat <0,001 km dan RMSE <10⁻⁵ menit). Ini verifikasi implementasi pada asumsi yang benar, bukan akurasi asal pada data nyata.

Masing-masing model kemudian diuji dengan **35 perturbasi waktu uniform independen ±2,5 menit**:

| Model | Median galat posisi asal | Galat maksimum |
|---|---:|---:|
| Lingkaran tanpa angin | 0,421 km | 0,823 km |
| Lingkaran + angin | 0,542 km | 1,064 km |
| Elips + angin | 0,537 km | 1,258 km |

Ini sensitivitas bersyarat pada jaringan simulasi yang baik, model benar, dan U diketahui. ±2,5 menit hanya pendekatan kuantisasi optimistis untuk data 5 menit; ketidakpastian onset nyata bisa lebih besar. Nilai ini **bukan confidence interval asal BAMNET**.

## 4. Geometri BAMNET dan mengapa penting

Koordinat diambil dari metadata NetCDF; koordinat tiap stasiun konsisten pada semua bulan yang tersedia. Proyeksi lokal equirectangular memakai pusat sekitar 100,32516°E, −0,53416° dan satuan km. Untuk wilayah ini pendekatan planar dipakai, bukan peta topografi.

Lima titik memanjang pada arah sekitar 58,6° terhadap timur, dengan rasio sebaran minor/mayor **0,108**. Konfigurasi memanjang mengurangi informasi arah tegak lurus jaringan. Elevasi berbeda dari **27 hingga 945 m ASL**, sedangkan model mengabaikan topografi. Karena itu hasil tidak boleh langsung dianggap jalur aliran permukaan sebenarnya.

Uji simulasi pada **koordinat lima BAMNET** menghasilkan banyak elips berbeda dengan waktu kedatangan yang sama hampir eksak; rank Jacobian maksimal 5 untuk 6 parameter. Gambar `ellipse_nonuniqueness.png` memperlihatkan asal yang berbeda meskipun fit waktunya sama baik. Lingkaran secara teori lebih terbatasi, tetapi data nyata tetap memerlukan asosiasi dan asumsi yang tepat.

## 5. Menentukan kandidat event bersama

197 deteksi bukan otomatis 197 cold pool fisik berbeda: satu cold pool bisa melewati beberapa stasiun. Untuk membentuk hipotesis, deteksi diurutkan menurut waktu dan dikelompokkan secara **greedy dengan rentang maksimum dari deteksi pertama**, bukan chaining tanpa batas. Rentang utama 180 menit; sensitivitas 60/90/120 menit ikut disimpan.

| Rentang | Kelompok temporal | Kelompok ≥3 stasiun | Kelompok ≥4 stasiun | Kelompok 5 stasiun |
|---|---:|---:|---:|---:|
| 60 menit | 155 | 9 | 0 | 0 |
| 90 menit | 150 | 10 | 1 | 0 |
| 120 menit | 146 | 11 | 1 | 0 |
| 180 menit | 143 | 11 | 2 | 0 |

Pada rentang utama: 103 kelompok hanya satu stasiun, 40 kelompok melibatkan setidaknya dua. Ada satu kelompok yang berisi pengulangan stasiun; kelompok ambigu seperti itu tidak dipakai untuk fit sumber. Semua 197 deteksi tetap tersimpan tepat sekali di tabel anggota kelompok.

**Kedekatan waktu tidak membuktikan satu kejadian fisik.** Kelompok dapat mencampur cold pool berbeda, mencerminkan pemicu konveksi yang berbeda, atau gagal memasangkan satu cold pool yang onsetnya jauh. Asosiasi perlu diverifikasi dengan meteorologi/spasial sebelum sumber diterima. Grafik suhu kasus contoh disertakan untuk inspeksi, bukan dinyatakan sebagai pembuktian asosiasi.

## 6. Skenario asal data nyata

Untuk 11 kelompok dengan ≥3 stasiun unik, fit utama menetapkan **c=5 m/s** dan mencari x,y,t₀. Skenario sensitivitas memakai **c=3 dan 8 m/s**. Nilai ini bukan pengukuran kecepatan propagasi; nilainya dipilih untuk eksperimen konseptual.

Proxy U adalah median antarstasiun dari rata-rata komponen u/v selama **t−40..t−20 menit**, dihitung per stasiun terhadap onsetnya sendiri. Mean membutuhkan setidaknya 4 dari 5 nilai 5 menit yang valid. U kelompok dipakai konstan; ini tidak menjamin merepresentasikan lingkungan prakejadian di semua tempat/waktu. Bila hanya sebagian stasiun memiliki data cukup, jumlahnya ditulis pada tabel.

Skenario elips menetapkan r=1,5 dan theta sejajar proxy U; ketika U<0,5 m/s, theta memakai sumbu utama jaringan. Theta+90° ikut diuji. **Rasio dan orientasi elips data nyata tidak diestimasi bebas**. Angka r=1,5 adalah skenario, bukan nilai yang dibuktikan dari data.

Batas inversi: x,y dalam ±100 km dari pusat proyeksi; waktu mulai antara 360 menit sebelum hingga sesaat sebelum deteksi pertama. c yang dicari, hanya pada kelompok empat stasiun, dibatasi 1–30 m/s dan harus cukup cepat agar adveksi subkritis pada metrik bentuk. Skenario fixed-c yang superkritis dilewati dengan status eksplisit. Batas ini adalah prior/aturan model; titik yang menyentuhnya perlu diperhatikan, bukan dianggap sumber terpercaya. Sensitivitas domain ±50/100/150 km ditampilkan untuk kasus contoh.

### Contoh 28 Juni 2025, kelompok G059

Deteksi: BAMNET01 05:10, BAMNET03 05:40, BAMNET04 06:05, BAMNET05 06:20 **UTC**. Proxy U=(0,646,1,772) m/s.

Dengan c tetap 5 m/s:

| Skenario | Estimasi longitude | Estimasi latitude | Waktu mulai UTC, dibulatkan | RMSE kedatangan |
|---|---:|---:|---|---:|
| Lingkaran tanpa angin | 100,179°E | −0,626° | 04:29 | 2,56 menit |
| Lingkaran + proxy angin | 100,102°E | −0,692° | 04:22 | 5,35 menit |
| Elips r=1,5 + proxy angin | 100,037°E | −0,756° | 03:57 | 5,35 menit |

Ini **tiga hipotesis bersyarat**, bukan tiga titik asal yang terverifikasi. Perubahan asumsi menggeser lokasi dan waktu. Kesalahan lebih kecil pada lingkaran tanpa angin dalam kasus ini tidak membuktikan angin sebenarnya nol; proxy angin/model/asosiasi dapat tidak sesuai.

Jika c juga dicari pada empat stasiun, RMSE dapat mendekati nol pada lebih dari satu model. Karena ada empat parameter dan hanya empat data, **RMSE nol bukan validasi**. Model lingkaran tanpa angin, misalnya, memberi c≈4,09 m/s pada kasus ini, tetapi itu tetap tidak merupakan kecepatan yang diverifikasi independen.

Kelompok empat stasiun lain, G043 pada 21 April, memberikan residual lebih besar dan beberapa solusi menyentuh batas waktu mulai −360 menit atau batas spasial. Itu peringatan bahwa model/asosiasi/prior mungkin tidak memadai. Tidak ada pemilihan model berdasarkan residual saja.

Peta `real_origin_loss_profile.png` menampilkan galat sebagai fungsi asal untuk circle c=5 m/s, dengan t₀ dioptimalkan pada setiap grid. Kontur putih adalah **RMSE grid terbaik +2,5 menit**, bukan 95% CI. Seluruh array grid tersedia di NPZ. Batasan kecepatan, bentuk, U, asosiasi, dan prior tidak ikut dirangkum sebagai satu confidence interval.

## 7. Animasi dan cara menjelaskannya

- `three_models.gif`: sumber simulasi sama, stasiun sama, tiga mekanisme berbeda. Titik hijau berarti **front model** sudah mencapai stasiun.
- `real_conditional_example.gif`: contoh G059 pada tiga asumsi; titik hijau berarti **onset telah teramati dalam katalog**, bukan tanda front model pasti sudah tiba. Perbedaan front dan waktu observasi dapat dilihat langsung.

Koordinat asal dalam CSV adalah hasil optimasi dengan digit komputasi; banyak digit bukan bukti ketelitian geografis/temporal. Waktu input, semua grouping, dan tabel kedatangan memakai UTC (WIB=UTC+7). Onset katalog berasal dari kriteria suhu, sehingga tidak otomatis menandai fase front yang persis sama di semua stasiun; respons lokal dan offset suhu–gust dapat membiasakan inversi. t₀ adalah waktu awal efektif hasil ekstrapolasi model, bukan waktu pembentukan konveksi yang terukur langsung.

## 8. Yang perlu dibahas pada meeting

1. Persamaan dan animasi sudah menggambarkan tiga tahap yang diminta.
2. Inversi berhasil pada simulasi dengan asumsi diketahui; ini bukan bukti akurasi data nyata.
3. Untuk BAMNET, verifikasi event yang sama antarstasiun merupakan prasyarat penting sebelum menetapkan asal.
4. Geometri memanjang, kurangnya kelompok ≥5 stasiun, asumsi c/U, dan topografi membatasi identifikasi sumber.
5. Elips bebas tidak dapat diidentifikasi dari data sekarang; gunakan bentuk yang ditentukan atau tambah observasi independen.
6. Berikutnya: verifikasi satelit/radar/struktur suhu–angin, pembatasan c dan U dari informasi independen, stasiun lintas sumbu jaringan, dan evaluasi topografi. Pusat model tidak harus sama persis dengan lokasi hujan atau inisiasi konveksi.

## Reproduksi

Python 3.12; dependensi pada `requirements.txt`. File input 197 dan koordinat sudah disertakan, tetapi NetCDF besar tetap berada pada data pengguna.

```powershell
python -m pip install -r requirements.txt
python -m unittest -v test_propagation
python run_study.py --data-root "C:\Users\bamnet\Downloads\TUTI\THESIS\OlahDataBamnet\Level2Data_New" --catalog "inputs\coldpool_197_5min.csv" --output "."
python finish_study.py --data-root "C:\Users\bamnet\Downloads\TUTI\THESIS\OlahDataBamnet\Level2Data_New"
```

Seed acak ditetapkan. Tujuh unit test lulus; input dan Eq. (3) diverifikasi; semua deteksi masuk satu kelompok; pemulihan simulasi, rank elips, prediksi kasus nyata, jumlah frame GIF, dan workbook diperiksa. Lihat `validation_report.json`, `manifest.json`, dan tabel multi-start untuk hasil lengkap. Semua sumber nyata tetap bersyarat, bukan lokasi asal yang sudah dikonfirmasi.
