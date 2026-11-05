# New monotonicity experiment: BAMNET 2025

This folder contains a separate experiment using **1-minute data**, changing **only the monotonicity constraint**. Earlier results are preserved.

## Files

- [ZIP containing all new results](BAMNET_New_Monotonic_All_Results.zip).
- [Full detection workbook](detection_1min/BAMNET_semua_event.xlsx).
- [English meeting PDF](Meeting_New_Monotonic_EN.pdf).
- [Indonesian meeting PDF](Meeting_New_Monotonic_ID.pdf).
- [Analysis workbook](Analisis_Meeting_BAMNET.xlsx).
- [Six-panel composite, 310 events](plots/composite_6panel_310.png).
- `coldpool_flexible.py`, `analyze.py`, `finalize_reports.py`: detection, analysis, and bilingual PDF generation.

On GitHub, open a file and choose **Download raw file**.

## Exact averaging convention

Every candidate onset remains on the **1-minute grid**. Construct five centered 5-minute means at **t−5, t, t+5, t+10, t+15**. Each mean uses five source temperatures from **center−2 through center+2 minutes**, requiring complete data. The four consecutive differences of these five means must be ≤0. Equal means are allowed.

The source blocks are t−7..t−3, t−2..t+2, t+3..t+7, t+8..t+12, and t+13..t+17. Thus mean support extends to **t−7..t+17**. Four increments require five means; five disjoint blocks of five values cannot be constructed from only the original 21 samples. This support extension is explicitly documented, rather than silently claiming all averaging occurs within t−5..t+15.

The input is `air_temperature_1min`, not the pre-existing 5-minute temperature variable. The averaging convention is centered, not a calendar-aligned resampling that restricts event times to multiples of five minutes.

## Everything else stays unchanged

- **Raw** T(t+15)−T(t−5) ≤−2°C.
- **Raw** T(t)−T(t−5) ≤−0.5°C.
- Wind input remains maximum raw wind within each **1-minute** bin.
- Eq. (2), conditional Eq. (3), thresholds, windows, event grouping, overlap-based categories, month-by-month processing, and legacy boundary behavior remain unchanged.

The means are used **only for monotonicity**, not for the two raw-temperature drop tests. Timing shifts remain resolution-relative.

## Detection results

| Station | Cold-pool candidates | Wind-only | Temperature-only, legacy |
|---|---:|---:|---:|
| BAMNET01 | 92 | 763 | 46 |
| BAMNET02 | 66 | 479 | 23 |
| BAMNET03 | 79 | 877 | 43 |
| BAMNET04 | 42 | 934 | 48 |
| BAMNET05 | 31 | 947 | 45 |
| **Total** | **310** | **4,000** | **205** |

There are **4,515 catalog entries**. The 205 thermal-only entries include **148 valid thermal episodes and 57 legacy missing-mask artifacts**. All entries remain visible and flagged; physical interpretation excludes the QA artifacts.

The previous strict 1-minute algorithm produced 187/4,117/148, while the original 5-minute algorithm produced 197/2,294/178. More detections do not automatically indicate better accuracy.

## Matching and wind analysis

One-to-one matching requires the same station. At ±20 minutes, the new catalog has **189 cold-pool pairs with the original 5-minute catalog** (61.0% of new candidates and 95.9% of the 5-minute candidates). There are **187 pairs with the previous 1-minute cold-pool catalog**; all earlier candidates have a nearby new candidate, without claiming identical timestamps.

The 148 valid thermal-only episodes show **81 increases, 61 decreases, 1 unchanged, and 5 insufficient-data cases** in mean speed. The median change among 143 evaluable cases is **+0.110 m/s**. Background is −40..−20 minutes and event conditions 0..20 minutes, with at least 17 of 21 observations per mean.

Vector-mean directions are reliable for 131 episodes when both pre/post resultant magnitudes are at least 0.5 m/s. Median absolute rotation is **96.2°**: 45 rotations below 45°, 50 from 45–135°, and 36 above 135°. These are changes in measured wind direction, not cold-pool propagation directions. Event-level component changes and change-vector orientations are also provided.

Background speed correlates negatively with mean-speed change (Spearman rho **−0.489**), but the peak-change relationship is weaker (**−0.117**). Baseline subtraction creates mathematical coupling and possible regression to the mean; these associations do not establish causation.

## Composite and normalization

The composite now contains **310 event traces**, aligned over −40..+120 minutes on a 1-minute grid. Rows show mean speed, gust, and gust minus speed; columns show absolute values and **X_i(t)−X_i(0)**. Subtract onset separately for each event before computing the ensemble mean and population standard deviation. Units remain m/s; this is neither a percentage nor a z-score.

Missing samples are not interpolated; per-minute sample counts are exported. A viewing zoom retains all observations in the statistics. Extreme-sample flags (speed>20 or gust>30 m/s) are review markers, not automatic rejection criteria.

Among 308 candidates with sufficient paired windows, median peak increases are **+2.265 m/s for speed** and **+3.888 m/s for gust**. Gust-selected samples cannot prove superior detection accuracy. Retain gust provisionally for method consistency and use speed/u/v for response and direction analysis; independent validation is still required.

## Validation and reproducibility

Three controlled tests passed, including AST checks confirming unchanged wind and grouping algorithms. All **458 valid thermal events** were independently checked against the NetCDF source: both raw-temperature drops still pass and all four mean increments are nonpositive. Composite arrays have shape 310×161. See `validation_report.json` for coverage and matching checks.

Use Python 3.12 and `requirements.txt`. Run the detector with `--resolution-minutes 1 --boundary-policy legacy`. The full Windows commands are in the Indonesian README. `catalogs/1min/all_events.csv` is the input snapshot for the analysis; refresh it if rerunning detection with altered settings. NetCDF source data are not duplicated in this output folder.

Limitations remain: no satellite verification; unequal observational coverage; sample-based rolling/grouping may cross gaps; inherited duration metrics are not validated physical cold-pool durations. No source data or earlier results were overwritten.
