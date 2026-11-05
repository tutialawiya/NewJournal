"""Deteksi cold pool dengan durasi fisik, bukan jumlah sampel hardcode.

Contoh: python coldpool_flexible.py --data-root DATA --output HASIL
        --resolution-minutes 5 --boundary-policy legacy

Mode legacy mempertahankan hasil kode pengguna (termasuk perilaku mask NaN).
Mode exclude membuang mask kosong dan memperbaiki kasus daftar kriteria kosong.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr


@dataclass(frozen=True)
class Timing:
    resolution_minutes: int = 5
    temperature_before_minutes: int = 5
    temperature_after_minutes: int = 15
    initial_drop_minutes: int = 5
    wind_before_minutes: int = 20
    wind_after_minutes: int = 40
    wind_reference_minutes: int = 120
    anomaly_reference_minutes: int = 1440
    fallback_lag_minutes: int = 80
    event_gap_minutes: int = 20

    def __post_init__(self):
        if self.resolution_minutes <= 0:
            raise ValueError('Resolusi harus positif.')
        for name, value in asdict(self).items():
            if name != 'resolution_minutes':
                self.steps(value)  # Jangan membulatkan durasi yang tidak habis dibagi.
        if self.initial_drop_minutes > self.temperature_before_minutes + self.temperature_after_minutes:
            raise ValueError('Jendela penurunan awal melebihi jendela suhu.')

    def steps(self, minutes: int) -> int:
        """Ubah menit ke jumlah langkah; tanda negatif tetap dipertahankan."""
        duration = pd.Timedelta(minutes=minutes)
        resolution = pd.Timedelta(minutes=self.resolution_minutes)
        if duration % resolution:
            raise ValueError(f'{minutes} menit tidak habis dibagi resolusi {self.resolution_minutes} menit.')
        return int(duration // resolution)

    @property
    def frequency(self):
        return f'{self.resolution_minutes}min'


def temperature_criterion(temperature: xr.DataArray, timing: Timing):
    before = timing.steps(timing.temperature_before_minutes)
    after = timing.steps(timing.temperature_after_minutes)
    initial = timing.steps(timing.initial_drop_minutes)
    # +1 memasukkan kedua ujung interval, bukan tambahan satu menit.
    window = temperature.rolling(time=before + after + 1).construct('window')
    total_drop = window.isel(window=-1) - window.isel(window=0) <= -2.0
    monotonic = (window.diff('window') <= 0).all('window')
    first_drop = window.isel(window=initial) - window.isel(window=0) <= -0.5
    # 15/5 = -3 sampel; 15/1 = -15 sampel.
    return (total_drop & monotonic & first_drop).shift(time=-after)


def wind_criterion(wind_gust: xr.DataArray, timing: Timing, n_w=3, n_miss=2):
    count_2h = timing.steps(timing.wind_reference_minutes) + 1
    # Tidak ditambah 1: mempertahankan 288 sampel pada kode asli 5 menit.
    count_24h = timing.steps(timing.anomaly_reference_minutes)
    window_count = timing.steps(timing.wind_before_minutes + timing.wind_after_minutes) + 1
    future_shift = -timing.steps(timing.wind_after_minutes)
    reference_shift = timing.steps(timing.fallback_lag_minutes)
    mean_2h = wind_gust.rolling(time=count_2h, center=True).mean(skipna=True)
    std_2h = wind_gust.rolling(time=count_2h, center=True).std(skipna=True, ddof=0)
    anomaly = wind_gust - mean_2h
    mean_24h = anomaly.rolling(time=count_24h, center=True).mean(skipna=True)
    std_24h = anomaly.rolling(time=count_24h, center=True).std(skipna=True, ddof=0)
    maximum_anomaly = (anomaly.rolling(time=window_count).construct('window')
                       .max('window', skipna=False).shift(time=future_shift))
    maximum_speed = (wind_gust.rolling(time=window_count).construct('window')
                     .max('window', skipna=True).shift(time=future_shift))
    eq2 = maximum_anomaly >= mean_24h + n_w * std_24h
    # Eq. 3 tetap khusus saat anomali hilang, mengikuti kode pengguna.
    threshold_3 = mean_2h.shift(time=reference_shift) + n_miss * std_2h.shift(time=reference_shift)
    eq3 = maximum_anomaly.isnull() & (maximum_speed >= threshold_3)
    return eq2 | eq3, eq2, eq3


def event_groups(mask: xr.DataArray, timing: Timing, boundary_policy='legacy'):
    values = mask.values
    if boundary_policy == 'exclude':
        values = mask.fillna(False).astype(bool).values
    # np.where(NaN) dipertahankan HANYA dalam mode legacy agar cocok dengan 178.
    positions = np.flatnonzero(values)
    if not len(positions):
        return np.empty((0, 2), dtype=int)
    groups = np.split(positions, np.flatnonzero(np.diff(positions) > timing.steps(timing.event_gap_minutes)) + 1)
    return np.asarray([(g[0], g[-1]) for g in groups], dtype=int)


def classify_frame(frame: pd.DataFrame, timing: Timing, boundary_policy='legacy'):
    if boundary_policy not in {'legacy', 'exclude'}:
        raise ValueError('boundary_policy harus legacy atau exclude.')
    if not frame.index.is_unique or not frame.index.is_monotonic_increasing:
        raise ValueError('Timestamp harus unik dan berurutan.')
    if (frame.index.as_unit('ns').asi8 % pd.Timedelta(timing.frequency).value).any():
        raise ValueError('Timestamp tidak sejajar dengan resolusi yang dipilih.')
    def array(name):
        return xr.DataArray(frame[name].to_numpy(), dims=['time'], coords={'time': frame.index.to_numpy()})
    temperature, gust, pressure = array('temperature'), array('gust'), array('pressure')
    thermal = temperature_criterion(temperature, timing)
    wind, eq2, eq3 = wind_criterion(gust, timing)
    ta = event_groups(thermal, timing, boundary_policy)
    wa = event_groups(wind, timing, boundary_policy)
    times = frame.index.to_numpy()
    if boundary_policy == 'exclude' and (not len(ta) or not len(wa)):
        cp = np.array([], dtype=times.dtype)
        temp_only = times[ta[:, 0]] if len(ta) else np.array([], dtype=times.dtype)
        wind_only = times[wa[:, 0]] if len(wa) else np.array([], dtype=times.dtype)
    else:
        cp, temp_only, wind_only = get_coldpool_event(ta, wa, times, times)
    return (cp, temp_only, wind_only), (temperature, gust, pressure), (thermal, wind, eq2, eq3)


def read_time(variable):
    unit, origin = variable.units.split(' since ', 1)
    factor = {'seconds': 1, 'minutes': 60}[unit]
    if getattr(variable, 'calendar', 'proleptic_gregorian') != 'proleptic_gregorian':
        raise ValueError('Kalender NetCDF tidak didukung.')
    return pd.Timestamp(origin).as_unit('ns').value // 1_000_000_000 + np.asarray(variable[:], dtype=np.int64) * factor


def load_month(path: Path, start, end, timing: Timing):
    """Suhu dari seri resolusi terpilih; gust maksimum raw wind per interval."""
    resolution_seconds = int(pd.Timedelta(timing.frequency).total_seconds())
    suffix = f'{timing.resolution_minutes}min'
    with netCDF4.Dataset(path) as ds:
        for name in [f'time_{suffix}', f'air_temperature_{suffix}', f'air_pressure_{suffix}']:
            if name not in ds.variables:
                raise ValueError(f'{path.name}: variabel {name} tidak tersedia.')
        if ds[f'air_temperature_{suffix}'].units != 'degree_Celsius' or ds['wind_speed'].units != 'm s-1':
            raise ValueError(f'{path.name}: unit suhu/angin tidak sesuai.')
        seconds = read_time(ds['time_1sec'])
        if (np.diff(seconds) <= 0).any():
            raise ValueError(f'{path.name}: timestamp raw tidak unik/berurutan.')
        lo, hi = np.searchsorted(seconds, [start.as_unit('ns').value // 1_000_000_000,
                                           end.as_unit('ns').value // 1_000_000_000])
        seconds = seconds[lo:hi]
        if not len(seconds):
            raise ValueError(f'{path.name}: tidak ada data pada bulan yang diminta.')
        speed = np.ma.filled(ds['wind_speed'][lo:hi], np.nan).astype(float)
        bins, first = np.unique(seconds // resolution_seconds, return_index=True)
        maximum = np.maximum.reduceat(np.where(np.isfinite(speed), speed, -np.inf), first)
        maximum[maximum == -np.inf] = np.nan
        gust = pd.Series(maximum, index=pd.to_datetime(bins * resolution_seconds, unit='s'))
        gust = gust.reindex(pd.date_range(gust.index[0], gust.index[-1], freq=timing.frequency))
        index = pd.to_datetime(read_time(ds[f'time_{suffix}']), unit='s')
        frame = pd.DataFrame({
            'temperature': np.ma.filled(ds[f'air_temperature_{suffix}'][:], np.nan),
            'pressure': np.ma.filled(ds[f'air_pressure_{suffix}'][:], np.nan),
        }, index=index)
        common = frame.index.intersection(gust.index).sort_values()
        frame = frame.loc[common].copy()
        frame['gust'] = gust.reindex(common)
    return frame


def get_coldpool_event(temperature_criterion_arr, wind_criterion_arr, temp_time_vector, wind_time_vector):
    """
    Mendeteksi cold pool event berdasarkan overlap antara temperature dan wind criterion.
    """
    if len(temperature_criterion_arr) == 0 or len(wind_criterion_arr) == 0:
        return (np.array([]), np.array([]), np.array([]))
    start_T = temperature_criterion_arr[:, 0][:, None]
    end_T = temperature_criterion_arr[:, 1][:, None]
    start_WSPD = wind_criterion_arr[:, 0][None, :]
    end_WSPD = wind_criterion_arr[:, 1][None, :]
    overlap = (start_T <= end_WSPD) & (start_WSPD <= end_T)
    mask_T_combined = overlap.any(axis=1)
    coldpool_time = temp_time_vector[temperature_criterion_arr[mask_T_combined, 0]]
    temp_only_time = temp_time_vector[temperature_criterion_arr[~mask_T_combined, 0]]
    mask_W_combined = overlap.any(axis=0)
    wind_only_time = wind_time_vector[wind_criterion_arr[~mask_W_combined, 0]]
    return (coldpool_time, temp_only_time, wind_only_time)

def compute_cp_anomalies(t0, temperature, wind_speed_gust, pressure, ref_window_before=40, ref_window_after=20, analysis_window_before=40, analysis_window_after=120):
    t0 = pd.Timestamp(t0)
    t_ref_start = t0 - pd.Timedelta(minutes=ref_window_before)
    t_ref_end = t0 - pd.Timedelta(minutes=ref_window_after)
    t_ana_start = t0 - pd.Timedelta(minutes=analysis_window_before)
    t_ana_end = t0 + pd.Timedelta(minutes=analysis_window_after)
    try:
        T_ref = float(temperature.sel(time=slice(t_ref_start, t_ref_end)).mean().values)
        T_ref_temp = float(temperature.sel(time=slice(t0 - pd.Timedelta(minutes=10), t0)).max().values)
        T_ana = temperature.sel(time=slice(t_ana_start, t_ana_end))
        U_ana = wind_speed_gust.sel(time=slice(t_ana_start, t_ana_end))
        P_ana = pressure.sel(time=slice(t_ana_start, t_ana_end))
        U_ref = float(wind_speed_gust.sel(time=slice(t_ref_start, t_ref_end)).mean().values)
        P_ref = float(pressure.sel(time=slice(t_ref_start, t_ref_end)).mean().values)
        delta_T = float((T_ana - T_ref_temp).min().values)
        delta_U = float((U_ana - U_ref).max().values)
        delta_P = float((P_ana - P_ref).max().values)
        T_after = temperature.sel(time=slice(t0, t_ana_end))
        recovered = T_after.where(np.abs(T_after - T_ref_temp) <= 0.5).dropna('time')
        if len(recovered) > 0:
            t_end = pd.Timestamp(recovered.time.values[0])
            duration_min = int((t_end - t0).total_seconds() / 60)
        else:
            duration_min = np.nan
    except Exception:
        delta_T, delta_U, delta_P, duration_min = (np.nan, np.nan, np.nan, np.nan)
    return (delta_T, delta_U, delta_P, duration_min)

def build_cp_catalog(coldpool_time, temp_only_time, wind_only_time, temperature, wind_speed_gust, pressure, station_id):
    records = []
    for t0 in coldpool_time:
        t0 = pd.Timestamp(t0)
        dT, dU, dP, dur = compute_cp_anomalies(t0, temperature, wind_speed_gust, pressure)
        records.append({'station_id': station_id, 't0': t0, 'cp_type': 'coldpool', 'delta_T': round(dT, 3) if not np.isnan(dT) else np.nan, 'delta_U': round(dU, 3) if not np.isnan(dU) else np.nan, 'delta_P': round(dP, 3) if not np.isnan(dP) else np.nan, 'duration_min': dur})
    for t0 in temp_only_time:
        t0 = pd.Timestamp(t0)
        dT, dU, dP, dur = compute_cp_anomalies(t0, temperature, wind_speed_gust, pressure)
        records.append({'station_id': station_id, 't0': t0, 'cp_type': 'temperature_only', 'delta_T': round(dT, 3) if not np.isnan(dT) else np.nan, 'delta_U': round(dU, 3) if not np.isnan(dU) else np.nan, 'delta_P': round(dP, 3) if not np.isnan(dP) else np.nan, 'duration_min': dur})
    for t0 in wind_only_time:
        t0 = pd.Timestamp(t0)
        dT, dU, dP, dur = compute_cp_anomalies(t0, temperature, wind_speed_gust, pressure)
        records.append({'station_id': station_id, 't0': t0, 'cp_type': 'wind_only', 'delta_T': round(dT, 3) if not np.isnan(dT) else np.nan, 'delta_U': round(dU, 3) if not np.isnan(dU) else np.nan, 'delta_P': round(dP, 3) if not np.isnan(dP) else np.nan, 'duration_min': dur})
    catalog = pd.DataFrame(records).sort_values('t0').reset_index(drop=True)
    return catalog


def run(root: Path, output: Path, timing: Timing, boundary_policy='legacy', year=2025, stations=None):
    output.mkdir(parents=True, exist_ok=True)
    (output / 'monthly').mkdir(exist_ok=True)
    catalogs, month_rows, qa_rows = [], [], []
    for station in stations or [f'BAMNET{i:02d}' for i in range(1, 6)]:
        paths = sorted((root / station).glob(f'*Level2monthly_{year}*.nc'))
        if not paths:
            raise FileNotFoundError(f'Tidak ada file {year} di {root / station}')
        for path in paths:
            month = int(path.stem[-2:])
            start = pd.Timestamp(year, month, 1)
            frame = load_month(path, start, start + pd.offsets.MonthBegin(1), timing)
            detections, arrays, masks = classify_frame(frame, timing, boundary_policy)
            cp, temp_only, wind_only = detections
            if any(len(a) for a in detections):
                catalog = build_cp_catalog(cp, temp_only, wind_only, *arrays, station)
            else:
                catalog = pd.DataFrame(columns=['station_id', 't0', 'cp_type', 'delta_T', 'delta_U', 'delta_P', 'duration_min'])
            catalog['source_month'] = f'{year}{month:02d}'
            thermal = masks[0]
            catalog['temperature_mask_missing_at_t0'] = [bool(pd.isna(thermal.sel(time=t).item())) for t in catalog.t0]
            catalog['t0_WIB'] = pd.to_datetime(catalog.t0) + pd.Timedelta(hours=7)
            catalog.to_csv(output / 'monthly' / f'{station}_{year}{month:02d}.csv', index=False)
            catalogs.append(catalog)
            counts = catalog.cp_type.value_counts().to_dict()
            row = {'station_id': station, 'source_month': f'{year}{month:02d}', **{c: int(counts.get(c, 0)) for c in ['coldpool', 'wind_only', 'temperature_only']}}
            month_rows.append(row)
            qa_rows.append({'station_id': station, 'source_month': f'{year}{month:02d}',
                            'missing_mask_samples': int(pd.isna(thermal.values).sum()),
                            'irregular_time_steps': int((frame.index.to_series().diff().dropna() != pd.Timedelta(timing.frequency)).sum())})
            print(row, flush=True)
    combined = pd.concat(catalogs, ignore_index=True).sort_values(['t0', 'station_id']).reset_index(drop=True)
    combined.to_csv(output / 'all_events.csv', index=False)
    for category in ['coldpool', 'wind_only', 'temperature_only']:
        combined.loc[combined.cp_type == category].to_csv(output / f'{category}.csv', index=False)
    summary = pd.crosstab(combined.station_id, combined.cp_type).reindex(columns=['coldpool', 'wind_only', 'temperature_only'], fill_value=0)
    summary.to_csv(output / 'summary.csv')
    pd.DataFrame(month_rows).to_csv(output / 'monthly_summary.csv', index=False)
    pd.DataFrame(qa_rows).to_csv(output / 'quality_summary.csv', index=False)
    manifest = {'timing': asdict(timing), 'boundary_policy': boundary_policy, 'year': year,
                'input_file_count': len(month_rows), 'counts': combined.cp_type.value_counts().to_dict(),
                'wind_input': f'maximum wind_speed in {timing.resolution_minutes}-minute bins',
                'time_zone': 'UTC; t0_WIB is UTC+7', 'satellite_verified': False}
    (output / 'run_manifest.json').write_text(json.dumps(manifest, indent=2))
    export_excel(combined, summary, output / 'BAMNET_semua_event.xlsx', manifest)
    print(json.dumps(manifest), flush=True)
    return combined


def export_excel(events, summary, path, manifest):
    from openpyxl.worksheet.table import Table, TableStyleInfo
    notes = [
        ('Resolusi', f"{manifest['timing']['resolution_minutes']} menit"),
        ('Kebijakan batas', manifest['boundary_policy']),
        ('Jumlah', json.dumps(manifest['counts'])),
        ('Waktu', 't0 adalah UTC. t0_WIB adalah UTC+7. Tahun dan bulan mengikuti UTC.'),
        ('QA', 'temperature_mask_missing_at_t0=True menandai event dari mask suhu kosong. Mode legacy mempertahankannya untuk kesesuaian dengan hasil asli 178.'),
        ('Angin', manifest['wind_input']),
        ('Metode', 'Eq. 2 OR Eq. 3; Eq. 3 hanya ketika anomali angin hilang. Hoeller 2024, DOI 10.1029/2023GL107308.'),
        ('Verifikasi', 'Belum diverifikasi satelit. Kelengkapan pengamatan berbeda antarstasiun.'),
        ('Durasi', 'duration_min mengikuti fungsi katalog asli; definisi pemulihan mulai t0 dapat menghasilkan 0 dan belum tervalidasi.'),
    ]
    with pd.ExcelWriter(path, engine='openpyxl', datetime_format='yyyy-mm-dd hh:mm') as writer:
        events.to_excel(writer, sheet_name='SEMUA_EVENT', index=False)
        for name, category in [('COLD_POOL', 'coldpool'), ('WIND_ONLY', 'wind_only'), ('TEMPERATURE_ONLY', 'temperature_only')]:
            events.loc[events.cp_type == category].to_excel(writer, sheet_name=name, index=False)
        summary.reset_index().to_excel(writer, sheet_name='RINGKASAN', index=False)
        events.loc[events.temperature_mask_missing_at_t0].to_excel(writer, sheet_name='QA_MASK_KOSONG', index=False)
        pd.DataFrame(notes, columns=['item', 'keterangan']).to_excel(writer, sheet_name='CATATAN', index=False)
        for i, sheet in enumerate(writer.book.worksheets):
            sheet.freeze_panes = 'A2'
            sheet.sheet_view.showGridLines = False
            if sheet.max_row > 1:
                table = Table(displayName=f'Tabel{i}', ref=sheet.dimensions)
                table.tableStyleInfo = TableStyleInfo(name='TableStyleMedium2', showRowStripes=True)
                sheet.add_table(table)
            for cells in sheet.columns:
                sheet.column_dimensions[cells[0].column_letter].width = min(35, max(18, len(str(cells[0].value)) + 3))
            if sheet.title == 'CATATAN':
                sheet.column_dimensions['B'].width = 110


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resolution-minutes', type=int, choices=[1, 5], default=5)
    parser.add_argument('--boundary-policy', choices=['legacy', 'exclude'], default='legacy')
    parser.add_argument('--year', type=int, default=2025)
    parser.add_argument('--stations', nargs='+')
    args = parser.parse_args()
    run(args.data_root, args.output, Timing(args.resolution_minutes), args.boundary_policy, args.year, args.stations)
