"""Validate the new monotonicity experiment and prepare bilingual meeting PDFs."""
from pathlib import Path
import ast,json,argparse
import numpy as np
import pandas as pd
import netCDF4
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
from openpyxl import load_workbook
from coldpool_flexible import read_time

P=Path(__file__).parent
parser=argparse.ArgumentParser()
parser.add_argument('--data-root',type=Path,default=Path('/workspace/NewJournal'))
ROOT=parser.parse_args().data_root
new=pd.read_csv(P/'detection_1min/all_events.csv',parse_dates=['t0'])
old=pd.read_csv(P/'catalog_old_1min.csv',parse_dates=['t0'])
s=json.loads((P/'statistics.json').read_text())
metrics=pd.read_csv(P/'event_wind_metrics.csv')
thermal=metrics.loc[(metrics.category=='temperature_only')&~metrics.mask_QA]
direction_table=thermal.loc[thermal.direction_reliable].groupby('direction_group').agg(n=('delta_speed','count'),mean_speed_change=('delta_speed','mean'),median_speed_change=('delta_speed','median'),increase_count=('delta_speed',lambda x:int((x>1e-8).sum())))
direction_table.to_csv(P/'direction_groups_speed_changes.csv')
vector=thermal.loc[thermal.change_vector_angle_to_background.notna()].copy()
vector['orientation_group']=pd.cut(vector.change_vector_angle_to_background,[0,45,135,180],include_lowest=True,labels=['parallel_0to45','cross_45to135','opposed_135to180'])
vector_table=vector.groupby('orientation_group',observed=True).agg(n=('delta_speed','count'),mean_speed_change=('delta_speed','mean'),median_speed_change=('delta_speed','median'))
vector_table.to_csv(P/'change_vector_orientation_summary.csv')
tree=ast.parse((P/'analyze.py').read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='match')
ns={'np':np,'pd':pd,'linear_sum_assignment':linear_sum_assignment};exec(compile(ast.Module(body=[fn],type_ignores=[]),'match','exec'),ns)
a=old.loc[~old.temperature_mask_missing_at_t0].copy();b=new.loc[~new.temperature_mask_missing_at_t0].copy()
a['event_id']=[f'old1m-{i}' for i in range(len(a))];b['event_id']=[f'new1m-{i}' for i in range(len(b))]
over=[]
for cat in ['coldpool','temperature_only','wind_only']:
 for tol in [0,5,10,20,30]:
  matches=ns['match'](a.loc[a.cp_type==cat],b.loc[b.cp_type==cat],tol)
  over.append({'category':cat,'tolerance_minutes':tol,'old_valid':int((a.cp_type==cat).sum()),'new_valid':int((b.cp_type==cat).sum()),'matches':len(matches)})
  if tol==20:matches.rename(columns=lambda x:x.replace('_1min','_old_1min').replace('_5min','_new_1min')).to_csv(P/f'old_new_pairs_{cat}.csv',index=False)
pd.DataFrame(over).to_csv(P/'old_new_overlap_sensitivity.csv',index=False)
comparison=pd.DataFrame({'old_1min':old.cp_type.value_counts(),'new_1min':new.cp_type.value_counts(),'original_5min':pd.read_csv(P/'catalogs/5min/all_events.csv').cp_type.value_counts()})
comparison.to_csv(P/'counts_three_versions.csv')

checked=0
for (station,month),g in new.loc[(new.cp_type!='wind_only')&~new.temperature_mask_missing_at_t0].groupby(['station_id','source_month']):
 with netCDF4.Dataset(ROOT/station/f'{station.lower()}_AWS_Level2monthly_{month}.nc') as ds:
  temp=pd.Series(np.ma.filled(ds['air_temperature_1min'][:],np.nan),index=pd.to_datetime(read_time(ds['time_1min']),unit='s'))
 for _,e in g.iterrows():
  raw=temp.reindex(pd.date_range(e.t0-pd.Timedelta('5min'),e.t0+pd.Timedelta('15min'),freq='1min')).to_numpy()
  means=np.array([temp.reindex(pd.date_range(e.t0+pd.Timedelta(minutes=c-2),e.t0+pd.Timedelta(minutes=c+2),freq='1min')).mean(skipna=False) for c in [-5,0,5,10,15]])
  assert np.isfinite(raw).all() and np.isfinite(means).all()
  assert raw[-1]-raw[0]<=-2+1e-9 and raw[5]-raw[0]<=-.5+1e-9
  assert (np.diff(means)<=1e-9).all()
  checked+=1
assert checked==458
z=np.load(P/'composite_310_arrays.npz');minute=z['minutes'];arrs=[z['speed'],z['gust'],z['gust_minus_speed']]
zero=int(np.flatnonzero(minute==0)[0]);coverage=[]
extreme=[];event_index=pd.read_csv(P/'composite_event_index.csv')
for i,e in event_index.iterrows():
 for j in np.flatnonzero((arrs[0][i]>20)|(arrs[1][i]>30)):
  extreme.append({'event_id':e.event_id,'station_id':e.station_id,'t0':e.t0,'minute':int(minute[j]),'speed':float(arrs[0][i,j]),'gust':float(arrs[1][i,j])})
pd.DataFrame(extreme).to_csv(P/'extreme_sample_review_flags.csv',index=False)
for arr in arrs:
 assert arr.shape==(310,161)
 delta=arr-arr[:,zero,None]
 assert np.allclose(delta[np.isfinite(delta[:,zero]),zero],0)
 coverage.append({'min_n':int(np.isfinite(arr).sum(axis=0).min()),'max_n':int(np.isfinite(arr).sum(axis=0).max()),'finite_at_onset':int(np.isfinite(arr[:,zero]).sum())})
fig,axes=plt.subplots(3,2,figsize=(13,11),sharex=True)
for i,(name,arr) in enumerate(zip(['Speed','Gust','Gust − speed'],arrs)):
 for j in range(2):
  x=arr if j==0 else arr-arr[:,zero,None];ax=axes[i,j]
  for line in x:ax.plot(minute,line,lw=.5,alpha=.11,color=['#4477aa','#cc6677','#228833'][i])
  mean=np.nanmean(x,axis=0);sd=np.nanstd(x,axis=0)
  ax.fill_between(minute,mean-sd,mean+sd,color='grey',alpha=.3);ax.plot(minute,mean,color='black',lw=2.3);ax.axvline(0,color='black',ls=':')
  lim=(-1,12) if j==0 else (-8,8)
  if i==2:lim=(-1,7) if j==0 else (-5,5)
  ax.set(ylim=lim,ylabel='m/s',title=f'{name}: '+('absolute' if j==0 else 'X(t)−X(0)')+f'; {int(((x<lim[0])|(x>lim[1])).sum())} points outside view')
for ax in axes[-1]:ax.set_xlabel('Minutes from onset')
fig.suptitle('310 events: viewing zoom only; full-data mean ±1 SD retained');fig.tight_layout();fig.savefig(P/'plots/composite_6panel_310_zoom.png',dpi=160);plt.close(fig)

for language in ['EN','ID']:
 fig,ax=plt.subplots(figsize=(11,5));ax.axis('off')
 text=('ONLY MONOTONICITY CHANGES\n\nInput remains 1-minute temperature and 1-minute gust.\nCenters of the five means: −5, 0, +5, +10, +15 minutes.\nEach mean uses five values: center−2 .. center+2 minutes.\nCompare the four consecutive differences of these five means.\n\nRaw T(t+15)−T(t−5) <= −2 C remains unchanged.\nRaw T(t)−T(t−5) <= −0.5 C remains unchanged.\nWind equations, thresholds, grouping and legacy mode are unchanged.\n\nImportant: mean support extends from t−7 to t+17 minutes.') if language=='EN' else ('HANYA SYARAT MONOTON DIUBAH\n\nInput tetap suhu 1 menit dan gust 1 menit.\nPusat lima mean: −5, 0, +5, +10, +15 menit.\nSetiap mean memakai lima nilai: pusat−2 .. pusat+2 menit.\nEmpat selisih mean berturut-turut harus <=0.\n\nRaw T(t+15)−T(t−5) <= −2 C tidak berubah.\nRaw T(t)−T(t−5) <= −0,5 C tidak berubah.\nPersamaan angin, ambang, grouping dan mode legacy tetap.\n\nPenting: dukungan mean meluas dari t−7 hingga t+17 menit.')
 ax.text(.02,.97,text,va='top',fontfamily='monospace',fontsize=12)
 fig.tight_layout();fig.savefig(P/f'plots/new_monotonic_definition_{language}.png',dpi=160);plt.close(fig)

ov=pd.read_csv(P/'overlap_sensitivity.csv');cpov=ov.loc[(ov['sample']=='QA_excluded')&(ov.category=='coldpool')&(ov.tolerance_minutes==20)].iloc[0]
oldcp=next(r['matches'] for r in over if r['category']=='coldpool' and r['tolerance_minutes']==20)
perf=s['speed_vs_gust'][0]
report={'status':'passed','controlled_tests':3,'only_monotonic_algorithm_changed':True,'independent_actual_event_checks':checked,'counts':new.cp_type.value_counts().to_dict(),'QA_count':int(new.temperature_mask_missing_at_t0.sum()),'composite_shape':[310,161],'composite_coverage':coverage,'extreme_samples_for_review':len(extreme),'extreme_events_for_review':len(set(e['event_id'] for e in extreme)),'matches_vs_original_5min_at_20min':int(cpov.matched),'matches_vs_old_1min_at_20min':oldcp,'definition':'five centered 5-minute means at t-5,t,t+5,t+10,t+15; means use center-2..center+2; four nonpositive increments; raw drop tests unchanged'}
(P/'validation_report.json').write_text(json.dumps(report,indent=2))
notes=[('Method / Metode',report['definition']),('Counts / Jumlah','310 coldpool; 4000 wind-only; 205 temperature-only (148 valid +57 legacy QA).'),('Normalization','X_i(t)−X_i(0), per event before computing mean and population SD; units m/s.'),('Baseline','Mean −40..−20 min; event mean 0..20 min; at least 17/21 finite values.'),('Interpretation','Direction is from vector-mean u/v; direction rotation is not cold-pool propagation.'),('Limitations','Gust-selected sample; no independent accuracy/satellite validation. Baseline-change correlations have mathematical coupling.')]
with pd.ExcelWriter(P/'Analisis_Meeting_BAMNET.xlsx',engine='openpyxl',mode='a',if_sheet_exists='replace') as w:
 pd.DataFrame(notes,columns=['item','definition']).to_excel(w,sheet_name='METODE',index=False)
 pd.DataFrame(over).to_excel(w,sheet_name='Old_new_overlap',index=False)
 comparison.reset_index(names='category').to_excel(w,sheet_name='Three_versions_counts',index=False)
 direction_table.reset_index().to_excel(w,sheet_name='Direction_rotation_groups',index=False)
 vector_table.reset_index().to_excel(w,sheet_name='Change_vector_groups',index=False)

def textpage(pdf,title,paras):
 fig=plt.figure(figsize=(11.7,8.3));fig.text(.05,.94,title,fontsize=20,weight='bold');y=.85
 for text in paras:
  fig.text(.06,y,text,fontsize=12,va='top',linespacing=1.45);y-=.07+.04*len(text.split('\n'))
 pdf.savefig(fig);plt.close(fig)
def imagepage(pdf,title,name,caption):
 fig=plt.figure(figsize=(11.7,8.3));fig.text(.04,.95,title,fontsize=17,weight='bold')
 ax=fig.add_axes([.04,.14,.92,.75]);ax.imshow(Image.open(P/'plots'/name));ax.axis('off')
 fig.text(.04,.075,caption,fontsize=10,va='top');pdf.savefig(fig);plt.close(fig)
for lang in ['EN','ID']:
 en=lang=='EN'
 with PdfPages(P/f'Meeting_New_Monotonic_{lang}.pdf') as pdf:
  textpage(pdf,'New monotonicity experiment: BAMNET 2025' if en else 'Eksperimen monoton baru: BAMNET 2025',[
   '1-minute input retained. Only the monotonicity test uses five 5-minute means.' if en else 'Input tetap 1 menit. Hanya uji monoton memakai lima mean 5 menit.',
   'New: 310 cold pools; 4,000 wind-only; 205 temperature-only.' if en else 'Baru: 310 cold pool; 4.000 angin saja; 205 suhu saja.',
   'Previous 1-minute result: 187 / 4,117 / 148. Original 5-minute: 197 / 2,294 / 178.' if en else '1 menit lama: 187 / 4.117 / 148. Data 5 menit asli: 197 / 2.294 / 178.',
   '205 thermal-only entries = 148 valid episodes +57 legacy missing-mask QA entries.' if en else '205 suhu saja = 148 episode valid +57 QA dari mask kosong mode legacy.',
   'All 57 available files were processed. No satellite verification is performed.' if en else 'Seluruh 57 file tersedia diproses. Belum ada verifikasi satelit.'
  ])
  imagepage(pdf,'Exact new monotonicity definition' if en else 'Definisi tepat syarat monoton baru',f'new_monotonic_definition_{lang}.png','No minute-scale monotonicity test remains; raw drop amplitudes remain unchanged.' if en else 'Uji monoton per menit diganti; amplitudo penurunan suhu raw tetap sama.')
  imagepage(pdf,'A real temperature example' if en else 'Contoh suhu nyata','monotonic_actual_example.png','Purple: five centered means from minute data; compare four increments.' if en else 'Ungu: lima mean terpusat dari data 1 menit; bandingkan empat selisihnya.')
  imagepage(pdf,'Overlap with original 5-minute detections' if en else 'Overlap dengan deteksi 5 menit asli','overlap_sensitivity.png',f'±20 min: {int(cpov.matched)} cold-pool pairs; {cpov.percent_1min_matched:.1f}% of new, {cpov.percent_5min_matched:.1f}% of 5-minute events.\nSame station, one-to-one matching; QA excluded.' if en else f'±20 menit: {int(cpov.matched)} pasangan cold pool; {cpov.percent_1min_matched:.1f}% dari versi baru, {cpov.percent_5min_matched:.1f}% dari 5 menit.\nStasiun sama, pasangan satu-ke-satu, QA dikeluarkan.')
  textpage(pdf,'Changes relative to the previous 1-minute algorithm' if en else 'Perubahan terhadap algoritma 1 menit sebelumnya',[
   f'Cold-pool pairs within ±20 minutes: {oldcp} (old 187 versus new 310).' if en else f'Pasangan cold pool dalam ±20 menit: {oldcp} (lama 187, baru 310).',
   'More candidates are not automatically more accurate detections.\nThe modified constraint is an explicit sensitivity experiment.' if en else 'Lebih banyak kandidat tidak otomatis lebih akurat.\nAturan baru adalah eksperimen sensitivitas yang dinyatakan eksplisit.',
   'Five centered means span 25 source minutes (t−7..t+17).\nIt is impossible to form five disjoint 5-minute blocks from only 21 samples.' if en else 'Lima mean terpusat memakai 25 menit sumber (t−7..t+17).\nLima blok 5 menit tidak dapat dibuat dari hanya 21 sampel.',
   'Raw drop tests still use t−5..t+15 and t−5..t.\nWind equations and event grouping are unchanged.' if en else 'Uji penurunan raw tetap t−5..t+15 dan t−5..t.\nPersamaan angin dan pengelompokan event tidak diubah.'
  ])
  imagepage(pdf,'Thermal-only wind response and baseline' if en else 'Respons angin dan baseline suhu saja','temperature_only_wind_analysis.png',f'Valid thermal-only: {s["speed_sign_valid"]}.\nBackground −40..−20; event 0..20 min; minimum 17/21 samples.' if en else '148 suhu saja valid: 81 naik, 61 turun, 1 tetap, 5 data kurang.\nBackground −40..−20, event 0..20 menit; minimal 17/21 sampel.')
  textpage(pdf,'Direction and baseline results' if en else 'Hasil arah dan baseline',[
   f'Median speed change: {s["delta_speed_distribution"]["50%"]:.3f} m/s; n=143.' if en else f'Median perubahan speed: {s["delta_speed_distribution"]["50%"]:.3f} m/s; n=143.',
   f'Reliable vector-mean directions: {s["direction_reliable_count"]}; median absolute rotation {s["absolute_direction_change_reliable"]["50%"]:.1f}°.' if en else f'Arah rata-rata vektor andal: {s["direction_reliable_count"]}; median rotasi absolut {s["absolute_direction_change_reliable"]["50%"]:.1f}°.',
   f'Baseline vs mean-speed change: rho={s["baseline_delta_speed_spearman"]["rho"]:.3f}.\nBaseline vs peak change: rho={s["baseline_peak_change_spearman"]["rho"]:.3f}.' if en else f'Baseline vs perubahan mean-speed: rho={s["baseline_delta_speed_spearman"]["rho"]:.3f}.\nBaseline vs perubahan puncak: rho={s["baseline_peak_change_spearman"]["rho"]:.3f}.',
   'These are descriptive associations. Subtracting baseline creates\nmathematical coupling; local wind direction is not propagation direction.' if en else 'Ini asosiasi deskriptif. Pengurangan baseline menimbulkan\nmathematical coupling; arah angin lokal bukan arah propagasi.'
  ])
  imagepage(pdf,'Vector geometry: interpretation limits' if en else 'Geometri vektor: batas interpretasi','vector_geometry.png','Schematic only, not inferred cold-pool motion.' if en else 'Hanya skema, bukan estimasi gerak cold pool.')
  imagepage(pdf,'Six-panel composite: 310 selected events' if en else 'Komposit enam panel: 310 event terpilih','composite_6panel_310.png','Thin event traces; thick black mean; shading ±1 population SD. Left absolute; right X(t)−X(0).\nWindow −40..+120 min. Missing values are not interpolated.' if en else 'Garis event tipis; mean hitam tebal; shading ±1 SD populasi. Kiri absolut, kanan X(t)−X(0).\nJendela −40..+120 menit. Nilai kosong tidak diinterpolasi.')
  imagepage(pdf,'Viewing zoom, not data exclusion' if en else 'Zoom tampilan, bukan pembuangan data','composite_6panel_310_zoom.png','All data are retained in the mean and standard deviation.' if en else 'Seluruh data tetap masuk perhitungan mean dan standar deviasi.')
  textpage(pdf,'Speed versus gust and next decisions' if en else 'Speed versus gust dan keputusan berikutnya',[
   f'Paired eligible events: {perf["n_paired"]}. Median peak increase:\nspeed {perf["median_peak_speed_change"]:.3f} m/s; gust {perf["median_peak_gust_change"]:.3f} m/s.' if en else f'Event berpasangan layak: {perf["n_paired"]}. Median kenaikan puncak:\nspeed {perf["median_peak_speed_change"]:.3f} m/s; gust {perf["median_peak_gust_change"]:.3f} m/s.',
   'Gust-selected samples cannot establish that gust has higher accuracy.\nUse independent labels and the same thermal candidates for a fair comparison.' if en else 'Sampel yang dipilih memakai gust tidak membuktikan gust lebih akurat.\nGunakan label independen dan kandidat termal sama untuk pembandingan adil.',
   'Normalization: subtract each event onset value, then calculate ensemble mean/SD.\nUnits remain m/s; not division, percentage change, or z-score.' if en else 'Normalisasi: kurangi nilai onset tiap event, lalu hitung mean/SD ensemble.\nSatuan tetap m/s; bukan pembagian, persentase, atau z-score.',
   'Next meeting: confirm the centered averaging convention, assess QA artifacts,\nvalidate sensitivity against independent evidence, and investigate direction.' if en else 'Meeting berikut: sepakati konvensi mean terpusat, periksa artefak QA,\nvalidasi sensitivitas dengan bukti independen, lalu lanjut analisis arah.'
  ])
print(json.dumps(report,indent=2))
