from pathlib import Path
import sys,json,argparse
import numpy as np
import pandas as pd
import netCDF4
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).parent))
from coldpool_flexible import Timing,load_month,read_time

parser=argparse.ArgumentParser(description='Analisis perbandingan BAMNET dan komposit angin')
parser.add_argument('--data-root',type=Path,default=Path('/workspace/NewJournal'))
parser.add_argument('--catalog-root',type=Path,default=Path(__file__).parent/'catalogs')
parser.add_argument('--output',type=Path,default=Path(__file__).parent)
args=parser.parse_args()
ROOT=args.data_root
OUT=args.output
OUT.mkdir(exist_ok=True)
PLOTS=OUT/'plots';PLOTS.mkdir(exist_ok=True)
cats=['coldpool','temperature_only','wind_only']
data={r:pd.read_csv(args.catalog_root/f'{r}min'/'all_events.csv',parse_dates=['t0']) for r in [1,5]}
for r,d in data.items():
 d['event_id']=[f'{r}m-{i:04d}' for i in range(len(d))]
valid={r:d.loc[~d.temperature_mask_missing_at_t0].copy() for r,d in data.items()}

def match(a,b,tolerance):
 records=[]
 for station in sorted(set(a.station_id)|set(b.station_id)):
  x=a.loc[a.station_id==station].reset_index(drop=True);y=b.loc[b.station_id==station].reset_index(drop=True)
  if x.empty or y.empty:continue
  difference=(y.t0.to_numpy()[None,:]-x.t0.to_numpy()[:,None])/np.timedelta64(1,'m')
  costs=np.where(abs(difference)<=tolerance,abs(difference),1e9)
  costs=np.concatenate([costs,np.full((len(x),len(x)),1e5)],axis=1)
  rows,cols=linear_sum_assignment(costs)
  for i,j in zip(rows,cols):
   if j>=len(y) or costs[i,j]>=1e5:continue
   records.append({'station_id':station,'event_id_1min':x.iloc[i].event_id,'t0_1min':x.iloc[i].t0,'category_1min':x.iloc[i].cp_type,
    'event_id_5min':y.iloc[j].event_id,'t0_5min':y.iloc[j].t0,'category_5min':y.iloc[j].cp_type,'offset_5minus1_minutes':float(difference[i,j])})
 return pd.DataFrame(records)

sensitivity=[];primary=[]
for mode,ds in [('legacy',data),('QA_excluded',valid)]:
 for category in cats:
  a=ds[1].loc[ds[1].cp_type==category];b=ds[5].loc[ds[5].cp_type==category]
  for tolerance in [0,5,10,20,30]:
   matches=match(a,b,tolerance);m=len(matches)
   sensitivity.append({'sample':mode,'category':category,'tolerance_minutes':tolerance,'n_1min':len(a),'n_5min':len(b),'matched':m,
    'percent_1min_matched':100*m/len(a) if len(a) else np.nan,'percent_5min_matched':100*m/len(b) if len(b) else np.nan,'jaccard':m/(len(a)+len(b)-m) if len(a)+len(b) else np.nan})
   if mode=='QA_excluded' and tolerance==20 and len(matches):primary.append(matches)
pd.DataFrame(sensitivity).to_csv(OUT/'overlap_sensitivity.csv',index=False)
pd.concat(primary).to_csv(OUT/'overlap_same_category_20min.csv',index=False)
cross=match(valid[1],valid[5],20);cross.to_csv(OUT/'overlap_cross_category_20min.csv',index=False)
confusion=pd.crosstab(cross.category_1min,cross.category_5min).reindex(index=cats,columns=cats,fill_value=0)
confusion.to_csv(OUT/'category_transitions_matched.csv')
print('Matching complete',flush=True)

# One-minute station series, retaining padding to cross month boundaries.
station_data={}
for station in sorted(valid[1].station_id.unique()):
 frames=[]
 for path in sorted((ROOT/station).glob('*Level2monthly_2025*.nc')):
  month=int(path.stem[-2:]);start=pd.Timestamp(2025,month,1);end=start+pd.offsets.MonthBegin(1)
  frame=load_month(path,start-pd.Timedelta('1D'),end+pd.Timedelta('1D'),Timing(1))
  with netCDF4.Dataset(path) as ds:
   time=pd.to_datetime(read_time(ds['time_1min']),unit='s')
   extra=pd.DataFrame({name:np.ma.filled(ds[source][:].astype(float),np.nan) for name,source in [('speed','wind_speed_1min'),('direction','wind_from_direction_1min'),('u','u_wind_1min'),('v','v_wind_1min')]},index=time)
  frames.append(frame.join(extra,how='left'))
 series=pd.concat(frames).sort_index()
 duplicates=series.loc[series.index.duplicated(keep=False)]
 if ((duplicates.groupby(level=0).nunique(dropna=False)>1).any(axis=1)).any():raise ValueError(f'Conflicting padded records {station}')
 station_data[station]=series.loc[~series.index.duplicated()]
 print('Loaded',station,flush=True)

relative=np.arange(-40,121)
cp=valid[1].loc[valid[1].cp_type=='coldpool'].reset_index(drop=True)
assert len(cp)==310
speed=[];gust=[]
for _,e in cp.iterrows():
 window=station_data[e.station_id].reindex(e.t0+pd.to_timedelta(relative,unit='min'))
 speed.append(window.speed.to_numpy());gust.append(window.gust.to_numpy())
speed=np.array(speed);gust=np.array(gust);gap=gust-speed
np.savez_compressed(OUT/'composite_310_arrays.npz',minutes=relative,speed=speed,gust=gust,gust_minus_speed=gap,event_ids=cp.event_id.to_numpy(dtype=str))
cp.to_csv(OUT/'composite_event_index.csv',index=False)
fig,axes=plt.subplots(3,2,figsize=(14,12),sharex=True)
summary=[]
for row,(name,arr,color) in enumerate([('Mean wind speed',speed,'#4477aa'),('Wind gust',gust,'#cc6677'),('Gust minus speed',gap,'#228833')]):
 for col in range(2):
  values=arr if col==0 else arr-arr[:,relative.tolist().index(0),None]
  mean=np.nanmean(values,axis=0);std=np.nanstd(values,axis=0);n=np.isfinite(values).sum(axis=0)
  ax=axes[row,col]
  for line in values:ax.plot(relative,line,color=color,alpha=.12,lw=.55)
  ax.fill_between(relative,mean-std,mean+std,color='grey',alpha=.3,label='±1 SD across events')
  ax.plot(relative,mean,color='black',lw=2.5,label='Mean across events')
  ax.axvline(0,color='black',ls='--',lw=.8);ax.axhline(0,color='grey',lw=.5)
  ax.set_title(f'{name}: '+('absolute' if col==0 else 'X(t) − X(0)')+f'\n310 event traces; available n={n.min()}–{n.max()}')
  ax.set_ylabel('m/s');ax.grid(alpha=.15)
  if row==0 and col==0:ax.legend(fontsize=8)
  summary.extend({'variable':name,'mode':'absolute' if col==0 else 'onset_centered','minute':int(t),'mean':float(m),'std':float(s),'n':int(k)} for t,m,s,k in zip(relative,mean,std,n))
for ax in axes[-1]:ax.set_xlabel('Minutes from cold-pool onset (UTC anchor)')
fig.suptitle('310 station-derived cold pools, 1-minute data; missing samples are not interpolated',fontsize=14)
fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(PLOTS/'composite_6panel_310.png',dpi=180);fig.savefig(PLOTS/'composite_6panel_310.pdf');plt.close(fig)
pd.DataFrame(summary).to_csv(OUT/'composite_mean_std_coverage.csv',index=False)

def mean_checked(values,minimum=17):
 a=np.asarray(values,dtype=float)
 return float(np.nanmean(a)) if np.isfinite(a).sum()>=minimum else np.nan

def direction(u,v):return float(np.degrees(np.arctan2(-u,-v))%360) if np.isfinite(u+v) else np.nan

metrics=[]
for _,e in data[1].loc[data[1].cp_type.isin(['temperature_only','coldpool'])].iterrows():
 s=station_data[e.station_id]
 before=s.reindex(e.t0+pd.to_timedelta(np.arange(-40,-19),unit='min'))
 after=s.reindex(e.t0+pd.to_timedelta(np.arange(0,21),unit='min'))
 r={'event_id':e.event_id,'station_id':e.station_id,'t0':e.t0,'category':e.cp_type,'mask_QA':e.temperature_mask_missing_at_t0}
 for key in ['speed','gust','u','v']:
  r[key+'_before']=mean_checked(before[key]);r[key+'_after']=mean_checked(after[key]);r['delta_'+key]=r[key+'_after']-r[key+'_before']
 for key in ['speed','gust']:
  a=after[key].to_numpy();r['peak_'+key+'_change']=float(np.nanmax(a))-r[key+'_before'] if np.isfinite(a).sum()>=17 else np.nan
 r['speed_change_sign']='missing' if not np.isfinite(r['delta_speed']) else ('increase' if r['delta_speed']>1e-8 else 'decrease' if r['delta_speed']<-1e-8 else 'unchanged')
 r['speed_change_01ms']='missing' if not np.isfinite(r['delta_speed']) else ('increase' if r['delta_speed']>.1 else 'decrease' if r['delta_speed']<-.1 else 'within_0.1ms')
 r['direction_before']=direction(r['u_before'],r['v_before']);r['direction_after']=direction(r['u_after'],r['v_after'])
 r['direction_change_signed']=((r['direction_after']-r['direction_before']+180)%360)-180
 r['direction_change_absolute']=abs(r['direction_change_signed'])
 r['baseline_vector_speed']=np.hypot(r['u_before'],r['v_before']);r['event_vector_speed']=np.hypot(r['u_after'],r['v_after'])
 r['direction_reliable']=bool(r['baseline_vector_speed']>=.5 and r['event_vector_speed']>=.5)
 angle=r['direction_change_absolute']
 r['direction_group']='missing' if not np.isfinite(angle) else 'aligned_0to45' if angle<45 else 'cross_45to135' if angle<135 else 'opposed_135to180'
 # Change-vector orientation is a diagnostic, not a measured propagation direction.
 b=np.array([r['u_before'],r['v_before']]);d=np.array([r['delta_u'],r['delta_v']])
 r['change_vector_magnitude']=float(np.linalg.norm(d))
 r['delta_parallel_to_background']=float(np.dot(b,d)/np.linalg.norm(b)) if np.linalg.norm(b)>=.5 else np.nan
 r['change_vector_angle_to_background']=float(np.degrees(np.arccos(np.clip(np.dot(b,d)/(np.linalg.norm(b)*np.linalg.norm(d)),-1,1)))) if np.linalg.norm(b)>=.5 and np.linalg.norm(d)>=.5 else np.nan
 metrics.append(r)
metrics=pd.DataFrame(metrics);metrics.to_csv(OUT/'event_wind_metrics.csv',index=False)
temp=metrics.loc[(metrics.category=='temperature_only') & ~metrics.mask_QA].copy()
temp_all=metrics.loc[metrics.category=='temperature_only']
temp.to_csv(OUT/'temperature_only_148_valid_metrics.csv',index=False)

def corr(frame,x,y):
 d=frame[[x,y]].dropna()
 if len(d)<3:return {'n':len(d),'rho':None,'p':None}
 stat=spearmanr(d[x],d[y]);return {'n':len(d),'rho':float(stat.statistic),'p':float(stat.pvalue)}

stats={'temperature_all_legacy':len(temp_all),'temperature_QA':int(temp_all.mask_QA.sum()),'temperature_valid':len(temp),
 'speed_sign_valid':temp.speed_change_sign.value_counts().to_dict(),'speed_sign_all205':temp_all.speed_change_sign.value_counts().to_dict(),
 'speed_change_deadband_01ms_valid':temp.speed_change_01ms.value_counts().to_dict(),
 'delta_speed_distribution':temp.delta_speed.describe(percentiles=[.1,.25,.5,.75,.9]).to_dict(),
 'direction_reliable_count':int(temp.direction_reliable.sum()),
 'direction_groups_reliable':temp.loc[temp.direction_reliable].direction_group.value_counts().to_dict(),
 'absolute_direction_change_reliable':temp.loc[temp.direction_reliable].direction_change_absolute.describe().to_dict(),
 'baseline_delta_speed_spearman':corr(temp,'speed_before','delta_speed'),
 'baseline_peak_change_spearman':corr(temp,'speed_before','peak_speed_change')}
station_corr=[]
for station,g in temp.groupby('station_id'):station_corr.append({'station_id':station,**corr(g,'speed_before','delta_speed')})
pd.DataFrame(station_corr).to_csv(OUT/'baseline_correlations_by_station.csv',index=False)
q=temp.dropna(subset=['speed_before','delta_speed']).copy();q['baseline_bin']=pd.qcut(q.speed_before,4,duplicates='drop')
baseline=q.groupby('baseline_bin',observed=True).agg(n=('delta_speed','count'),baseline_mean=('speed_before','mean'),delta_mean=('delta_speed','mean'),delta_median=('delta_speed','median'),increase_fraction=('delta_speed',lambda x:float((x>1e-8).mean())))
baseline.to_csv(OUT/'baseline_quartile_summary.csv')
stats['baseline_quartiles']=baseline.reset_index().astype({'baseline_bin':str}).to_dict('records')

fig,axes=plt.subplots(2,2,figsize=(12,9))
axes[0,0].hist(temp.delta_speed.dropna(),bins=18,color='#4477aa');axes[0,0].axvline(0,color='black');axes[0,0].set(xlabel='Post-event mean − background mean speed (m/s)',ylabel='Valid thermal-only events',title='148 real thermal-only episodes; available windows only')
reliable=temp.loc[temp.direction_reliable]
axes[0,1].hist(reliable.direction_change_signed.dropna(),bins=np.arange(-180,181,30),color='#228833');axes[0,1].set(xlabel='Signed wind-direction rotation (degrees)',ylabel='Events',title='Vector mean directions; pre/post resultant ≥0.5 m/s')
for station,g in temp.groupby('station_id'):axes[1,0].scatter(g.speed_before,g.delta_speed,s=23,alpha=.7,label=station)
axes[1,0].axhline(0,color='black',lw=.7);axes[1,0].set(xlabel='Background mean speed (m/s)',ylabel='Change of mean speed (m/s)',title='Baseline effect (association, not causation)');axes[1,0].legend(fontsize=7)
axes[1,1].scatter(reliable.direction_change_absolute,reliable.delta_speed,c=reliable.speed_before,cmap='viridis',s=30)
axes[1,1].axhline(0,color='black',lw=.7);axes[1,1].set(xlabel='Absolute direction rotation (degrees)',ylabel='Change of mean speed (m/s)',title='Observed rotation ≠ cold-pool propagation direction')
fig.tight_layout();fig.savefig(PLOTS/'temperature_only_wind_analysis.png',dpi=180);fig.savefig(PLOTS/'temperature_only_wind_analysis.pdf');plt.close(fig)

# Paired performance diagnostics; detection was selected using gust, so these
# describe the selected sample and cannot prove superiority or specificity.
perf=[]
for category in ['coldpool','temperature_only']:
 g=metrics.loc[(metrics.category==category)&~metrics.mask_QA]
 paired=g.dropna(subset=['peak_speed_change','peak_gust_change','delta_speed','delta_gust'])
 perf.append({'category':category,'n_catalog':len(g),'n_paired':len(paired),'mean_speed_increase_fraction':float((paired.delta_speed>0).mean()),'mean_gust_increase_fraction':float((paired.delta_gust>0).mean()),'median_mean_speed_change':float(paired.delta_speed.median()),'median_mean_gust_change':float(paired.delta_gust.median()),'median_peak_speed_change':float(paired.peak_speed_change.median()),'median_peak_gust_change':float(paired.peak_gust_change.median())})
pd.DataFrame(perf).to_csv(OUT/'speed_vs_gust_diagnostics.csv',index=False);stats['speed_vs_gust']=perf

# Actual event where coarse monotonicity hides minute-scale fluctuations.
examples=[]
for _,e in valid[5].loc[valid[5].cp_type=='coldpool'].iterrows():
 series=station_data[e.station_id].temperature.reindex(pd.date_range(e.t0-pd.Timedelta('5min'),e.t0+pd.Timedelta('15min'),freq='1min'))
 if series.notna().all() and (series.diff().dropna()>0).any():examples.append((e,series))
if examples:
 e,series=examples[0];path=ROOT/e.station_id/f'{e.station_id.lower()}_AWS_Level2monthly_{e.source_month}.nc'
 with netCDF4.Dataset(path) as ds:
  five=pd.Series(np.ma.filled(ds['air_temperature_5min'][:],np.nan),index=pd.to_datetime(read_time(ds['time_5min']),unit='s')).reindex(pd.date_range(e.t0-pd.Timedelta('5min'),e.t0+pd.Timedelta('15min'),freq='5min'))
 fig,ax=plt.subplots(figsize=(10,5));ax.plot((series.index-e.t0).total_seconds()/60,series,'o-',label='Raw 1-min temperature; raw drop tests unchanged')
 ax.plot((five.index-e.t0).total_seconds()/60,five,'s--',lw=2,label='Original 5-min aggregated temperature')
 centers=np.array([-5,0,5,10,15])
 means=np.array([station_data[e.station_id].temperature.reindex(pd.date_range(e.t0+pd.Timedelta(minutes=int(c)-2),e.t0+pd.Timedelta(minutes=int(c)+2),freq='1min')).mean(skipna=False) for c in centers])
 ax.plot(centers,means,'D-',lw=2,color='purple',label='New rule: five centered means from 1-min data')
 ax.axvline(0,color='black',ls=':');ax.set(xlabel='Minutes from 5-min onset',ylabel='Temperature (°C)',title=f'New monotonicity: four increments between five centered 5-minute means\n{e.station_id}, {e.t0} UTC')
 ax.legend();fig.tight_layout();fig.savefig(PLOTS/'monotonic_actual_example.png',dpi=180);plt.close(fig)
 pd.DataFrame({'minute':(series.index-e.t0).total_seconds()/60,'temperature_1min':series.to_numpy(),'increment_1min':series.diff().to_numpy()}).to_csv(OUT/'monotonic_example_values.csv',index=False)
 stats['monotonic_example']={'station':e.station_id,'t0':str(e.t0),'positive_1min_increments':int((series.diff()>0).sum())}

sens=pd.DataFrame(sensitivity);fig,ax=plt.subplots(figsize=(8,5))
for c in cats:
 d=sens.loc[(sens['sample']=='QA_excluded')&(sens.category==c)];ax.plot(d.tolerance_minutes,d.percent_1min_matched,'o-',label=c)
ax.set(xlabel='Allowed onset separation (minutes)',ylabel='% of 1-min events matched one-to-one',title='Same station and category; QA artifacts excluded');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(PLOTS/'overlap_sensitivity.png',dpi=180);plt.close(fig)

with pd.ExcelWriter(OUT/'Analisis_Meeting_BAMNET.xlsx',engine='openpyxl') as w:
 for name,df in [('Overlap_sensitivity',sens),('Pairs_same_category',pd.concat(primary)),('Pairs_cross_category',cross),('Category_transitions',confusion.reset_index()),('Temp_all_205',temp_all),('Temp_valid_148',temp),('All_wind_metrics',metrics),('Baseline_quartiles',baseline.reset_index().astype({'baseline_bin':str})),('Speed_vs_gust',pd.DataFrame(perf)),('Composite_summary',pd.DataFrame(summary))]:df.to_excel(w,sheet_name=name,index=False)
 for sheet in w.book.worksheets:sheet.freeze_panes='A2';sheet.auto_filter.ref=sheet.dimensions
(OUT/'statistics.json').write_text(json.dumps(stats,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
print(json.dumps(stats,indent=2,default=str),flush=True)
