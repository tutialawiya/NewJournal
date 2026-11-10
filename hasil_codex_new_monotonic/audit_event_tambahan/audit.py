"""Audit the 121 new-1-minute cold-pool entries unmatched to 5-minute CPs."""
from pathlib import Path
import sys,json,argparse
import numpy as np
import pandas as pd
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data-root',type=Path,default=Path('/workspace/NewJournal'))
parser.add_argument('--results-root',type=Path,default=Path(__file__).resolve().parents[1])
args=parser.parse_args()
sys.path.insert(0,str(args.results_root.resolve()))
import coldpool_flexible as new
import reference_old_detector as old

P=Path(__file__).parent;SOURCE=args.results_root
a=pd.read_csv(SOURCE/'catalogs/5min/all_events.csv',parse_dates=['t0'])
b=pd.read_csv(SOURCE/'detection_1min/all_events.csv',parse_dates=['t0'])
same=pd.read_csv(SOURCE/'overlap_same_category_20min.csv',parse_dates=['t0_1min','t0_5min'])
same=same.loc[same.category_1min=='coldpool'];keys=set(zip(same.station_id,same.t0_1min))
extra=b.loc[(b.cp_type=='coldpool')&np.array([ (r.station_id,r.t0) not in keys for r in b.itertuples()],dtype=bool)].copy()
assert len(extra)==121
cross=pd.read_csv(SOURCE/'overlap_cross_category_20min.csv',parse_dates=['t0_1min','t0_5min'])
cross=cross.loc[cross.category_1min=='coldpool']
extra=extra.merge(cross[['station_id','t0_1min','category_5min','t0_5min','offset_5minus1_minutes']],left_on=['station_id','t0'],right_on=['station_id','t0_1min'],how='left')
extra['previous_label']=extra.category_5min.fillna('no_matched_5min_event_within20')
extra=extra.sort_values(['t0','station_id']).reset_index(drop=True)

def wind_numbers(frame,minutes):
 u=frame.gust
 mean=u.rolling(120//minutes+1,center=True).mean()
 std=u.rolling(120//minutes+1,center=True).std(ddof=0)
 anomaly=u-mean
 mean24=anomaly.rolling(1440//minutes,center=True).mean()
 std24=anomaly.rolling(1440//minutes,center=True).std(ddof=0)
 def maximum(series,skip):
  matrix=pd.concat([series.shift(-k) for k in range(-20//minutes,40//minutes+1)],axis=1)
  return matrix.max(axis=1,skipna=skip)
 maximum_anomaly=maximum(anomaly,False)
 maximum_speed=maximum(u,True)
 return pd.DataFrame({'max_anomaly':maximum_anomaly,'threshold_eq2':mean24+3*std24,'max_speed':maximum_speed,'threshold_eq3':mean.shift(80//minutes)+2*std.shift(80//minutes)},index=frame.index)

rows=[];cache={};audit_windows=[]
for i,e in extra.iterrows():
 key=(e.station_id,e.t0.month)
 if key not in cache:
  start=pd.Timestamp(2025,e.t0.month,1);end=start+pd.offsets.MonthBegin(1)
  path=args.data_root/e.station_id/f'{e.station_id.lower()}_AWS_Level2monthly_{e.t0:%Y%m}.nc'
  f=new.load_month(path,start,end,new.Timing(1));_,_,m=new.classify_frame(f,new.Timing(1),'legacy')
  f5=old.load_month(path,start,end,old.Timing(5));_,_,m5=old.classify_frame(f5,old.Timing(5),'legacy')
  cache[key]=(f,m,f5,m5,wind_numbers(f,1),wind_numbers(f5,5))
 f,m,f5,m5,w,w5=cache[key]
 tg=new.event_groups(m[0],new.Timing(1),'legacy')
 group=next((g for g in tg if f.index[g[0]]==e.t0),None)
 assert group is not None
 start,end=group
 thermal=m[0].fillna(False).astype(bool)
 confirmed=thermal.isel(time=slice(start,end+1))&m[1].isel(time=slice(start,end+1))
 confirmation=pd.Timestamp(confirmed.time.values[np.flatnonzero(confirmed.values)[0]]) if confirmed.any().item() else None
 row={'audit_id':f'NEW-{i+1:03d}','station_id':e.station_id,'new_onset_UTC':e.t0,'new_onset_WIB':e.t0+pd.Timedelta('7h'),'new_thermal_episode_end_UTC':f.index[end],
 'previous_5min_label_paired':e.previous_label,'paired_5min_onset_UTC':e.t0_5min,'paired_offset_5minus1_min':e.offset_5minus1_minutes,
 'first_joint_confirmation_UTC':confirmation,'new_wind_pass_at_onset':bool(m[1].sel(time=e.t0).item())}
 # A nearby event is context only; it is not a second one-to-one assignment.
 old_same=a.loc[a.station_id==e.station_id].copy();old_same['offset']=(old_same.t0-e.t0).dt.total_seconds()/60
 near=old_same.loc[old_same.offset.abs().idxmin()]
 row.update(nearest_5min_label_context=near.cp_type,nearest_5min_onset_context=near.t0,nearest_5min_offset_context_min=float(near.offset))
 ref=e.t0_5min if pd.notna(e.t0_5min) else e.t0.round('5min')
 row['diagnostic_5min_anchor_UTC']=ref
 if ref not in f5.index:ref=f5.index[f5.index.get_indexer([ref],method='nearest')[0]]
 vals=f5.temperature.reindex(pd.date_range(ref-pd.Timedelta('5min'),ref+pd.Timedelta('15min'),freq='5min'))
 row.update(old_anchor_T20_C=float(vals.iloc[-1]-vals.iloc[0]),old_anchor_T5_C=float(vals.iloc[1]-vals.iloc[0]),old_anchor_monotonic=bool(vals.notna().all() and (np.diff(vals)<=0).all()),
            old_anchor_temperature_pass=bool(m5[0].fillna(False).astype(bool).sel(time=ref).item()),old_anchor_wind_pass=bool(m5[1].sel(time=ref).item()))
 for name,value in w5.loc[ref].items():row['old_anchor_wind_'+name]=value
 if confirmation is not None:
  t=confirmation;raw=f.temperature.reindex(pd.date_range(t-pd.Timedelta('5min'),t+pd.Timedelta('15min'),freq='1min'))
  means=np.array([f.temperature.reindex(pd.date_range(t+pd.Timedelta(minutes=c-2),t+pd.Timedelta(minutes=c+2),freq='1min')).mean(skipna=False) for c in [-5,0,5,10,15]])
  row.update(new_confirm_T20_C=float(raw.iloc[-1]-raw.iloc[0]),new_confirm_T5_C=float(raw.iloc[5]-raw.iloc[0]),new_confirm_max_mean_increment_C=float(np.diff(means).max()),
   new_confirmation_equation='eq2' if m[2].sel(time=t).item() else 'eq3',new_temperature_and_wind_pass=True)
  for name,value in w.loc[t].items():row['new_confirm_wind_'+name]=value
  assert raw.iloc[-1]-raw.iloc[0]<=-2+1e-9 and raw.iloc[5]-raw.iloc[0]<=-.5+1e-9 and (np.diff(means)<=1e-9).all()
 else:row.update(new_confirmation_equation='overlap_of_group_intervals_only',new_temperature_and_wind_pass=False)
 # Examine the original 5-minute candidate masks around the new onset.
 left=e.t0-pd.Timedelta('20min');right=e.t0+pd.Timedelta('20min')
 tmask=m5[0].fillna(False).astype(bool).sel(time=slice(left,right));wmask=m5[1].sel(time=slice(left,right))
 row['old_thermal_pass_anchors_within20']=int(tmask.sum().item());row['old_joint_pass_anchors_within20']=int((tmask&wmask).sum().item())
 if e.previous_label=='temperature_only':reason='Previously temperature-only: 1-min wind/episode alignment now confirms thermal episode'
 elif e.previous_label=='wind_only':reason='Previously wind-only: new raw-temperature/mean-monotonic candidate now overlaps wind episode'
 elif row['old_joint_pass_anchors_within20']>0:reason='5-min joint candidate exists nearby; event onset/grouping or one-to-one assignment differs'
 elif row['old_thermal_pass_anchors_within20']>0:reason='Local 5-min thermal candidate exists but lacks simultaneous wind confirmation; inspect group timing'
 else:reason='No 5-min thermal candidate within20; temperature sampling/drop/monotonic or missing-data difference'
 row['diagnostic_reason']=reason
 rows.append(row)
 for stamp in pd.date_range(e.t0-pd.Timedelta('25min'),e.t0+pd.Timedelta('25min'),freq='1min'):
  if stamp in f.index:audit_windows.append({'audit_id':row['audit_id'],'time_UTC':stamp,'relative_minute':(stamp-e.t0).total_seconds()/60,'temperature_pass_1min':bool(thermal.sel(time=stamp).item()),'wind_pass_1min':bool(m[1].sel(time=stamp).item())})
 print(row['audit_id'],e.station_id,e.t0,e.previous_label,flush=True)
result=pd.DataFrame(rows);result.to_csv(P/'new_121_event_audit.csv',index=False)
pd.DataFrame(audit_windows).to_csv(P/'criterion_flags_near_new_onsets.csv',index=False)
summary=result.groupby(['station_id','previous_5min_label_paired']).size().unstack(fill_value=0)
assert len(result)==121 and result.audit_id.is_unique
assert result.previous_5min_label_paired.value_counts().to_dict()=={'no_matched_5min_event_within20':72,'temperature_only':44,'wind_only':5}
notes=pd.DataFrame([
 ('Definition','Additional = new-1min cold pool unmatched to 5min coldpool at same station within ±20 min, using one-to-one pairing. Not proof of a new physical cold pool.'),
 ('Previous label','Assigned using the existing cross-category one-to-one matching, same station ±20 min. No match is not the same as no meteorological event.'),
 ('Diagnostic anchor','Paired old onset when available; otherwise nearest 5min grid. This is not a controlled same-time counterfactual.'),
 ('Onset vs confirmation','Catalog onset is start of the thermal episode. first_joint_confirmation_UTC is first actual anchor in that episode passing both criteria.'),
 ('Raw temperature','T20 <=−2 C, T5 <=−0.5 C. Only monotonicity uses five centered means; four increments must be <=0.'),
 ('Wind','Eq2 uses anomaly; Eq3 is conditional on missing anomaly. Threshold values are recalculated at each anchor and depend on gust resolution.'),
 ('Caution','Diagnostics show criterion differences and temporal associations, not independently validated physical identities or causal attribution of a single changed component.'),
 ('Time','UTC and WIB=UTC+7. Thresholds and differences are in C or m/s.'),
],columns=['item','definition'])
with pd.ExcelWriter(P/'New_121_ColdPool_Event_Audit.xlsx',engine='openpyxl',datetime_format='yyyy-mm-dd hh:mm') as writer:
 result.to_excel(writer,sheet_name='ALL_121_ADDITIONAL',index=False)
 for label,sheet in [('temperature_only','PREV_TEMPERATURE_44'),('wind_only','PREV_WIND_5'),('no_matched_5min_event_within20','NO_PAIRED_EVENT_72')]:result.loc[result.previous_5min_label_paired==label].to_excel(writer,sheet_name=sheet,index=False)
 summary.reset_index().to_excel(writer,sheet_name='SUMMARY',index=False)
 notes.to_excel(writer,sheet_name='METHOD',index=False)
 for ws in writer.book.worksheets:ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
report={'additional':len(result),'previous_labels':result.previous_5min_label_paired.value_counts().to_dict(),'by_station':result.station_id.value_counts().to_dict(),'confirmed_joint_anchors':int(result.new_temperature_and_wind_pass.sum()),'reasons':result.diagnostic_reason.value_counts().to_dict(),'new_wind_equations':result.new_confirmation_equation.value_counts().to_dict()}
(P/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
