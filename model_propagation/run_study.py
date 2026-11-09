"""Run conceptual propagation experiments and exploratory BAMNET source fits."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from propagation import arrival_times,front_points,invert

parser=argparse.ArgumentParser()
parser.add_argument('--data-root',type=Path,default=Path('/workspace/NewJournal'))
parser.add_argument('--output',type=Path,default=Path(__file__).parent)
parser.add_argument('--catalog',type=Path,default=None)
args=parser.parse_args();ROOT=args.data_root;OUT=args.output
for p in [OUT,OUT/'inputs',OUT/'tables',OUT/'plots',OUT/'animations']:p.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(20261109)
catalog=args.catalog or ROOT/'hasil_codex_new_monotonic/catalogs/5min/all_events.csv'
all_events=pd.read_csv(catalog,parse_dates=['t0']);events=all_events.loc[all_events.cp_type=='coldpool'].sort_values(['t0','station_id']).reset_index(drop=True)
events=events.drop(columns=['event_id'],errors='ignore')
assert len(events)==197,'This study must use the requested 197-event 5-minute catalog'
events.insert(0,'event_id',[f'CP5-{i:03d}' for i in range(1,len(events)+1)])
events.to_csv(OUT/'inputs/coldpool_197_5min.csv',index=False)
metadata=[];station_series={}
for station in sorted(events.station_id.unique()):
 frames=[]
 for f in sorted((ROOT/station).glob('*Level2monthly_2025*.nc')):
  with netCDF4.Dataset(f) as ds:
   if not any(x['station_id']==station for x in metadata):metadata.append({'station_id':station,'longitude':float(ds.getncattr('station_longitude')),'latitude':float(ds.getncattr('station_latitude')),'altitude_m_ASL':float(ds.getncattr('station_altitudeASL'))})
   t=ds['time_5min'];unit,origin=t.units.split(' since ');assert unit=='minutes'
   index=pd.Timestamp(origin)+pd.to_timedelta(np.asarray(t[:]),unit='min')
   frames.append(pd.DataFrame({key:np.ma.filled(ds[var][:].astype(float),np.nan) for key,var in [('u','u_wind_5min'),('v','v_wind_5min'),('T','air_temperature_5min')]},index=index))
 s=pd.concat(frames).sort_index();station_series[station]=s.loc[~s.index.duplicated()]
stations=pd.DataFrame(metadata);lon0=stations.longitude.mean();lat0=stations.latitude.mean()
stations['x_km']=(stations.longitude-lon0)*111.195*np.cos(np.deg2rad(lat0));stations['y_km']=(stations.latitude-lat0)*111.195
stations.to_csv(OUT/'inputs/station_coordinates.csv',index=False)
xy=stations[['x_km','y_km']].to_numpy();_,sv,vt=np.linalg.svd(xy-xy.mean(axis=0),full_matrices=False)
axis_angle=float(np.degrees(np.arctan2(vt[0,1],vt[0,0]))%180)
geometry={'longitude_origin':lon0,'latitude_origin':lat0,'projection':'local equirectangular, km; planar approximation','principal_spreads_km':sv.tolist(),'minor_major_spread_ratio':float(sv[1]/sv[0]),'principal_axis_degrees_CCW_from_east':axis_angle,'altitude_range_m':stations.altitude_m_ASL.tolist()}
(OUT/'tables/geometry.json').write_text(json.dumps(geometry,indent=2))
dist=np.linalg.norm(xy[:,None]-xy[None,:],axis=2)
pd.DataFrame(dist,index=stations.station_id,columns=stations.station_id).to_csv(OUT/'tables/station_distances_km.csv')
fig,ax=plt.subplots(figsize=(8,8));sc=ax.scatter(xy[:,0],xy[:,1],c=stations.altitude_m_ASL,s=90,cmap='terrain')
for _,r in stations.iterrows():ax.annotate(r.station_id,(r.x_km,r.y_km),xytext=(5,5),textcoords='offset points')
fig.colorbar(sc,ax=ax,label='Altitude (m ASL)');ax.set(xlabel='East (km)',ylabel='North (km)',aspect='equal',title=f'BAMNET geometry: minor/major spread = {sv[1]/sv[0]:.3f}')
ax.grid(alpha=.2);fig.tight_layout();fig.savefig(OUT/'plots/station_geometry.png',dpi=170);plt.close(fig)

def associate(window):
 groups=[];current=[];first=None
 for _,e in events.iterrows():
  if first is None or (e.t0-first).total_seconds()/60>window:
   if current:groups.append(pd.DataFrame(current))
   current=[e];first=e.t0
  else:current.append(e)
 if current:groups.append(pd.DataFrame(current))
 return groups

sens=[]
for w in [60,90,120,180]:
 groups=associate(w)
 sens.append({'window_minutes':w,'groups':len(groups),'single_station':sum(g.station_id.nunique()==1 for g in groups),'at_least_2':sum(g.station_id.nunique()>=2 for g in groups),'at_least_3':sum(g.station_id.nunique()>=3 for g in groups),'at_least_4':sum(g.station_id.nunique()>=4 for g in groups),'at_least_5':sum(g.station_id.nunique()>=5 for g in groups),'repeated_station_groups':sum(g.station_id.duplicated().any() for g in groups)})
pd.DataFrame(sens).to_csv(OUT/'tables/association_sensitivity.csv',index=False)
groups=associate(180);members=[];group_summary=[];real_fits=[];real_solutions=[];group_info={}
for j,g in enumerate(groups,1):
 gid=f'G{j:03d}';g=g.copy();g['group_id']=gid;members.append(g)
 repeated=g.station_id.duplicated().any()
 row={'group_id':gid,'start_UTC':g.t0.min(),'end_UTC':g.t0.max(),'span_minutes':(g.t0.max()-g.t0.min()).total_seconds()/60,'detections':len(g),'unique_stations':g.station_id.nunique(),'repeated_station':bool(repeated),'association_status':'temporal candidate, not confirmed same physical cold pool'}
 wind_rows=[]
 for _,e in g.iterrows():
  pre=station_series[e.station_id].reindex(pd.date_range(e.t0-pd.Timedelta('40min'),e.t0-pd.Timedelta('20min'),freq='5min'))
  valid=np.isfinite(pre[['u','v']]).all(axis=1)
  wind_rows.append({'event_id':e.event_id,'station_id':e.station_id,'u_before':float(pre.loc[valid,'u'].mean()) if valid.sum()>=4 else np.nan,'v_before':float(pre.loc[valid,'v'].mean()) if valid.sum()>=4 else np.nan,'n_background':int(valid.sum())})
 background=pd.DataFrame(wind_rows)
 wind=np.array([background.u_before.median(),background.v_before.median()]);wind_ok=np.isfinite(wind).all()
 row.update(background_u_ms=wind[0],background_v_ms=wind[1],background_stations=int(background.u_before.notna().sum()))
 group_summary.append(row)
 if len(g)<3 or repeated:continue
 positions=stations.set_index('station_id').loc[g.station_id,['x_km','y_km']].to_numpy()
 times=(g.t0-g.t0.min()).dt.total_seconds().to_numpy()/60
 theta=float(np.degrees(np.arctan2(wind[1],wind[0]))) if wind_ok and np.linalg.norm(wind)>=.5 else axis_angle
 group_info[gid]=(g,positions,times,wind,theta)
 background.assign(group_id=gid).to_csv(OUT/'tables'/f'{gid}_background.csv',index=False)
 scenarios=[('circle_no_wind',np.zeros(2),1.,0.)]
 if wind_ok:scenarios.extend([('circle_background',wind,1.,0.),('ellipse_fixed_shape',wind,1.5,theta),('ellipse_cross_orientation',wind,1.5,theta+90)])
 for label,w,r,angle in scenarios:
  for speed in [3.,5.,8.]:
   best,solutions=invert(positions,times,wind=w,fixed_speed=speed,ratio=r,angle=angle,starts=14,seed=j)
   if best is None:
    real_fits.append({'group_id':gid,'model':label,'fixed_speed_ms':speed,'status':'skipped: prescribed wind is supercritical for this shape/speed'})
    continue
   origin=g.t0.min()+pd.Timedelta(minutes=best['t_start'])
   result={'group_id':gid,'model':label,'fixed_speed_ms':speed,'n_stations':len(g),'degrees_of_freedom':len(g)-3,'status':'conditional exploratory fit; not validated origin','origin_time_UTC':origin,
    'origin_longitude':lon0+best['x']/(111.195*np.cos(np.deg2rad(lat0))),'origin_latitude':lat0+best['y']/111.195,
    'wind_u_ms':w[0],'wind_v_ms':w[1],'domain_halfwidth_km':100.,**best}
   real_fits.append(result)
   if speed==5:
    for k,sol in enumerate(solutions):real_solutions.append({'group_id':gid,'model':label,'start_id':k,**sol})
 if len(g)>=4:
  for label,w,r,angle in scenarios[:3]:
   best,solutions=invert(positions,times,wind=w,ratio=r,angle=angle,starts=24,seed=j)
   if best:
    real_fits.append({'group_id':gid,'model':label+'_unknown_speed','fixed_speed_ms':np.nan,'n_stations':len(g),'degrees_of_freedom':len(g)-4,'status':'exploratory, zero residual degrees of freedom; speed also fitted','origin_time_UTC':g.t0.min()+pd.Timedelta(minutes=best['t_start']),
     'origin_longitude':lon0+best['x']/(111.195*np.cos(np.deg2rad(lat0))),'origin_latitude':lat0+best['y']/111.195,'wind_u_ms':w[0],'wind_v_ms':w[1],'domain_halfwidth_km':100.,**best})
 print('Fitted candidate',gid,len(g),flush=True)
pd.concat(members).to_csv(OUT/'tables/candidate_group_members.csv',index=False)
pd.DataFrame(group_summary).to_csv(OUT/'tables/candidate_groups.csv',index=False)
real=pd.DataFrame(real_fits);real.to_csv(OUT/'tables/real_conditional_fits.csv',index=False)
pd.DataFrame(real_solutions).to_csv(OUT/'tables/real_multistart_solutions.csv',index=False)

# Controlled synthetic experiment on a spatially distributed nine-station array.
synthetic_xy=np.array([[-30,-25],[0,-30],[30,-25],[-35,0],[0,0],[35,0],[-30,25],[0,30],[30,25]],float)
truth={'source':np.array([-12.,-8.]),'t_start':-60.,'speed':6.}
models=[('circle_no_wind',(0.,0.),1.,0.,False),('circle_background',(2.,1.),1.,0.,False),('ellipse_background',(2.,1.),1.6,30.,True)]
synthetic_results=[];synthetic_times=[];mc=[];synthetic_best={}
for k,(name,wind,ratio,angle,free) in enumerate(models):
 exact=arrival_times(synthetic_xy,truth['source'],truth['t_start'],truth['speed'],wind,ratio,angle)
 best,_=invert(synthetic_xy,exact,wind=wind,free_shape=free,starts=30,seed=100+k)
 best['source_error_km']=float(np.linalg.norm(np.array([best['x'],best['y']])-truth['source']))
 synthetic_results.append({'model':name,'truth_x':-12.,'truth_y':-8.,'truth_t_start':-60.,'truth_speed_ms':6.,'truth_ratio':ratio,'truth_angle':angle,**best})
 synthetic_best[name]=best
 for i,t in enumerate(exact):synthetic_times.append({'model':name,'station':f'S{i+1}','x':synthetic_xy[i,0],'y':synthetic_xy[i,1],'arrival_minutes':t})
 for run in range(35):
  observed=exact+rng.uniform(-2.5,2.5,len(exact))
  estimate,_=invert(synthetic_xy,observed,wind=wind,free_shape=free,starts=2,initial=best['parameters_vector'],seed=run)
  mc.append({'model':name,'run':run,'source_error_km':float(np.linalg.norm(np.array([estimate['x'],estimate['y']])-truth['source'])),**estimate})
 print('Synthetic model',name,'source error',best['source_error_km'],flush=True)
pd.DataFrame(synthetic_results).to_csv(OUT/'tables/synthetic_recovery.csv',index=False)
pd.DataFrame(synthetic_times).to_csv(OUT/'tables/synthetic_arrivals.csv',index=False)
pd.DataFrame(mc).to_csv(OUT/'tables/synthetic_timing_sensitivity.csv',index=False)
mc_summary=pd.DataFrame(mc).groupby('model').source_error_km.agg(['count','median','max'])
mc_summary.to_csv(OUT/'tables/synthetic_uncertainty_summary.csv')

fig,axes=plt.subplots(1,3,figsize=(16,5.5),sharex=True,sharey=True)
for ax,(name,wind,ratio,angle,free) in zip(axes,models):
 ax.scatter(synthetic_xy[:,0],synthetic_xy[:,1],color='black',s=24,label='Point stations')
 ax.scatter(*truth['source'],marker='*',s=160,color='red',label='True origin')
 best=synthetic_best[name];ax.scatter(best['x'],best['y'],marker='x',s=75,color='blue',label='Inverted origin')
 for elapsed in [25,50,75,100,125]:
  front=front_points(truth['source'],-60,-60+elapsed,6,wind,ratio,angle)
  ax.plot(front[:,0],front[:,1],alpha=.6,label=f'{elapsed} min' if elapsed==25 else None)
 ax.set(title=name.replace('_',' '),xlabel='East (km)',aspect='equal',xlim=(-65,75),ylim=(-65,65));ax.grid(alpha=.15)
axes[0].set_ylabel('North (km)');axes[0].legend(fontsize=7);fig.suptitle('Synthetic demonstration: source known, wind prescribed; ellipse shape fitted only on 9 stations')
fig.tight_layout();fig.savefig(OUT/'plots/conceptual_three_models.png',dpi=160);plt.close(fig)

fig,axes=plt.subplots(1,3,figsize=(15,5),sharex=True,sharey=True)
for ax,(name,*_) in zip(axes,models):
 d=pd.DataFrame(mc).loc[lambda x:x.model==name]
 ax.scatter(d.x,d.y,alpha=.5,s=15);ax.scatter(-12,-8,marker='*',s=150,color='red')
 ax.set(title=name.replace('_',' '),xlabel='Estimated origin east (km)',aspect='equal');ax.grid(alpha=.15)
axes[0].set_ylabel('Estimated origin north (km)');fig.suptitle('Timing sensitivity: 35 uniform ±2.5-minute perturbations; conditional on correct model and known wind')
fig.tight_layout();fig.savefig(OUT/'plots/synthetic_timing_sensitivity.png',dpi=160);plt.close(fig)

# Identifiability with the actual geometry; no claim of ellipse recovery from five stations.
actual_synthetic=[]
for name,wind,ratio,angle,free in models:
 exact=arrival_times(xy,truth['source'],-60,6,wind,ratio,angle)
 best,solutions=invert(xy,exact,wind=wind,free_shape=free,starts=40,seed=30)
 for i,sol in enumerate(solutions):actual_synthetic.append({'model':name,'start_id':i,'true_x':-12.,'true_y':-8.,**sol})
pd.DataFrame(actual_synthetic).to_csv(OUT/'tables/bamnet_geometry_synthetic_multisolutions.csv',index=False)

# Real-data origin profile for a four-station temporal candidate, c=5 m/s.
four=[gid for gid,info in group_info.items() if len(info[0])>=4]
profile_info={}
if four:
 example=four[-1];g,pos,ts,wind,theta=group_info[example]
 x=np.linspace(-100,100,161);y=np.linspace(-100,100,161);xx,yy=np.meshgrid(x,y)
 distance=np.sqrt((xx[:,:,None]-pos[:,0])**2+(yy[:,:,None]-pos[:,1])**2)
 travel=distance/(5*.06);t0=np.clip(np.mean(ts-travel,axis=2),ts.min()-360,ts.min()-1e-5)
 loss=np.sqrt(np.mean((t0[:,:,None]+travel-ts)**2,axis=2))
 fig,ax=plt.subplots(figsize=(8,7));image=ax.pcolormesh(xx,yy,np.minimum(loss,20),shading='auto',cmap='viridis_r');fig.colorbar(image,ax=ax,label='Arrival RMSE (min); scale capped at 20')
 if loss.min()+2.5<loss.max():ax.contour(xx,yy,loss,levels=[loss.min()+2.5],colors='white')
 ax.scatter(xy[:,0],xy[:,1],color='black',s=30)
 for _,r in stations.iterrows():ax.text(r.x_km,r.y_km,r.station_id,fontsize=7)
 d=real.loc[(real.group_id==example)&(real.model=='circle_no_wind')&(real.fixed_speed_ms==5)]
 for _,r in d.iterrows():ax.scatter(r.x,r.y,marker='*',s=150,color='red')
 ax.set(xlabel='Origin east (km)',ylabel='Origin north (km)',aspect='equal',title=f'{example}: conditional source loss, circle with fixed c=5 m/s\nWhite = best-grid RMSE +2.5 min, NOT a confidence contour')
 fig.tight_layout();fig.savefig(OUT/'plots/real_origin_loss_profile.png',dpi=170);plt.close(fig)
 np.savez_compressed(OUT/'tables/real_origin_profile.npz',x=x,y=y,rmse=loss)
 profile_info={'group_id':example,'start_UTC':str(g.t0.min()),'stations':g.station_id.tolist(),'grid_best_rmse_minutes':float(loss.min()),'near_optimal_grid_cells':int((loss<=loss.min()+2.5).sum())}
 domain=[]
 for extent in [50.,100.,150.]:
  best,_=invert(pos,ts,fixed_speed=5,extent=extent,starts=25)
  domain.append({'group_id':example,'extent_km':extent,**best})
 pd.DataFrame(domain).to_csv(OUT/'tables/domain_sensitivity.csv',index=False)

# Three forward-model animations from the same synthetic source.
fig,axes=plt.subplots(1,3,figsize=(13,4.6));lines=[];dots=[]
for ax,(name,wind,ratio,angle,free) in zip(axes,models):
 ax.scatter(synthetic_xy[:,0],synthetic_xy[:,1],c='grey',s=24);ax.scatter(*truth['source'],marker='*',s=100,color='red')
 line,=ax.plot([],[],lw=2);dot=ax.scatter([],[],color='green',s=36);lines.append(line);dots.append(dot)
 ax.set(xlim=(-65,75),ylim=(-65,65),aspect='equal',xlabel='East (km)',title=name.replace('_',' '))
axes[0].set_ylabel('North (km)');title=fig.suptitle('')
def update(elapsed):
 for i,(name,wind,ratio,angle,free) in enumerate(models):
  front=front_points(truth['source'],-60,-60+elapsed,6,wind,ratio,angle)
  lines[i].set_data(front[:,0],front[:,1])
  arrived=arrival_times(synthetic_xy,truth['source'],-60,6,wind,ratio,angle)<=-60+elapsed
  dots[i].set_offsets(synthetic_xy[arrived])
 title.set_text(f'Forward model: {elapsed:.0f} minutes after release; green = front has reached station')
 return lines+dots
animation=FuncAnimation(fig,update,frames=np.linspace(0,150,61),interval=100)
animation.save(OUT/'animations/three_models.gif',writer=PillowWriter(fps=10));plt.close(fig)

# A static summary of observed arrival times for multi-station temporal candidates.
fig,ax=plt.subplots(figsize=(11,6))
selected=[(gid,info) for gid,info in group_info.items()]
for j,(gid,(g,pos,ts,w,theta)) in enumerate(selected):
 ax.scatter(ts,np.full(len(g),j),c=[int(s[-2:]) for s in g.station_id],cmap='viridis',vmin=1,vmax=5,s=55)
 for t,station in zip(ts,g.station_id):ax.annotate(station,(t,j),xytext=(3,4),textcoords='offset points',fontsize=6)
ax.set(yticks=np.arange(len(selected)),yticklabels=[f'{gid} {info[0].t0.min():%Y-%m-%d}' for gid,info in selected],xlabel='Minutes from earliest detection in group',title='Candidate associations: time proximity alone does not prove a common cold pool')
fig.tight_layout();fig.savefig(OUT/'plots/candidate_arrivals.png',dpi=160);plt.close(fig)

parameter_table=pd.DataFrame([
 {'model':'circle, no wind','unknowns':'x0,y0,t0,c','p':4,'minimum_arrivals_for_local_count':4,'note':'4 observations leave zero residual df; nonlinear/global uniqueness not guaranteed'},
 {'model':'circle, prescribed wind','unknowns':'x0,y0,t0,c','p':4,'minimum_arrivals_for_local_count':4,'note':'wind prescribed, not inferred; unknown wind would add 2 parameters'},
 {'model':'ellipse, prescribed wind','unknowns':'x0,y0,t0,c,axis_ratio,orientation','p':6,'minimum_arrivals_for_local_count':6,'note':'cannot freely identify from at most 5 single-pass BAMNET arrivals'},
 {'model':'real-data conditional scenarios','unknowns':'x0,y0,t0','p':3,'minimum_arrivals_for_local_count':3,'note':'speed, wind, shape are fixed assumptions; not independent estimates'}])
parameter_table.to_csv(OUT/'tables/identifiability.csv',index=False)
with pd.ExcelWriter(OUT/'Propagation_Results.xlsx',engine='openpyxl') as writer:
 for name,frame in [('Input_197',events),('Stations',stations),('Association_sensitivity',pd.DataFrame(sens)),('Groups',pd.DataFrame(group_summary)),('Group_members',pd.concat(members)),('Conditional_real_fits',real.drop(columns=['parameters_vector'],errors='ignore')),('Real_multistart',pd.DataFrame(real_solutions).drop(columns=['parameters_vector'],errors='ignore')),('Synthetic_recovery',pd.DataFrame(synthetic_results).drop(columns=['parameters_vector'],errors='ignore')),('Synthetic_timing',pd.DataFrame(mc).drop(columns=['parameters_vector'],errors='ignore')),('Identifiability',parameter_table)]:
  frame.to_excel(writer,sheet_name=name,index=False)
 for sheet in writer.book.worksheets:sheet.freeze_panes='A2';sheet.auto_filter.ref=sheet.dimensions
manifest={'source_catalog':str(catalog.relative_to(ROOT)) if catalog.is_relative_to(ROOT) else str(catalog),'source_catalog_sha256':hashlib.sha256(catalog.read_bytes()).hexdigest(),'coldpool_count':197,'source_resolution_minutes':5,'reference_projection':geometry,'association_window_minutes':180,'association_sensitivity':sens,'candidate_groups_with_at_least_3':len(group_info),'candidate_groups_with_4':len(four),'free_ellipse_real_data':'not fitted: at most four associated detections, fewer than six free parameters','synthetic_MC':'35 uniform independent ±2.5-minute timing perturbations per model; assumes correct model and fixed wind','real_profile':profile_info,'seed':20261109,'units':'km, minutes, m/s; angle degrees CCW from east','real_fits':'hypotheses conditional on association, fixed speed/shape/wind and domain; not validated origins'}
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2,default=str))
print(json.dumps(manifest,indent=2,default=str),flush=True)
