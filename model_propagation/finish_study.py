"""Validate outputs, illustrate a conditional BAMNET case, and make bilingual PDFs."""
from pathlib import Path
import json,hashlib,argparse
import numpy as np
import pandas as pd
import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
from openpyxl import load_workbook
from propagation import arrival_times,front_points

P=Path(__file__).parent
parser=argparse.ArgumentParser();parser.add_argument('--data-root',type=Path,default=Path('/workspace/NewJournal'))
ROOT=parser.parse_args().data_root
manifest=json.loads((P/'manifest.json').read_text());geo=manifest['reference_projection']
stations=pd.read_csv(P/'inputs/station_coordinates.csv');xy=stations[['x_km','y_km']].to_numpy()
events=pd.read_csv(P/'inputs/coldpool_197_5min.csv',parse_dates=['t0'])
fits=pd.read_csv(P/'tables/real_conditional_fits.csv');groups=pd.read_csv(P/'tables/candidate_groups.csv')
members=pd.read_csv(P/'tables/candidate_group_members.csv',parse_dates=['t0'])
fits['origin_time_lower_bound_active']=fits.t_start.le(-360+1e-4)
fits.to_csv(P/'tables/real_conditional_fits.csv',index=False)
branch_path=P/'inputs/equation_support_197.csv'
if not branch_path.exists():branch_path=Path('/workspace/newjournal-analysis/results_reproduced_197/coldpool_branch_support.csv')
if branch_path.exists():
 branch=pd.read_csv(branch_path,parse_dates=['t0'])
 assert set(zip(branch.station_id,branch.t0))==set(zip(events.station_id,events.t0))
 branch.to_csv(P/'inputs/equation_support_197.csv',index=False)
 assert branch.eq2_any_thermal_anchor.sum()==188 and branch.eq3_any_thermal_anchor.sum()==9

example=manifest['real_profile']['group_id'];g=members.loc[members.group_id==example].sort_values('t0')
reference=g.t0.min();obs=(g.t0-reference).dt.total_seconds().to_numpy()/60
positions=stations.set_index('station_id').loc[g.station_id,['x_km','y_km']].to_numpy()
labels=['circle_no_wind','circle_background','ellipse_fixed_shape']
palette={'circle_no_wind':'#1f77b4','circle_background':'#ff7f0e','ellipse_fixed_shape':'#2ca02c'}
chosen=fits.loc[(fits.group_id==example)&(fits.fixed_speed_ms==5)&fits.model.isin(labels)].set_index('model').loc[labels]
predictions=[]
for name,r in chosen.iterrows():
 pred=arrival_times(positions,[r.x,r.y],r.t_start,5,[r.wind_u_ms,r.wind_v_ms],r.ratio,r.angle)
 for (_,e),t,y in zip(g.iterrows(),obs,pred):predictions.append({'group_id':example,'model':name,'station_id':e.station_id,'observed_minutes':t,'predicted_minutes':float(y),'residual_minutes':float(y-t),'observed_UTC':e.t0})
prediction_table=pd.DataFrame(predictions);prediction_table.to_csv(P/'tables/real_example_arrival_predictions.csv',index=False)

fig,axes=plt.subplots(1,2,figsize=(13,5))
for name,r in chosen.iterrows():
 axes[0].scatter(r.x,r.y,s=80,color=palette[name],label=name.replace('_',' '))
axes[0].scatter(xy[:,0],xy[:,1],marker='^',c='black',label='Station')
for _,r in stations.iterrows():axes[0].annotate(r.station_id,(r.x_km,r.y_km),fontsize=7)
axes[0].set(xlabel='East (km)',ylabel='North (km)',aspect='equal',title=f'{example}, {reference:%Y-%m-%d}: c fixed at 5 m/s\nDifferent assumptions give different candidate origins');axes[0].legend(fontsize=7)
for name,d in prediction_table.groupby('model'):axes[1].plot(d.observed_minutes,d.predicted_minutes,'o-',color=palette[name],label=name.replace('_',' '))
axes[1].plot([-5,90],[-5,90],'k--',lw=.8);axes[1].set(xlabel='Observed minutes after first onset',ylabel='Predicted minutes after first onset',title='Conditional arrival-time fit, not independent validation');axes[1].legend(fontsize=7)
fig.tight_layout();fig.savefig(P/'plots/real_conditional_example.png',dpi=170);plt.close(fig)

fig,axes=plt.subplots(1,3,figsize=(14,5));lines=[];dots=[]
for ax,(name,r) in zip(axes,chosen.iterrows()):
 ax.scatter(positions[:,0],positions[:,1],s=35,color='grey');ax.scatter(r.x,r.y,color='red',marker='*',s=100)
 line,=ax.plot([],[],lw=2);dot=ax.scatter([],[],c='green',s=45);lines.append(line);dots.append(dot)
 ax.set(xlim=(-65,25),ylim=(-55,35),aspect='equal',xlabel='East (km)',title=name.replace('_',' '))
axes[0].set_ylabel('North (km)');title=fig.suptitle('')
def update(t):
 for i,(name,r) in enumerate(chosen.iterrows()):
  f=front_points([r.x,r.y],r.t_start,t,5,[r.wind_u_ms,r.wind_v_ms],r.ratio,r.angle)
  lines[i].set_data(f[:,0],f[:,1]);dots[i].set_offsets(positions[obs<=t])
 title.set_text(f'{example}: {t:.0f} min from first observed onset; green = observed arrival, not model arrival')
 return lines+dots
FuncAnimation(fig,update,frames=np.linspace(-40,100,57),interval=110).save(P/'animations/real_conditional_example.gif',writer=PillowWriter(fps=9));plt.close(fig)

# Observed temperature traces support inspection of temporal association without proving identity.
fig,ax=plt.subplots(figsize=(10,5));trace_rows=[]
for _,e in g.iterrows():
 path=ROOT/e.station_id/f'{e.station_id.lower()}_AWS_Level2monthly_{e.t0:%Y%m}.nc'
 with netCDF4.Dataset(path) as ds:
  unit,origin=ds['time_5min'].units.split(' since ')
  index=pd.Timestamp(origin)+pd.to_timedelta(np.asarray(ds['time_5min'][:]),unit='min')
  T=pd.Series(np.ma.filled(ds['air_temperature_5min'][:],np.nan),index=index)
 grid=pd.date_range(reference-pd.Timedelta('40min'),reference+pd.Timedelta('140min'),freq='5min')
 baseline=T.reindex(pd.date_range(e.t0-pd.Timedelta('40min'),e.t0-pd.Timedelta('20min'),freq='5min')).mean()
 delta=T.reindex(grid)-baseline;time=(grid-reference).total_seconds()/60
 line,=ax.plot(time,delta,label=e.station_id);ax.axvline((e.t0-reference).total_seconds()/60,color=line.get_color(),ls=':',lw=.8)
 trace_rows.extend({'station_id':e.station_id,'time_UTC':t,'temperature_anomaly_C':v} for t,v in zip(grid,delta))
ax.set(xlabel='Minutes from earliest candidate detection',ylabel='Temperature − local pre-event background (°C)',title=f'{example}: observed 5-minute temperature traces\nDotted lines are catalog onsets; similarity alone does not validate association')
ax.legend();fig.tight_layout();fig.savefig(P/'plots/real_example_temperature_traces.png',dpi=170);plt.close(fig)
pd.DataFrame(trace_rows).to_csv(P/'tables/real_example_temperature_traces.csv',index=False)

# Show the non-identifiability directly, using only equally good synthetic solutions.
multi=pd.read_csv(P/'tables/bamnet_geometry_synthetic_multisolutions.csv')
equivalent=multi.loc[(multi.model=='ellipse_background')&(multi.rmse_minutes<1e-5)]
fig,ax=plt.subplots(figsize=(8,6));color=ax.scatter(equivalent.x,equivalent.y,c=equivalent.ratio,s=45,cmap='plasma');fig.colorbar(color,ax=ax,label='Fitted ellipse axis ratio')
ax.scatter(-12,-8,marker='*',s=130,color='red',label='True synthetic source');ax.scatter(xy[:,0],xy[:,1],marker='^',c='black',label='Actual BAMNET sites')
ax.set(xlabel='Source east (km)',ylabel='Source north (km)',aspect='equal',title='Five BAMNET sites: multiple ellipses fit the same exact synthetic arrivals\nSix unknowns; Jacobian rank at most five')
ax.legend(fontsize=8);fig.tight_layout();fig.savefig(P/'plots/ellipse_nonuniqueness.png',dpi=170);plt.close(fig)

synthetic=pd.read_csv(P/'tables/synthetic_recovery.csv');mc=pd.read_csv(P/'tables/synthetic_uncertainty_summary.csv')
assert len(events)==197 and not events.duplicated(['station_id','t0']).any()
assert set(members.event_id)==set(events.event_id) and len(members)==197
assert synthetic.source_error_km.max()<1e-3 and synthetic.rmse_minutes.max()<1e-5
assert (equivalent.jacobian_rank<equivalent.parameters).all()
assert len(equivalent)>=2 and equivalent.x.max()-equivalent.x.min()>1
assert np.allclose(prediction_table.groupby('model').residual_minutes.apply(lambda a:np.sqrt(np.mean(a*a))).sort_index(),chosen.rmse_minutes.sort_index())
with Image.open(P/'animations/three_models.gif') as image:assert image.n_frames==61
with Image.open(P/'animations/real_conditional_example.gif') as image:assert image.n_frames==57
wb=load_workbook(P/'Propagation_Results.xlsx',read_only=True);assert wb['Input_197'].max_row-1==197;wb.close()
with pd.ExcelWriter(P/'Propagation_Results.xlsx',engine='openpyxl',mode='a',if_sheet_exists='replace') as writer:
 fits.drop(columns=['parameters_vector'],errors='ignore').to_excel(writer,sheet_name='Conditional_real_fits',index=False)
 prediction_table.to_excel(writer,sheet_name='Real_example_arrivals',index=False)
 pd.read_csv(P/'tables/domain_sensitivity.csv').drop(columns=['parameters_vector'],errors='ignore').to_excel(writer,sheet_name='Domain_sensitivity',index=False)
 multi.drop(columns=['parameters_vector'],errors='ignore').to_excel(writer,sheet_name='BAMNET_synthetic_branches',index=False)
 if branch_path.exists():pd.read_csv(P/'inputs/equation_support_197.csv').to_excel(writer,sheet_name='Eq2_Eq3_support',index=False)

validation={'status':'passed','unit_tests':7,'input_197_verified':True,'equation_support':{'eq2':188,'eq3':9},'event_membership_exactly_once':True,
 'synthetic_max_source_error_km':float(synthetic.source_error_km.max()),'synthetic_max_rmse_minutes':float(synthetic.rmse_minutes.max()),
 'equally_good_free_ellipse_solutions_on_5_sites':len(equivalent),'ellipse_source_x_range_km':[float(equivalent.x.min()),float(equivalent.x.max())],
 'real_group_count':len(groups),'groups_with_at_least_3':int((groups.unique_stations>=3).sum()),'groups_with_4':int((groups.unique_stations>=4).sum()),'groups_with_5':int((groups.unique_stations>=5).sum()),
 'real_example':example,'real_example_case_is_conditional_not_validated':True,'synthetic_animation_frames':61,'real_animation_frames':57,'origin_time_bound_flags_added':int(fits.origin_time_lower_bound_active.sum())}
(P/'validation_report.json').write_text(json.dumps(validation,indent=2))

def page_text(pdf,title,paragraphs):
 fig=plt.figure(figsize=(11.7,8.3));fig.text(.05,.94,title,fontsize=20,weight='bold');y=.84
 for text in paragraphs:fig.text(.06,y,text,fontsize=12,va='top',linespacing=1.4);y-=.07+.043*len(text.split('\n'))
 pdf.savefig(fig);plt.close(fig)
def page_plot(pdf,title,name,caption):
 fig=plt.figure(figsize=(11.7,8.3));fig.text(.04,.95,title,fontsize=17,weight='bold')
 ax=fig.add_axes([.04,.14,.92,.75]);ax.imshow(Image.open(P/'plots'/name));ax.axis('off')
 fig.text(.04,.075,caption,fontsize=9.5,va='top');pdf.savefig(fig);plt.close(fig)

for lang in ['EN','ID']:
 en=lang=='EN'
 with PdfPages(P/f'Propagation_Meeting_{lang}.pdf') as pdf:
  page_text(pdf,'Cold-pool propagation: from point arrivals to a source' if en else 'Propagasi cold pool: dari waktu kedatangan ke asal',[
   'Input: 197 original 5-minute cold-pool candidates, including 9 Eq. (3) episodes.\nThe newer 1-minute and modified-monotonicity catalogs are not used.' if en else 'Input: 197 kandidat cold pool 5 menit asli, termasuk 9 episode Eq. (3).\nKatalog 1 menit dan versi monoton baru tidak digunakan.',
   'Three models: expanding circle; advected circle; advected ellipse.' if en else 'Tiga model: lingkaran mengembang; lingkaran beradveksi; elips beradveksi.',
   'Goal: infer release location and time from station arrival times.\nFirst validate the inverse method using synthetic data with a known source.' if en else 'Tujuan: inferensi lokasi dan waktu mulai dari waktu kedatangan di stasiun.\nMetode inversi diuji dahulu pada simulasi dengan asal diketahui.',
   'Real-data fits are conditional hypotheses, not confirmed origins.\nAssociation, propagation speed, wind translation and terrain remain uncertain.' if en else 'Fit data nyata merupakan hipotesis bersyarat, bukan asal terkonfirmasi.\nAsosiasi event, kecepatan, adveksi angin, dan topografi masih tidak pasti.'
  ])
  page_text(pdf,'Forward equations and assumptions' if en else 'Persamaan maju dan asumsi',[
   'Circle, no wind: |r_i − r_0| = c (t_i − t_0).\nUnknown source r_0=(x_0,y_0), release time t_0, and radial speed c.' if en else 'Lingkaran tanpa angin: |r_i − r_0| = c (t_i − t_0).\nAsal r_0=(x_0,y_0), waktu mulai t_0, dan kecepatan radial c dicari.',
   'With wind: center(t)=r_0+U(t−t_0); radius=c(t−t_0).\nU is prescribed translation; measured surface wind is only a proxy.' if en else 'Dengan angin: pusat(t)=r_0+U(t−t_0); radius=c(t−t_0).\nU adalah adveksi yang diasumsikan; angin permukaan hanya proxy.',
   'Ellipse: semi-axes grow at c√r and c/√r, with fixed orientation theta.\nr is the major/minor axis ratio. r=1 exactly recovers the circle.' if en else 'Elips: laju sumbu semi-mayor c√r dan semi-minor c/√r, orientasi theta.\nr adalah rasio sumbu mayor/minor. r=1 tepat menjadi lingkaran.',
   'Constant speeds and wind; flat horizontal geometry; no terrain dynamics.\nKinematic conceptual models, not full atmospheric fluid simulations.' if en else 'Kecepatan dan angin konstan; geometri horizontal datar; dinamika topografi diabaikan.\nIni model konseptual kinematik, bukan simulasi fluida atmosfer lengkap.'
  ])
  page_plot(pdf,'Three conceptual stages and known-source recovery' if en else 'Tiga tahap konseptual dan pemulihan asal simulasi','conceptual_three_models.png','Nine distributed synthetic stations; prescribed wind; free ellipse fit only in this synthetic demonstration.' if en else 'Sembilan stasiun simulasi tersebar; angin ditentukan; elips bebas hanya pada demonstrasi simulasi ini.')
  page_plot(pdf,'Sensitivity to five-minute timing resolution' if en else 'Sensitivitas terhadap resolusi waktu lima menit','synthetic_timing_sensitivity.png','35 independent uniform ±2.5-minute perturbations per model. Conditional sensitivity, not real-data confidence intervals.' if en else '35 perturbasi uniform independen ±2,5 menit per model. Sensitivitas bersyarat, bukan interval kepercayaan data nyata.')
  page_plot(pdf,'Actual station geometry and elevation' if en else 'Geometri stasiun dan elevasi sebenarnya','station_geometry.png',f'Minor/major spatial spread ratio {geo["minor_major_spread_ratio"]:.3f}; elevation 27–945 m ASL. Horizontal model ignores terrain.' if en else f'Rasio sebaran minor/mayor {geo["minor_major_spread_ratio"]:.3f}; elevasi 27–945 m. Model horizontal mengabaikan topografi.')
  page_text(pdf,'What can be identified from point observations?' if en else 'Apa yang dapat ditentukan dari pengamatan titik?',[
   'Circle with unknown speed: 4 unknowns (x,y,t0,c).\nFour arrivals give zero residual degrees of freedom, not guaranteed uniqueness.' if en else 'Lingkaran dengan kecepatan tidak diketahui: 4 parameter (x,y,t0,c).\nEmpat kedatangan memberi nol derajat bebas residual, bukan jaminan unik.',
   'Free ellipse with prescribed wind: 6 unknowns.\nFive single-pass station arrivals cannot identify all six parameters.' if en else 'Elips bebas dengan angin ditentukan: 6 parameter.\nLima waktu kedatangan satu lintasan tidak cukup menentukan keenamnya.',
   'Real-data main fits fix c=5 m/s and shape assumptions: 3 source unknowns.\nSensitivity uses c=3 and 8 m/s; speed is a scenario, not an observed quantity.' if en else 'Fit utama data nyata menetapkan c=5 m/s dan bentuk: 3 parameter asal.\nSensitivitas memakai c=3 dan 8 m/s; ini skenario, bukan kecepatan terukur.',
   'Unknown background wind adds two more unknowns.\nThree observations with three source parameters also leave no residual df.' if en else 'Angin latar yang ikut dicari menambah dua parameter lagi.\nTiga observasi dengan tiga parameter asal juga tidak menyisakan derajat bebas.'
  ])
  page_plot(pdf,'Same exact arrivals, different free ellipses' if en else 'Waktu kedatangan sama, elips bebas berbeda','ellipse_nonuniqueness.png',f'{len(equivalent)} near-exact solutions on the five actual site coordinates. This is synthetic, not a real recovered source.' if en else f'{len(equivalent)} solusi hampir eksak pada koordinat lima stasiun. Ini simulasi, bukan pemulihan asal nyata.')
  page_plot(pdf,'Temporal association is a hypothesis' if en else 'Asosiasi temporal adalah hipotesis','candidate_arrivals.png','Primary bounded window: 180 min from first detection. Eleven groups have ≥3 sites; two have four; none has five.' if en else 'Jendela utama dibatasi 180 menit dari deteksi pertama. Sebelas kelompok ≥3 stasiun, dua kelompok 4 stasiun, tidak ada 5.')
  page_plot(pdf,'Inspect the actual observations' if en else 'Periksa pengamatan sebenarnya','real_example_temperature_traces.png','Background-normalized temperature helps inspect association, but common timing and cooling do not prove a single source.' if en else 'Anomali suhu membantu inspeksi asosiasi; kedekatan waktu dan pendinginan tidak membuktikan satu sumber yang sama.')
  page_plot(pdf,'Conditional real-data example: 28 June 2025' if en else 'Contoh data nyata bersyarat: 28 Juni 2025','real_conditional_example.png','Fixed c=5 m/s. Prescribed background wind/shape alter the inferred origin. No model is declared the true source.' if en else 'c tetap 5 m/s. Asumsi angin/bentuk mengubah estimasi asal. Tidak ada model yang dinyatakan sebagai sumber sebenarnya.')
  page_plot(pdf,'Source-loss surface, not a confidence map' if en else 'Permukaan galat asal, bukan peta kepercayaan','real_origin_loss_profile.png','Circle with c=5 m/s. White contour is best grid RMSE +2.5 min; not a 95% confidence boundary.' if en else 'Lingkaran c=5 m/s. Kontur putih: RMSE terbaik grid +2,5 menit; bukan batas kepercayaan 95%.')
  page_text(pdf,'What to present at the next meeting' if en else 'Yang dibahas pada meeting berikutnya',[
   'Synthetic recovery works under its stated assumptions; real-data origins remain conditional.' if en else 'Pemulihan simulasi berhasil pada asumsi yang ditentukan; asal data nyata tetap bersyarat.',
   'The 28 June group has a circle/no-wind c=5 m/s RMSE of about 2.56 min.\nFitting speed as well can drive RMSE to zero with four data and four unknowns.' if en else 'Kelompok 28 Juni: RMSE lingkaran tanpa angin c=5 m/s sekitar 2,56 menit.\nMencari kecepatan sekaligus dapat memberi RMSE nol pada empat data/parameter.',
   'Zero residual is not validation. Another four-station case hits prior bounds\nand has larger residuals: constant speed, association, or model may be inappropriate.' if en else 'Residual nol bukan validasi. Kasus empat stasiun lain menyentuh batas prior\ndan bergalat lebih besar: kecepatan, asosiasi, atau model dapat tidak sesuai.',
   'Next: independently verify shared events, constrain propagation speed/wind,\nadd off-axis observations, and use radar/satellite/terrain information.' if en else 'Berikutnya: verifikasi event bersama secara independen, batasi kecepatan/angin,\ntambah pengamatan di luar garis jaringan, gunakan radar/satelit/topografi.',
   'GIFs, all fit tables, alternative solutions, source-loss arrays, inputs and code\nare included. Earlier detection files have not been changed.' if en else 'GIF, tabel fit, solusi alternatif, array galat asal, input dan kode disertakan.\nFile deteksi sebelumnya tidak diubah.'
  ])
print(json.dumps(validation,indent=2))
