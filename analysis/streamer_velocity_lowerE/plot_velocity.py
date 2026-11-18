"""Leading-front axial velocity; exclude intervals after front resolution is lost."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
b=Path(__file__).resolve().parent
f=np.genfromtxt(b/'fronts.csv',delimiter=',',names=True)
t=np.genfromtxt(b/'times.csv',delimiter=',',names=True)['time_s']
assert len(f)==len(t)==52 and np.all(np.diff(t)>0)
z=np.column_stack([f[n] for n in ['z_1e4_m','z_1e3_m','z_1e2_m']])
v=np.diff(z,axis=0)/np.diff(t)[:,None]; tm=(t[1:]+t[:-1])/2
# Transparent quality screen: threshold spread <= 1 mm and local dz <= 1 mm.
good=np.all(np.isfinite(z),axis=1)&(np.ptp(z,axis=1)<=.001)&(f['dz_m']<=.001)
valid=good[:-1]&good[1:]
vg=np.where(valid,v[:,1],np.nan)
np.savetxt(b/'head_positions.csv',np.column_stack([f['frame'],t,z,f['dz_m'],good]),delimiter=',',header='frame,time_s,z_sigma_ge_1e-4_m,z_sigma_ge_1e-3_m,z_sigma_ge_1e-2_m,cell_dz_m,front_resolved',comments='')
np.savetxt(b/'head_velocity.csv',np.column_stack([f['frame'][:-1],f['frame'][1:],tm,v,valid,vg]),delimiter=',',header='frame_start,frame_end,midpoint_time_s,raw_front_velocity_1e-4_m_s,raw_front_velocity_1e-3_m_s,raw_front_velocity_1e-2_m_s,reliable_interval,estimated_head_axial_velocity_m_s',comments='')
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,(ax,ap)=plt.subplots(2,1,figsize=(9,7),sharex=True,gridspec_kw={'height_ratios':[2,1]})
for k,c in enumerate(['10⁻⁴','10⁻³','10⁻²']):
 ax.plot(tm*1e9,np.where(valid,v[:,k],np.nan)/1e6,lw=1.5,ls=['--','-',':'][k],label=f'σ ≥ {c} S/m')
 ap.plot(t*1e9,z[:,k],lw=1.3,ls=['--','-',':'][k])
last=np.flatnonzero(good)[-1]; cutoff=t[last]*1e9
for a in [ax,ap]:
 a.axvspan(cutoff,t[-1]*1e9,color='0.9',zorder=0)
 a.axvline(cutoff,color='0.5',ls='--',lw=1)
 a.grid(alpha=.25)
ax.text(.73,.62,'Head velocity unresolved\nCoarse mesh and fading\nconductivity front',transform=ax.transAxes,ha='center',color='0.35')
ax.set(title='LowerE 3D streamer: estimated axial head velocity, frames 0–51',ylabel='Axial velocity (10⁶ m/s)',ylim=(0,None))
ax.legend(frameon=False,loc='upper right')
ap.set(xlabel='Time (ns)',ylabel='Conductivity front z (m)',xlim=(0,t[-1]*1e9))
fig.text(.11,.015,'Velocity = Δz / Δt over saved-frame intervals (9 ns). No smoothing.\nFrame 0 has no channel; after frame 26 the conductivity front cannot reliably track a head.',fontsize=9,color='0.3')
fig.tight_layout(rect=(0,.065,1,1))
for ext in ['png','pdf']:fig.savefig(b/f'streamer_head_velocity.{ext}',dpi=200)
print('Last resolved frame:',last,'time ns:',cutoff)
print('Reliable velocity range m/s:',np.nanmin(vg),np.nanmax(vg))
print('Intervals:',valid.sum())
