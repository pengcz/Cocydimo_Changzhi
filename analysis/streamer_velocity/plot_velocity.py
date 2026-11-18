"""Estimate leading-channel axial velocity from Silo conductivity fronts.

Run extract output > fronts.csv, read_times output > times.csv, then this script.
Front = greatest non-ghost cell-center z with sigma_tot >= threshold.
This is a leading-front estimate, not individual head trajectories or 3D speed.
Frame zero has no conductivity and remains undefined. No smoothing is applied.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base = Path(__file__).resolve().parent
fronts = np.genfromtxt(base/'fronts.csv', delimiter=',', names=True)
times = np.genfromtxt(base/'times.csv', delimiter=',', names=True)
assert len(fronts)==len(times)==52
assert np.array_equal(fronts['frame'], times['frame'])
t = times['time_s']
assert np.all(np.diff(t)>0)
z = np.column_stack([fronts[n] for n in ['z_1e4_m','z_1e3_m','z_1e2_m']])
velocity = np.diff(z,axis=0)/np.diff(t)[:,None]
midpoint = .5*(t[1:]+t[:-1])
np.savetxt(base/'head_positions.csv', np.column_stack([fronts['frame'],t,z,fronts['dz_m']]), delimiter=',',header='frame,time_s,z_sigma_ge_1e-4_m,z_sigma_ge_1e-3_m,z_sigma_ge_1e-2_m,cell_dz_at_1e-3_front_m',comments='')
np.savetxt(base/'head_velocity.csv', np.column_stack([fronts['frame'][:-1],fronts['frame'][1:],midpoint,velocity]),delimiter=',',header='frame_start,frame_end,midpoint_time_s,axial_velocity_1e-4_m_s,axial_velocity_1e-3_m_s,axial_velocity_1e-2_m_s',comments='')
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,(ax,ap)=plt.subplots(2,1,figsize=(9,7),sharex=True,gridspec_kw={'height_ratios':[2,1]})
for k,cut in enumerate(['10⁻⁴','10⁻³','10⁻²']):
    ax.plot(midpoint*1e9,velocity[:,k]/1e6,lw=1.4 if k==1 else .9,alpha=1 if k==1 else .6,label=f'σ ≥ {cut} S/m',linestyle=['--','-',':'][k])
ap.plot(t*1e9,z[:,1],color='tab:blue',lw=1.5)
ax.set(title='3D streamer: leading-channel axial velocity, frames 0–51',ylabel='Axial velocity (10⁶ m/s)')
ax.legend(frameon=False)
ap.set(xlabel='Time (ns)',ylabel='Front position z (m)',xlim=(0,t[-1]*1e9))
for a in [ax,ap]: a.grid(alpha=.25)
fig.text(.11,.015,'Front tracked from conductivity cell centers; interval velocities, no smoothing.\nFrame 0 has no channel; its velocity is undefined. Branch identities are not tracked.',fontsize=9,color='0.3')
fig.tight_layout(rect=(0,.065,1,1))
for ext in ['png','pdf']: fig.savefig(base/f'streamer_head_velocity.{ext}',dpi=200)
valid = np.isfinite(velocity[:,1])
print(f'Time range: {t[0]*1e9:.3f}–{t[-1]*1e9:.3f} ns')
print(f'Valid velocity intervals: {valid.sum()}')
print(f'Axial velocity range: {velocity[valid,1].min():.6g}–{velocity[valid,1].max():.6g} m/s')
print(f'Mean front velocity: {(z[-1,1]-z[1,1])/(t[-1]-t[1]):.6g} m/s')
print(f'Max threshold position spread: {np.nanmax(np.ptp(z,axis=1)):.6g} m')
