import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import scipy.interpolate as si
import matplotlib.gridspec as gridspec
import pathlib
import pdb

#this just suppresses a deprecation warning from netCDF4
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

main_sim_path = '../sims_main/'

#list of simulations that will go into figure 2
simnames =[
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_16hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_18hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_20hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21-5hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-5hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_23hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_25hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_20hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_21hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_21-5hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_21-75hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_22hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_23hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_24hr',

           ]

rotper = np.array([16,18,20,21,21.5,22,22.25,22.5,23,24,25,20,21,21.5,21.75,22,23,24])
lines = np.array([0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1])

#create the arrays
Ts = np.zeros(len(simnames))
q500 = np.zeros(len(simnames))
icwp = np.zeros_like(Ts)
lcwp = np.zeros_like(Ts)
panom = np.zeros_like(Ts)
lhflx = np.zeros_like(Ts)
qflx = np.zeros_like(Ts)
tmq = np.zeros_like(Ts)
T500 = np.zeros_like(Ts)
ASR = np.zeros_like(Ts)
prect = np.zeros_like(Ts)
OLR = np.zeros_like(Ts)
ASR_trop = np.zeros_like(Ts)
imag_dps = np.zeros_like(Ts)

for i in np.arange(len(simnames)):
  print(simnames[i])

  #Global mean file----------------------
  simpath = main_sim_path + simnames[i]
  file = simpath + "/merged_hist/" + simnames[i] + ".cam.h0.globmean_0031_0060.nc"

  data = nc.Dataset(file,'r')
#  import pdb; pdb.set_trace()
  #read in data from global mean
  p = data['lev'][:].squeeze()
  Ts[i] = data['TS'][:].squeeze()
  q = data['Q'][:].squeeze()
  icwp[i] = data['TGCLDIWP'][:].squeeze()/1000.0
  lcwp[i] = data['TGCLDLWP'][:].squeeze()/1000.0
  lhflx[i] = data['LHFLX'][:].squeeze()
  qflx[i] = data['QFLX'][:].squeeze()
  tmq[i] = data['TMQ'][:].squeeze()
  T = data['T'][:].squeeze()
  fus = data['FUS'][:].squeeze()
  fds = data['FDS'][:].squeeze()
  prect[i] = data['PRECT'][:].squeeze()
  ful = data['FUL'][:].squeeze()
  fdl = data['FDL'][:].squeeze()

#  pdb.set_trace()
  #interpolate these to 500 hPa
  q500[i] = si.interp1d(p,q)(500.0)
  T500[i] = si.interp1d(p,T)(500.0)
  ASR[i] = fds[0] - fus[0]
  OLR[i] = ful[0] - fdl[0]

  ASR_trop[i] = si.interp1d(data['ilev'][:].squeeze(),fds-fus)(200.0)

  #pressure anomaly
  if simnames[i] == 'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_23hr' or simnames[i] == 'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr':
    file = simpath + "_branch2/merged_hist/" + simnames[i] + "_branch2_p_anom_save.npz"
  else:
    file = simpath + "_branch/merged_hist/" + simnames[i] + "_branch_p_anom_save.npz"
  arc = np.load(file)
  if 'field_anom_mean' in arc:
      panom[i] = np.max(np.abs(arc['field_anom_mean']))
  else:
      # in some of the files I used the name ps_anom_mean instead of field_anom_mean
      panom[i] = np.max(np.abs(arc['ps_anom_mean']))
  imag_dps[i] = -1*arc['clm'][1,2,2]*np.sqrt(2*15/8)

print(panom)
print(rotper)
print(ASR)
print(ASR_trop)

#this will be the synthesis plot
plt.rcParams.update({'font.size':6})
fig = plt.figure(figsize=(7.5,5))

#Set up the plot arrangement
#outer_grid = gridspec.GridSpec(2,1,wspace=0.1,hspace=0.3,left=0.065,right=0.96,bottom=0.09,top=0.96,height_ratios=(2,2))
#row1 = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_grid[0],wspace=0.5,hspace=0.1,width_ratios=(1,1,1,1))
#row2 = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_grid[1],wspace=0.5,hspace=0.1,width_ratios=(1,1,1,1))
#row3 = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_grid[2],wspace=0.5,hspace=0.1,width_ratios=(1,1,1,1))
#row4 = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_grid[3],wspace=0.5,hspace=0.1,width_ratios=(1,1,1,1))
#row5 = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_grid[4],wspace=0.5,hspace=0.1,width_ratios=(1,1,1,1))


bottom = 0.7
height = 0.27
width = 0.26
#top row
#ax = fig.add_subplot(row1[0])
ax = fig.add_axes([0.05,bottom,width,height])
ax.plot(rotper[lines==0],Ts[lines==0],color='k',marker='.',linestyle='-',label='30 ppm CH$_4$',lw=2)
ax.plot(rotper[lines==1],Ts[lines==1],color='0.5',marker='.',linestyle='-',label='3000 ppm CH$_4$',lw=2)
#ax.plot(rotper[6:],Ts[6:],color='k',marker='.',linestyle='-',label='Low CH$_4$')
ax.set(ylabel='Surface temperature (K)')
ax.tick_params(direction='in')
ax.set_xticks([16,18,20,22,24])
ax.legend(loc='upper left',fontsize=5)

#ax = fig.add_subplot(row1[1])
ax = fig.add_axes([0.37,bottom,width,height])
ax.plot(rotper[lines==0],panom[lines==0],color='k',marker='.',linestyle='-',lw=2)
ax.plot(rotper[lines==1],panom[lines==1],color='0.5',marker='.',linestyle='-',lw=2)
ax.set(ylabel='Peak pressure anomaly (Pa)')
ax.tick_params(direction='in')
ax.set_xticks([16,18,20,22,24])

ax = fig.add_axes([0.72,bottom,width,height])
ax.plot(rotper[lines==0],imag_dps[lines==0],color='k',marker='.',linestyle='-',lw=2)
ax.plot(rotper[lines==1],imag_dps[lines==1],color='0.5',marker='.',linestyle='-',lw=2)
ax.set(ylabel='Semidiurnal pressure anomaly (Pa)\n(imaginary component)')
ax.tick_params(direction='in')
ax.set_xticks([16,18,20,22,24])

#next row
bottom = 0.38
width= 0.23

dy = 1.2*np.max((np.max(q500[lines==0]) - np.min(q500[lines==0]),
             np.max(q500[lines==1]) - np.min(q500[lines==1])))*1000
#ax = fig.add_subplot(row1[2])
ax = fig.add_axes([0.07,bottom,width,height])
ax.plot(rotper[lines==0],q500[lines==0]*1000,color='k',marker='.',linestyle='-',lw=2)
ax2 = ax.twinx()
ax2.plot(rotper[lines==1],q500[lines==1]*1000,color='0.5',marker='.',linestyle='-',lw=2)
#ax.plot(rotper[6:],q500[6:],color='k',marker='.',linestyle='-')
ax.set(ylabel='Specific Humidity \nat 500 hPa (10$^{-3}$ kg kg$^{-1}$)')
ax.tick_params(direction='in')
#ax.set_xticks([1e-8,1e-6,1e-4,1e-2])
ax.set_xticks([16,18,20,22,24])
ycen = 0.5*(ax.get_ylim()[1] + ax.get_ylim()[0])
ax.set_ylim(ycen-0.5*dy,ycen+0.5*dy)
ycen = 0.5*(ax2.get_ylim()[1] + ax2.get_ylim()[0])
ax2.set_ylim(ycen-0.5*dy,ycen+0.5*dy)
ax2.tick_params(axis='y', labelcolor='0.4')

dy = 1.2*np.max((np.max(tmq[lines==0]) - np.min(tmq[lines==0]),
             np.max(tmq[lines==1]) - np.min(tmq[lines==1])))
#ax = fig.add_subplot(row1[3])
ax = fig.add_axes([0.41,bottom,width,height])
ax.plot(rotper[lines==0],tmq[lines==0],color='k',marker='.',linestyle='-',lw=2)
ax2 = ax.twinx()
ax2.plot(rotper[lines==1],tmq[lines==1],color='0.5',marker='.',linestyle='-',lw=2)
ax.set(ylabel='Total precipitable water\n(kg m$^{-2}$)')
ax.tick_params(direction='in')
ax.set_xticks([16,18,20,22,24])
ycen = 0.5*(ax.get_ylim()[1] + ax.get_ylim()[0])
ax.set_ylim(ycen-0.5*dy,ycen+0.5*dy)
ycen = 0.5*(ax2.get_ylim()[1] + ax2.get_ylim()[0])
ax2.set_ylim(ycen-0.5*dy,ycen+0.5*dy)
ax2.tick_params(axis='y', labelcolor='0.4')

dy = 1.2*np.max((np.max(lhflx[lines==0]) - np.min(lhflx[lines==0]),
             np.max(lhflx[lines==1]) - np.min(lhflx[lines==1])))
#ax = fig.add_subplot(row2[2])
ax = fig.add_axes([0.73,bottom,width,height])
ax.plot(rotper[lines==0],lhflx[lines==0],color='k',marker='.',linestyle='-',lw=2)
ax2=ax.twinx()
ax2.plot(rotper[lines==1],lhflx[lines==1],color='0.5',marker='.',linestyle='-',lw=2)
ax.set(ylabel='Surface latent heat\nflux (W m$^{-2}$)') #,xlabel='Length of day (hours)')
ax.tick_params(direction='in')
ax.set_xticks([16,18,20,22,24])
ycen = 0.5*(ax.get_ylim()[1] + ax.get_ylim()[0])
ax.set_ylim(ycen-0.5*dy,ycen+0.5*dy)
ycen = 0.5*(ax2.get_ylim()[1] + ax2.get_ylim()[0])
ax2.set_ylim(ycen-0.5*dy,ycen+0.5*dy)
ax2.tick_params(axis='y', labelcolor='0.4')


#3rd row
dy_cwp = 1.2*np.max((np.abs(np.max(icwp[lines==0])-np.min(icwp[lines==0])),np.abs(np.max(lcwp[lines==0])-np.min(lcwp[lines==0]))))

bottom = 0.06
width= 0.23
#ax = fig.add_subplot(row2[0])
ax = fig.add_axes([0.07,bottom,width,height])
ax.plot(rotper[lines==0],lcwp[lines==0],color='k',marker='.',linestyle='-',lw=2)
ax2 = ax.twinx()
ax2.plot(rotper[lines==1],lcwp[lines==1],color='0.5',marker='.',linestyle='-',lw=2)
#ax.plot(rotper[6:],lcwp[6:],color='k',marker='.',linestyle='-')
ax.set(ylabel='Liquid water path (kg m$^{-2}$)',xlabel='Length of day (hours)')
ax.tick_params(direction='in')
#ax.set_xticks([1e-8,1e-6,1e-4,1e-2])
ylims = ax.get_ylim()
ax.set_ylim((0.5*(ylims[0]+ylims[1]-dy_cwp),0.5*(ylims[0]+ylims[1]+dy_cwp)))
ylims = ax2.get_ylim()
ax2.set_ylim((0.5*(ylims[0]+ylims[1]-dy_cwp),0.5*(ylims[0]+ylims[1]+dy_cwp)))
ax.set_xticks([16,18,20,22,24])
ax2.tick_params(axis='y', labelcolor='0.4')

#ax = fig.add_subplot(row2[1])
ax = fig.add_axes([0.42,bottom,width,height])
ax.plot(rotper[lines==0],icwp[lines==0],color='k',marker='.',linestyle='-',lw=2)
ax.plot(rotper[lines==1],icwp[lines==1],color='0.5',marker='.',linestyle='-',lw=2)
#ax.plot(rotper[6:],icwp[6:],color='k',marker='.',linestyle='-')
ax.set(ylabel='Ice water path (kg m$^{-2}$)',xlabel='Length of day (hours)')
ax.tick_params(direction='in')
#ax.set_xticks([1e-8,1e-6,1e-4,1e-2])
ylims = ax.get_ylim()
ax.set_ylim((0.5*(ylims[0]+ylims[1]-dy_cwp),0.5*(ylims[0]+ylims[1]+dy_cwp)))
ax.set_xticks([16,18,20,22,24])

dy = 1.2*np.max((np.max(ASR[lines==0]) - np.min(ASR[lines==0]),
             np.max(ASR[lines==1]) - np.min(ASR[lines==1])))
#ax = fig.add_subplot(row2[3])
ax = fig.add_axes([0.72,bottom,width,height])
#ax.plot(rotper,qflx*1e5,color='k',marker='.',linestyle='-',lw=2)
#ax.set(ylabel='Surface water flux\n(10$^{-5}$ kg m$^{-2}$ s$^{-1}$)',xlabel='Rotation period (hours)')
ax.plot(rotper[lines==0],ASR[lines==0],color='k',marker='.',linestyle='-',lw=2)
ax2 = ax.twinx()
ax2.plot(rotper[lines==1],ASR[lines==1],color='0.5',marker='.',linestyle='-',lw=2)
ax.set(ylabel='Absorbed Solar radiation\n(W m$^{-2}$)',xlabel='Length of day (hours)')
ax.tick_params(direction='in')
ax.set_xticks([16,18,20,22,24])
ycen = 0.5*(ax.get_ylim()[1] + ax.get_ylim()[0])
ax.set_ylim(ycen-0.5*dy,ycen+0.5*dy)
ycen = 0.5*(ax2.get_ylim()[1] + ax2.get_ylim()[0])
ax2.set_ylim(ycen-0.5*dy,ycen+0.5*dy)
ax2.tick_params(axis='y', labelcolor='0.4')

plt.savefig('figure_1.pdf')
plt.close()
