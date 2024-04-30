import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pdb
import pathlib
import warnings
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import scipy.interpolate as sint
#from scipy.integrate import simps
#import pyshtools as sh
#pdb.set_trace()
warnings.filterwarnings("ignore", category=DeprecationWarning)

parent = '../sims_main/'


simnames = [
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_16hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_18hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_20hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_23hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_25hr_branch',
           ]

rotpers = np.array([0.6667, 0.75, 0.8333,0.875, 0.8958, 0.9167, 0.9375, 1.0, 1.04167])
label = ['16 hour', '18 hour', '20 hour', '21 hour', '21.5 hour', '22 hour', '22.5 hour', '23 hour', '24 hour', '25 hour']
xlabels=['Deep convective mass\nflux (10$^{-3}$ kg m$^{-2}$ s$^{-1}$)', 'Deep convective heating\nrate (K day$^{-1}$)', 'Horizontal divergence\n(s$^{-1}$)']
#colors = ['lightskyblue','skyblue','blue','mediumblue','navy','darkred','red','coral','lightcoral']
blues = plt.cm.Blues
reds = plt.cm.Reds
colors = [blues(0.1),blues(0.25),blues(0.5),blues(0.75),blues(1.0),'k',reds(1.),reds(0.75),reds(0.5),reds(0.25)]
xmax = [6.6,3.5,2.5]
xmin = [-0.5,-6.5,-2]
scale = [1e-3,1,1e-6]

ilons = [54, 0, 18]
ilat = 22
titles = ['Sunrise', 'Noon', 'Sunset']
#cmap = plt.cm.RdBu_r

fig, axes = plt.subplots(ncols=3,nrows=3,figsize=(7.5,6))

for i in np.arange(len(simnames)):
  print(simnames[i])
  ifiles = [
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_CMFMCDZM_save.npz',
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_ZMDT_save.npz',
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_DIVV_save.npz',
]

  for ifile in np.arange(len(ifiles)):
    arc = np.load(ifiles[ifile])
    lon = arc['lon']
    lat = arc['lat']
    p = arc['p']
    rotrate = arc['rotrate']
    field_anom_mean = arc['field_anom_mean']
    field_mean = arc['field_mean']
    clm_mean = arc['clm_mean']
    clm_anom = arc['clm_anom']
    lmax = arc['lmax']
#    if ifile == 2:
#      pdb.set_trace()
    if len(field_mean) > len(p):
      field_mean = 0.5*(field_mean[:-1] + field_mean[1:])

    if ifile in (1,):
      field_mean *= 86400

    for jlon in np.arange(len(ilons)):
#      if ifile == 1:
#        axes[jlon][0].plot(field_mean[:,ilat,ilon],p,c=colors[i],linestyle='--')

      ax = axes[jlon][ifile]

      ilon = ilons[jlon]
#      ax.semilogy(field_mean[:,ilat,ilon],p,c=cmap((rotpers[i]-rotpers[0])/(rotpers[-1]-rotpers[0])),label=label[i])
#      pdb.set_trace()
      ax.plot(field_mean[:,ilat,ilon]/scale[ifile],p,c=colors[i],label=label[i])

for jlon in np.arange(len(ilons)):
  for ifile in np.arange(len(ifiles)):
    ax = axes[jlon][ifile]
    #if ifile < 2:
    ax.set(xlim=(xmin[ifile],xmax[ifile]))
    if ifile == 0:
      ax.set(ylabel='Pressure (hPa)')
    ax.set_ylim(80,1000)

    #else:
    #  ax.set(ylabel='Pressure (hPa)') #,xlim=(200,310))
    #  ax.set_ylim(1,1000)

    ax.invert_yaxis()
    ax.set_title(titles[jlon],fontsize=8)
    if jlon == 2:
      ax.legend(loc='best',ncols=1,fontsize=5)
      ax.set(xlabel=xlabels[ifile])

plt.tight_layout()
plt.savefig('figure_11.pdf')
plt.close()

