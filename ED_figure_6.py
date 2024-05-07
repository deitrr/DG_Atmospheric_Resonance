import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pathlib
import warnings
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

warnings.filterwarnings("ignore", category=DeprecationWarning)

#change this path to match the location of the sim_main download from the
#data repository
parent = '../sims_main/'


simnames = [
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_16hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_18hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_20hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_23hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_25hr_branch',
           ]

rotpers = np.array([0.6667, 0.75, 0.8333,0.875, 0.8958, 0.9167, 0.9271,
                    0.9375, 1.0, 1.04167])
label = ['16 hour', '18 hour', '20 hour', '21 hour', '21.5 hour', '22 hour',
                '22.25 hour', '22.5 hour', '23 hour', '24 hour', '25 hour']
xlabels=['Deep convective mass\nflux (10$^{-3}$ kg m$^{-2}$ s$^{-1}$)',
        'Deep convective heating\nrate (K day$^{-1}$)', 'Horizontal divergence\n(s$^{-1}$)']
blues = plt.cm.Blues
reds = plt.cm.Reds
colors = [blues(0.1),blues(0.25),blues(0.5),blues(0.75),blues(0.88),blues(1.0),
                        'k',reds(1.),reds(0.75),reds(0.5),reds(0.25)]
xmax = [6.6,3.5,2.5]
xmin = [-0.5,-6.5,-2]
scale = [1e-3,1,1e-6]

ilons = [54, 0, 18] #longitudes near terminators and noon
ilat = 22  #equatorial latitude
titles = ['Sunrise', 'Noon', 'Sunset']

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

    if len(field_mean) > len(p):
      field_mean = 0.5*(field_mean[:-1] + field_mean[1:])

    if ifile in (1,):
      field_mean *= 86400

    for jlon in np.arange(len(ilons)):

      ax = axes[jlon][ifile]

      ilon = ilons[jlon]
      ax.plot(field_mean[:,ilat,ilon]/scale[ifile],p,c=colors[i],label=label[i])

for jlon in np.arange(len(ilons)):
  for ifile in np.arange(len(ifiles)):
    ax = axes[jlon][ifile]
    #if ifile < 2:
    ax.set(xlim=(xmin[ifile],xmax[ifile]))
    if ifile == 0:
      ax.set(ylabel='Pressure (hPa)')
    ax.set_ylim(80,1000)

    ax.invert_yaxis()
    ax.set_title(titles[jlon],fontsize=8)
    if jlon == 2:
      ax.legend(loc='best',ncols=1,fontsize=4)
      ax.set(xlabel=xlabels[ifile])

plt.tight_layout()
plt.savefig('figures/ED_figure_6.eps')
plt.close()
