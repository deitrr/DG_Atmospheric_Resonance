import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pdb
import pathlib
import warnings
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.interpolate as sint
import pyshtools as sh
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore", category=DeprecationWarning)

#change this path to match the location of the sim_main download from the
#data repository
main_path = '../sims_main/'

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
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_25hr_branch',
]

label = ['16 hour', '18 hour', '20 hour', '21 hour', '21.5 hour', '22 hour',
            '22.25 hour', '22.5 hour', '23 hour', '24 hour', '25 hour']
rotpers = np.array([0.6667, 0.75, 0.8333, 0.875, 0.8958, 0.9167, 0.9271, 0.9375,
                    0.9583, 1.0, 1.04167])

clabels = ['Surface pressure anomaly (Pa)',
            'Surface pressure anomaly (Pa)\ndiurnal mode',
            'Surface pressure anomaly (Pa)\nsemidiurnal mode']
cmap = plt.cm.RdBu_r
cscale = [400,50,400]

fig = plt.figure(figsize=(7.5,9.5))

outer_grid = gridspec.GridSpec(1,1,wspace=0.2,hspace=0.1,left=0.05,right=0.98,
                                        bottom=0.03,top=0.94,height_ratios=(1,))
inner_grid = gridspec.GridSpecFromSubplotSpec(12,3,subplot_spec=outer_grid[0],
                wspace=0.05,hspace=0.15,height_ratios=(1,6,6,6,6,6,6,6,6,6,6,6))

for i in np.arange(len(simnames)):
    p_anom_file = main_path + simnames[i]+'/merged_hist/'+simnames[i]+'_p_anom_save.npz'

    arc = np.load(p_anom_file)
    lon = arc['lon']
    lat = arc['lat']
    ps = arc['ps']
    rotrate = arc['rotrate']
    if 'field_anom_mean' in arc:
        ps_anom_mean = arc['field_anom_mean']
    else:
        ps_anom_mean = arc['ps_anom_mean']
    lmax = arc['lmax']
    clm = arc['clm']
    lon2d, lat2d = np.meshgrid(lon,lat)

  ax = fig.add_subplot(inner_grid[3*(i+1)])
  m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection='cea')
  m.drawparallels([-45,0,45],labels = [True,False,False,False], fontsize=8)
  m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)

  c = m.pcolormesh(lon2d, lat2d, ps_anom_mean, cmap=cmap,rasterized=True,
                    latlon='True',vmax=cscale[0],vmin=-1*cscale[0])
  ax.yaxis.set_label_coords(-0.1,0.5)
  ax.text(0.02,0.85,label[i],rotation=0,transform=ax.transAxes,fontsize=10,
                    color='k',fontweight='bold')

  if i == len(simnames)-1:
    xlim = ax.get_xlim()
    dxlim = xlim[1] - xlim[0]
    ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
    ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'])
    ax.tick_params(direction='in')

  if i == 0:
    cax = fig.add_subplot(inner_grid[0])
    cbar = plt.colorbar(c,cax=cax,orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cax.xaxis.set_label_position('top')
    cbar.set_label(clabels[0])

  for imode in [1,2]:
    clm_copy = clm.copy()
    clm_copy[:,:imode,:] = 0.0
    clm_copy[:,imode+1:,:] = 0.0
    clm_copy[:,:,:imode] = 0.0
    clm_copy[:,:,imode+1:] = 0.0
    shmap = np.real(sh.expand.MakeGridDH(clm_copy,sampling=2))
    #stuff for setting up new grid
    n = 2*lmax+2
    shlat = np.arange(-90+180/n,90+180/n,180/n)
    shlon = np.arange(0,360,180/n)
    shlon2d, shlat2d = np.meshgrid(shlon,shlat)
    power = sh.spectralanalysis.spectrum(clm)

    ax = fig.add_subplot(inner_grid[3*(i+1)+imode])
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection='cea')
    m.drawparallels([-45,0,45],labels = [False,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(shlon2d, shlat2d, shmap, cmap=cmap,latlon='True',
                    rasterized=True,vmax=cscale[imode],vmin=-1*cscale[imode])

    if i == len(simnames)-1:
      xlim = ax.get_xlim()
      dxlim = xlim[1] - xlim[0]
      ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
      ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'])
      ax.tick_params(direction='in')

    if i == 0:
      cax = fig.add_subplot(inner_grid[imode])
      cbar = plt.colorbar(c,cax=cax,orientation='horizontal')
      cax.xaxis.set_ticks_position('top')
      cax.xaxis.set_label_position('top')
      cbar.set_label(clabels[imode],fontsize=8)

plt.savefig('figure_6.pdf')
plt.close()
