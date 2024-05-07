import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pathlib
import warnings
from mpl_toolkits.basemap import Basemap
import pyshtools as sh
import matplotlib.gridspec as gridspec


warnings.filterwarnings("ignore", category=DeprecationWarning)

#change this path to match the location of the sim_main download from the
#data repository
parent_path = '../sims_main/'

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

rotpers = np.array([
  0.6667, 0.75, 0.8333, 0.875, 0.8958, 0.9167, 0.9271, 0.9375, 0.9583, 1.0, 1.04167
])
label = ['16 hour', '18 hour', '20 hour', '21 hour', '21.5 hour', '22 hour',
            '22.25 hour', '22.5 hour', '23 hour', '24 hour', '25 hour']
mean_lim = [3e-9, 11e-8]
mean_scale = [1e-7, '10$^{-7}$']

anom_lim = [-2.9e-8, 2.9e-8]
anom_scale = [1e-8, '10$^{-8}$']

lm_lim = np.array([[-6e-9,6e-9],[-8e-9,8e-9]])
lm_scale = [[1e-9, '10$^{-9}$'],[1e-9, '10$^{-9}$']]

field = 'PRECT'
clabel= 'Total precip\nrate (m s$^{-1}$)'
savelabel = field

cmap = plt.cm.RdBu_r
fig = plt.figure(figsize=(10,7.5))

outer_grid = gridspec.GridSpec(1,1,wspace=0.2,hspace=0.1,left=0.03,right=0.98,
                                    bottom=0.05,top=0.93,height_ratios=(1,))
inner_grid = gridspec.GridSpecFromSubplotSpec(12,4,subplot_spec=outer_grid[0]
                ,wspace=0.15,hspace=0.1,height_ratios=(1,6,6,6,6,6,6,6,6,6,6,6))

def smoothing_lon(field,ntimes):
   field_tmp = field.copy()
   field_cp = np.zeros_like(field)
   for i in np.arange(ntimes):
     field_cp[:,1:-1] = 0.25*field_tmp[:,:-2] + 0.5*field_tmp[:,1:-1] + 0.25*field_tmp[:,2:]
     field_cp[:,0] = 0.25*field_tmp[:,-1] + 0.5*field_tmp[:,0] + 0.25*field_tmp[:,1]
     field_cp[:,-1] = 0.25*field_tmp[:,-2] + 0.5*field_tmp[:,-1] + 0.25*field_tmp[:,0]
     field_tmp = field_cp.copy()

   return field_cp

for i in np.arange(len(simnames)):
    out_file = parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'_'+savelabel+'_save.npz'
    print(simnames[i])


    arc = np.load(out_file)
    lon = arc['lon']
    lat = arc['lat']
    rotrate = arc['rotrate']
    field_anom_mean = arc['field_anom_mean']
    field_mean = arc['field_mean']
    clm_mean = arc['clm_mean']
    clm_anom = arc['clm_anom']
    lmax = arc['lmax']
    lon2d, lat2d = np.meshgrid(lon,lat)

    field_mean = smoothing_lon(field_mean,10)
    field_anom_mean = smoothing_lon(field_anom_mean,10)

    #stuff for setting up new grid
    n = 2*lmax+2
    shlat = np.arange(-90+180/n,90+180/n,180/n)
    shlon = np.arange(0,360,180/n)
    shlon2d, shlat2d = np.meshgrid(shlon,shlat)

    ax = fig.add_subplot(inner_grid[4*(i+1)])
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection='cea')
    m.drawparallels([-45,0,45],labels = [True,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(lon2d, lat2d, field_mean/mean_scale[0], cmap='plasma',
                    rasterized=True,latlon='True',vmax=mean_lim[1]/mean_scale[0],
                    vmin=mean_lim[0]/mean_scale[0])
    ax.yaxis.set_label_coords(-0.1,0.5)
    xlim = ax.get_xlim()
    dxlim = xlim[1] - xlim[0]
    ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
    ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'])
    ax.tick_params(direction='in')
    ax.text(0.02,0.8,label[i],rotation=0,transform=ax.transAxes,fontsize=10,
                            color='w',fontweight='bold')
    if i == 0:
      cax = fig.add_subplot(inner_grid[0])
      cbar = plt.colorbar(c,cax=cax,orientation='horizontal')
      cax.xaxis.set_ticks_position('top')
      cax.xaxis.set_label_position('top')
      cbar.set_label('Total precipitation rate (%s m s$^{-1}$)'%mean_scale[1])

    ax = fig.add_subplot(inner_grid[4*(i+1)+1])
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection='cea')
    m.drawparallels([-45,0,45],labels = [True,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(lon2d, lat2d, field_anom_mean/anom_scale[0], cmap=cmap,
                    rasterized=True,latlon='True',vmax=anom_lim[1]/anom_scale[0],
                    vmin=anom_lim[0]/anom_scale[0])
    ax.yaxis.set_label_coords(-0.1,0.5)
    xlim = ax.get_xlim()
    dxlim = xlim[1] - xlim[0]
    ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
    ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'])
    ax.tick_params(direction='in')
    if i == 0:
      cax = fig.add_subplot(inner_grid[1])
      cbar = plt.colorbar(c,cax=cax,orientation='horizontal')
      cax.xaxis.set_ticks_position('top')
      cax.xaxis.set_label_position('top')
      cbar.set_label('Precip rate anomaly (%s m s$^{-1}$)'%anom_scale[1])

    mode_label = ['Diurnal', 'Semidiurnal']
    for imode in np.array([1,2]):
      clm_mean_cp = clm_mean.copy()
      clm_mean_cp[:,:imode,:] = 0.0
      clm_mean_cp[:,imode+1:,:] = 0.0
      clm_mean_cp[:,:,:imode] = 0.0
      clm_mean_cp[:,:,imode+1:] = 0.0
      shmap_mean = np.real(sh.expand.MakeGridDH(clm_mean_cp,sampling=2))

      clm_anom_cp = clm_mean.copy()
      clm_anom_cp[:,:imode,:] = 0.0
      clm_anom_cp[:,imode+1:,:] = 0.0
      clm_anom_cp[:,:,:imode] = 0.0
      clm_anom_cp[:,:,imode+1:] = 0.0
      shmap_anom = np.real(sh.expand.MakeGridDH(clm_anom_cp,sampling=2))

      ax = fig.add_subplot(inner_grid[4*(i+1) + imode+1])
      m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection='cea')
      m.drawparallels([-45,0,45],labels = [True,False,False,False], fontsize=6)
      m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
      c = m.pcolormesh(shlon2d, shlat2d, shmap_anom/lm_scale[imode-1][0], cmap=cmap,
            latlon='True',rasterized=True,vmax=lm_lim[imode-1][1]/lm_scale[imode-1][0],
            vmin=lm_lim[imode-1][0]/lm_scale[imode-1][0])
      xlim = ax.get_xlim()
      dxlim = xlim[1] - xlim[0]
      ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
      ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'])
      ax.tick_params(direction='in')
      if i == 0:
        cax = fig.add_subplot(inner_grid[imode+1])
        cbar = plt.colorbar(c,cax=cax,orientation='horizontal')
        cax.xaxis.set_ticks_position('top')
        cax.xaxis.set_label_position('top')
        cbar.set_label(mode_label[imode-1]+' mode (%s m s$^{-1}$)'%lm_scale[imode-1][1])

plt.savefig('figures/ED_figure_5.eps')
plt.close()
