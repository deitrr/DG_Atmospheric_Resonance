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

label = ['16 hour', '18 hour', '20 hour', '21 hour', '21.5 hour', '22 hour', '22.25 hour', '22.5 hour', '23 hour', '24 hour', '25 hour']
#rotpers = [0.8333, 0.8542, 0.875, 0.8958, 0.9167, 1.0]
rotpers = np.array([0.6667, 0.75, 0.8333, 0.875, 0.8958, 0.9167, 0.9271, 0.9375, 0.9583, 1.0, 1.04167])
#offset = ((365*60)*360/rotpers-180)%360

clabels = ['Surface pressure anomaly (Pa)', 'Surface pressure anomaly (Pa)\ndiurnal mode', 'Surface pressure anomaly (Pa)\nsemidiurnal mode']
cmap = plt.cm.RdBu_r
cscale = [400,50,400]

#fig, axes = plt.subplots(ncols=2,nrows=2,figsize=(7.5,4))
#fig, axes = plt.subplots(ncols=2,nrows=10,figsize=(7.5,11))
fig = plt.figure(figsize=(7.5,9.5))

outer_grid = gridspec.GridSpec(1,1,wspace=0.2,hspace=0.1,left=0.05,right=0.98,bottom=0.03,top=0.94,height_ratios=(1,))
inner_grid = gridspec.GridSpecFromSubplotSpec(12,3,subplot_spec=outer_grid[0],wspace=0.05,hspace=0.15,height_ratios=(1,6,6,6,6,6,6,6,6,6,6,6))

for i in np.arange(len(simnames)):
  #if simnames[i] == 'lr_cam4_4x5_modern_20hr_branch':
  if i < 0:
    data = nc.Dataset(main_path+simnames[i]+'/merged_hist/'+simnames[i]+'.cam.h1.merged_0061_0061.nc','r')

    ps = data['PS'][:]
    lat = data['lat'][:]
    lon = data['lon'][:]

    #need to compute averages over rotation periods
    t0 = data['time'][0]
    #try using annual average
    ps_mean = np.mean(ps,axis=0)

    #pressure anomaly as function of time, lat, lon
    ps_anom  = ps - ps_mean[None,:,:] #using annual mean instead of diurnal

    # rotation rate in rad/Earth solar day
    rotrate = 2*np.pi/rotpers[i]

    #longitude of sun as a function of time
    offset = np.mean(lon[np.argmax(data['FSDS'][:][0,6:-6,:],axis=1)])
    print(offset)
    lon_of_sun = (data['lon'][:][None,:] + (data['time'][:][:,None]-t0)*rotrate*180/np.pi-offset)%360

    ftime = len(data['time'][:])

    ps_anom_reorder = np.zeros_like(ps_anom)
    for itime in np.arange(ftime):
      interp = sint.interp2d(lon_of_sun[itime],lat,ps_anom[itime,:,:])
      ps_anom_reorder[itime,:,:] = interp(lon,lat)

    ps_anom_mean = np.mean(ps_anom_reorder[:ftime,:,:], axis=0)

    print('p anomaly peaks = ',np.max(ps_anom_mean),np.min(ps_anom_mean))

    lon2d, lat2d = np.meshgrid(lon,lat)

    lmax = np.floor(np.sqrt(len(lat)*len(lon))/2)
    clm, chi2 = sh.expand.SHExpandLSQ(ps_anom_mean.ravel(),lat2d.ravel(),lon2d.ravel(),lmax)

  else:
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


#  ax = axes[i][0]
  ax = fig.add_subplot(inner_grid[3*(i+1)])
  m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection='cea')
  m.drawparallels([-45,0,45],labels = [True,False,False,False], fontsize=8)
  m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)

  c = m.pcolormesh(lon2d, lat2d, ps_anom_mean, cmap=cmap,rasterized=True,latlon='True',vmax=cscale[0],vmin=-1*cscale[0])
#  ax.set(title=label[i],ylabel='Latitude')
  ax.yaxis.set_label_coords(-0.1,0.5)
  #if i == len(simnames)-1:
    #ax.set_xlabel('Solar longitude')
    #ax.xaxis.set_label_coords(0.5,-0.08)
  ax.text(0.02,0.85,label[i],rotation=0,transform=ax.transAxes,fontsize=10,color='k',fontweight='bold')

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

  #cax = make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05)
  #cbar = plt.colorbar(c,cax=cax)

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

    print('SH peaks = ',np.max(shmap),np.min(shmap))
    power = sh.spectralanalysis.spectrum(clm)

#  ax = axes[i][1]
    ax = fig.add_subplot(inner_grid[3*(i+1)+imode])
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection='cea')
    m.drawparallels([-45,0,45],labels = [False,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(shlon2d, shlat2d, shmap, cmap=cmap,latlon='True',rasterized=True,vmax=cscale[imode],vmin=-1*cscale[imode])
#  ax.set(title=label[i]+' (l = 2, m = 2 mode)')
#  if i == len(simnames)-1:
#    ax.set_xlabel('Solar longitude')
#    ax.xaxis.set_label_coords(0.5,-0.1)
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

 # cax = make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05)
 # cbar = plt.colorbar(c,cax=cax)
 # cbar.set_label('Pressure anomaly\n(Pa)',fontsize=8)

#plt.tight_layout()
plt.savefig('figure_6.pdf')
plt.close()
