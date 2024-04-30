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

warnings.filterwarnings("ignore", category=DeprecationWarning)

parent_path = '/home/rdeitrick/waccm_rot_sims/'

simnames = [
            #'lambres_4x5_modern_20hr_branch', #0
            #'lambres_4x5_modern_20-5hr_branch2',
            #'lambres_4x5_modern_21hr_branch2',
            #'lambres_4x5_modern_21-5hr_branch2',
            #'lr_cam4_4x5_modern_16hr_branch', #4
            #'lr_cam4_4x5_modern_18hr_branch',
            #'lr_cam4_4x5_modern_20hr_branch',
            #'lr_cam4_4x5_modern_21hr_branch',
            #'lr_cam4_4x5_modern_22hr_branch',
            #'lr_cam4_4x5_modern_23hr_branch',
            #'lr_cam4_4x5_modern_24hr_branch',
            #'lr_cam4_4x5_modern_25hr_branch',
            #'lr_exocam_4x5_lowCH4_16hr_branch', #12
            #'lr_exocam_4x5_lowCH4_18hr_branch',
            #'lr_exocam_4x5_lowCH4_20hr_branch',
            #'lr_exocam_4x5_lowCH4_21hr_branch',
            #'lr_exocam_4x5_lowCH4_22hr_branch',
            #'lr_exocam_4x5_lowCH4_23hr_branch',
            #'lr_exocam_4x5_lowCH4_24hr_branch2',
            #'lr_exocam_4x5_hiCH4_16hr_branch', #19
            #'lr_exocam_4x5_hiCH4_18hr_branch',
            #'lr_exocam_4x5_hiCH4_20hr_branch',
            #'lr_exocam_4x5_hiCH4_21hr_branch',
            #'lr_exocam_4x5_hiCH4_22hr_branch2',
            #'lr_exocam_4x5_hiCH4_23hr_branch',
            #'lr_exocam_4x5_hiCH4_24hr_branch'
#            'nomg1_lr_exocam_4x5_hiCH4_16hr_branch2', #12
#            'nomg1_lr_exocam_4x5_hiCH4_18hr_branch2',
#            'nomg1_lr_exocam_4x5_hiCH4_20hr_branch2',
#            'nomg1_lr_exocam_4x5_hiCH4_21hr_branch2',
#            'nomg1_lr_exocam_4x5_hiCH4_21-5hr_branch2',
#            'nomg1_lr_exocam_4x5_hiCH4_22hr_branch2',
#            'nomg1_lr_exocam_4x5_hiCH4_23hr_branch2',
#            'nomg1_lr_exocam_4x5_hiCH4_24hr_branch2',
#            'nomg1_lr_exocam_4x5_hiCH4_25hr_branch2',
            #'nomg1_lr_exocam_4x5_lowCH4_16hr_branch', #21
            #'nomg1_lr_exocam_4x5_lowCH4_18hr_branch',
            #'nomg1_lr_exocam_4x5_lowCH4_20hr_branch',
            #'nomg1_lr_exocam_4x5_lowCH4_21hr_branch',
            #'nomg1_lr_exocam_4x5_lowCH4_22hr_branch',
            #'nomg1_lr_exocam_4x5_lowCH4_23hr_branch',
            #'nomg1_lr_exocam_4x5_lowCH4_24hr_branch',
            #'nomg1_lr_exocam_4x5_lowCH4_25hr_branch',
            #'wcsc_1850_pop2_4x5_branch' #29
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

#rotpers = np.array([0.8333, 0.8542, 0.875, 0.8958,
#                    0.6667, 0.75, 0.8333, 0.875, 0.9167, 0.9583, 1.0, 1.04167,
#                    0.6667, 0.75, 0.8333, 0.875, 0.8958, 0.9167, 0.9583, 1.0, 1.04167,
#                    0.6667, 0.75, 0.8333, 0.875, 0.9167, 0.9583, 1.0, 1.04167,
#                    1.0  ])
rotpers = np.array([0.6667, 0.75, 0.8333,0.875, 0.8958, 0.9167,0.9375, 0.9583, 1.0, 1.04167])
label = ['16 hour', '18 hour', '20 hour', '21 hour', '21.5 hour', '22 hour', '22.5 hour', '23 hour', '24 hour', '25 hour']
mean_lim = [-1.3e-6, 1.8e-6]  #for 900 hPa level
anom_lim = [-8e-7, 8e-7]
lm_lim = np.array([[-8e-8,8e-8],[-5e-7,5e-7]])
#mean_lim = [-1.3e-6, 1.8e-6]  #for 250 hPa level
#anom_lim = [-8.5e-7, 8.5e-7]
#lm_lim = np.array([[-3e-8,3e-8],[-6e-7,6e-7]])
field = 'DIVV'
clabel= 'Wind divergence\nat 200 hPa'
savelabel = field + '200'
#ilev = 38 # roughly pressure = 900 hPa
#ilev = 28  #roughly 266 hPa
ilev = 26 #200 hPa
#ilev = 22 #100 hPa
#ilev = 32 #500 hPa
#ilev = 29 #300 hPa
ilev = 39 #surface

cmap = plt.cm.RdBu_r
fig, axes = plt.subplots(ncols=4,nrows=10,figsize=(16,12))

for i in np.arange(len(simnames)):
  out_file = parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'_'+field+'_save.npz'
#  out_file_u = parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'_U_save.npz'
#  out_file_v = parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'_V_save.npz'
  print(simnames[i])

  if not pathlib.Path(out_file).exists():
    if simnames[i] == 'lr_exocam_4x5_hiCH4_22hr_branch2':
      data = nc.Dataset(simnames[i]+'/merged_hist/'+simnames[i]+'.cam.h1.merged_0065_0065.nc','r')
    else:
      data = nc.Dataset(parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'.cam.h1.merged_0061_0061.nc','r')

    field0 = data[field][:,32,:,:]
    lat = data['lat'][:]
    lon = data['lon'][:]

    #need to compute averages over rotation periods
    t0 = data['time'][0]
    #try using annual average
    field0_mean = np.mean(field0,axis=0)

    #pressure anomaly as function of time, lat, lon
    field_anom  = field0 - field0_mean[None,:,:] #using annual mean instead of diurnal

    # rotation rate in rad/Earth solar day
    rotrate = 2*np.pi/rotpers[i]

    #longitude of sun as a function of time
    #for ilat in np.arange(len(lat)):
    offset = np.mean(lon[np.argmax(data['FSDS'][:][0,6:-6,:],axis=1)])
    #print(offset)
    #offset = 0.0 #180.0 + ((t0+1) * rotrate * 180/np.pi)%360
    #additional correction from fsds file
    fsds_file = parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'_fsds_phase_save.npz'
    arc2 = np.load(fsds_file)
    fsds_mean = arc2['ps_anom_mean'] #forgot to change the name in this file!
    clm_fsds = arc2['clm']
    offset2 = np.arctan2(clm_fsds[1,1,1],clm_fsds[0,1,1])*180/np.pi

    lon_of_sun = (data['lon'][:][None,:] + (data['time'][:][:,None]-t0)*rotrate*180/np.pi-offset-offset2)%360

    ftime = len(data['time'][:])

    field_reorder = np.zeros_like(field0)
    field_anom_reorder = np.zeros_like(field_anom)
    for itime in np.arange(ftime): #1000):
      interp = sint.interp2d(lon_of_sun[itime],lat,field0[itime,:,:])
      field_reorder[itime,:,:] = interp(lon,lat)
      interp = sint.interp2d(lon_of_sun[itime],lat,field_anom[itime,:,:])
      field_anom_reorder[itime,:,:] = interp(lon,lat)

    lon2d, lat2d = np.meshgrid(lon,lat)
    field_mean = np.mean(field_reorder[:ftime,:,:], axis=0)
    field_anom_mean = np.mean(field_anom_reorder[:ftime,:,:], axis=0)

    lmax = np.floor(np.sqrt(len(lat)*len(lon))/2)
    clm_mean, chi2 = sh.expand.SHExpandLSQ(field_mean.ravel(),lat2d.ravel(),lon2d.ravel(),lmax)
    clm_anom, chi2 = sh.expand.SHExpandLSQ(field_anom_mean.ravel(),lat2d.ravel(),lon2d.ravel(),lmax)

    np.savez(out_file,lon=lon,lat=lat,rotrate=rotrate,field_anom_mean=field_anom_mean,field_mean=field_mean, clm_mean=clm_mean,clm_anom=clm_anom,lmax=lmax)

  else:

    arc = np.load(out_file)
    lon = arc['lon']
    lat = arc['lat']
    #ps = arc['ps']
    rotrate = arc['rotrate']
    field_anom_mean = arc['field_anom_mean'][ilev]
    field_mean = arc['field_mean'][ilev]
    clm_mean = arc['clm_mean'][ilev]
    clm_anom = arc['clm_anom'][ilev]
    lmax = arc['lmax']
    lon2d, lat2d = np.meshgrid(lon,lat)

    #arc_u = np.load(out_file_u)
    #u = arc_u['field_mean'][ilev]
    #uanom = arc_u['field_anom_mean'][ilev]

#    arc_v = np.load(out_file_v)
#    v = arc_v['field_mean'][ilev]
#    vanom = arc_v['field_anom_mean'][ilev]

#    pdb.set_trace()

  print('field peaks = ',np.max(field_mean),np.min(field_mean))
  print('field anom peaks = ',np.max(field_anom_mean),np.min(field_anom_mean))

#    lmax = arc['lmax']
#    clm = arc['clm']

  #stuff for setting up new grid
  n = 2*lmax+2
  shlat = np.arange(-90+180/n,90+180/n,180/n)
  shlon = np.arange(0,360,180/n)
  shlon2d, shlat2d = np.meshgrid(shlon,shlat)

#  dp22[i] = np.sqrt(clm[0,2,2]**2 + clm[1,2,2]**2)
#  p_anom_amp[i] = np.max(shmap)


  ax = axes[i][0]
  m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False)
  m.drawparallels([-45,0,45],labels = [True,False,False,False], fontsize=6)
  m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
  c = m.pcolormesh(lon2d, lat2d, field_mean, cmap='plasma',rasterized=True,latlon='True',vmax=mean_lim[1],vmin=mean_lim[0])
  ax.set(title=label[i]+' (full)',ylabel='Latitude')
  ax.yaxis.set_label_coords(-0.1,0.5)
  ax.xaxis.set_ticks([-90,0,90])
  ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'])
  ax.tick_params(direction='in')

#  lon[lon>180] -= 360
#  ax.quiver(lon[::5], lat[::4], u[::4,::5], v[::4,::5],color='w')

#  if i == len(simnames)-1:
#    ax.set_xlabel('Solar longitude')
#    ax.xaxis.set_label_coords(0.5,-0.08)
  cax = make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05)
  cbar = plt.colorbar(c,cax=cax)

#  ax = axes[i][1]
#  m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False)
#  m.drawparallels([-45,0,45],labels = [True,False,False,False], fontsize=6)
#  m.drawmeridians([-90,0,90],labels = [False,False,False,True], fontsize=6)
#  c = m.pcolormesh(shlon2d, shlat2d, shmap_mean, cmap=cmap,latlon='True',rasterized=True,vmax=lm_lim[1],vmin=lm_lim[0])
#  ax.set(title=label[i]+' (l = 2, m = 2 mode)')
#  if i == len(simnames)-1:
#    ax.set_xlabel('Solar longitude')
#    ax.xaxis.set_label_coords(0.5,-0.1)
#  cax = make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05)
#  cbar = plt.colorbar(c,cax=cax)

  ax = axes[i][1]
  m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False)
  m.drawparallels([-45,0,45],labels = [True,False,False,False], fontsize=6)
  m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
  c = m.pcolormesh(lon2d, lat2d, field_anom_mean, cmap=cmap,rasterized=True,latlon='True',vmax=anom_lim[1],vmin=anom_lim[0])
  ax.set(title=label[i]+' (anomaly)')
  ax.yaxis.set_label_coords(-0.1,0.5)
  ax.xaxis.set_ticks([-90,0,90])
  ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'])
  ax.tick_params(direction='in')

#  ax.quiver(lon[::5], lat[::4], uanom[::4,::5], vanom[::4,::5])

#  if i == len(simnames)-1:
#    ax.set_xlabel('Solar longitude')
#    ax.xaxis.set_label_coords(0.5,-0.08)
  cax = make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05)
  cbar = plt.colorbar(c,cax=cax)

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

    print('SH peaks (mean) = ',np.max(shmap_mean),np.min(shmap_mean))
    print('SH peaks (anomaly) = ',np.max(shmap_anom),np.min(shmap_anom))

    ax = axes[i][imode+1]
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False)
    m.drawparallels([-45,0,45],labels = [True,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(shlon2d, shlat2d, shmap_anom, cmap=cmap,latlon='True',rasterized=True,vmax=lm_lim[imode-1][1],vmin=lm_lim[imode-1][0])
    ax.set(title=label[i]+' (l = %d, m = %d mode)'%(imode,imode))
    ax.xaxis.set_ticks([-90,0,90])
    ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'])
    ax.tick_params(direction='in')

   # if i == len(simnames)-1:
   #   ax.set_xlabel('Solar longitude')
   #   ax.xaxis.set_label_coords(0.5,-0.1)
    cax = make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05)
    cbar = plt.colorbar(c,cax=cax)

  cbar.set_label(clabel)

plt.tight_layout()
plt.savefig(parent_path + savelabel+'_dasl_s0p9_ch4-30.pdf')
plt.close()


#fig, axes = plt.subplots(ncols=1,nrows=1,figsize=(6,5))

#for iset in np.arange(len(sets)-1):
  #axes.plot(rotpers[sets[iset]:sets[iset+1]]*24, dp22[sets[iset]:sets[iset+1]], '^', color=colors[iset], linestyle='-')
#  axes.plot(rotpers[sets[iset]:sets[iset+1]]*24, p_anom_amp[sets[iset]:sets[iset+1]], 's', color=colors[iset], linestyle='-', label=labels[iset])
#axes.set_xlabel('Rotation period (hours)')
#axes.set_ylabel('Pressure anomaly (Pa)')
#axes.legend(loc='best')

#plt.tight_layout()
#plt.savefig('PanomVsRot.pdf')
#plt.close()
