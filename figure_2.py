import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pathlib
from mpl_toolkits.basemap import Basemap
import scipy.interpolate as sint
import pyshtools as sh
import matplotlib.gridspec as gridspec

#this just suppresses a deprecation warning from netCDF4
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

#change this path to match the location of the sim_main download from the
#data repository
parent_path = '../sims_main/'

simnames = [
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch2',
           ]

label = ['22.25 hr', '24 hr']

fig = plt.figure(figsize=(7.5,3))
proj = 'cea'  #cylindrical equal area projection from basemap

outer_grid = gridspec.GridSpec(1,3,wspace=0.1,hspace=0.15,left=0.06,right=0.91,
                                bottom=0.06,top=0.93,width_ratios=(15,15,16))

left_grid1 = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=outer_grid[0,0],
                                                        wspace=0.1,hspace=0.1)
cen_grid1 = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=outer_grid[0,1],
                                                        wspace=0.1,hspace=0.1)
right_grid1 = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec=outer_grid[0,2],
                                      wspace=0.1,hspace=0.1,width_ratios=(15,1))

cscale = [1, 1, 1]
crange = np.array([ [1, 3], [-470, 470],
                    [0, 1.5], [-500, 500],
                    [0, 1.5], [-500,500] ])
crange_sh = np.array([ [-470,470], [-470,470], [-0.1,0.1] ])

clabel = ['Pressure anomaly (Pa)']
tlabels = ['Average anomaly', 'Diurnal mode', 'Semidiurnal mode']

#smoothing in longitude (necessary to remove spikes generate by storms)
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
    #plot trop pressure
    out_file = parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'_p_anom_save.npz'

    arc = np.load(out_file)
    lon = arc['lon']
    lat = arc['lat']
    rotrate = arc['rotrate']
    lmax = arc['lmax']
    lon2d, lat2d = np.meshgrid(lon,lat)

    #stuff for setting up new grid
    n = 2*lmax+2
    shlat = np.arange(-90+180/n,90+180/n,180/n)
    shlon = np.arange(0,360,180/n)
    shlon2d, shlat2d = np.meshgrid(shlon,shlat)

    if 'field_anom_mean' in arc:
        field_anom_mean = arc['field_anom_mean']
    else:
        # in some of the files I used the name ps_anom_mean instead of field_anom_mean
        field_anom_mean = arc['ps_anom_mean']

    clm_mean = arc['clm']
    print(np.max(np.abs(field_anom_mean/cscale[0])))

    field_anom_mean = smoothing_lon(field_anom_mean,10)

    #anomaly
    ax = fig.add_subplot(left_grid1[i])
    if i == 0:
    ax.set_title(tlabels[0],fontsize=8)
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
    m.drawparallels([-60,-30,0,30,60],labels = [True,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(lon2d, lat2d, field_anom_mean/cscale[0], cmap='RdBu_r',rasterized=True,latlon='True',vmax=crange[1][1],vmin=crange[1][0])
    ax.set(ylabel='Latitude')
    ax.yaxis.set_label_coords(-0.15,0.5)
    xlim = ax.get_xlim()
    dxlim = xlim[1] - xlim[0]
    ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
    if i == 1:
    ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=8)
    else:
    ax.xaxis.set_ticklabels([])
    ax.tick_params(direction='in')
    ax.text(0.03,0.9,label[i],rotation=0,transform=ax.transAxes,fontsize=9,color='k',fontweight='bold')

    #spharm
    for imode in [1,2]:
    clm_anom_cp = clm_mean.copy()
    clm_anom_cp[:,:imode,:] = 0.0
    clm_anom_cp[:,imode+1:,:] = 0.0
    clm_anom_cp[:,:,:imode] = 0.0
    clm_anom_cp[:,:,imode+1:] = 0.0
    shmap_anom = np.real(sh.expand.MakeGridDH(clm_anom_cp,sampling=2))
    if imode == 1:
      ax = fig.add_subplot(cen_grid1[i])
    else:
      ax = fig.add_subplot(right_grid1[2*i])
    if i == 0:
      ax.set_title(tlabels[imode],fontsize=8)
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
    m.drawparallels([-60,-30,0,30,60],labels = [False,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(shlon2d, shlat2d, shmap_anom/cscale[0], cmap='RdBu_r',rasterized=True,latlon='True',vmax=crange_sh[imode-1][1],vmin=crange_sh[imode-1][0])
    print(np.max(np.abs(shmap_anom/cscale[0])))

    ax.yaxis.set_label_coords(-0.1,0.5)
    xlim = ax.get_xlim()
    dxlim = xlim[1] - xlim[0]
    ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
    if i == 1:
      ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=8)
    else:
      ax.xaxis.set_ticklabels([])
    ax.tick_params(direction='in')

    if imode == 2:
      cax = fig.add_subplot(right_grid1[2*i+1])
      cbar = plt.colorbar(c,cax=cax)
      cax.tick_params(axis='y',labelsize=6)

    cbar.set_label(clabel[0],fontsize=8)

plt.savefig('figure_2.pdf')
plt.close()
