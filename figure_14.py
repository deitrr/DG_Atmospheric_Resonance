import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pdb
import pathlib
import warnings
#from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import scipy.interpolate as sint
#from scipy.integrate import simps
#import pyshtools as sh
#pdb.set_trace()
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore", category=DeprecationWarning)

parent = '../sims_main/'


simnames = [
#            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_16hr_branch',
#            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_18hr_branch',
#            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_20hr_branch',
#            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21hr_branch',
#            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21-5hr_branch',
#            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr_branch2',
#            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-5hr_branch',
#            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_23hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch',
#            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_25hr_branch',
           ]

#rotpers = np.array([0.6667, 0.75, 0.8333,0.875, 0.8958, 0.9167, 0.9375, 1.0, 1.04167])

#label = ['16 hour', '18 hour', '20 hour', '21 hour', '21.5 hour', '22 hour', '22.5 hour', '23 hour', '24 hour', '25 hour']
#label = ['22 hour', '24 hour']
label = ['22.25 hour', '24 hour']
clabel = ['Mass flux anomaly (shallow)\n(10$^{-3}$ kg m$^{-2}$ s$^{-1}$)', 'Mass flux anomaly (deep)\n(10$^{-3}$ kg m$^{-2}$ s$^{-1}$)']

#ilons = [54, 0, 18]
ilat = 22

#clevs = [np.linspace(-0.18,0.18,19)*1e-5,np.linspace(-1.2,1.2,25),np.linspace(-5,5,21),np.linspace(-2.,2,21),
#         np.linspace(-0.012,0.012,13),
#         np.linspace(1.5,8,21)*1e-8]
clevs = [np.linspace(-3,3,25)*1e-3,np.linspace(-3,3,25)*1e-3,np.linspace(1.5,8,21)*1e-8]
#np.linspace(-2.75,2.75,21)*1e-8]
cscale = [ 1e-3, 1e-3, 1e-8]

#fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7.5,7))
fig = plt.figure(figsize=(7.5,6))

outer_grid = gridspec.GridSpec(1,1,wspace=0.2,hspace=0.1,left=0.1,right=0.87,bottom=0.05,top=0.96,height_ratios=(1,))
inner_grid = gridspec.GridSpecFromSubplotSpec(3,3,subplot_spec=outer_grid[0],wspace=0.05,hspace=0.15,height_ratios=(1,1,1),width_ratios=(15,15,1))

def smoothing_lon(field,ntimes):
#   print("smoothing")
   field_tmp = field.copy()
   field_cp = np.zeros_like(field)

   if len(np.shape(field)) == 2:
     for i in np.arange(ntimes):
       field_cp[:,1:-1] = 0.25*field_tmp[:,:-2] + 0.5*field_tmp[:,1:-1] + 0.25*field_tmp[:,2:]
       field_cp[:,0] = 0.25*field_tmp[:,-1] + 0.5*field_tmp[:,0] + 0.25*field_tmp[:,1]
       field_cp[:,-1] = 0.25*field_tmp[:,-2] + 0.5*field_tmp[:,-1] + 0.25*field_tmp[:,0]
       field_tmp = field_cp.copy()
   elif len(np.shape(field)) == 3:
     for i in np.arange(ntimes):
       field_cp[:,:,1:-1] = 0.25*field_tmp[:,:,:-2] + 0.5*field_tmp[:,:,1:-1] + 0.25*field_tmp[:,:,2:]
       field_cp[:,:,0] = 0.25*field_tmp[:,:,-1] + 0.5*field_tmp[:,:,0] + 0.25*field_tmp[:,:,1]
       field_cp[:,:,-1] = 0.25*field_tmp[:,:,-2] + 0.5*field_tmp[:,:,-1] + 0.25*field_tmp[:,:,0]
       field_tmp = field_cp.copy()
   
   return field_cp

def smoothing_lon7(field,ntimes):
#   print("smoothing")
   field_tmp = field.copy()
   field_cp = np.zeros_like(field)

   if len(np.shape(field)) == 2:
     for i in np.arange(ntimes):
       field_cp[:,0] = 1.0/64*(field_tmp[:,-3]+6*field_tmp[:,-2]+15*field_tmp[:,-1] 
                                  + 20*field_tmp[:,0] + 15*field_tmp[:,1] + 6*field_tmp[:,2] + field_tmp[:,3])
       field_cp[:,1] = 1.0/64*(field_tmp[:,-2]+6*field_tmp[:,-1]+15*field_tmp[:,0] 
                                  + 20*field_tmp[:,1] + 15*field_tmp[:,2] + 6*field_tmp[:,3] + field_tmp[:,4])
       field_cp[:,2] = 1.0/64*(field_tmp[:,-1]+6*field_tmp[:,0]+15*field_tmp[:,1] 
                                  + 20*field_tmp[:,2] + 15*field_tmp[:,3] + 6*field_tmp[:,4] + field_tmp[:,5])
       field_cp[:,3:-3] = 1.0/64*(field_tmp[:,:-6]+6*field_tmp[:,1:-5]+15*field_tmp[:,2:-4] 
			          + 20*field_tmp[:,3:-3] + 15*field_tmp[:,4:-2] + 6*field_tmp[:,5:-1] + field_tmp[:,6:])
       field_cp[:,-3] = 1.0/64*(field_tmp[:,-6]+6*field_tmp[:,-5]+15*field_tmp[:,-4] 
                                  + 20*field_tmp[:,-3] + 15*field_tmp[:,-2] + 6*field_tmp[:,-1] + field_tmp[:,0])
       field_cp[:,-2] = 1.0/64*(field_tmp[:,-5]+6*field_tmp[:,-4]+15*field_tmp[:,-3] 
                                  + 20*field_tmp[:,-2] + 15*field_tmp[:,-1] + 6*field_tmp[:,0] + field_tmp[:,1])
       field_cp[:,-1] = 1.0/64*(field_tmp[:,-4]+6*field_tmp[:,-3]+15*field_tmp[:,-2] 
                                  + 20*field_tmp[:,-1] + 15*field_tmp[:,0] + 6*field_tmp[:,1] + field_tmp[:,2])
       field_tmp = field_cp.copy()
   elif len(np.shape(field)) == 3:
     for i in np.arange(ntimes):
       field_cp[:,:,0] = 1.0/64*(field_tmp[:,:,-3]+6*field_tmp[:,:,-2]+15*field_tmp[:,:,-1] 
                                  + 20*field_tmp[:,:,0] + 15*field_tmp[:,:,1] + 6*field_tmp[:,:,2] + field_tmp[:,:,3])
       field_cp[:,:,1] = 1.0/64*(field_tmp[:,:,-2]+6*field_tmp[:,:,-1]+15*field_tmp[:,:,0] 
                                  + 20*field_tmp[:,:,1] + 15*field_tmp[:,:,2] + 6*field_tmp[:,:,3] + field_tmp[:,:,4])
       field_cp[:,:,2] = 1.0/64*(field_tmp[:,:,-1]+6*field_tmp[:,:,0]+15*field_tmp[:,:,1] 
                                  + 20*field_tmp[:,:,2] + 15*field_tmp[:,:,3] + 6*field_tmp[:,:,4] + field_tmp[:,:,5])
       field_cp[:,:,3:-3] = 1.0/64*(field_tmp[:,:,:-6]+6*field_tmp[:,:,1:-5]+15*field_tmp[:,:,2:-4] 
                                  + 20*field_tmp[:,:,3:-3] + 15*field_tmp[:,:,4:-2] + 6*field_tmp[:,:,5:-1] + field_tmp[:,:,6:])
       field_cp[:,:,-3] = 1.0/64*(field_tmp[:,:,-6]+6*field_tmp[:,:,-5]+15*field_tmp[:,:,-4] 
                                  + 20*field_tmp[:,:,-3] + 15*field_tmp[:,:,-2] + 6*field_tmp[:,:,-1] + field_tmp[:,:,0])
       field_cp[:,:,-2] = 1.0/64*(field_tmp[:,:,-5]+6*field_tmp[:,:,-4]+15*field_tmp[:,:,-3] 
                                  + 20*field_tmp[:,:,-2] + 15*field_tmp[:,:,-1] + 6*field_tmp[:,:,0] + field_tmp[:,:,1])
       field_cp[:,:,-1] = 1.0/64*(field_tmp[:,:,-4]+6*field_tmp[:,:,-3]+15*field_tmp[:,:,-2] 
                                  + 20*field_tmp[:,:,-1] + 15*field_tmp[:,:,0] + 6*field_tmp[:,:,1] + field_tmp[:,:,2])
       field_tmp = field_cp.copy()
   
   return field_cp

for i in np.arange(len(simnames)):
  print(simnames[i])
  ifiles = [
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_CMFSH_save.npz',
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_CMFMCDZM_save.npz',
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_PRECT_save.npz',
  ]

  for ifile in np.arange(len(ifiles)):
    arc = np.load(ifiles[ifile])
    lon = arc['lon']
    lat = arc['lat']
    rotrate = arc['rotrate']
    field_anom_mean = arc['field_anom_mean']
    field_mean = arc['field_mean']
    clm_mean = arc['clm_mean']
    clm_anom = arc['clm_anom']
    lmax = arc['lmax']
   
    field_mean = smoothing_lon(field_mean,3)
    field_anom_mean = smoothing_lon(field_anom_mean,3)


    if 'PRECT' in ifiles[ifile]:
      land = parent+simnames[i]+'/merged_hist/'+simnames[i]+'_land_PRECT_save.npz'
      arcl = np.load(land)
      prec_anom_land = arcl['field_anom_mean']
      prec_land = arcl['field_mean']
      ocean = parent+simnames[i]+'/merged_hist/'+simnames[i]+'_ocean_PRECT_save.npz'
      arco = np.load(ocean)
      prec_anom_ocean = arco['field_anom_mean']
      prec_ocean = arco['field_mean']
 
      prec_land = smoothing_lon(prec_land,4)
      prec_ocean = smoothing_lon(prec_ocean,4)
#      prec_land = smoothing_lon7(prec_land,1)
#      prec_ocean = smoothing_lon7(prec_ocean,)

#    if 'CMFMC' in ifiles[ifile]:
#       deep = parent+simnames[i]+'/merged_hist/'+simnames[i]+'_CMFMCDZM_save.npz'
#       arcd = np.load(deep)
#       field_mean_deep = arc['field_mean']

    tropics = np.logical_and(lat>=-15,lat<=15)

    print(lat[ilat])
    #ax = axes[ifile][i]
    ax = fig.add_subplot(inner_grid[3*ifile+i])

#    pdb.set_trace()
    f = np.zeros_like(field_anom_mean)
    if ifile < len(ifiles) - 1:
      if ifile == 0:
        f[:,:,lon<180] = field_anom_mean[:,:,lon>=180] # - field_mean_deep[:,:,lon>=180]
        f[:,:,lon>=180] = field_anom_mean[:,:,lon<180] #- field_mean_deep[:,:,lon<180]
      else:
        f[:,:,lon<180] = field_anom_mean[:,:,lon>=180]
        f[:,:,lon>=180] = field_anom_mean[:,:,lon<180]

      p = arc['p']
      if len(f) > len(p):  #fluxes are defined on interfaces rather than midpoints
        f = 0.5*(f[:-1] + f[1:])

      lon = np.hstack([lon[lon>=180]-360,lon[lon<180]])
      fplot = np.mean(f[:,tropics,:]*np.cos(lat[tropics]*np.pi/180)[None,:,None],axis=1)

      c = ax.contourf(lon,p,fplot/cscale[ifile],clevs[ifile]/cscale[ifile],cmap='RdBu_r')
      for cc in c.collections:
        cc.set_edgecolor("face")
      ax.contour(lon,p,fplot/cscale[ifile],[0],colors='k')
      ax.invert_yaxis()
      ax.set(ylim=(1000,50))
      ax.xaxis.set_ticks([-90,0,90])
#      ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'])
      ax.xaxis.set_ticklabels([])
      ax.vlines([-90,0,90],1000,50,'k',linestyles=':',linewidths=1)
      #cax = make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05)
      #cbar = plt.colorbar(c,cax=cax)
      if i == 1:
        cax = fig.add_subplot(inner_grid[3*ifile+2])
        cbar = plt.colorbar(c,cax=cax)
        cbar.set_label(clabel[ifile],fontsize=8)
      if i == 0:
        ax.set(ylabel = 'Pressure (hPa)')
    else:
      f[:,lon<180] = field_mean[:,lon>=180]
      f[:,lon>=180] = field_mean[:,lon<180]
      fplot = np.mean(f[tropics,:]*np.cos(lat[tropics]*np.pi/180)[:,None],axis=0)
      #pdb.set_trace()

      f[:,lon<180] = prec_land[:,lon>=180]
      f[:,lon>=180] = prec_land[:,lon<180]
      fland = np.nanmean(f[tropics,:]*np.cos(lat[tropics]*np.pi/180)[:,None],axis=0)

      f[:,lon<180] = prec_ocean[:,lon>=180]
      f[:,lon>=180] = prec_ocean[:,lon<180]
      focean = np.nanmean(f[tropics,:]*np.cos(lat[tropics]*np.pi/180)[:,None],axis=0)

      lon = np.hstack([lon[lon>=180]-360,lon[lon<180]])
      ax.plot(lon,fplot/cscale[ifile],label='land +\nocean')
      ax.plot(lon,fland/cscale[ifile],c='orange',label='land')
      ax.plot(lon,focean/cscale[ifile],c='darkblue',label='ocean')

      ax.set(ylim=(clevs[ifile][0]/cscale[ifile],clevs[ifile][-1]/cscale[ifile]),xlim=(-180,176))
#      ax.set(xlim=(-180,180))
      ax.xaxis.set_ticks([-90,0,90])
      ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'])
      ax.vlines([-90,0,90],clevs[ifile][0]/cscale[ifile],clevs[ifile][-1]/cscale[ifile],'k',linestyles=':',linewidths=1)
      if i == 0:
        ax.set(ylabel = 'Precipation rate\n(10$^{-8}$ m s$^{-1}$)')

    if i == 1:
      ax.yaxis.set_ticklabels([])
      if ifile == len(ifiles)-1:
        ax.legend(loc = 'lower left',bbox_to_anchor=(1.02,0),fontsize=7)
    if ifile == 0:
      ax.set_title(label[i])

#plt.tight_layout()
plt.savefig('figure_14.pdf')
plt.close()
