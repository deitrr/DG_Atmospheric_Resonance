import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pathlib
import warnings
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore", category=DeprecationWarning)

#change this path to match the location of the sim_main download from the
#data repository
parent = '../sims_main/'

simnames = [
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch2',
           ]

label = ['22.25 hour', '24 hour']
clabel = ['Horizontal divergence\n(10$^{-6}$ s$^{-1}$)', 'Temperature anomaly\n(K)',
          'Relative humidity\nanomaly (%)', '$\log$[Cloud water (kg m$^{-3}$)]\nanomaly ',
          '$\log$[Humidity (kg kg$^{-1}$)]\nanomaly ']

clevs = [np.linspace(-0.18,0.18,19)*1e-5,np.linspace(-1.2,1.2,25),
         np.linspace(-5,5,21),np.linspace(-2.,2,21),
         np.linspace(-0.012,0.012,13), np.linspace(1.5,8,21)*1e-8]

cscale = [1e-6, 1, 1, 1, 1, 1e-8]

cm = 1./2.54
fig = plt.figure(figsize=(18*cm,21.6*cm))

outer_grid = gridspec.GridSpec(1,1,wspace=0.2,hspace=0.1,left=0.1,right=0.87,
                                    bottom=0.05,top=0.97,height_ratios=(1,))
inner_grid = gridspec.GridSpecFromSubplotSpec(6,3,subplot_spec=outer_grid[0],
    wspace=0.05,hspace=0.15,height_ratios=(1,1,1,1,1,1),width_ratios=(15,15,1))

def smoothing_lon(field,ntimes):
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

for i in np.arange(len(simnames)):
  print(simnames[i])
  ifiles = [
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_DIVV_save.npz',
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_T_save.npz',
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_RELHUM_save.npz',
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_CWC_save.npz',
            parent+simnames[i]+'/merged_hist/'+simnames[i]+'_log_Q_save.npz',
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

    tropics = np.logical_and(lat>=-15,lat<=15)

    ax = fig.add_subplot(inner_grid[3*ifile+i])

    f = np.zeros_like(field_anom_mean)
    if ifile < len(ifiles) - 1:
      if ifile == 0:
        f[:,:,lon<180] = field_mean[:,:,lon>=180]
        f[:,:,lon>=180] = field_mean[:,:,lon<180]
      else:
        f[:,:,lon<180] = field_anom_mean[:,:,lon>=180]
        f[:,:,lon>=180] = field_anom_mean[:,:,lon<180]
      lon = np.hstack([lon[lon>=180]-360,lon[lon<180]])
      fplot = np.mean(f[:,tropics,:]*np.cos(lat[tropics]*np.pi/180)[None,:,None],axis=1)
      p = arc['p']
      c = ax.contourf(lon,p,fplot/cscale[ifile],clevs[ifile]/cscale[ifile],cmap='RdBu_r')
      for cc in c.collections:
        cc.set_edgecolor("face")
      ax.contour(lon,p,fplot/cscale[ifile],[0],colors='k')
      ax.invert_yaxis()
      ax.set(ylim=(1000,50))
      ax.xaxis.set_ticks([-90,0,90])
      ax.xaxis.set_ticklabels([])
      ax.vlines([-90,0,90],1000,50,'k',linestyles=':',linewidths=1)
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

plt.savefig('figures/figure_5.pdf')
plt.close()
