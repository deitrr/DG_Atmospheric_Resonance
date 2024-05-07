import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pathlib
import warnings
from mpl_toolkits.basemap import Basemap
import pyshtools as sh
import matplotlib.gridspec as gridspec
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':7})

warnings.filterwarnings("ignore", category=DeprecationWarning)

#change this path to match the location of the sim_main download from the
#data repository
parent_path = '../DG_Atmospheric_Resonance_archive/sims_main/'

simnames = [
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch2'
           ]

label = ['22.25 hr', '24 hr']

cm = 1./2.54
fig = plt.figure(figsize=(18*cm,21.6*cm))
proj = 'cea'  #cylindrical equal area projection

outer_grid = gridspec.GridSpec(4,2,wspace=0.15,hspace=0.15,left=0.06,right=0.88,
                                    bottom=0.03,top=0.97,width_ratios=(16,31))

left_grid1 = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec=outer_grid[0,0],
                                    wspace=0.1,hspace=0.15,width_ratios=(15,1))
right_grid1 = gridspec.GridSpecFromSubplotSpec(2,3,subplot_spec=outer_grid[0,1],
                                    wspace=0.1,hspace=0.15,width_ratios=(15,15,1))
left_grid2 = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec=outer_grid[1,0],
                                    wspace=0.1,hspace=0.15,width_ratios=(15,1))
right_grid2 = gridspec.GridSpecFromSubplotSpec(2,3,subplot_spec=outer_grid[1,1],
                                    wspace=0.1,hspace=0.15,width_ratios=(15,15,1))
left_grid3 = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec=outer_grid[2,0],
                                    wspace=0.1,hspace=0.15,width_ratios=(15,1))
right_grid3 = gridspec.GridSpecFromSubplotSpec(2,3,subplot_spec=outer_grid[2,1],
                                    wspace=0.1,hspace=0.15,width_ratios=(15,15,1))
left_grid4 = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec=outer_grid[3,0],
                                    wspace=0.1,hspace=0.15,width_ratios=(15,1))
right_grid4 = gridspec.GridSpecFromSubplotSpec(2,3,subplot_spec=outer_grid[3,1],
                                    wspace=0.1,hspace=0.15,width_ratios=(15,15,1))

levs = [26,39]
imode = 2

cscale = [1e-6, 1e-6, 1000., 1e-7]
crange = np.array([ [-2.1, 2.1], [-0.9, 0.9],
                    [-2, 2], [-0.7, 0.7],
                    [0.07, 0.2], [-0.065,0.065],
                    [0.1, 1.1], [-0.29, 0.29]
 ])
clabel=[r'$\nabla\cdot\vec{v}$'+' (10$^{-6}$ s$^{-1}$)\nat 200 hPa',
        r'$\nabla\cdot\vec{v}$'+' (10$^{-6}$ s$^{-1}$)\nat 1000 hPa',
        'Cloud water path\n(kg m$^{-2}$)',
        'Precipitation\n(10$^{-7}$ m s$^{-1}$)']

overplot = False
overc = [ [ [-0.6,-0.3,0,0.3,0.6], [-0.4,-0.2,0.,0.2,0.4] ],
          [ [-0.3,-0.15,0,0.15,0.3], [-0.1,-0.05,0,0.05,0.1] ],
          [ [-0.5,-0.25,0,0.25,0.5], [-0.3,-0.15,0,0.15,0.3] ],
          [ [-0.3,-0.15,0,0.15,0.3], [-0.1,-0.05,0,0.05,0.1] ],
          [ [-0.04,-0.02,0,0.02,0.04], [-0.02,-0.01,0,0.01,0.02] ],
          [ [-0.02,-0.01,0,0.01,0.02], [-0.0024,-0.0012,0,0.0012,0.0024] ],
          [ [-0.2,-0.1,0,0.1,0.2], [-0.05,-0.025,0,0.025,0.05] ],
          [ [-0.08,-0.04,0,0.04,0.08], [-0.03,-0.015,0,0.015,0.03] ] ]
overfmt = [ [ '%.1f', '%.1f' ],
            [ '%.2f', '%.2f' ],
            [ '%.2f', '%.2f' ],
            [ '%.2f', '%.2f' ],
            [ '%.2f', '%.2f' ],
            [ '%.2f', '%.4f' ],
            [ '%.2f', '%.3f' ],
            [ '%.2f', '%.3f' ] ]

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

  #plot divergence
  out_file = parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'_DIVV_save.npz'

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

  for j in np.arange(len(levs)):
    ilev = levs[j]
    field_anom_mean = arc['field_anom_mean'][ilev]
    field_mean = arc['field_mean'][ilev]
    clm_mean = arc['clm_mean'][ilev]
    clm_anom = arc['clm_anom'][ilev]

    field_mean = smoothing_lon(field_mean,10)
    field_anom_mean = smoothing_lon(field_anom_mean,10)

    #full field
    if j == 0:
      ax = fig.add_subplot(left_grid1[2*i])
      if i == 0:
        ax.set_title('Average field',fontsize=7)
    else:
      ax = fig.add_subplot(left_grid2[2*i])
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
    m.drawparallels([-60,-30,0,30,60],labels = [True,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(lon2d, lat2d, field_mean/cscale[j], cmap='PiYG',rasterized=True,latlon='True',vmax=crange[2*j][1],vmin=crange[2*j][0])
    ax.set(ylabel='Latitude')
    ax.yaxis.set_label_coords(-0.15,0.5)
    xlim = ax.get_xlim()
    dxlim = xlim[1] - xlim[0]
    ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
    if i == 1:
      ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=7)
    else:
      ax.xaxis.set_ticklabels([])
    ax.tick_params(direction='in')
    ax.text(0.01,1.01,label[i],rotation=0,transform=ax.transAxes,fontsize=7,color='k',fontweight='bold')
    if j == 0:
      cax = fig.add_subplot(left_grid1[2*i+1])
    else:
      cax = fig.add_subplot(left_grid2[2*i+1])
    cbar = plt.colorbar(c,cax=cax)

    #anomaly
    if j == 0:
      ax = fig.add_subplot(right_grid1[3*i])
      if i == 0:
        ax.set_title('Average anomaly',fontsize=7)
    else:
      ax = fig.add_subplot(right_grid2[3*i])
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
    m.drawparallels([-60,-30,0,30,60],labels = [False,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(lon2d, lat2d, field_anom_mean/cscale[j], cmap='RdBu_r',rasterized=True,latlon='True',vmax=crange[2*j+1][1],vmin=crange[2*j+1][0])
    if overplot:
      #overplot contours
      #OMFG this is so annoying
      x, y = m(lon2d,lat2d)
      spliti = np.argmin(x[0])
      xtmp = np.hstack([x[:,spliti:],x[:,:spliti]])
      field_anom_tmp = np.hstack([field_anom_mean[:,spliti:],field_anom_mean[:,:spliti]])
      cline = m.contour(xtmp, y, field_anom_tmp/cscale[j], overc[2*j+i][0], colors='k', linewidths=1)
      ax.clabel(cline,inline=True,fontsize=6,fmt=overfmt[2*j+i][0])
    ax.yaxis.set_label_coords(-0.1,0.5)
    xlim = ax.get_xlim()
    dxlim = xlim[1] - xlim[0]
    ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
    if i == 1:
      ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=7)
    else:
      ax.xaxis.set_ticklabels([])
    ax.tick_params(direction='in')


    #spharm
    clm_anom_cp = clm_mean.copy()
    clm_anom_cp[:,:imode,:] = 0.0
    clm_anom_cp[:,imode+1:,:] = 0.0
    clm_anom_cp[:,:,:imode] = 0.0
    clm_anom_cp[:,:,imode+1:] = 0.0
    shmap_anom = np.real(sh.expand.MakeGridDH(clm_anom_cp,sampling=2))
    if j == 0:
      ax = fig.add_subplot(right_grid1[3*i+1])
      if i == 0:
        ax.set_title('Semidiurnal mode',fontsize=7)
    else:
      ax = fig.add_subplot(right_grid2[3*i+1])
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
    m.drawparallels([-60,-30,0,30,60],labels = [False,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(shlon2d, shlat2d, shmap_anom/cscale[j], cmap='RdBu_r',rasterized=True,latlon='True',vmax=crange[2*j+1][1],vmin=crange[2*j+1][0])
    if overplot:
      x, y = m(shlon2d,shlat2d)
      spliti = np.argmin(x[0])
      xtmp = np.hstack([x[:,spliti:],x[:,:spliti]])
      shmap_tmp = np.hstack([shmap_anom[:,spliti:],shmap_anom[:,:spliti]])
      cline = m.contour(xtmp, y, shmap_tmp/cscale[j], overc[2*j+i][1], colors='k', linewidths=1)
      ax.clabel(cline,inline=True,fontsize=6,fmt=overfmt[2*j+i][1])
    ax.yaxis.set_label_coords(-0.1,0.5)
    xlim = ax.get_xlim()
    dxlim = xlim[1] - xlim[0]
    ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
    if i == 1:
      ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=7)
    else:
      ax.xaxis.set_ticklabels([])
    ax.tick_params(direction='in')
    if j == 0:
      cax = fig.add_subplot(right_grid1[3*i+2])
    else:
      cax = fig.add_subplot(right_grid2[3*i+2])
    cbar = plt.colorbar(c,cax=cax)
    cbar.set_label(clabel[j],fontsize=7)


  #plot CWP
  out_file = parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'_CWP_save.npz'

  arc = np.load(out_file)
  field_anom_mean = arc['field_anom_mean']
  field_mean = arc['field_mean']
  clm_mean = arc['clm_mean']
  clm_anom = arc['clm_anom']

  field_mean = smoothing_lon(field_mean,10)
  field_anom_mean = smoothing_lon(field_anom_mean,10)

  ax = fig.add_subplot(left_grid3[2*i])
  m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
  m.drawparallels([-60,-30,0,30,60],labels = [True,False,False,False], fontsize=6)
  m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
  c = m.pcolormesh(lon2d, lat2d, field_mean/cscale[2], cmap='plasma_r',rasterized=True,latlon='True',vmax=crange[4][1],vmin=crange[4][0])
  ax.set(ylabel='Latitude')
  ax.yaxis.set_label_coords(-0.15,0.5)
  xlim = ax.get_xlim()
  dxlim = xlim[1] - xlim[0]
  ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
  if i == 1:
    ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=7)
  else:
    ax.xaxis.set_ticklabels([])
  ax.tick_params(direction='in')
  ax.text(0.01,1.01,label[i],rotation=0,transform=ax.transAxes,fontsize=7,color='k',fontweight='bold')
  cax = fig.add_subplot(left_grid3[2*i+1])
  cbar = plt.colorbar(c,cax=cax)


  ax = fig.add_subplot(right_grid3[3*i])
  m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
  m.drawparallels([-60,-30,0,30,60],labels = [False,False,False,False], fontsize=6)
  m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
  c = m.pcolormesh(lon2d, lat2d, field_anom_mean/cscale[2], cmap='RdBu_r',rasterized=True,latlon='True',vmax=crange[5][1],vmin=crange[5][0])
  #ax.set(ylabel='Latitude')
  if overplot:
    #overplot contours
    #OMFG this is so annoying
    x, y = m(lon2d,lat2d)
    spliti = np.argmin(x[0])
    xtmp = np.hstack([x[:,spliti:],x[:,:spliti]])
    field_anom_tmp = np.hstack([field_anom_mean[:,spliti:],field_anom_mean[:,:spliti]])
    cline = m.contour(xtmp, y, field_anom_tmp/cscale[2], overc[4+i][0], colors='k', linewidths=1)
    ax.clabel(cline,inline=True,fontsize=6,fmt=overfmt[4+i][0])
  ax.yaxis.set_label_coords(-0.1,0.5)
  xlim = ax.get_xlim()
  dxlim = xlim[1] - xlim[0]
  ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
  if i == 1:
    ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=7)
  else:
    ax.xaxis.set_ticklabels([])
  ax.tick_params(direction='in')

  #spharm
  clm_anom_cp = clm_mean.copy()
  clm_anom_cp[:,:imode,:] = 0.0
  clm_anom_cp[:,imode+1:,:] = 0.0
  clm_anom_cp[:,:,:imode] = 0.0
  clm_anom_cp[:,:,imode+1:] = 0.0
  shmap_anom = np.real(sh.expand.MakeGridDH(clm_anom_cp,sampling=2))
  ax = fig.add_subplot(right_grid3[3*i+1])
  m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
  m.drawparallels([-60,-30,0,30,60],labels = [False,False,False,False], fontsize=6)
  m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
  c = m.pcolormesh(shlon2d, shlat2d, shmap_anom/cscale[2], cmap='RdBu_r',rasterized=True,latlon='True',vmax=crange[5][1],vmin=crange[5][0])
  if overplot:
    x, y = m(shlon2d,shlat2d)
    spliti = np.argmin(x[0])
    xtmp = np.hstack([x[:,spliti:],x[:,:spliti]])
    shmap_tmp = np.hstack([shmap_anom[:,spliti:],shmap_anom[:,:spliti]])
    cline = m.contour(xtmp, y, shmap_tmp/cscale[2], overc[4+i][1], colors='k', linewidths=1)
    ax.clabel(cline,inline=True,fontsize=6,fmt=overfmt[4+i][1])
  ax.yaxis.set_label_coords(-0.1,0.5)
  ax.xaxis.set_ticks([-90,0,90])
  xlim = ax.get_xlim()
  dxlim = xlim[1] - xlim[0]
  ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
  if i == 1:
    ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=7)
  else:
    ax.xaxis.set_ticklabels([])
  ax.tick_params(direction='in')
  cax = fig.add_subplot(right_grid3[3*i+2])
  cbar = plt.colorbar(c,cax=cax)
  cbar.set_label(clabel[2],fontsize=7)


  #plot precip
  out_file = parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'_PRECT_save.npz'

  arc = np.load(out_file)
  field_anom_mean = arc['field_anom_mean']
  field_mean = arc['field_mean']
  clm_mean = arc['clm_mean']
  clm_anom = arc['clm_anom']

  field_mean = smoothing_lon(field_mean,10)
  field_anom_mean = smoothing_lon(field_anom_mean,10)

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

  ax = fig.add_subplot(left_grid4[2*i])
  m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
  m.drawparallels([-60,-30,0,30,60],labels = [True,False,False,False], fontsize=6)
  m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
  c = m.pcolormesh(lon2d, lat2d, field_mean/cscale[3], cmap='plasma_r',rasterized=True,latlon='True',vmax=crange[6][1],vmin=crange[6][0])
  ax.set(ylabel='Latitude')
  ax.yaxis.set_label_coords(-0.15,0.5)
  xlim = ax.get_xlim()
  dxlim = xlim[1] - xlim[0]
  ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
  if i == 1:
    ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=7)
  else:
    ax.xaxis.set_ticklabels([])
  ax.tick_params(direction='in')
  ax.text(0.01,1.01,label[i],rotation=0,transform=ax.transAxes,fontsize=7,color='k',fontweight='bold')
  cax = fig.add_subplot(left_grid4[2*i+1])
  cbar = plt.colorbar(c,cax=cax)
  cax.tick_params(axis='y')

  ax = fig.add_subplot(right_grid4[3*i])
  m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
  m.drawparallels([-60,-30,0,30,60],labels = [False,False,False,False], fontsize=6)
  m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
  c = m.pcolormesh(lon2d, lat2d, field_anom_mean/cscale[3], cmap='RdBu_r',rasterized=True,latlon='True',vmax=crange[7][1],vmin=crange[7][0])
  if overplot:
    #overplot contours
    #OMFG this is so annoying
    x, y = m(lon2d,lat2d)
    spliti = np.argmin(x[0])
    xtmp = np.hstack([x[:,spliti:],x[:,:spliti]])
    field_anom_tmp = np.hstack([field_anom_mean[:,spliti:],field_anom_mean[:,:spliti]])
    cline = m.contour(xtmp, y, field_anom_tmp/cscale[3], overc[6+i][0], colors='k', linewidths=1)
    ax.clabel(cline,inline=True,fontsize=6,fmt=overfmt[6+i][0])
  ax.yaxis.set_label_coords(-0.1,0.5)
  xlim = ax.get_xlim()
  dxlim = xlim[1] - xlim[0]
  ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
  if i == 1:
    ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=7)
  else:
    ax.xaxis.set_ticklabels([])
  ax.tick_params(direction='in')


  #spharm
  for imode in [2,]:
    clm_anom_cp = clm_mean.copy()
    clm_anom_cp[:,:imode,:] = 0.0
    clm_anom_cp[:,imode+1:,:] = 0.0
    clm_anom_cp[:,:,:imode] = 0.0
    clm_anom_cp[:,:,imode+1:] = 0.0
    shmap_anom = np.real(sh.expand.MakeGridDH(clm_anom_cp,sampling=2))
    ax = fig.add_subplot(right_grid4[3*i+1])
    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection=proj)
    m.drawparallels([-60,-30,0,30,60],labels = [False,False,False,False], fontsize=6)
    m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)
    c = m.pcolormesh(shlon2d, shlat2d, shmap_anom/cscale[3], cmap='RdBu_r',rasterized=True,latlon='True',vmax=crange[7][1],vmin=crange[7][0])
    if overplot:
      #overplot contours
      #OMFG this is so annoying
      x, y = m(shlon2d,shlat2d)
      spliti = np.argmin(x[0])
      xtmp = np.hstack([x[:,spliti:],x[:,:spliti]])
      shmap_tmp = np.hstack([shmap_anom[:,spliti:],shmap_anom[:,:spliti]])
      cline = m.contour(xtmp, y, shmap_tmp/cscale[3], overc[6+i][1], colors='k', linewidths=1)
      ax.clabel(cline,inline=True,fontsize=6,fmt=overfmt[6+i][1])
    ax.yaxis.set_label_coords(-0.1,0.5)
    ax.xaxis.set_ticks([-90,0,90])
    xlim = ax.get_xlim()
    dxlim = xlim[1] - xlim[0]
    ax.xaxis.set_ticks([xlim[0]+0.25*dxlim,xlim[0]+0.5*dxlim,xlim[0]+0.75*dxlim])
    if i == 1:
      ax.xaxis.set_ticklabels(['Sunrise','Noon','Sunset'],fontsize=7)
    else:
      ax.xaxis.set_ticklabels([])
    ax.tick_params(direction='in')
  cax = fig.add_subplot(right_grid4[3*i+2])
  cbar = plt.colorbar(c,cax=cax)
  cbar.set_label(clabel[3],fontsize=7)
  cax.tick_params(axis='y')

plt.savefig('figures/figure_4.pdf')
plt.close()
