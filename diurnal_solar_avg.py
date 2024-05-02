import numpy as np
import netCDF4 as nc
import pdb
import pathlib
import warnings
import os
import scipy.interpolate as sint
import pyshtools as sh
from windspharm.standard import VectorWind
#suppress deprecation warnings coming from netCDF4
warnings.filterwarnings("ignore", category=DeprecationWarning)

#Change this as needed to point to the 'main' simulations
parent_path = '../sims_main/'

#Create a separate folder so that old output isn't overwritten
output_path = '../sims_main_new/'

#list of simulations here, with rotperiod (units of modern day), and heat capacity
simlist = [
#           {'name': 'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr_branch2',
#               'rotper': 0.9167, 'cp': 1034.93},
           {'name': 'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch2',
                'rotper': 1.0, 'cp': 1034.93},
          ]

#These lists tell what quantities to use in averaging process
#[name in CAM (some exceptions to this),
# desired output name (need to match plot code),
# dimensions (2 or 3),
# whether to recenter solar longitude using FSDS output (True or False),
# geography setting ('all, 'ocean', or 'land')]
# NOTE that 3-D fields use a lot of memory! (20 GB)
# ALSO, running all of these at once will take a very long time
fields = [
#          ['FSDS','fsds_phase',2, False, 'all'],
#          ['PS', 'p_anom', 2, True, 'all'],
#          ['TGCLDCWP', 'CWP', 2, True, 'all'],
#          ['PRECT', 'PRECT', 2, True, 'all'],
#          ['PRECT', 'PRECT', 2, True, 'land'],
#          ['PRECT', 'PRECT', 2, True, 'ocean'],
#          ['DIVV', 'DIVV', 3, True, 'all'], #special case, need to calc from U and V
          ['QRS', 'QRS', 3, True, 'all'],   #special case, convert units
          ['CMFMC', 'CMFMC', 3, True, 'all'],
          ['CMFMCDZM', 'CMFMCDZM', 3, True, 'all'],
          ['ZMDT', 'ZMDT', 3, True, 'all'],
          ['T', 'T', 3, True, 'all'],
          ['RELHUM', 'RELHUM', 3, True, 'all'],
          ['CWC', 'CWC', 3, True, 'all'],   #special case, use IWC + LWC and take log
          ['Q', 'log_Q', 3, True, 'all']    #special case, take log of Q
          ]


def dasl(sim, cam_field_name, out_field_name, ndim, recenter=True, geo='all'):
    #diurnal averaging about solar longitude for 2D field
    #need to have SW flux file already calculated to center the sun at long = 0

    simname = sim['name']
    rotper = sim['rotper']

    #find the merged file
    path = parent_path + simname + '/merged_hist/'
    files_list = [filename for filename in os.listdir(path) if 'merged' in filename]
    files_list.sort()
    path_merge = pathlib.Path(path) / files_list[0]

    print('Reading file %s...'%path_merge)
    print('Processing %s...'%out_field_name)
    data = nc.Dataset(path_merge,'r')

    #make sure output path exists
    out_path_full = pathlib.Path(output_path + simname + '/merged_hist')
    if not out_path_full.exists():
        out_path_full.mkdir(parents=True)

    if geo == 'land' or geo == 'ocean':
        out_file = np.str(out_path_full / (simname + '_' + geo + '_' + out_field_name + '_save.npz'))
        #need to get land fraction, just take it from any merged 60yr history file
        landfrac = nc.Dataset(parent_path +
                    'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr/merged_hist/'+
                    'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr.cam.h0.merged_0001_0060.nc','r')['LANDFRAC'][:].squeeze()
    else:
        out_file = np.str(out_path_full / (simname + '_' + out_field_name + '_save.npz'))

    #now get the necessary fields
    lat = data['lat'][:]
    lon = data['lon'][:]
    p = data['lev'][:]

    #some quantities need to be computed or modified in some why, handle these here
    if cam_field_name == 'DIVV':
        ftime = len(data['time'][:])

        field0 = np.zeros_like(data['U'][:])
        for itime in np.arange(ftime):
            u = data['U'][itime]
            v = data['V'][itime]
            u1 = np.swapaxes(np.swapaxes(u,0,1),1,2)[::-1,:,:]
            v1 = np.swapaxes(np.swapaxes(v,0,1),1,2)[::-1,:,:]
            V = VectorWind(u1,v1)
            divV = V.divergence()
            divV1 = np.swapaxes(np.swapaxes(divV,1,2),0,1)[:,::-1,:]
            field0[itime] = divV1

        del u, v, u1, V, divV, divV1

    elif cam_field_name == 'QRS':
        #heating rate, converted to W / kg (CAM/WACCM) or J/day/kg (ExoCAM)
        #code for figures 12 and 13 will correct to W/kg for ExoCAM
        field0 = data[cam_field_name][:] * sim['cp']

    elif cam_field_name == 'CWC':
        #total cloud water content, in log space
        field0 = np.log10(data['IWC'][:] + data['LWC'][:])

    elif cam_field_name == 'Q':
        field0 = np.log10(data['Q'][:])

    else:
        field0 = data[cam_field_name][:]

    #need to compute averages over rotation periods
    t0 = data['time'][0]
    #annual average
    field0_mean = np.mean(field0,axis=0)

    #pressure anomaly as function of time, lat, lon
    if ndim == 2:
        field_anom  = field0 - field0_mean[None,:,:] #using annual mean instead of diurnal
    elif ndim == 3:
        field_anom  = field0 - field0_mean[None,:,:,:] #using annual mean instead of diurnal
    del field0_mean  #just trying to beat down memory usage

    # rotation rate in rad/Earth solar day
    rotrate = 2*np.pi/rotper

    #longitude of sun as a function of time
    offset = np.mean(lon[np.argmax(data['FSDS'][:][0,6:-6,:],axis=1)])

    if recenter == True:
        #additional correction from fsds file
        fsds_file = np.str(out_path_full / (simname + '_fsds_phase_save.npz'))
        arc2 = np.load(fsds_file)
        fsds_mean = arc2['field_anom_mean']
        clm_fsds = arc2['clm_anom']
        offset2 = np.arctan2(clm_fsds[1,1,1],clm_fsds[0,1,1])*180/np.pi
    else:
        offset2 = 0.0

    lon_of_sun = (data['lon'][:][None,:] + (data['time'][:][:,None]-t0)*rotrate*180/np.pi-offset-offset2)%360

    ftime = len(data['time'][:])

    field_reorder = np.zeros_like(field0)
    field_anom_reorder = np.zeros_like(field_anom)
    if geo == 'land' or geo == 'ocean':
        landfrac_reorder = np.zeros_like(field0)

    for itime in np.arange(ftime):
      if ndim == 3:
          for ilev in np.arange(len(p)):
              interp = sint.interp2d(lon_of_sun[itime],lat,field0[itime,ilev,:,:])
              field_reorder[itime,ilev,:,:] = interp(lon,lat)
              interp = sint.interp2d(lon_of_sun[itime],lat,field_anom[itime,ilev,:,:])
              field_anom_reorder[itime,ilev,:,:] = interp(lon,lat)
      else:
          interp = sint.interp2d(lon_of_sun[itime],lat,field0[itime,:,:])
          field_reorder[itime,:,:] = interp(lon,lat)
          interp = sint.interp2d(lon_of_sun[itime],lat,field_anom[itime,:,:])
          field_anom_reorder[itime,:,:] = interp(lon,lat)
      if geo == 'land' or geo == 'ocean':
          interp = sint.interp2d(lon_of_sun[itime],lat,landfrac[0,:,:])
          landfrac_reorder[itime,:,:] = interp(lon,lat)

    del field0, field_anom   #keeping memory down
    lon2d, lat2d = np.meshgrid(lon,lat)
    lmax = np.floor(np.sqrt(len(lat)*len(lon))/2)

    #below geo options not written to handle 3-D fields
    if geo == 'land':
        #quick and dirty way to mask ocean areas
        field_mean = field_reorder * landfrac_reorder
        field_mean[landfrac_reorder==0] = np.nan
        field_anom_mean = field_anom_reorder * landfrac_reorder
        field_anom_mean[landfrac_reorder==0] = np.nan
        clm_mean = np.zeros((2,np.int(lmax+1),np.int(lmax+1)))
        clm_anom = np.zeros((2,np.int(lmax+1),np.int(lmax+1)))
    elif geo == 'ocean':
        #quick and dirty way to mask land areas
        field_mean = field_reorder * (1.0 - landfrac_reorder)
        field_mean[landfrac_reorder==1] = np.nan
        field_anom_mean = field_anom_reorder * (1.0 - landfrac_reorder)
        field_anom_mean[landfrac_reorder==1] = np.nan
        clm_mean = np.zeros((2,np.int(lmax+1),np.int(lmax+1)))
        clm_anom = np.zeros((2,np.int(lmax+1),np.int(lmax+1)))
    else:
        if ndim == 3:
            field_mean = np.mean(field_reorder[:ftime,:,:,:], axis=0)
            field_anom_mean = np.mean(field_anom_reorder[:ftime,:,:,:], axis=0)
            del field_reorder, field_anom_reorder
            clm_mean = np.zeros((len(p),2,np.int(lmax+1),np.int(lmax+1)))
            clm_anom = np.zeros((len(p),2,np.int(lmax+1),np.int(lmax+1)))
            for ilev in np.arange(len(p)):
              clm_mean[ilev], chi2 = sh.expand.SHExpandLSQ(field_mean[ilev,:,:].ravel(),lat2d.ravel(),lon2d.ravel(),lmax)
              clm_anom[ilev], chi2 = sh.expand.SHExpandLSQ(field_anom_mean[ilev,:,:].ravel(),lat2d.ravel(),lon2d.ravel(),lmax)
        else:
            field_mean = np.nanmean(field_reorder[:ftime,:,:], axis=0)
            field_anom_mean = np.nanmean(field_anom_reorder[:ftime,:,:], axis=0)
            del field_reorder, field_anom_reorder
            clm_mean, chi2 = sh.expand.SHExpandLSQ(field_mean.ravel(),lat2d.ravel(),lon2d.ravel(),lmax)
            clm_anom, chi2 = sh.expand.SHExpandLSQ(field_anom_mean.ravel(),lat2d.ravel(),lon2d.ravel(),lmax)

    np.savez(out_file,lon=lon,lat=lat,rotrate=rotrate,field_anom_mean=field_anom_mean,field_mean=field_mean, clm_mean=clm_mean,clm_anom=clm_anom,lmax=lmax,p=p)


for ifield in np.arange(len(fields)):
    for isim in np.arange(len(simlist)):
        dasl(simlist[isim], fields[ifield][0], fields[ifield][1],
              fields[ifield][2], recenter = fields[ifield][3], geo = fields[ifield][4])
