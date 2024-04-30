import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pdb
import pathlib
import warnings
import os
# from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.interpolate as sint
import pyshtools as sh
warnings.filterwarnings("ignore", category=DeprecationWarning)

#Change this as needed to point to the 'main' simulations
parent_path = '../sims_main/'

#Create a separate folder so that old output isn't overwritten
output_path = '../sims_main_new/'

#prob turn this into a dictionary with the rotation rate built in
# simnames =  [
#             'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr_branch2',
#             'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch2'
#             ]
#
# rotpers = np.array([0.9167, 1.0])

simlist = {
           'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr_branch2': 0.9167,
           'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch2': 1.0
           }


fields = [
          ['TGCLDCWP', 'CWP', 2],
          ]


def dasl_2D(simname, rotper, cam_field_name, out_field_name):
    #diurnal averaging about solar longitude for 2D field
    #need to have SW flux file already calculated to center the sun at long = 0

    #find the merged file
    path = parent_path + simname + '/merged_hist/'
    files_list = [filename for filename in os.listdir(path) if 'merged' in filename]
    files_list.sort()
    path_merge = pathlib.Path(path) / files_list[0]

    print('Reading file %s...'%path_merge)
    data = nc.Dataset(path_merge,'r')

    #make sure output path exists
    out_path_full = pathlib.Path(output_path + simname + '/merged_hist')
    if not out_path_full.exists():
        out_path_full.mkdir(parents=True)

    out_file = np.str(out_path_full / (simname + '_' + out_field_name + '_save.npz'))

    #now get the necessary fields
#    z3 = data['Z3'][:]
    lat = data['lat'][:]
    lon = data['lon'][:]

    field0 = data[cam_field_name][:]
#    field0 = simps(-1.0*liwc,z3,axis=1)

    #need to compute averages over rotation periods
    t0 = data['time'][0]
    #try using annual average
    field0_mean = np.mean(field0,axis=0)

    #pressure anomaly as function of time, lat, lon
    field_anom  = field0 - field0_mean[None,:,:] #using annual mean instead of diurnal

    # rotation rate in rad/Earth solar day
    rotrate = 2*np.pi/rotper

    #longitude of sun as a function of time
    #for ilat in np.arange(len(lat)):
    offset = np.mean(lon[np.argmax(data['FSDS'][:][0,6:-6,:],axis=1)])
    #print(offset)
    #offset = 0.0 #180.0 + ((t0+1) * rotrate * 180/np.pi)%360
    #additional correction from fsds file
    fsds_file = path + simname + '_fsds_phase_save.npz'
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


for ifield in np.arange(len(fields)):
    for simname in simlist.keys():
        if fields[ifield][2] == 2:
            dasl_2D(simname, simlist[simname], fields[ifield][0], fields[ifield][1])
