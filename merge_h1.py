from cdo import *
cdo = Cdo()
import argparse
import pathlib

#This script merges the original CAM history files (hourly data) into one file.
#It is included here only for record keeping as the merged history files are
#those preserved in the repository. File paths are specific to RD's machines and
#directories.

parser = argparse.ArgumentParser()
parser.add_argument('sim_name', nargs='*', help='name of simulation')
parser.add_argument('-i','--init_year', nargs=1, default=[61], type=int, help='first year of averaging')
parser.add_argument('-l','--last_year', nargs=1, default=[61], type=int, help='last year of averaging')
parser.add_argument('-w','--overwrite', action='store_true', help = 'force overwrite of output files')

args = parser.parse_args()

start_year = args.init_year[0]
end_year = args.last_year[0]

archive_path = pathlib.Path('/kenes/data/rdeitrick/ccsm/archive')
atm_hist_files = archive_path / args.sim_name[0] / 'atm' / 'hist'

#build up list of files to merge
list_files = ''
for i in range(start_year,end_year+1):
  for j in range(1,13):
    for day in range(1,32):
      day_file = str(atm_hist_files / args.sim_name[0]) + '.cam.h1.%04d-%02d-%02d-03600.nc'%(i,j,day)
      if pathlib.Path(day_file).exists():
        list_files = list_files + day_file + ' '

output_path = pathlib.Path('/home/rdeitrick/waccm_rot_sims/'+args.sim_name[0]) / 'merged_hist'
if not output_path.exists():
  output_path.mkdir(parents=True,exist_ok=True)

#merge files
output_merged = output_path / (args.sim_name[0] + '.cam.h1.merged_%04d_%04d.nc'%(start_year,end_year))
if not output_merged.exists() or args.overwrite:
  cdo.mergetime(input=list_files,output=str(output_merged))

#global mean time evolution
output_timeevol = output_path / (args.sim_name[0] + '.cam.h1.timeevol_%04d_%04d.nc'%(start_year,end_year))
if not output_timeevol.exists() or args.overwrite:
  cdo.fldmean(input=str(output_merged),output=str(output_timeevol))

#compute daily averages
output_daymean = output_path / (args.sim_name[0] + '.cam.h1.daymean_%04d_%04d.nc'%(start_year,end_year))
if not output_daymean.exists() or args.overwrite:
  cdo.ydaymean(input=str(output_merged),output=str(output_daymean))
