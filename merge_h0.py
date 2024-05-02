from cdo import *
cdo = Cdo()
import argparse
import pathlib

#This script merges the original CAM history files (monthly averages) into one file.
#It is included here only for record keeping as the merged history files are
#those preserved in the repository. File paths are specific to RD's machines and
#directories.

parser = argparse.ArgumentParser()
parser.add_argument('sim_name', nargs='*', help='name of simulation')
parser.add_argument('-i','--init_year', nargs=1, default=[1], type=int, help='first year of averaging')
parser.add_argument('-l','--last_year', nargs=1, default=[5], type=int, help='last year of averaging')
parser.add_argument('-w','--overwrite', action='store_true', help = 'force overwrite of output files')

args = parser.parse_args()

start_year = args.init_year[0]
end_year = args.last_year[0]

#path to history files in ccsm/archive
archive_path = pathlib.Path('/kenes/data/rdeitrick/ccsm/archive')

atm_hist_files = archive_path / args.sim_name[0] / 'atm' / 'hist'

#build up list of files to merge
list_files = ''
for i in range(start_year,end_year+1):
  for j in range(1,13):
    list_files = list_files + str(atm_hist_files / args.sim_name[0]) + '.cam.h0.%04d-%02d.nc'%(i,j) + ' '

#check paths
output_path = pathlib.Path('/home/rdeitrick/waccm_rot_sims/'+args.sim_name[0]) / 'merged_hist'
if not output_path.exists():
  output_path.mkdir(parents=True,exist_ok=True)

#merge the history files
output_merged = output_path / (args.sim_name[0] + '.cam.h0.merged_%04d_%04d.nc'%(start_year,end_year))
if not output_merged.exists() or args.overwrite:
  cdo.mergetime(input=list_files,output=str(output_merged))

#global mean time evolution
output_timeevol = output_path / (args.sim_name[0] + '.cam.h0.timeevol_%04d_%04d.nc'%(start_year,end_year))
if not output_timeevol.exists() or args.overwrite:
  cdo.fldmean(input='-select,name=T,TS,TBOT,FLUT,FSNT,FLNT,ozone '+str(output_merged),output=str(output_timeevol))

#annual average
output_annavg = output_path / (args.sim_name[0] + '.cam.h0.annavg_%04d_%04d.nc'%(start_year,end_year))
if not output_annavg.exists() or args.overwrite:
  cdo.yearmonmean(input='-select,name=T,TS,TBOT,FLUT,FSNT,FLNT '+str(output_timeevol),output=str(output_annavg))
