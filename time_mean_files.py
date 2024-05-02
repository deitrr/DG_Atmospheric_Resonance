from cdo import *
cdo = Cdo()
import argparse
import pathlib

#This script post-processes some average quantities from the merged history
#files containing monthly data.
#It is included here only for record keeping. File paths are specific to RD's machines and
#directories.

parser = argparse.ArgumentParser()
parser.add_argument('sim_name', nargs='*', help='name of simulation')
parser.add_argument('-i','--init_year', nargs=1, default=[31], type=int, help='first year of averaging')
parser.add_argument('-l','--last_year', nargs=1, default=[60], type=int, help='last year of averaging')
parser.add_argument('-w','--overwrite', action='store_true', help = 'force overwrite of output files')

args = parser.parse_args()

start_year = args.init_year[0]
end_year = args.last_year[0]

archive_path = pathlib.Path('/kenes/user/rdeitrick/waccm_rot_sims/ccsm/archive')
atm_hist_files = archive_path / args.sim_name[0] / 'atm' / 'hist'

output_path = pathlib.Path('/home/rdeitrick/waccm_rot_sims/' + args.sim_name[0]) / 'merged_hist'
if not output_path.exists():
  output_path.mkdir(parents=True,exist_ok=True)

#set path of merged file years 1 to end_year
output_merged = output_path / (args.sim_name[0] + '.cam.h0.merged_%04d_%04d.nc'%(1,end_year))
if not output_merged.exists():
  output_merged = output_path / (args.sim_name[0] + '.cam.h0.merged_%04d_%04d.nc'%(start_year,end_year))

#calc time mean from merged file
selyear = '-selyear,%02d/%02d '%(start_year,end_year)
output_timmean = output_path / (args.sim_name[0] + '.cam.h0.timmean_%04d_%04d.nc'%(start_year,end_year))
print(str(output_timmean))
if not output_timmean.exists() or args.overwrite:
  cdo.timmean(input=selyear+str(output_merged),output=str(output_timmean))

#calc zonal mean from time mean file
output_zonmean = output_path / (args.sim_name[0] + '.cam.h0.zonmean_%04d_%04d.nc'%(start_year,end_year))
print(str(output_zonmean))
if not output_zonmean.exists() or args.overwrite:
  cdo.zonmean(input='-select,name=CLOUD,TS,T,Q,OMEGA,CLDICE,CLDLIQ,PCLDTOP,PCLDBOT,IWC,LWC,RELHUM,LWCF,SWCF,TMQ,FLNS,FLNSC,FLUT,FLUTC,FSDS,FSDSC,FSNT,FSNTC,FSNTOA,FLNT,U,FLNTC,V '+str(output_timmean),output=str(output_zonmean))

#calc global means?
output_globmean = output_path / (args.sim_name[0] + '.cam.h0.globmean_%04d_%04d.nc'%(start_year,end_year))
print(str(output_globmean))
if not output_globmean.exists() or args.overwrite:
  cdo.fldmean(input=str(output_timmean),output=str(output_globmean))

#tropical mean
output_tropmean = output_path / (args.sim_name[0] + '.cam.h0.tropmean_%04d_%04d.nc'%(start_year,end_year))
print(str(output_tropmean))
if not output_tropmean.exists() or args.overwrite:
  cdo.fldmean(input='-select,name=CLOUD,TS,T,Q,OMEGA,CLDICE,CLDLIQ,PCLDTOP,PCLDBOT,IWC,LWC,RELHUM,LWCF,SWCF,TMQ,FLNS,FLNSC,FLUT,FLUTC,FSDS,FSDSC,FSNT,FSNTC,FSNTOA,FLNT,U,FLNTC,V -sellonlatbox,0,360,-15,15 '+str(output_timmean),output=str(output_tropmean))

#south polar mean
output_spolmean = output_path / (args.sim_name[0] + '.cam.h0.spolmean_%04d_%04d.nc'%(start_year,end_year))
print(str(output_spolmean))
if not output_spolmean.exists() or args.overwrite:
  cdo.fldmean(input='-select,name=CLOUD,TS,T,Q,OMEGA,CLDICE,CLDLIQ,PCLDTOP,PCLDBOT,IWC,LWC,RELHUM,LWCF,SWCF,TMQ,FLNS,FLNSC,FLUT,FLUTC,FSDS,FSDSC,FSNT,FSNTC,FSNTOA,FLNT,U,FLNTC,V -sellonlatbox,0,360,-90,-60 '+str(output_timmean),output=str(output_spolmean))

#north polar mean
output_npolmean = output_path / (args.sim_name[0] + '.cam.h0.npolmean_%04d_%04d.nc'%(start_year,end_year))
print(str(output_npolmean))
if not output_npolmean.exists() or args.overwrite:
  cdo.fldmean(input='-select,name=CLOUD,TS,T,Q,OMEGA,CLDICE,CLDLIQ,PCLDTOP,PCLDBOT,IWC,LWC,RELHUM,LWCF,SWCF,TMQ,FLNS,FLNSC,FLUT,FLUTC,FSDS,FSDSC,FSNT,FSNTC,FSNTOA,FLNT,U,FLNTC,V -sellonlatbox,0,360,60,90 '+str(output_timmean),output=str(output_npolmean))
