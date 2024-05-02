from cdo import *
cdo = Cdo()
import argparse
import pathlib
import pdb

#pdb.set_trace()

parser = argparse.ArgumentParser()
parser.add_argument('sim_name', nargs='*', help='name of simulation')
parser.add_argument('-i','--init_year', nargs=1, default=[1], type=int, help='first year of averaging')
parser.add_argument('-l','--last_year', nargs=1, default=[5], type=int, help='last year of averaging')
parser.add_argument('-w','--overwrite', action='store_true', help = 'force overwrite of output files')

args = parser.parse_args()

start_year = args.init_year[0]
end_year = args.last_year[0]

archive_path = pathlib.Path('/kenes/data/rdeitrick/ccsm/archive')
#archive_path = pathlib.Path('~/data-rdeitrick/ccsm/archive')
atm_hist_files = archive_path / args.sim_name[0] / 'atm' / 'hist'

#build up list of files to merge
list_files = ''
for i in range(start_year,end_year+1):
  for j in range(1,13):
    list_files = list_files + str(atm_hist_files / args.sim_name[0]) + '.cam.h0.%04d-%02d.nc'%(i,j) + ' '

#pdb.set_trace()
output_path = pathlib.Path('/home/rdeitrick/waccm_rot_sims/'+args.sim_name[0]) / 'merged_hist'
if not output_path.exists():
  output_path.mkdir(parents=True,exist_ok=True)

#merge files
output_merged = output_path / (args.sim_name[0] + '.cam.h0.merged_%04d_%04d.nc'%(start_year,end_year))

if not output_merged.exists() or args.overwrite:
  cdo.mergetime(input=list_files,output=str(output_merged))

#calc time mean from merged file
#selyear = '-selyear,%02d/%02d '%(start_year,end_year)
#output_timmean = output_path / (args.sim_name[0] + '.cam.h0.timmean_%04d_%04d.nc'%(start_year,end_year))

#if not output_timmean.exists() or args.overwrite:
#  cdo.timmean(input=selyear+str(output_merged),output=str(output_timmean))

#lts and eis files---------------------------------
#output_merged_lts = output_path / (args.sim_name[0] + '.cam.h0.merged_%04d_%04d_lts.nc'%(1,end_year))
#output_timmean_lts = output_path / (args.sim_name[0] + '.cam.h0.timmean_%04d_%04d_lts.nc'%(start_year,end_year))
#if not output_timmean.exists() or args.overwrite:
#cdo.timmean(input=selyear+str(output_merged_lts),output=str(output_timmean_lts))
#--------
#output_merged_eis = output_path / (args.sim_name[0] + '.cam.h0.merged_%04d_%04d_eis.nc'%(1,end_year))
#output_timmean_eis = output_path / (args.sim_name[0] + '.cam.h0.timmean_%04d_%04d_eis.nc'%(start_year,end_year))
#if not output_timmean.exists() or args.overwrite:
#cdo.timmean(input=selyear+str(output_merged_eis),output=str(output_timmean_eis))

#calc seasonal means from merged file
#output_djfmean = output_path / (args.sim_name[0] + '.cam.h0.djfmean_%04d_%04d.nc'%(start_year,end_year))
#output_jjamean = output_path / (args.sim_name[0] + '.cam.h0.jjamean_%04d_%04d.nc'%(start_year,end_year))

#if not output_djfmean.exists() or args.overwrite: #use JFM instead of DJF because of date stamp of avg occuring on 1st of next month
#  cdo.timmean(input=selyear+'-select,season=JFM '+str(output_merged),output=str(output_djfmean))

#if not output_jjamean.exists() or args.overwrite: #use JAS instead of JJA (see above)
#  cdo.timmean(input=selyear+'-select,season=JAS '+str(output_merged),output=str(output_jjamean))

#calc zonal mean from time mean file
#output_zonmean = output_path / (args.sim_name[0] + '.cam.h0.zonmean_%04d_%04d.nc'%(start_year,end_year))
#if not output_zonmean.exists() or args.overwrite:
#  cdo.zonmean(input='-select,name=CLOUD,TS,T,Q,OMEGA,CLDICE,CLDLIQ,PCLDTOP,PCLDBOT,IWC,LWC,RELHUM,LWCF,SWCF,TMQ,FLNS,FLNSC,FLUT,FLUTC,FSDS,FSDSC,FSNT,FSNTC,FSNTOA,FLNT,U,FLNTC,V '+str(output_timmean),output=str(output_zonmean))

#calc global means?
#output_globmean = output_path / (args.sim_name[0] + '.cam.h0.globmean_%04d_%04d.nc'%(start_year,end_year))
#if not output_globmean.exists() or args.overwrite:
#cdo.fldmean(input='-select,name=CLDICE,CLDLIQ,TCLDT,ICLDT,TCLDP,ICLDP,LCLDT,LCLDP,ICLDTWP,PMID,TS,TGCLDCWP,TGCLDIWP,TGCLDLWP,T,Q,LWCF,SWCF,ozone,cb_ozone_c,IWC,LWC,PS,AREI,AREL,CLOUD,CLDHGH,CLDMED,CLDLOW,RELHUM,FSNS,FSDS,CLDTOT '+str(output_timmean),output=str(output_globmean))

#output_globmean_lts = output_path / (args.sim_name[0] + '.cam.h0.globmean_%04d_%04d_lts.nc'%(start_year,end_year))
#if not output_globmean.exists() or args.overwrite:
#cdo.fldmean(input=str(output_timmean_lts),output=str(output_globmean_lts))

#output_globmean_eis = output_path / (args.sim_name[0] + '.cam.h0.globmean_%04d_%04d_eis.nc'%(start_year,end_year))
#if not output_globmean.exists() or args.overwrite:
#cdo.fldmean(input=str(output_timmean_eis),output=str(output_globmean_eis))

#need to do some annoying things to deal with cloud top T
#output_timmean_cldt = output_path / (args.sim_name[0] + '.cam.h0.timmean_%04d_%04d_cldt.nc'%(start_year,end_year))
#cdo.setmissval('nan',input='-select,name=ICLDT,LCLDT,TCLDT,LCLDP,ICLDP,TCLDP '+str(output_timmean),output=str(output_timmean_cldt))
#output_globmean_cldt = output_path / (args.sim_name[0] + '.cam.h0.globmean_%04d_%04d_cldt.nc'%(start_year,end_year))
#cdo.fldmean(input=str(output_timmean_cldt),output=str(output_globmean_cldt))

#tropical mean
#output_tropmean = output_path / (args.sim_name[0] + '.cam.h0.tropmean_%04d_%04d.nc'%(start_year,end_year))
#if not output_tropmean.exists() or args.overwrite:
#  cdo.fldmean(input='-select,name=TS,T,Q,ozone,IWC,LWC,PS,AREI,AREL,CLOUD,CLDHGH,CLDMED,CLDLOW,RELHUM,FSNS,FSDS -sellonlatbox,0,360,-15,15 '+str(output_timmean),output=str(output_tropmean))

#south polar mean
#output_spolmean = output_path / (args.sim_name[0] + '.cam.h0.spolmean_%04d_%04d.nc'%(start_year,end_year))
#if not output_spolmean.exists() or args.overwrite:
#  cdo.fldmean(input='-select,name=TS,T,Q,ozone,IWC,LWC,PS,AREI,AREL,CLOUD,CLDHGH,CLDMED,CLDLOW,RELHUM,FSNS,FSDS -sellonlatbox,0,360,-90,-60 '+str(output_timmean),output=str(output_spolmean))

#north polar mean
#output_npolmean = output_path / (args.sim_name[0] + '.cam.h0.npolmean_%04d_%04d.nc'%(start_year,end_year))
#if not output_npolmean.exists() or args.overwrite:
#  cdo.fldmean(input='-select,name=TS,T,Q,ozone,IWC,LWC,PS,AREI,AREL,CLOUD,CLDHGH,CLDMED,CLDLOW,RELHUM,FSNS,FSDS -sellonlatbox,0,360,60,90 '+str(output_timmean),output=str(output_npolmean))

#global mean time evolution3
output_timeevol = output_path / (args.sim_name[0] + '.cam.h0.timeevol_%04d_%04d.nc'%(start_year,end_year))
if not output_timeevol.exists() or args.overwrite:
  cdo.fldmean(input='-select,name=T,TS,TBOT,FLUT,FSNT,FLNT,ozone '+str(output_merged),output=str(output_timeevol))

#annual average
output_annavg = output_path / (args.sim_name[0] + '.cam.h0.annavg_%04d_%04d.nc'%(start_year,end_year))
if not output_annavg.exists() or args.overwrite:
  cdo.yearmonmean(input='-select,name=T,TS,TBOT,FLUT,FSNT,FLNT '+str(output_timeevol),output=str(output_annavg))

#output_flut = output_path / (args.sim_name[0] + '.cam.h0.flut_%04d_%04d.nc'%(1,end_year))
#if not output_timeevol.exists() or args.overwrite:
#  cdo.select(input='name=FLUT '+str(output_merged),output=str(output_flut))

#vertical integrals
#output_vertsum = output_path / (args.sim_name[0] + '.cam.h0.vertsum_%04d_%04d.nc'%(start_year,end_year))
#if not output_vertsum.exists() or args.overwrite:
#  cdo.vertsum(input='-select,name=LWC,IWC '+str(output_timmean),output=str(output_vertsum))
