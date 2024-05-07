import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pathlib
import warnings
import pyshtools as sh
import matplotlib.gridspec as gridspec
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':7})

warnings.filterwarnings("ignore", category=DeprecationWarning)

#change these paths to match the location of the sim_main and sim_suppl
# download from the data repository
parent_path = '../DG_Atmospheric_Resonance_archive/sims_main/'
suppl_path = '../DG_Atmospheric_Resonance_archive/sims_suppl/'

#list of simulations to put on this plot (30 year averages)
#yes, this is a big mess
simnamesT = [
            'lr_cam4_4x5_modern_16hr', #0
            'lr_cam4_4x5_modern_18hr',
            'lr_cam4_4x5_modern_20hr',
            'lr_cam4_4x5_modern_21hr',
            'lr_cam4_4x5_modern_22hr',
            'lr_cam4_4x5_modern_23hr',
            'lr_cam4_4x5_modern_24hr',
            'lr_cam4_4x5_modern_25hr',
            'nspladj_lr_wcsc_4x5_modern_21-5hr_ZM2', #8
            'nspladj_lr_wcsc_4x5_modern_22hr_ZM2',
            'lambres_4x5_modern_22-5hr_ZM2',
            'nspladj_lr_wcsc_4x5_modern_24hr_ZM2',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-369_21-5hr_ZM2', #12
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-369_22hr_ZM2',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-369_22-5hr_ZM2',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-369_24hr_ZM2',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-1476_21-5hr_ZM2', #16
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-1476_22hr_ZM2',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-1476_22-5hr_ZM2',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-1476_24hr_ZM2',
            'wcsc_1850_pop2_4x5', #20
            'nomg1_lr_exocam_4x5_lowCH4_16hr', #21
            'nomg1_lr_exocam_4x5_lowCH4_18hr',
            'nomg1_lr_exocam_4x5_lowCH4_20hr',
            'nomg1_lr_exocam_4x5_lowCH4_21hr',
            'nomg1_lr_exocam_4x5_lowCH4_22hr',
            'nomg1_lr_exocam_4x5_lowCH4_22-2hr',
            'nomg1_lr_exocam_4x5_lowCH4_22-3hr',
            'nomg1_lr_exocam_4x5_lowCH4_23hr',
            'nomg1_lr_exocam_4x5_lowCH4_24hr',
            'nomg1_lr_exocam_4x5_lowCH4_25hr',
            'nomg1_lr_exocam_4x5_hiCH4_16hr', #31
            'nomg1_lr_exocam_4x5_hiCH4_18hr',
            'nomg1_lr_exocam_4x5_hiCH4_20hr',
            'nomg1_lr_exocam_4x5_hiCH4_21hr',
            'nomg1_lr_exocam_4x5_hiCH4_21-5hr',
            'nomg1_lr_exocam_4x5_hiCH4_21-6hr',
            'nomg1_lr_exocam_4x5_hiCH4_21-7hr',
            'nomg1_lr_exocam_4x5_hiCH4_22hr',
            'nomg1_lr_exocam_4x5_hiCH4_23hr',
            'nomg1_lr_exocam_4x5_hiCH4_24hr',
            'nomg1_lr_exocam_4x5_hiCH4_25hr',
            'solar0p9_lr_exocam_4x5_ch4-0_co2-5250_22hr', #42
            'solar0p9_lr_exocam_4x5_ch4-0_co2-5250_24hr',
            'solar0p9_lr_exocam_4x5_ch4-10_co2-5250_21hr', #44
            'solar0p9_lr_exocam_4x5_ch4-10_co2-5250_21-5hr',
            'solar0p9_lr_exocam_4x5_ch4-10_co2-5250_22hr',
            'solar0p9_lr_exocam_4x5_ch4-10_co2-5250_24hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_16hr', #48
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_18hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_20hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21-5hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-5hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_23hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_25hr',
            'solar0p9_lr_exocam_4x5_ch4-100_co2-5250_21hr', #59
            'solar0p9_lr_exocam_4x5_ch4-100_co2-5250_21-5hr',
            'solar0p9_lr_exocam_4x5_ch4-100_co2-5250_22hr',
            'solar0p9_lr_exocam_4x5_ch4-100_co2-5250_24hr',
            'solar0p9_lr_exocam_4x5_ch4-1000_co2-5250_21-5hr', #63
            'solar0p9_lr_exocam_4x5_ch4-1000_co2-5250_24hr',
            'solar0p9_lr_exocam_4x5_ch4-2000_co2-5250_21-5hr', #65
            'solar0p9_lr_exocam_4x5_ch4-2000_co2-5250_24hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_20hr', #67
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_21hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_21-5hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_21-75hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_22hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_22hr',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_24hr',
            'solar0p9_lr_4x5_ch4-30_co2-2000_21-5hr', #74
            'solar0p9_lr_4x5_ch4-30_co2-2000_22hr',
            'solar0p9_lr_4x5_ch4-30_co2-2000_22-5hr',
            'solar0p9_lr_4x5_ch4-30_co2-2000_23hr',
            'solar0p9_lr_4x5_ch4-30_co2-2000_24hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_16hr', #79
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_18hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_20hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_20-5hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_21hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_21-5hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_22hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_24hr',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr', #87
            ] #88

#list of simulations to put on this plot (diurnally averaged around solar longitude)
simnames = [
            'lr_cam4_4x5_modern_16hr_branch', #0
            'lr_cam4_4x5_modern_18hr_branch',
            'lr_cam4_4x5_modern_20hr_branch',
            'lr_cam4_4x5_modern_21hr_branch',
            'lr_cam4_4x5_modern_22hr_branch',
            'lr_cam4_4x5_modern_23hr_branch',
            'lr_cam4_4x5_modern_24hr_branch',
            'lr_cam4_4x5_modern_25hr_branch',
            'nspladj_lr_wcsc_4x5_modern_21-5hr_ZM2_branch', #8
            'nspladj_lr_wcsc_4x5_modern_22hr_ZM2_branch',
            'lambres_4x5_modern_22-5hr_ZM2_branch',
            'nspladj_lr_wcsc_4x5_modern_24hr_ZM2_branch',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-369_21-5hr_ZM2_branch', #12
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-369_22hr_ZM2_branch',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-369_22-5hr_ZM2_branch',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-369_24hr_ZM2_branch',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-1476_21-5hr_ZM2_branch', #16
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-1476_22hr_ZM2_branch',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-1476_22-5hr_ZM2_branch',
            'nspladj_lr_wcsc_4x5_o2-1perc_co2-1476_24hr_ZM2_branch',
            'wcsc_1850_pop2_4x5_branch', #20
            'nomg1_lr_exocam_4x5_lowCH4_16hr_branch', #21
            'nomg1_lr_exocam_4x5_lowCH4_18hr_branch',
            'nomg1_lr_exocam_4x5_lowCH4_20hr_branch',
            'nomg1_lr_exocam_4x5_lowCH4_21hr_branch',
            'nomg1_lr_exocam_4x5_lowCH4_22hr_branch',
            'nomg1_lr_exocam_4x5_lowCH4_22-2hr_branch2',
            'nomg1_lr_exocam_4x5_lowCH4_22-3hr_branch2',
            'nomg1_lr_exocam_4x5_lowCH4_23hr_branch',
            'nomg1_lr_exocam_4x5_lowCH4_24hr_branch',
            'nomg1_lr_exocam_4x5_lowCH4_25hr_branch',
            'nomg1_lr_exocam_4x5_hiCH4_16hr_branch2', #31
            'nomg1_lr_exocam_4x5_hiCH4_18hr_branch2',
            'nomg1_lr_exocam_4x5_hiCH4_20hr_branch2',
            'nomg1_lr_exocam_4x5_hiCH4_21hr_branch2',
            'nomg1_lr_exocam_4x5_hiCH4_21-5hr_branch2',
            'nomg1_lr_exocam_4x5_hiCH4_21-6hr_branch2',
            'nomg1_lr_exocam_4x5_hiCH4_21-7hr_branch2',
            'nomg1_lr_exocam_4x5_hiCH4_22hr_branch2',
            'nomg1_lr_exocam_4x5_hiCH4_23hr_branch2',
            'nomg1_lr_exocam_4x5_hiCH4_24hr_branch2',
            'nomg1_lr_exocam_4x5_hiCH4_25hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-0_co2-5250_22hr_branch', #42
            'solar0p9_lr_exocam_4x5_ch4-0_co2-5250_24hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-10_co2-5250_21hr_branch', #44
            'solar0p9_lr_exocam_4x5_ch4-10_co2-5250_21-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-10_co2-5250_22hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-10_co2-5250_24hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_16hr_branch', #48
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_18hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_20hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_21-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_23hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_25hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-100_co2-5250_21hr_branch', #59
            'solar0p9_lr_exocam_4x5_ch4-100_co2-5250_21-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-100_co2-5250_22hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-100_co2-5250_24hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-1000_co2-5250_21-5hr_branch', #63
            'solar0p9_lr_exocam_4x5_ch4-1000_co2-5250_24hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-2000_co2-5250_21-5hr_branch', #65
            'solar0p9_lr_exocam_4x5_ch4-2000_co2-5250_24hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_20hr_branch', #67
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_21hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_21-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_21-75hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_22hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_23hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-3000_co2-5250_24hr_branch',
            'solar0p9_lr_4x5_ch4-30_co2-2000_21-5hr_branch', #74
            'solar0p9_lr_4x5_ch4-30_co2-2000_22hr_branch',
            'solar0p9_lr_4x5_ch4-30_co2-2000_22-5hr_branch',
            'solar0p9_lr_4x5_ch4-30_co2-2000_23hr_branch',
            'solar0p9_lr_4x5_ch4-30_co2-2000_24hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_16hr_branch', #79
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_18hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_20hr_branch2',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_20-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_21hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_21-5hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_22hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-50000_24hr_branch',
            'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr_branch_div2' #87
            ]

#list of rotation periods
rotpers = np.array([
                    0.6667, 0.75, 0.8333, 0.875, 0.9167, 0.9583, 1.0, 1.04167,
                    0.8958, 0.9167, 0.9375, 1.0,
                    0.8958, 0.9167, 0.9375, 1.0,
                    0.8958, 0.9167, 0.9375, 1.0,
                    1.0,
                    0.6667, 0.75, 0.8333, 0.875, 0.9167, 0.925, 0.9292, 0.9583, 1.0, 1.04167,
                    0.6667, 0.75, 0.8333, 0.875, 0.8958, 0.9, 0.9042, 0.9167, 0.9583, 1.0, 1.04167,
                    0.9167, 1.0,
                    0.875, 0.8958, 0.9167, 1.0,
                    0.6667, 0.75, 0.8333, 0.875, 0.8958, 0.9167, 0.9271, 0.9375, 0.9583, 1.0, 1.04167,
                    0.875, 0.8958, 0.9167, 1.0,
                    0.8958, 1.0,
                    0.8958, 1.0,
                    0.8333, 0.875, 0.8958, 0.9063, 0.9167, 0.9583, 1.0,
                    0.8958, 0.9167, 0.9375, 0.9583, 1.0,
                    0.6667, 0.75, 0.8333, 0.8542, 0.875, 0.8958, 0.9167, 1.0,
                    0.9271
  ])

cmap = plt.cm.seismic

#these lists determine how the simulations get sorted and placed on the plot
sets = [0, 8, 12, 16, 20, 21, 31, 42, 44, 48, 59, 63, 65, 67, 74, 79, 87, 88]
rows = [0, 0, 0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 1]

#these lists determine how each set of simulations is displayed
colors = [cmap(0.0), cmap(0.1), cmap(0.3), cmap(0.4), cmap(0.55),cmap(0.7),
          cmap(0.9), cmap(0.2), cmap(0.3), 'k', cmap(0.4), cmap(0.55),
          cmap(0.65),'0.3',cmap(0.75),cmap(0.9), 'gold']
labels = ['CAM4', 'WACCM', 'WACCM, low O$_3$', 'WACCM, low O$_3$, 4x CO$_2$',
            'WACCM, POP', 'ExoCAM, 0.8 ppm CH$_4$','ExoCAM, 3000 ppm CH$_4$',
            'ExoCAM, 0 ppm CH$_4$', 'ExoCAM, 10 ppm CH$_4$',
            'ExoCAM, 30 pmm CH$_4$', 'ExoCAM, 100 ppm CH$_4$',
            'ExoCAM, 1000 ppm CH$_4$','ExoCAM, 2000 ppm CH$_4$',
            'ExoCAM, 3000 ppm CH$_4$', 'ExoCAM, 30 ppm CH$_4$,\n 2000 ppm CO$_2$',
            'ExoCAM, 30 ppm CH$_4$,\n 50000 ppm CO$_2$',
            'ExoCAM, 30 pmm CH$_4$,\n2nd-order div.-damping']
markers = ['*', '*', '*', '*', '*', 's', 's', 's', 's', 's', 's', 's', 's',
            '^', 's', 's','*']
ls = [':', ':', ':', ':', 'None', ':', ':', 'None', ':', '-', ':', 'None',
            'None', '--', ':', ':','None']

dp22 = np.zeros_like(rotpers)
p_anom_amp = np.zeros_like(rotpers)
Ts = np.zeros_like(rotpers)
reference_index = np.where(rotpers==1.0)[0]

for i in np.arange(len(simnames)):
  if pathlib.Path(parent_path+simnames[i]).exists():
    p_anom_file = parent_path + simnames[i]+'/merged_hist/'+simnames[i]+'_p_anom_save.npz'
  else:
    p_anom_file = suppl_path + simnames[i]+'/merged_hist/'+simnames[i]+'_p_anom_save.npz'

  arc = np.load(p_anom_file)
  lon = arc['lon']
  lat = arc['lat']
  rotrate = arc['rotrate']
  if 'field_anom_mean' in arc:
      ps_anom_mean = arc['field_anom_mean']
  else:
      ps_anom_mean = arc['ps_anom_mean']
  lmax = arc['lmax']
  if 'clm_anom' in arc:
      clm = arc['clm_anom']
  else:
      clm = arc['clm']


  clm[:,:2,:] = 0.0
  clm[:,3:,:] = 0.0
  clm[:,:,:2] = 0.0
  clm[:,:,3:] = 0.0
  shmap = np.real(sh.expand.MakeGridDH(clm,sampling=2))
  #stuff for setting up new grid
  n = 2*lmax+2
  shlat = np.arange(-90+180/n,90+180/n,180/n)
  shlon = np.arange(0,360,180/n)
  shlon2d, shlat2d = np.meshgrid(shlon,shlat)

  dp22[i] = np.sqrt(clm[0,2,2]**2 + clm[1,2,2]**2)
  p_anom_amp[i] = np.max(shmap)

  #get glob mean surf T
  if pathlib.Path(parent_path+simnamesT[i]).exists():
    file = parent_path + simnamesT[i] + "/merged_hist/" + simnamesT[i] + ".cam.h0.globmean_0031_0060.nc"
  else:
    file = suppl_path + simnamesT[i] + "/merged_hist/" + simnamesT[i] + ".cam.h0.globmean_0031_0060.nc"

  if not pathlib.Path(file).exists():
    print(file+" is missing")
  else:
    data = nc.Dataset(file,'r')
    Ts[i] = data['TS'][:].squeeze()

cm = 1./2.54
fig = plt.figure(figsize=(18*cm,14.4*cm))

outer_grid = gridspec.GridSpec(2,1,wspace=0.2,hspace=0.3,left=0.08,right=0.8,
                                bottom=0.1,top=0.98,height_ratios=(2,2))
upper_grid = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=outer_grid[0],
                                    wspace=0.17,hspace=0.1,width_ratios=(1,1))
lower_grid = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=outer_grid[1],
                                    wspace=0.17,hspace=0.1,width_ratios=(1,1))

ax1 = fig.add_subplot(upper_grid[0])
ax2 = fig.add_subplot(upper_grid[1])
ax3 = fig.add_subplot(lower_grid[0])
ax4 = fig.add_subplot(lower_grid[1])

lines2 = []
lines4 = []
labels2 = []
labels4 = []
for iset in np.arange(len(sets)-1):

  if sets[iset] == 48 or sets[iset] == 67:
    lw = 2
    zorder = 100
  else:
    lw = 1
    zorder = -1
  if rows[iset] == 0:
    ax = ax1
  else:
    ax = ax3

  if markers[iset] == '*':
    ms = 4
  else:
    ms = 3

  line, =  ax.plot(rotpers[sets[iset]:sets[iset+1]]*24, p_anom_amp[sets[iset]:sets[iset+1]], marker=markers[iset], color=colors[iset], linestyle=ls[iset], label=labels[iset],lw=lw,zorder=zorder,ms=ms)
  if rows[iset] == 0:
    lines2.append(line)
    labels2.append(labels[iset])
  else:
    lines4.append(line)
    labels4.append(labels[iset])

  if iset != len(sets)-2:
    if rows[iset] == 0:
      ax = ax2
    else:
      ax = ax4
    ax.plot(p_anom_amp[sets[iset]:sets[iset+1]], (Ts[sets[iset]:sets[iset+1]]-Ts[reference_index[iset]]),marker=markers[iset], color=colors[iset], label=labels[iset],ls='None',ms = ms)

ax1.set_xlabel('Length of day (hours)')
ax1.set_ylabel('Semidiurnal pressure amplitude (Pa)')
ax1.set(ylim = (0,400))
ax2.set_ylabel('$T_{\mathrm{surf}} - T_{\mathrm{surf,24hr}}$ (K)')
ax2.set_xlabel('Semidiurnal pressure amplitude (Pa)')
ax2.hlines(0.0,0,400,colors='k',linestyles=':')
ax2.set_xlim(0,400)
ax2.set(ylim = (-0.7,4.1))
ax2.text(410, 3.5,'S = S$_0$\n simulations',fontsize=7)
ax2.legend(lines2,labels2,loc='best',bbox_to_anchor=(1.,0.1,0.2,0.5),fontsize=6,handlelength=3,ncols=1)


ax3.set_xlabel('Length of day (hours)')
ax3.set_ylabel('Semidiurnal pressure amplitude (Pa)')
ax3.set(ylim = (0,400))

ax4.set_ylabel('$T_{\mathrm{surf}} - T_{\mathrm{surf,24hr}}$ (K)')
ax4.set_xlabel('Semidiurnal pressure amplitude (Pa)')
ax4.hlines(0.0,0,400,colors='k',linestyles=':')
ax4.set_xlim(0,400)
ax4.set(ylim = (-0.7,4.1))
ax4.legend(lines4,labels4,loc='best',bbox_to_anchor=(1.,0.3,0.2,0.5),fontsize=6,handlelength=3,ncols=1)
ax4.text(410, 3.5,'S = 0.9 S$_0$\n simulations',fontsize=7)

plt.savefig('figures/figure_3.pdf')
plt.close()
