import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import scipy.interpolate as sint
import pyshtools as sh

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':7})

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

#change these paths to match the location of the sim_main and sim_suppl
# download from the data repository
sims_main = '../DG_Atmospheric_Resonance_archive/sims_main/'
sims_suppl = '../DG_Atmospheric_Resonance_archive/sims_suppl/'

waccm_sim1 = 'lambres_4x5_modern_24hr_ZM2'
waccm_sim2 = 'lambres_4x5_modern_22-5hr_ZM2'

exo_sim1 = 'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_24hr'
exo_sim2 = 'solar0p9_lr_exocam_4x5_ch4-30_co2-5250_22-25hr'

globfiles1 = [sims_suppl + waccm_sim1 + '/merged_hist/' + waccm_sim1 + '.cam.h0.globmean_0031_0060.nc',
              sims_main + exo_sim1 + '/merged_hist/' + exo_sim1 + '.cam.h0.globmean_0031_0060.nc']
swfiles1 = [sims_suppl + waccm_sim1 + '_branch2/merged_hist/' + waccm_sim1 +'_branch2_QRS_save.npz',
            sims_main + exo_sim1 + '_branch2/merged_hist/' + exo_sim1 +'_branch2_QRS_save.npz']
mffiles1 = [sims_suppl + waccm_sim1 + '_branch2/merged_hist/' + waccm_sim1 + '_branch2_CMFMC_save.npz',
             sims_main + exo_sim1 + '_branch2/merged_hist/' + exo_sim1 +'_branch2_CMFMC_save.npz']


globfiles2 = [sims_suppl + waccm_sim2 + '/merged_hist/' + waccm_sim2 + '.cam.h0.globmean_0031_0060.nc',
              sims_main + exo_sim2 + '/merged_hist/' + exo_sim2 + '.cam.h0.globmean_0031_0060.nc']
swfiles2 = [sims_suppl + waccm_sim2 + '_branch2/merged_hist/' + waccm_sim2 +'_branch2_QRS_save.npz',
            sims_main + exo_sim2 + '_branch2/merged_hist/' + exo_sim2 + '_branch2_QRS_save.npz']
mffiles2 = [sims_suppl + waccm_sim2 + '_branch2/merged_hist/' + waccm_sim2 +'_branch2_CMFMC_save.npz',
            sims_main + exo_sim2 + '_branch2/merged_hist/' + exo_sim2 + '_branch2_CMFMC_save.npz']

globfiles = [globfiles1, globfiles2]
swfiles = [swfiles1,swfiles2]
mffiles = [mffiles1,mffiles2]

swfactor = [1.0, 1.0/86400]
cps = [1005., 1035.]
model = ['WACCM', 'ExoCAM']

def solve_vse(h,x2,J2_full,Hf):
    #x2, J2_full, Hf should have same shape (J2_full can be complex)
    a = np.zeros_like(J2_full)
    b = np.zeros_like(J2_full)
    dHfdx = np.gradient(Hf,x2)
    A = 1
    C = 1

    y2 = np.zeros((len(h),len(x2)),dtype=complex)

    J2 = J2_full.copy()
    dx = x2[1] - x2[0]  #assumes uniform x grid

    #loop over h
    for ih in np.arange(len(h)):
      a[0] = 1.0/ (1-(Hf[0]/h[ih] - 0.5)*dx)
      b[0] = 0    # no gravitational tide
      kx2 = 0.25*(4*kappa*Hf/h[ih] - 1)
      lam = np.zeros_like(kx2,dtype=complex)

      lam[kx2>=0] = np.sqrt(kx2[kx2>=0])
      lam[kx2<0] = 1j*np.sqrt(-kx2[kx2<0])

      for ix in np.arange(1,len(x2)):
          B = - ( 2 + dx**2/4*(1-4/h[ih]*(kappa*Hf[ix] + dHfdx[ix])) )
          D = dx**2 * kappa * J2[ix] / (gamma*g*h[ih]) * np.exp(-x2[ix]/2)
          a[ix] = - A / (B + a[ix-1]*C)
          b[ix] = (D - b[ix-1]*C) / (B + a[ix-1]*C)

      #top level y
      y2[ih,-1] = (b[-3] + b[-2]*(2j*lam[-2]*dx + a[-3])) / (1 - a[-2]*(2j*lam[-2]*dx + a[-3]))
      for ix in np.arange(len(x2)-2,-1,-1):
          y2[ih,ix] = a[ix]*y2[ih,ix+1] + b[ix]

    dps2 = 1j*gamma/sigma * y2[:,0] * ps
    return dps2

#constants
g = 9.8
kappa = 2./7
gamma = 1.4
R = 287.
Rp = 6371e3
Lam0 = 11.1

def h2prot(h):
    rotrate = 1/(2*Rp)*np.sqrt(h*Lam0*g)
    prot = 2*np.pi/rotrate/3600
    return prot

def prot2h(prot):
    rotrate = 2*np.pi/prot/3600
    h = (2*Rp*rotrate)**2/Lam0/g/1000
    return h

#equivalent depths to calculate (meters)
h = np.linspace(7,12,500)*1000
#rotation parameters corresponding to h
rotrate = 1/(2*Rp)*np.sqrt(h*Lam0*g)
sigma = 2*rotrate
prot = 2*np.pi/rotrate/3600

cm = 1./2.54
fig, axes = plt.subplots(ncols=3,nrows=2,figsize=(18*cm,12*cm))

colors = ['k','orange','blue','0.5']
labels = [['Total','Stratosphere only','Troposphere only','Radiation only'],[None,None,None,None]]
mfc = [colors,colors] #['none','none','none','none']]
symbols = ['-',':']
setlabels = ['Forcing from 24hr case','Forcing from resonant case']

for iset in np.arange(len(globfiles)):
    for i in np.arange(len(globfiles[iset])):
        #read in climatology and flip it around so bottom is at index 0
        data = nc.Dataset(globfiles[iset][i],'r')
        T = (data['T'][:].squeeze()[::-1]).data
        p = data['lev'][:].squeeze()[::-1]*100  #pressure in Pa (roughly)
        ps = data['PS'][:].squeeze() #surface pressure
        z = data['Z3'][:].squeeze()[::-1]
        data.close()

        H = R*T/g

        #read in sw file, multiple coefficients by conversion factor and normalization factors
        arc_sw = np.load(swfiles[iset][i])
        clm_sw = arc_sw['clm_mean'][::-1]*swfactor[i]*np.sqrt(2*15/8)
        #forcing for semidiurnal component
        Jsw = clm_sw[:,0,2,2] + 1j*clm_sw[:,1,2,2]

        #read in mf file, multiply by normalization and convert to heating per kg
        arc_mf = np.load(mffiles[iset][i])
        clm_mf = arc_mf['clm_mean'][::-1]*np.sqrt(2*15/8)
        dTdz = np.gradient(T,z)
        rho = p / (R*T)
        Jmf = cps[i]/rho*(clm_mf[:,0,2,2]+1j*clm_mf[:,1,2,2])*(dTdz + g/cps[i])

        x = (-np.log(p/np.max(p)).squeeze()).data

        #interpolate to uniform grid
        x2 = np.linspace(x[0],x[-1],300)
        p2 = np.max(p)*np.exp(-x2)
        Hf = sint.interp1d(x,H,fill_value='extrapolate')(x2)
        J2_full = sint.interp1d(x,Jsw+Jmf,fill_value='extrapolate')(x2)
        J2_rad = sint.interp1d(x,Jsw,fill_value='extrapolate')(x2)

        dps2 = solve_vse(h,x2,J2_full,Hf)
        res_per = prot[np.argmax(dps2)]

        ax = axes[i][0]
        ax.plot(prot,np.absolute(dps2),symbols[iset],c=colors[0],ms=1,mfc=mfc[iset][0])
        ax.set(xlim=(np.min(prot),np.max(prot)),ylim=(0,3000))
        axtmp = ax.secondary_xaxis('top',functions=(prot2h,h2prot))
        axtmp.set_xlabel('Equivalent depth (km)')

        ax = axes[i][1]
        ax.plot(prot,np.absolute(dps2),symbols[iset],c=colors[0],ms=1,mfc=mfc[iset][0])
        ax.set(xlim=(res_per-0.5,res_per+0.5),ylim=(0,3000))
        axtmp = ax.secondary_xaxis('top',functions=(prot2h,h2prot))
        axtmp.set_xlabel('Equivalent depth (km)')

        ax = axes[i][2]
        if i == 0:
            ax.plot(prot,np.absolute(dps2),symbols[iset],c=colors[0],label=setlabels[iset],ms=1,mfc=mfc[iset][0])
        else:
            ax.plot(prot,np.absolute(dps2),symbols[iset],c=colors[0],label=labels[iset][0],ms=1,mfc=mfc[iset][0])
        ax.set(xlim=(23.5,24.5),ylim=(0,200))
        axtmp = ax.secondary_xaxis('top',functions=(prot2h,h2prot))
        axtmp.set_xlabel('Equivalent depth (km)')

        J2_stratos = J2_full.copy()
        J2_stratos[p2>100e2] = 0.0
        dps2_s = solve_vse(h,x2,J2_stratos,Hf)
        ax = axes[i][0]
        ax.plot(prot,np.absolute(dps2_s),symbols[iset],c=colors[1],ms=1,mfc=mfc[iset][1])
        ax.set(ylim=(0,3000))
        ax.set_xlabel('Length of day (hours)')

        ax = axes[i][1]
        ax.plot(prot,np.absolute(dps2_s),symbols[iset],c=colors[1],ms=1,mfc=mfc[iset][1])
        ax.set(xlim=(res_per-0.5,res_per+0.5),ylim=(0,3000))
        ax.set_xlabel('Length of day (hours)')

        ax = axes[i][2]
        if i == 0:
            ax.plot(prot,np.absolute(dps2_s),symbols[iset],c=colors[1],ms=1,mfc=mfc[iset][1])
        else:
            ax.plot(prot,np.absolute(dps2_s),symbols[iset],c=colors[1],label=labels[iset][1],ms=1,mfc=mfc[iset][1])
        ax.set(xlim=(23.5,24.5),ylim=(0,200))
        ax.set_xlabel('Length of day (hours)')

        J2_tropos = J2_full.copy()
        J2_tropos[p2<100e2] = 0.0
        dps2_t = solve_vse(h,x2,J2_tropos,Hf)
        ax = axes[i][0]
        ax.plot(prot,np.absolute(dps2_t),symbols[iset],c=colors[2],ms=1,mfc=mfc[iset][2])
        ax.set(ylim=(0,3000),ylabel='Pressure anomaly (Pa)\n(%s forcing/profile)'%model[i])
        ax.set_xlabel('Length of day (hours)')

        ax = axes[i][1]
        ax.plot(prot,np.absolute(dps2_t),symbols[iset],c=colors[2],ms=1,mfc=mfc[iset][2])
        ax.set(xlim=(res_per-0.5,res_per+0.5),ylim=(0,3000))
        ax.set_xlabel('Length of day (hours)')

        ax = axes[i][2]
        if i == 0:
            ax.plot(prot,np.absolute(dps2_t),symbols[iset],c=colors[2],ms=1,mfc=mfc[iset][2])
        else:
            ax.plot(prot,np.absolute(dps2_t),symbols[iset],c=colors[2],label=labels[iset][2],ms=1,mfc=mfc[iset][2])
        ax.set(xlim=(23.5,24.5),ylim=(0,200))
        ax.set_xlabel('Length of day (hours)')

        dps2_r = solve_vse(h,x2,J2_rad,Hf)
        ax = axes[i][0]
        ax.plot(prot,np.absolute(dps2_r),symbols[iset],c=colors[3],ms=1,mfc=mfc[iset][3])
        ax.set(ylim=(0,3000))
        ax.set_xlabel('Length of day (hours)')

        ax = axes[i][1]
        ax.plot(prot,np.absolute(dps2_r),symbols[iset],c=colors[3],ms=1,mfc=mfc[iset][3])
        ax.set(xlim=(res_per-0.5,res_per+0.5),ylim=(0,3000))
        ax.set_xlabel('Length of day (hours)')

        ax = axes[i][2]
        if i == 0:
            ax.plot(prot,np.absolute(dps2_r),symbols[iset],c=colors[3],ms=1,mfc=mfc[iset][3])
        else:
            ax.plot(prot,np.absolute(dps2_r),symbols[iset],c=colors[3],label=labels[iset][3],ms=1,mfc=mfc[iset][3])
        ax.set(xlim=(23.5,24.5),ylim=(0,200))
        ax.set_xlabel('Length of day (hours)')

ax.legend(loc='best',fontsize=6)
axes[0][2].legend(loc='best',fontsize=6)

plt.tight_layout()
plt.savefig('figures/ED_figure_2.eps')
plt.close()
