#%matplotlib inline
from ppclass import pp
import numpy as np
import matplotlib.pyplot as mpl
from matplotlib.cm import get_cmap
import ppplot
from scipy import *

ppplot.changefont(16)
#mpl.rcParams['font.family'] = "times" # for good log major/minor ticks
mpl.rcParams['lines.linewidth'] = 2.0
mpl.rcParams['lines.markersize'] = 4

#############################################
dafile = "/planeto/dbardet/results/simulation_referente/ref_61lvls/ccat_precast_318k.nc"
#############################################
colorm = "magma" #"brg" #"magma"
pal = get_cmap(name=colorm)
coltab = [pal(i) for i in np.linspace(0,0.9,14)]
#############################################
def roughindex(ls,year):
    ## this is the year_th simulated year
    return int(((year-1)*49.) + (ls*49./360.))
#############################################
## Conversion planetographic to planetocentric 
def convert_to_planetocentric(lat_planetographic):
    e = (60268./54364.)**2.
    return np.arctan(np.tan(lat_planetographic*3.14159/180.)/e)*180./3.14159
#############################################
pressure = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/pressure.txt")
#############################################
Lstarget=int(input("Ls target? (2005-2006==Ls310, 2014-2015==Ls65) \n"))
if Lstarget == 310:
    Earthyear = 2005
    ## CIRS data
    dataCIRS = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/Temperature_Limb2005-2006_avg.txt", dtype="float")
    CIRSlimb = dataCIRS[1:,:]
    lat_limb = convert_to_planetocentric(dataCIRS[0,:])
if Lstarget == 65:
    Earthyear = 2015
    ## CIRS data
    dataCIRS = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/Temperature_Limb2015.txt", dtype="float")
    CIRSlimb = dataCIRS[1:,:]
    lat_limb = convert_to_planetocentric(dataCIRS[0,:])
#############################################
## Average over average over all latitudes except the equatorial region within +/-20$^\circ$
global_mean=np.zeros(len(pressure))
for i in range(len(pressure)):
    length=0
    for j in range(len(lat_limb)):
        if lat_limb[j] <= -20 or lat_limb[j] >= 20:
            global_mean[i]+=CIRSlimb[i][j]
            length+=1
    global_mean[i]=global_mean[i]/length

#############################################
zefig = ppplot.figuref(x=8,y=14)
mer = ppplot.plot1d(fig=zefig,linestyle="", logy=True, invert=True, xlabel="Temperature (K) at Ls={}$^\circ$".format(Lstarget),ylabel="Pressure (mbar)",xmin=80,xmax=170,nxticks=10,ymin=1.e-3,ymax=1.e2,nyticks=8)
#mer = ppplot.plot1d(fig=zefig,linestyle="", logy=True, xlabel="Temperature (K) at Ls={}$^\circ$".format(Lstarget),ylabel="Pressure (mbar)",xmin=80,xmax=170,nxticks=10,ymin=2.e2,ymax=1.e-3,nyticks=8)

############ Observations ###################
mer.f, mer.x, mer.legend, mer.color, mer.marker, mer.linestyle = pressure*1.e3, global_mean, "CIRS {} Limb".format(Earthyear),"b","","dotted"; mer.make()

############ Rad-Conv #######################
ff,xx,yy,zz,tt = pp(file="/planeto/dbardet/results/simulation_referente/ref_61lvls/ccat_precast_radiatif_300k.nc", var="temp",x=999,useindex="1000",t=roughindex(Lstarget,9)).getfd()
ff_mean=np.zeros(len(zz))
for i in range(len(zz)):
    length=0
    for j in range(len(yy)):
        if yy[j] <= -20 or yy[j] >= 20:
            ff_mean[i]+=ff[i][j]
            length+=1
    ff_mean[i]=ff_mean[i]/length
mer.f, mer.x, mer.legend, mer.color, mer.marker, mer.linestyle,mer.invert = zz*1.e-2, ff_mean, "Rad-conv model", "g","","--", False ; mer.make()

############ GCM outputs ####################

## NO-RING GCM OUTPUTS
ff,xx,yy,zz,tt = pp(file="/planeto/dbardet/results/simulation_referente/ref_61lvls/ccat_precast_318k_noring.nc", var="temperature",x=999,useindex="1000",t=roughindex(Lstarget,4)).getfd()
ff_mean=np.zeros(len(zz))
for i in range(len(zz)):
    length=0
    for j in range(len(yy)):
        if yy[j] <= -20 or yy[j] >= 20:
            ff_mean[i]+=ff[i][j]
            length+=1
    ff_mean[i]=ff_mean[i]/length
mer.f, mer.x, mer.legend, mer.color, mer.marker, mer.linestyle = zz*1.e-2, ff_mean, "GCM, No-Ring", "k","", None ; mer.make()

## RINGS GCM OUTPUTS
for year in np.arange(7,14):
    ff,xx,yy,zz,tt = pp(file=dafile,var="temperature",x=999,useindex="1000",t=roughindex(Lstarget,year)).getfd()
    ff_mean=np.zeros(len(zz))
    for i in range(len(zz)):
        length=0
        for j in range(len(yy)):
            if yy[j] <= -20 or yy[j] >= 20:
                ff_mean[i]+=ff[i][j]
                length+=1
        ff_mean[i]=ff_mean[i]/length
    mer.f, mer.x, mer.legend = zz*1.e-2, ff_mean, "GCM, Year %i" % (year)
    mer.marker = ""
    mer.linestyle = None
    mer.color = coltab[year]
    mer.make()
#
#if ptarget <= 100.:
#   mpl.legend(loc='lower center',ncol=3)
#else:
#   mpl.legend(loc='upper center',ncol=3)
ppplot.save(filename="global_vert-struct_temp_comparison_CIRS-GCM_Ls{}".format(Lstarget),mode="png")
mpl.show()


