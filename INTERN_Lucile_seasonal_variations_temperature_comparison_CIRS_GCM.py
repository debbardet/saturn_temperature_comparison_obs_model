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
pressure = pressure*1.e5 #(conversion in Pascal)
##
ptarget=float(input("Pressure target? (in Pa)\n"))  # choix de la pression
if ptarget >=3.e-1 and ptarget <= 3.e3:
    print("CIRS Limb (2005-2006,2010,2015)")
if ptarget >= 10. and ptarget <= 1.e3:
    print("TEXES (2017)")
if ptarget >= 8.e3 and ptarget <= 3.e4 or (ptarget >= 50. and ptarget <= 5.e2):
    print ("CIRS Nadir (2017, 2015, 2014, 2005)")
#############################################
Lstarget=int(input("Ls target? (2005-2006==Ls310, 2010==Ls, 2014-2015==Ls65, 2017==Ls85) \n"))
if Lstarget == 310:
    Earthyear = 2005
    ## CIRS data
    CIRSlimb = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/Temperature_Limb2005-2006_indiv.txt")
    lat_limb = convert_to_planetocentric(CIRSlimb[0,:])
    ## Fletcher data 
    fletCIRS = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/data_fletcher_2005_format_guerlet.txt")
    flet_Latitude_pc = convert_to_planetocentric(fletCIRS[0,:])
    pres_flet = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/pressure_fletcher_2005.txt")
    pres_flet = pres_flet*1.e5 #(conversion in Pascal)
#if Lstarget == ?:
    #Earthyear = 2010
    ## CIRS data
    #CIRSlimb = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/Temperature_Limb2010.txt")
    #lat_limb = convert_to_planetocentric(CIRSlimb[0,:])
if Lstarget == 65:
    Earthyear = 2015
    ## CIRS data
    CIRSlimb = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/Temperature_Limb2015.txt")
    lat_limb = convert_to_planetocentric(CIRSlimb[0,:])
    CIRSnadir1 = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/Temperature_NadirMay2015.txt")
    lat_nadir1 = convert_to_planetocentric(CIRSnadir1[0,:])
    CIRSnadir2 = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/Temperature_NadirOct2014.txt")
    lat_nadir2 = convert_to_planetocentric(CIRSnadir2[0,:])
if Lstarget == 85:
    Earthyear = 2017
    ## CIRS data
    CIRSnadir1 = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/Temperature_NadirFeb2017.txt")
    lat_nadir1 = convert_to_planetocentric(CIRSnadir1[0,:])
    CIRSnadir2 = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/Temperature_NadirJul2017_HemSud.txt")
    lat_nadir2 = convert_to_planetocentric(CIRSnadir2[0,:])
    ## TEXES data
    daTEXES = np.loadtxt("/planeto/dbardet/results/simulation_referente/stage_Lucile/Retrieved_Temperature_Profiles/Temperature_NadirTEXES_Mar2017.txt")
    lat_TEXES = convert_to_planetocentric(daTEXES[0,:])
#############################################

diff = np.abs(pressure - ptarget); w=np.where(diff==np.min(diff))
if Lstarget == 310 or Lstarget == 65:
    CIRSlimb = np.squeeze(CIRSlimb[w[0],:])
if Lstarget == 65 or Lstarget == 85:
    CIRSnadir1 = np.squeeze(CIRSnadir1[w[0],:])
    CIRSnadir2 = np.squeeze(CIRSnadir2[w[0],:])
if Lstarget == 85:
    daTEXES = np.squeeze(daTEXES[w[0],:])
if Lstarget == 310:
    diff = np.abs(pres_flet - ptarget); w_flet=np.where(diff==np.min(diff))
    fletCIRS = np.squeeze(fletCIRS[w_flet[0],:])

#############################################
zefig = ppplot.figuref(x=12,y=8)
mer = ppplot.plot1d(fig=zefig,linestyle="",xlabel="Latitude (deg N)",ylabel="Ls={}$^\circ$: Temperature (K) at {} mbar".format(Lstarget,ptarget*1.e-2),ymin=100,ymax=170,nyticks=15)

############ Observations ###################
if Lstarget == 310 or Lstarget == 65:
    mer.f, mer.x, mer.legend, mer.color, mer.marker = CIRSlimb, lat_limb,"CIRS {} Limb".format(Earthyear),"b","s"; mer.make()
if Lstarget == 65:
    mer.f, mer.x, mer.legend, mer.color, mer.marker = CIRSnadir1, lat_nadir1,"CIRS May {} Nadir".format(Earthyear),"b","x"; mer.make() 
    mer.f, mer.x, mer.legend, mer.color, mer.marker = CIRSnadir2, lat_nadir2,"CIRS Oct {} Nadir".format(Earthyear-1),"r","x"; mer.make()
if Lstarget == 310:
    mer.f, mer.x, mer.legend, mer.color, mer.marker = fletCIRS, flet_Latitude_pc,"CIRS {} Nadir".format(Earthyear),"b","x"; mer.make() 
if Lstarget == 85:
    mer.f, mer.x, mer.legend, mer.color, mer.marker = daTEXES, lat_TEXES,"TEXES {} Nadir".format(Earthyear),"r","o"; mer.make() 
    mer.f, mer.x, mer.legend, mer.color, mer.marker = CIRSnadir1, lat_nadir1,"CIRS Feb {} Nadir".format(Earthyear),"b","x"; mer.make() 
    mer.f, mer.x, mer.legend, mer.color, mer.marker = CIRSnadir2, lat_nadir2,"CIRS Jul {} Nadir".format(Earthyear),"r","x"; mer.make() 

############ Rad-Conv #######################
ff,xx,yy,zz,tt = pp(file="/planeto/dbardet/results/simulation_referente/ref_61lvls/ccat_precast_radiatif_300k.nc", var="temp",x=999,z=ptarget,useindex="1000",t=roughindex(Lstarget,9)).getfd()
mer.f, mer.x, mer.legend, mer.color, mer.marker, mer.linestyle = ff, yy, "Radiative-convective model", "g","","--" ; mer.make()

############ GCM outputs ####################
## NO-RING GCM OUTPUTS
ff,xx,yy,zz,tt = pp(file="/planeto/dbardet/results/simulation_referente/ref_61lvls/ccat_precast_318k_noring.nc", var="temperature",x=999,z=ptarget,useindex="1000",t=roughindex(Lstarget,4)).getfd()
mer.f, mer.x, mer.legend, mer.color, mer.marker, mer.linestyle = ff, yy, "GCM, No-Ring", "k","", None ; mer.make()
## RINGS GCM OUTPUTS
for year in np.arange(7,14):
    ff,xx,yy,zz,tt = pp(file=dafile,var="temperature",x=999,z=ptarget,useindex="1000",t=roughindex(Lstarget,year)).getfd()
    mer.f, mer.x, mer.legend = ff, yy, "GCM, Year %i" % (year)
    mer.marker = ""
    mer.linestyle = None
    mer.color = coltab[year]
    mer.make()

if ptarget <= 100.:
   mpl.legend(loc='lower center',ncol=3)
else:
   mpl.legend(loc='upper center',ncol=3)
ppplot.save(filename="temperature_comparison_CIRS-GCM-96lvl_{}mbar_Ls{}".format(ptarget*1.e-2,Lstarget),mode="png")
mpl.show()

