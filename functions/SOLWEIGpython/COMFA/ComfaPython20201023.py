# -*- coding: utf-8 -*-
"""
This code is a translation of radiation related function in the COMFA model

Translated on Fri Oct 23 
@author: xlinfr
"""

from radiationfunctionsCOMFA import CNRRabs_Total
from radiationfunctionsCOMFA import Rad_Total_solweig
import numpy as np


# Radiation components in w m-2
Lin = 200.
Lup = 300.
Kin = 800.
Kup = 200.

# Time variables used in solweig version
yyyy = 2020
month = 6
day = 23
hour = 12
minu = 0

# yyyy = 2020
# month = 12
# day = 20
# hour = 12
# minu = 0

# Time vaiables used in original
t=hour + minu/60. # decimal time of day (up to 24)
if (yyyy % 4) == 0:
    if (yyyy % 100) == 0:
        if (yyyy % 400) == 0:
            leapyear = 1
        else:
            leapyear = 0
    else:
        leapyear = 1
else:
    leapyear = 0
if leapyear == 1:
    dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
else:
    dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

doy = np.sum(dayspermonth[0:month - 1]) + day # day of year (Jan 1 = 1)

# Lubbock 33.5779,-101.8552

lat = 33.5779 # latitude
lon = -101.8552 # longitude
alt = 992. # Altitude (m asl)
utc = -6 # standard time lubbock

Atr=0.7 # atmospheric transmittance
alpha = 0.37 # albedo of the cylinder
emis = 0.95 # emissivity of the cylinder
L=0.1 # length of cylinder (cm)
D=0.01 # Diameter of cylinder (cm)

Aeff=0.78 # Effective area of body. 0.78 for standing from Campebll and Normal (1998) (0.70 for sitting)


Total_CNRRabs_solweig, zensolweig, Kdsolweig = Rad_Total_solweig(L,D,Lin,Lup,Kin,Kup,yyyy,month,day,hour,minu,lat,lon,alt, emis,alpha,Aeff,utc, Kd=None, Ta=None, RH=None)

Total_CNRRabs, zenCOMFA, kdCOMFA = CNRRabs_Total(alpha, L, D, Lin, Lup, Kin, Kup, Atr, doy, t, np.array([lat]), alt, emis, Aeff)  #np.array([lat]) to deal with timeseries?

print('Total_CNRRabs: ' + str(Total_CNRRabs))
print('altComfa: ' + str(90-zenCOMFA))
print('Total_CNRRabs_solweig: ' + str(Total_CNRRabs_solweig))
print('ZenSolweig: ' + str(90-zensolweig))
