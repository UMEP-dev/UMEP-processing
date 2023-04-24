# -*- coding: utf-8 -*-
"""
This code is a translation of radiation related function in the COMFA model

Translated on Fri Oct 23 
@author: xlinfr
"""

import numpy as np
from ....util.SEBESOLWEIGCommonFiles import Solweig_v2015_metdata_noload as metload
from ....util.SEBESOLWEIGCommonFiles.diffusefraction import diffusefraction


def CNRRabs_Total(alpha,L,D,Lin,Lup,Kin,Kup,Atr,d,t,lat,A,emis,Aeff):
    # VERTICAL CYLINDER MODEL 
    # Used to calculate the total radiation absorbed by a cylinder, such as the CRT, with inputs
    # alpha (albedo of cylinder), Kin (incoming shortwave measured by CNR net
    # radiometer, W/m2), Kup (reflected shortwave measured by CNR, W/m2), Lin
    # (CNR incoming longwave radiation, W/m2), Lup (CNR1 outgoing longwave
    # radition, W/m2), L (length of cylinder, m), D (diamter of cylinder, m).
    # 
    #Also need solar zenith, which is calculated from: lat (latitude, degrees), d (days, Jan 1=1), t (decimal
    # time, 24h decimal time).
        
    #Atr - atmospheric transmittance - can calculate for the day/time of day or
    #estimate based on sky conditions.
        
        # LOAD data ....
        
        #set constants
        # D=0.01
    # CNRRabs_Total.m:20
        
        # L=0.1
    # CNRRabs_Total.m:21
        
        # alpha=0.37
    # CNRRabs_Total.m:22
        
        #     E.g., 
    #     Kin = metdata(:,10); #put in column for Kin from metdata file
    #     Kup = metdata(:,11); 
    #     Lin  = metdata(:,20);
    #     Lup = metdata(:,21); 
    #     d = metdata(:,24);
    #     Atr = 0.70; 
    #     lat = 33.75;
    #     A = 992;
    
    # Aeff=0.78 # Effective area of body. 0.78 for standing from Campebll and Normal (1998) (0.70 for sitting) Moved to main function

    Kin_abs, zen, Kd = CNR_Kinabs_meas(alpha, Kin, L, D, A, lat, d, t, Atr)
    Kup_abs = CNR_Kup(alpha, Kup, L, D)
    Lin_abs = LinMeas_abs(L, D, Lin, emis)
    Lup_abs = LupMeas_abs(L, D, Lup, emis)
    Acyl = CRT_Acyl(L,D)
    Total_CNRRabs = np.dot(Aeff,((Kin_abs + Kup_abs + Lin_abs + Lup_abs) / (Acyl)))

    return Total_CNRRabs, zen, Kd

def COMFA_RAD_SPATIAL_TC(L, D, Lin, Lup, Kin, Kup, emis, alpha, Aeff, Kd, metdata, location, utc):
    """
    Same as CNRRabs_Total but using a SPN1 sensor to obtain Kd or
    Same as CNRRabs_Total but using reindl et al for Kd and sun_position.py
    """

    deg2rad = np.pi/180.

    # location = {'longitude': lon, 'latitude': lat, 'altitude': alt}
    YYYY, altitude, azimuth, zen, jday, leafon, dectime, altmax = metload.Solweig_2015a_metdata_noload(metdata, location, utc)

    Kb = Kin - Kd
    # Not perpendicular but horrisontal surface in COMFA Kb = (Kin - Kd)/(np.sin(altitude*deg2rad))

    Acs = CRT_Acs(L,D)

    if altitude > 0:
        # Kdirect perpendicular to the beam
        Kp = Kb/np.cos(zen)
        Kb_abs = (1 - alpha) * (Kp * np.sin(zen)) * Acs
    else:
        Kb_abs = 0
        Kd = 0
    # Kb_abs = np.dot(np.dot((1 - alpha),(np.dot(Kb,np.sin(zen)))),Acs)
    Acyl = CRT_Acyl(L,D)
    # Kd_abs = np.dot(np.dot(np.dot((1 - alpha),Kd),0.5),Acyl)
    Kd_abs = (1 - alpha) * Kd * 0.5 * Acyl

    Kin_abs = (Kb_abs + Kd_abs)

    # same as original
    Kup_abs = CNR_Kup(alpha, Kup, L, D)
    Lin_abs = LinMeas_abs(L, D, Lin, emis)
    Lup_abs = LupMeas_abs(L, D, Lup, emis)
    Acyl = CRT_Acyl(L,D)
    # Total_CNRRabs_solweig = np.dot(Aeff,((Kin_abs + Kup_abs + Lin_abs + Lup_abs) / (Acyl)))
    Total_CNRRabs_solweig = Aeff * ((Kin_abs + Kup_abs + Lin_abs + Lup_abs) / Acyl)

    return Total_CNRRabs_solweig

def Rad_Total_solweig(L,D,Lin,Lup,Kin,Kup,yyyy,hour,minu,doy,lat,lon,alt,emis,alpha,Aeff,utc, Kd=None, Ta=None, RH=None):
    """
    Same as CNRRabs_Total but using a SPN1 sensor to obtain Kd or
    Same as CNRRabs_Total but using reindl et al for Kd and sun_position.py
    """

    deg2rad = np.pi/180.
    metdata = np.zeros((1, 24)) - 999. #TODO: Move out of function
    # doy = day_of_year(yyyy, month, day)
    metdata[0, 0] = yyyy
    metdata[0, 1] = doy
    metdata[0, 2] = hour
    metdata[0, 3] = minu
    # metdata[0, 11] = Ta
    # metdata[0, 10] = RH
    # metdata[0, 14] = Kin
    # metdata[0, 21] = Kd

    location = {'longitude': lon, 'latitude': lat, 'altitude': alt}
    YYYY, altitude, azimuth, zen, jday, leafon, dectime, altmax = metload.Solweig_2015a_metdata_noload(metdata, location, utc)

    if Kd == None: # Kd and Kb is estimated with Reindl et al. 1990
        Itoa = 1370.0  # Effective solar constant
        D_ = sun_distance(doy)
        I0et = Itoa * np.cos(zen) * D_  # extra terrestial solar radiation
        Kt = Kin / I0et
        _, Kd = diffusefraction(Kin, altitude, Kt, Ta, RH)
        Kb = Kin - Kd
    else:
        Kb = Kin - Kd
        # Not perpendicular but horrisontal surface in COMFA Kb = (Kin - Kd)/(np.sin(altitude*deg2rad))

    Acs = CRT_Acs(L,D)

    if altitude > 0:
        # Kdirect perpendicular to the beam
        Kp = Kb/np.cos(zen)
        Kb_abs = (1 - alpha) * (Kp * np.sin(zen)) * Acs
    else:
        Kb_abs = 0
        Kd = 0
    # Kb_abs = np.dot(np.dot((1 - alpha),(np.dot(Kb,np.sin(zen)))),Acs)
    Acyl = CRT_Acyl(L,D)
    # Kd_abs = np.dot(np.dot(np.dot((1 - alpha),Kd),0.5),Acyl)
    Kd_abs = (1 - alpha) * Kd * 0.5 * Acyl

    Kin_abs = (Kb_abs + Kd_abs)

    # same as original
    Kup_abs = CNR_Kup(alpha, Kup, L, D)
    Lin_abs = LinMeas_abs(L, D, Lin, emis)
    Lup_abs = LupMeas_abs(L, D, Lup, emis)
    Acyl = CRT_Acyl(L,D)
    # Total_CNRRabs_solweig = np.dot(Aeff,((Kin_abs + Kup_abs + Lin_abs + Lup_abs) / (Acyl)))
    Total_CNRRabs_solweig = Aeff * ((Kin_abs + Kup_abs + Lin_abs + Lup_abs) / Acyl)

    return Total_CNRRabs_solweig, zen/(np.pi/180), Kd


def CNR_Kinabs_meas(alpha,Kin,L,D,A,lat,d,t,Atr):
    #VERTICAL CYLINDER MODEL USING THE PERPENDICULAR DIRECT IRRADIANCE AND
    #DIFFUSE BEAMS SEPARATELY
    #Used to determine the total incoming solar radiation (K) absorbed by the
    #human using the measured CNR1 net radiometer readings (W/m2), with inputs alpha (albedo of the cylinder), Kin (incoming solar
    #radiation measured by the net radiometer), L (length of the cylinder, cm),
    #D (diameter of the cylinder, cm), A is altitude in m, lat (latitude,
    #degrees), d (day, Jan 1 =1), and t (time, in decimal 24 hours)
    
    deg2rad = np.pi/180.
    Kb=Ratio_Kb(Kin,A,lat,d,t,Atr)
    Kd= Kin - Kb
    zen=solar_zenith(lat,d,t)
    Acs=CRT_Acs(L,D)
    Kb_abs=np.dot(np.dot((1 - alpha),(np.dot(Kb,np.sin(zen*deg2rad)))),Acs)
    Acyl=CRT_Acyl(L,D)
    Kd_abs=np.dot(np.dot(np.dot((1 - alpha),Kd),0.5),Acyl)
    Kin_abs=(Kb_abs + Kd_abs)

    return Kin_abs, zen, Kd
    
    # Written August, 2006 by Natasha Kenny
    # Revised January, 2007 to break up the diffuse and beam radiation
    # updated Jan 2016 Aaron Hardin & Jenni Vanos with explanations from Terry
    # Gillespie for ratio model

def Ratio_Kb(Kin,A,lat,d,t,Atr):
    # function is used to estimate incoming shortwave diffuse radiation under
    # clear sky conditions with inputs A(alititude, m), lat(latitude, degrees),
    # d (days, Jan 1 =1), t (time, 24 hours)
    
    m=opt_m(A,lat,d,t)
    Kb=Kin / ((1 + (np.dot(0.3,(1 - Atr ** m)))) / (Atr ** m))

    return Kb

def opt_m(A,lat,d,t):
    #optical air mass number used for input into the equations which estimate
    #incoming shortwave radiation under clear sky conditions
    deg2rad = np.pi/180.
    Pa=atm_P(A)
    zen=solar_zenith(lat,d,t)
    m=Pa / (np.dot(101.3,np.cos(zen*deg2rad)))

    # written January/07 by Natasha Kenny

    return m


def atm_P(A):

    # Written by Jenni Vanos, Aug 2009 for use in COMFA CNR to find Rabs.
    
    Po=101.3 #TODO move out? Pressure? 
    Pa=np.dot(Po,np.exp(- A / 8200))
    
    return Pa


def solar_zenith(lat,d,t):
    #function to calculate the cosine of the zenith angle based on the time of day, the
    #latitude of the site and the time of year. t is in hours ranging from 0 to
    #24.  d is the calender day with J = 1 at January 1. LC is the longitude
    #correction (+7.5 to -7.5 either side of standard meridian).
        
    #Solar declination ranges from +23.45 at summer soltice to -23.45 at winter
    #s olstice and is calculated as:
    deg2rad = np.pi/180.
    #make sure to change latitudinal correction for whatever city you are in.
    
    LC=-49/60. #TODO Remove since is should vary with location. What us this?

    ET=solar_ET(d)
    to=12 - LC - ET
    dec=solar_dec(d)
    
    cos_z=(np.dot(np.sin(lat*deg2rad),np.sin(dec*deg2rad))) + (np.dot(np.dot(np.cos(lat*deg2rad),np.cos(dec*deg2rad)),np.cos(np.dot(15.0,(t - to))*deg2rad)))
    zen=np.arccos(cos_z)*(180/np.pi)
    zen[zen>90]=90.

    # bad=np.where(zen > 90)
    # zen[bad]=90

    #Written January 2007 by Natasha Kenny

    return zen


def solar_ET(d):
    # used to calculate the equation of time for use in calcuating the zenith
    # angle of the sun.
    deg2rad = np.pi/180.
    f=279.575 + np.dot(0.9856,d)
    ET=(np.dot(- 104.7,np.sin(f*deg2rad)) + np.dot(596.2,np.sin(np.dot(2.0,f)*deg2rad)) + np.dot(4.3,np.sin(np.dot(3.0,f)*deg2rad)) - np.dot(12.7,np.sin(np.dot(4.0,f)*deg2rad)) - np.dot(429.3,np.cos(f*deg2rad)) - np.dot(2.0,np.cos(np.dot(2.0,f)*deg2rad)) + np.dot(19.3,np.cos(np.dot(3.0,f)*deg2rad))) / 3600.

    return ET


def solar_dec(d):
    # written to calcuate the solar declination with ranges from +23.45 at
    # summer solstice to -23.45 at the winter solstice.  d is the calendar day
    # with d = 1 being January 1st.
    
    #sin_dec = 0.39785 .* sind(278.97 + 0.9856 .* d + 1.9165 .* sind(356.6 +
    #0.9856 .* d)); from Campbell and norman pg. 168
    
    #dec = asind(sin_dec);
    deg2rad = np.pi/180.
    dec=np.dot(- 23.4,np.cos((np.dot(360,(d + 10)) / 365)*deg2rad))

    return dec


def CRT_Acs(L,D):
    # vertical cross sectional area of the cylinder in m
    
    # Written by Jenni Vanos, Aug 2009, to find the XC area of a cylinder,
    # where L = length (m) and D = diameter (m)- #for CRT L = ~0.10m and D =
    # 0.01m - ratio is most important part to keep constant for changing
    # cylidner size. 

    #Used in Kb_abs = (1-alpha).* (Kb .* sind(zen)).* Acs;
    
    Acs=np.dot(D,L)

    return Acs


def CNR_Kup(alpha,Kup,L,D):
    # Used to calculate the total reflected solar radiation (K) absorbed by the
    # a cylinder with inputs alpha (cylinder albedo), Kup (measured reflected solar
    # radiation by the CNR net radiometer (W/m2), Acyl (area of the cylinder, m2)
    
    Acyl=CRT_Acyl(L,D)
    Kup_abs=np.dot(np.dot(np.dot((1 - alpha),Kup),0.5),Acyl)

    return Kup_abs


def CRT_Acyl(L,D):
    #Function written August 2009 by Jenni Vanos for calculating the Area of
    #the CRT for use in the CNR Comfa model for finidng Rabs.
    
    #L = length (m) - #for CRT this is ~0.10m
    #D = diamater (m) - #for CRT this is ~0.01m
    
    # pi=3.14
    Acyl=np.dot(np.dot(2,np.pi),(D / 2) ** 2) + np.dot(np.dot(np.dot(2.0,np.pi),(D / 2)),L)

    return Acyl


def LinMeas_abs(L,D,Lin,e):
    #Used to Calculate the incoming longwave radiation absorbed by the cylinder with
    #inputs epsilon (emissivity of the cylinder, 0.95), Ta (air temp, degreesC)
    #L (length of cylinder, cm), D (diameter of cylinder, m)
    
    epsilon= e #0.95 moved to main function
    Acyl=CRT_Acyl(L,D)
    LinMeas_abs=np.dot(np.dot(np.dot(epsilon,Lin),0.5),Acyl)

    return LinMeas_abs


def LupMeas_abs(L, D, Lup, e):

    #Used to calculate the estimated outgoing terrestrial longwave absorbed by the
    #cylinder  for use in the Rabs (Estimation) Model with inputs Ta (air temperature),
    #L (length of cylinder, m), D (diameter of cylinder, cm)
    
    Acyl = CRT_Acyl(L,D)
    epsilon = e # 0.95 #mover to main function
    LupMeas_abs = np.dot(np.dot(np.dot(epsilon,Lup),0.5),Acyl)

    return LupMeas_abs

def day_of_year(yyyy, month, day):
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

        doy = np.sum(dayspermonth[0:month - 1]) + day

        return doy

def sun_distance(jday):
    """
    #% Calculatesrelative earth sun distance
    #% with day of year as input.
    #% Partridge and Platt, 1975
    """
    b = 2.*np.pi*jday/365.
    D = np.sqrt((1.00011+np.dot(0.034221, np.cos(b))+np.dot(0.001280, np.sin(b))+np.dot(0.000719,
                                        np.cos((2.*b)))+np.dot(0.000077, np.sin((2.*b)))))
    return D