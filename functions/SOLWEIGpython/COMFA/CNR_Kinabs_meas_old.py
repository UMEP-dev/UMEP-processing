# Generated with SMOP  0.41
from libsmop import *
# CNR_Kinabs_meas.m

    
@function
def CNR_Kinabs_meas(alpha=None,Kin=None,L=None,D=None,A=None,lat=None,d=None,t=None,Atr=None,*args,**kwargs):
    varargin = CNR_Kinabs_meas.varargin
    nargin = CNR_Kinabs_meas.nargin

    #VERTICAL CYLINDER MODEL USING THE PERPENDICULAR DIRECT IRRADIANCE AND
#DIFFUSE BEAMS SEPARATELY
#Used to determine the total incoming solar radiation (K) absorbed by the
#human using the measured CNR1 net radiometer readings (W/m2), with inputs alpha (albedo of the cylinder), Kin (incoming solar
#radiation measured by the net radiometer), L (length of the cylinder, cm),
#D (diameter of the cylinder, cm), A is altitude in m, lat (latitude,
#degrees), d (day, Jan 1 =1), and t (time, in decimal 24 hours)
    
    Kb=Ratio_Kb(Kin,A,lat,d,t,Atr)
# CNR_Kinabs_meas.m:11
    
    Kd=Kin - Kb
# CNR_Kinabs_meas.m:13
    
    zen=solar_zenith(lat,d,t)
# CNR_Kinabs_meas.m:15
    
    Acs=CRT_Acs(L,D)
# CNR_Kinabs_meas.m:17
    
    Kb_abs=multiply(multiply((1 - alpha),(multiply(Kb,sind(zen)))),Acs)
# CNR_Kinabs_meas.m:19
    
    Acyl=CRT_Acyl(L,D)
# CNR_Kinabs_meas.m:21
    
    Kd_abs=multiply(multiply(multiply((1 - alpha),Kd),0.5),Acyl)
# CNR_Kinabs_meas.m:23
    
    Kin_abs=(Kb_abs + Kd_abs)
# CNR_Kinabs_meas.m:25
    
    # Written August, 2006 by Natasha Kenny
# Revised January, 2007 to break up the diffuse and beam radiation
# updated Jan 2016 Aaron Hardin & Jenni Vanos with explanations from Terry
# Gillespie for ratio model