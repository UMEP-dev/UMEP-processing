# Generated with SMOP  0.41
from libsmop import *
# LinMeas_abs.m

    
@function
def LinMeas_abs(L=None,D=None,Lin=None,*args,**kwargs):
    varargin = LinMeas_abs.varargin
    nargin = LinMeas_abs.nargin

    #Used to Calculate the incoming longwave radiation absorbed by the cylinder with
#inputs epsilon (emissivity of the cylinder, 0.95), Ta (air temp, degreesC)
#L (length of cylinder, cm), D (diameter of cylinder, m)
    
    epsilon=0.95
# LinMeas_abs.m:7
    Acyl=CRT_Acyl(L,D)
# LinMeas_abs.m:9
    LinMeas_abs=multiply(multiply(multiply(epsilon,Lin),0.5),Acyl)
# LinMeas_abs.m:11
    
    #Written August, 2006 by Natasha Kenny
#Updated/Commented Dec 2015 Jenni Vanos