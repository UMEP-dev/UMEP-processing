# Generated with SMOP  0.41
from libsmop import *
# LupMeas_abs.m

    
@function
def LupMeas_abs(L=None,D=None,Lup=None,*args,**kwargs):
    varargin = LupMeas_abs.varargin
    nargin = LupMeas_abs.nargin

    #Used to calculate the estimated outgoing terrestrial longwave absorbed by the
#cylinder  for use in the Rabs (Estimation) Model with inputs Ta (air temperature),
#L (length of cylinder, m), D (diameter of cylinder, cm)
    
    Acyl=CRT_Acyl(L,D)
# LupMeas_abs.m:7
    
    epsilon=0.95
# LupMeas_abs.m:9
    
    
    LupMeas_abs=multiply(multiply(multiply(epsilon,Lup),0.5),Acyl)
# LupMeas_abs.m:11
    #Written August, 2006 by Natasha Kenny
#Updated/Commented Dec 2015 Jenni Vanos