# Generated with SMOP  0.41
from libsmop import *
# atm_P.m

    
@function
def atm_P(A=None,*args,**kwargs):
    varargin = atm_P.varargin
    nargin = atm_P.nargin

    # Written by Jenni Vanos, Aug 2009 for use in COMFA CNR to find Rabs.
    
    Po=101.3
# atm_P.m:6
    
    Pa=multiply(Po,exp(- A / 8200))
# atm_P.m:9
    
    #another equation that gives similar results found from website below:
    
    #Pa=Po.*(1-.0000225577.*A).^5.25588;
    
    # A = altitude, units don't matter because it cancels with Ao (which is 0 at surface). 
#Therefore, it will be exp^-1, and Po is surface pressure, which is 1bar, or 100kPa.
    
    # equation from: http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html