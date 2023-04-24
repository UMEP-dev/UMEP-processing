# Generated with SMOP  0.41
from libsmop import *
# opt_m.m

    
@function
def opt_m(A=None,lat=None,d=None,t=None,*args,**kwargs):
    varargin = opt_m.varargin
    nargin = opt_m.nargin

    #optical air mass number used for input into the equations which estimate
#incoming shortwave radiation under clear sky conditions
    
    Pa=atm_P(A)
# opt_m.m:6
    zen=solar_zenith(lat,d,t)
# opt_m.m:8
    m=Pa / (multiply(101.3,cosd(zen)))
# opt_m.m:10
    # written January/07 by Natasha Kenny