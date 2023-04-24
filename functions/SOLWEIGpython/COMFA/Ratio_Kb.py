# Generated with SMOP  0.41
from libsmop import *
# Ratio_Kb.m

    
@function
def Ratio_Kb(Kin=None,A=None,lat=None,d=None,t=None,Atr=None,*args,**kwargs):
    varargin = Ratio_Kb.varargin
    nargin = Ratio_Kb.nargin

    # function is used to estimate incoming shortwave diffuse radiation under
# clear sky conditions with inputs A(alititude, m), lat(latitude, degrees),
# d (days, Jan 1 =1), t (time, 24 hours)
    
    m=opt_m(A,lat,d,t)
# Ratio_Kb.m:8
    Kb=Kin / ((1 + (multiply(0.3,(1 - Atr ** m)))) / (Atr ** m))
# Ratio_Kb.m:10