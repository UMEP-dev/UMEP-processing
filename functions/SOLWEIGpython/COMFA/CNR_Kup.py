# Generated with SMOP  0.41
from libsmop import *
# CNR_Kup.m

    
@function
def CNR_Kup(alpha=None,Kup=None,L=None,D=None,*args,**kwargs):
    varargin = CNR_Kup.varargin
    nargin = CNR_Kup.nargin

    # Used to calculate the total reflected solar radiation (K) absorbed by the
# a cylinder with inputs alpha (cylinder albedo), Kup (measured reflected solar
# radiation by the CNR net radiometer (W/m2), Acyl (area of the cylinder,
# m2)
    
    Acyl=CRT_Acyl(L,D)
# CNR_Kup.m:8
    
    Kup_abs=multiply(multiply(multiply((1 - alpha),Kup),0.5),Acyl)
# CNR_Kup.m:10
    
    # Written August, 2006 by Natasha Kenny
# Updated Comments Jan 2016 Jenni Vanos