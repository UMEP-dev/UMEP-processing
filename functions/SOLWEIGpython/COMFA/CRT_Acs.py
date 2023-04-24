# Generated with SMOP  0.41
from libsmop import *
# CRT_Acs.m

    
@function
def CRT_Acs(L=None,D=None,*args,**kwargs):
    varargin = CRT_Acs.varargin
    nargin = CRT_Acs.nargin

    # vertical cross sectional area of the cylinder in m
    
    # Written by Jenni Vanos, Aug 2009, to find the XC area of a cylinder,
# where L = length (m) and D = diameter (m)- #for CRT L = ~0.10m and D =
# 0.01m - ratio is most important part to keep constant for changing
# cylidner size. 
# 
#Used in Kb_abs = (1-alpha).* (Kb .* sind(zen)).* Acs;
    
    Acs=multiply(D,L)
# CRT_Acs.m:13