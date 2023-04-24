# Generated with SMOP  0.41
from libsmop import *
# CRT_Acyl.m

    
@function
def CRT_Acyl(L=None,D=None,*args,**kwargs):
    varargin = CRT_Acyl.varargin
    nargin = CRT_Acyl.nargin

    #Function written August 2009 by Jenni Vanos for calculating the Area of
#the CRT for use in the CNR Comfa model for finidng Rabs.
    
    #L = length (m) - #for CRT this is ~0.10m
#D = diamater (m) - #for CRT this is ~0.01m
    
    pi=3.14
# CRT_Acyl.m:10
    Acyl=multiply(multiply(2,pi),(D / 2) ** 2) + multiply(multiply(dot(2.0,pi),(D / 2)),L)
# CRT_Acyl.m:12