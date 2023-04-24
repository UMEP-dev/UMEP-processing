# Generated with SMOP  0.41
from libsmop import *
# solar_ET.m

    
@function
def solar_ET(d=None,*args,**kwargs):
    varargin = solar_ET.varargin
    nargin = solar_ET.nargin

    # used to calculate the equation of time for use in calcuating the zenith
# angle of the sun.
    
    f=279.575 + multiply(0.9856,d)
# solar_ET.m:7
    ET=(multiply(- 104.7,sind(f)) + multiply(596.2,sind(dot(2.0,f))) + multiply(4.3,sind(dot(3.0,f))) - multiply(12.7,sind(dot(4.0,f))) - multiply(429.3,cosd(f)) - multiply(2.0,cosd(dot(2.0,f))) + multiply(19.3,cosd(dot(3.0,f)))) / 3600
# solar_ET.m:9