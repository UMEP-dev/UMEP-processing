# Generated with SMOP  0.41
from libsmop import *
# solar_zenith.m

    
@function
def solar_zenith(lat=None,d=None,t=None,*args,**kwargs):
    varargin = solar_zenith.varargin
    nargin = solar_zenith.nargin

    #function to calculate the cosine of the zenith angle based on the time of day, the
#latitude of the site and the time of year. t is in hours ranging from 0 to
#24.  d is the calender day with J = 1 at January 1. LC is the longitude
#correction (+7.5 to -7.5 either side of standard meridian).
    
    #Solar declination ranges from +23.45 at summer soltice to -23.45 at winter
#solstice and is calculated as:
    
    #make sure to change latitudinal correction for whatever city you are in.
    LC=- 49 / 60
# solar_zenith.m:12
    
    
    ET=solar_ET(d)
# solar_zenith.m:14
    to=12 - LC - ET
# solar_zenith.m:16
    dec=solar_dec(d)
# solar_zenith.m:18
    
    cos_z=(multiply(sind(lat),sind(dec))) + (multiply(multiply(cosd(lat),cosd(dec)),cosd(dot(15.0,(t - to)))))
# solar_zenith.m:20
    zen=acosd(cos_z)
# solar_zenith.m:22
    bad=find(zen > 90)
# solar_zenith.m:24
    zen[bad]=90
# solar_zenith.m:25
    #Written January 2007 by Natasha Kenny
    