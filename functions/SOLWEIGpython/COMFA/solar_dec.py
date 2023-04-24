# Generated with SMOP  0.41
from libsmop import *
# solar_dec.m

    
@function
def solar_dec(d=None,*args,**kwargs):
    varargin = solar_dec.varargin
    nargin = solar_dec.nargin

    # written to calcuate the solar declination with ranges from +23.45 at
# summer solstice to -23.45 at the winter solstice.  d is the calendar day
# with d = 1 being January 1st.
    
    #sin_dec = 0.39785 .* sind(278.97 + 0.9856 .* d + 1.9165 .* sind(356.6 +
#0.9856 .* d)); from Campbell and norman pg. 168
    
    #dec = asind(sin_dec);
    
    dec=dot(- 23.4,cosd(dot(360,(d + 10)) / 365))
# solar_dec.m:13
    
    # written January, 2007 by Natasha Kenny