# Generated with SMOP  0.41
# from libsmop import *
# CNRRabs_Total.m

    
# @function
def CNRRabs_Total(alpha=None,L=None,D=None,Lin=None,Lup=None,Kin=None,Kup=None,Atr=None,d=None,t=None,lat=None,A=None,*args,**kwargs):
    varargin = CNRRabs_Total.varargin
    nargin = CNRRabs_Total.nargin

    # VERTICAL CYLINDER MODEL 
    # Used to calculate the total radiation absorbed by a cylinder, such as the CRT, with inputs
    # alpha (albedo of cylinder), Kin (incoming shortwave measured by CNR net
    # radiometer, W/m2), Kup (reflected shortwave measured by CNR, W/m2), Lin
    # (CNR incoming longwave radiation, W/m2), Lup (CNR1 outgoing longwave
    # radition, W/m2), L (length of cylinder, m), D (diamter of cylinder, m).
    # 
    #Also need solar zenith, which is calculated from: lat (latitude, degrees), d (days, Jan 1=1), t (decimal
    # time, 24h decimal time).
        
    #Atr - atmospheric transmittance - can calculate for the day/time of day or
    #estimate based on sky conditions.
        
    # LOAD data ....
    
    #set constants
    D=0.01
    # CNRRabs_Total.m:20
    
    L=0.1
    # CNRRabs_Total.m:21
    
    alpha=0.37
    # CNRRabs_Total.m:22
    
    #     E.g., 
    #     Kin = metdata(:,10); #put in column for Kin from metdata file
    #     Kup = metdata(:,11); 
    #     Lin  = metdata(:,20);
    #     Lup = metdata(:,21); 
    #     d = metdata(:,24);
    #     Atr = 0.70; 
    #     lat = 33.75;
    #     A = 992;
    
    Aeff=0.78
    # CNRRabs_Total.m:33
    
    Kin_abs=CNR_Kinabs_meas(alpha,Kin,L,D,A,lat,d,t,Atr)
    # CNRRabs_Total.m:35
    Kup_abs=CNR_Kup(alpha,Kup,L,D)
    # CNRRabs_Total.m:37
    Lin_abs=LinMeas_abs(L,D,Lin)
    # CNRRabs_Total.m:39
    Lup_abs=LupMeas_abs(L,D,Lup)
    # CNRRabs_Total.m:41
    Acyl=CRT_Acyl(L,D)
    # CNRRabs_Total.m:43
    Total_CNRRabs=multiply(Aeff,((Kin_abs + Kup_abs + Lin_abs + Lup_abs) / (Acyl)))
    # CNRRabs_Total.m:46