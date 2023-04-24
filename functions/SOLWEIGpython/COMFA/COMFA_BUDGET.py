import numpy as np

def COMFA_Mact(weight, height, sex, age, activity, activity_unit):
    # Weight of  kid (kg)
    # weight = 20
    # Height of  kid (cm)
    # height = 114
    # Sex/gender of kid, 1 = boy, 2 = girl
    # sex = 1
    # Age of kid
    # age = 6
    # activity (MET), activity if kid
    # activity = 3.8
    # activity time, time of activity in minutes
    # activity_time = 30

    ### METABOLIC HEAT CALCULATION FOR CHILDREN ###
    # Body surface area of child or adult (different for infant according to Haylock et al. (1978))
    BSA = weight**0.5378 * height**0.3964 * 0.0242
    
    if activity_unit == 'MET':

        if sex == 2:
            if ((age >= 3) and (age <= 10)):
                # Resting metabolic rate of a 3-10 year old girl (kJ/day) (Cheng & Brown, 2020)
                RMR_J = 16.97*weight + 1.618*height + 371.2
            elif ((age > 10) and (age <= 18)):
                # Resting metabolic rate of a 10-18 year old girl (kJ/day) (Cheng & Brown, 2020)
                RMR_J = 8.365*weight + 4.65*height + 200
            else:
                RMR_J = 8.365*weight + 4.65*height + 200
        elif sex == 1:
            if ((age >= 3) and (age <= 10)):
                # Resting metabolic rate of a 3-10 year old boy (kJ/day) (Cheng & Brown, 2020)
                RMR_J = 19.6*weight + 1.033*height + 414.9
            elif ((age > 10) and (age <= 18)):
                # Resting metabolic rate of a 10-18 year old boy (kJ/day) (Cheng & Brown, 2020)
                RMR_J = 16.25*weight + 1.372*height + 515.5
            else:
                RMR_J = 8.365*weight + 4.65*height + 200
        # Convert to per hour. Original (kJ/day) = 24 hours.
        RMR_J = RMR_J / 24
        # Resting metabolic rate in Wm-2 for person
        RMR = RMR_J/BSA
        # Total metabolic energy by a person
        Mact = RMR * activity # * activity_time
        # Total metabolic energy by a person for PET
        Mact_PET = Mact * BSA
    elif activity_unit == 'W':
        Mact = activity / BSA
        Mact_PET = activity

    return Mact, Mact_PET

def COMFA_e(Ta, RH):
    # Tetens equation to calculate saturated vapour pressure (kPa) ata given air temperaure, 
    # with inputs Ta (air temperature, degrees C).
    # From Natasha Kenny and re-written by Jenni Vanos, April 2009, translated from MATLAB to Python by Nils Wallenberg 2022.
    # Saturated Vapour Pressure (kPa)
    es = 0.61354 * np.exp((17.502 * Ta) / (240.97 + Ta))
    # Ambient Vapour Pressure (kPa)
    ea = (RH / 100) * es

    # Saturation vapour pressure
    # es = 6.112 * np.exp( ( 17.67 * Ta ) / ( Ta + 243.5 ) )
    # Ambient vapour pressure
    # ea = (RH / 100) * es

    return ea

def COMFA_MET(Mact, Ta, RH):    
    # Ambient vapour pressure
    ea = COMFA_e(Ta, RH)    
    # Heat consumed when breathing
    f = 0.150 - 0.0173 * ea - 0.0014 * Ta

    MET = (1 - f) * Mact

    return MET

def COMFA_rc(va, rco):
    # Calculates the clothing resistance for 
    # input for the COMFA energy balance equation with inputs rco (insulation
    # value of clothing, s/m), va (activity wind velocity in m/s), vw
    # (windspeed, m/s). Note that rco must be known based on clothing worn. 

    if va == 0:
        rc = rco
    else:
        rc = rco * (-0.37 * (1 - np.exp(-va / 0.72)) + 1)
        # rc = rco * (-0.37 * (1 - np.exp(-va / 0.72)) + 1)

    # This is new equation from Kenny et al. (2009b) - IJB.

    # Alternative method for standing individual when clothing isn't known: -
    # UTCI regression equation (Havenith 2012 IJB). 
    # 
    # Calculates the clothing resistance for input for the COMFA energy balance 
    # equation with inputs of Tair (degrees C) to calculate the Icl (clo)based 
    # on air temperature from the UTCI equation. See Vanos, 2012b. 
    # For use on standing individuals only. Can be slowly walking or standing.

    # Written Oct 2009 by Jenni Vanos
    # Because we are assuming little to no activity speed (va) then rc = rco.

    # Icl = 1.372 - (0.01866*Tair) - (0.0004849*(Tair.^2)) - (0.000009333*(Tair.^3)); %Icl here is in clo.
    # Icl equation is from UTCI, and more info can be found in Vanos et al.
    # (2012b). 
    # rc = Icl*186.6; %186.6 is the conversion to put a clo into a s/m.

    return rc

def COMFA_vr(vw, va):
    # Written Oct 2009 by Jenni Vanos
    # This equation calculates the relative windspeed (vr) in m/s using the
    # activity windspeed (va) and windspeed (vw).

    vr = np.sqrt(va**2 + vw**2)

    # To acount for wind direction and activity direction, see Vanos et al.
    # (2012b) using angular equation from ISO for Vr. 

    return vr

def COMFA_ra(vw, va):
    # COMFA equation input under Convective Heat Loss or Gain (CONV)

    # Calculates the resistance of the boundary air layer around the body (s/m) with
    # inputs vr (Relative Wind velocity, m/s)
    # vw (wind speed, m/s), va (activity speed, m/s)

    # Thermal diffusivity of air (m2/s)
    k = 22e-6
    # kinematic viscosity of air (m2/s)
    v = 1.5e-5
    # diameter of cylinder as a model of a human (m)(Campbell, 1977)
    D = 0.17
    # Prandlt number
    Pr = 0.71

    vr = COMFA_vr(vw, va)
        
    Re =  vr * D / v

    if Re < 4000:
        A = 0.683
        n = 0.466
    elif ((Re >= 4000) and (Re < 40000)):
        A = 0.193
        n = 0.618
    else: # Re >= 40000
        A = 0.0266
        n = 0.805

    ra = 0.17 / (A * Re**n * Pr**0.33 * k)

    # A = Re
    # n = Re
    
    # L = Re < 4000
    # M = ((Re >= 4000) and (Re < 40000))
    # H = Re >= 40000
        
    # A[L] = 0.683  
    # n[L] = 0.466 

    # A[M] = 0.193
    # n[M] = 0.618
    
    # A[H] = 0.0266
    # n[H] = 0.805
    
    # Written by N. Kenny April 2006. Translated from MATLAB to Python by Nils Wallenberg 2022.

    return ra

def COMFA_Es(Mact, Ta, RH, age, kid):
    # Calculates sensible evaporative heat losses through Perspiration (W/m2)
    # with inputs Tair, Mact, RH

    MET = COMFA_MET(Mact, Ta, RH)

    if kid:
        Es = 0.42 * (MET - 58) * 0.5
    else:
        Es = 0.42 * (MET - 58)

    # Equation from Fanger (1970) for "comfortable" persons. 
    # Es = 0 if Met < 58.   

    # Written by Natasha Kenny, June 2006 for input into the COMFA energy
    # balance model under Evaporative Heat Loss (EVAP)
    # Now used by J. Vanos for Master's project, 2009. 

    # Based on data presented by Fanger (1972), Brown and
    # Gillespie (1986) suggest that body tissue resistance
    # decreases with Mact.

    return Es

def COMFA_rsk(Mact, Ta, RH, age, kid):
    # Same as rt (or skin tissue resistance) in papers. 
    # New equation added to COMFA model, from Kenny et al. (2009b in IJB).
    # Adapted April 2009 by Jenni Vanos for input into the COMFA equation under
    # Convective Heat Loss or Gain (CONV)

    # calculates the resistance (s/m) of the skin (body tissue) to heat transfer in
    # inputs: rho = density of air = ~1.16kg/m3; Cp = heat capacity of air =
    # 1000J/kg K), Es = sensible evaporative heat losses.

    # Evaporatve heat loss through perspiration. 
    Es = COMFA_Es(Mact, Ta, RH, age, kid) 
    # 1212 is volumetric heat capacity (or rho.*Cp)
    rsk = (1212) / (0.13 * Es + 15)

    # See Kenny et al. (2009b) pg 434 for derivation of rsk based on Kerslake 1972. 

    return rsk

def COMFA_Tc(Mact, Ta, RH):
    # Calculates core body temperature with inputs Mact (Metabolic Activity, W/m2),
    # Tair (air temperature in C), RH (relative humdity, %)

    MET = COMFA_MET(Mact, Ta, RH)

    Tc = 36.5 + (0.0043 * MET)

    # Written April 2009 by Jenni Vanos from Natasha Kenny, June 2006, for input into the COMFA energy balance 
    # model under Convective Heat Gain or Loss (CONV). See Kenny et al. (2009b)
    # in IJB for details. 

    # note this is not created for high intensity activity. 

    return Tc

def COMFA_Tsk(Mact, Ta, RH, vw, va, rco, age, kid):
    # Calculates skin temperature for input into the COMFA energy balance
    # equation with inputs Mact (metabolic activity, W/m2), Tair (air temp,
    # degrees C), rco (clothing resistance s/m), Vr (relative air
    # velocity m/s)

    # s m-1
    rsk = COMFA_rsk(Mact, Ta, RH, age, kid); 
    # s m-1
    ra = COMFA_ra(vw, va)
    # degC
    Tc = COMFA_Tc(Mact, Ta, RH)
    # s m-1
    rc = COMFA_rc(va, rco)
    # degC
    Tsk = (( Tc - Ta ) / (rsk + rc + ra)) * (ra + rc) + Ta

    # all resisteances are in s/m.
    
    # Equation for Tsk adapted from Kenny et al. (2009b) in IJB. Modified April 2009 by
    # Jenni Vanos. 

    return Tsk

def COMFA_CONV(Mact, Ta, RH, vw, va, rco, weight, height, age, kid):
    # Calculates the total heat loss (W/m2) by convection
    # with inputs Tair (air temperature, degrees C), Mact (metabolic activity,
    # W/m2), rco is the static clothing resistance (s/m),
    # vr (realtive wind velocity, m/s) (vr = vw when  va = 0) where va =
    # activity velocity, and vw = wind velocity

    # clothing resistance (Input should be va, not Ta as in original code?)
    rc = COMFA_rc(va, rco) 
    # resistance of boundary layer
    ra = COMFA_ra(vw, va)
    # Skin temperature
    Tsk = COMFA_Tsk(Mact, Ta, RH, vw, va, rco, age, kid) 
    # Final convection value 
    if kid:
        # Body to surface area of kid
        BSA_k = weight**0.5378 * height**0.3964 * 0.0242
        BMI_k = weight/(height/100)**2
        # Body-to-surface area of adult (weight = 65 kg and height = 176 cm, according to Cheng & Brown, 2020)
        adult_weight = 65
        adult_height = 176
        BSA_a = adult_weight**0.5378 * adult_height**0.3964 * 0.0242
        BMI_a = adult_weight/(adult_height/100)**2

        CONV = (1212 * (Tsk - Ta) / (rc + ra)) * ((BSA_k/weight)/(BSA_a/adult_weight))
    else:
        CONV = 1212 * (Tsk - Ta) / (rc + ra)

    # 1212 is just an  pre-calculated volumetric heat capapcity of air at 20oC ("rhoCp") of air (J m?3 K?1), 
    # assuming air density (rho) of ~1.2kg m-3 and Cp of 1000 J kg-1 K-1 at ~20oC,
    # which can be adapted based on air temperature 
    
    # Initally modified April, 2007 by Natasha Kenny to incorporate revisions of heat flow from surface
    # temp to air temp, rather than from core temp to air temp. See Kenny et al.
    # (2009a) IJB. 

    # Adapted Aug 2009 by Jenni Vanos in order to incorporate new Tsk
    # equation. Additional notes on rhoCp added by JV 2020.

    return CONV

def COMFA_rcv(vw, va, rcvo):
    # Written by Jenni Vanos April, 2009 for use in the updated COMFA model,
    # according to findings and changes made by Kenny et al. (2009). This is
    # equation 6 from Kenny et al. (2009b), and has inputs rcvo (or the static 
    # clothing vapour resistance, which is defined for standardized conditions 
    # (static body, wind still, i.e. speed < 0.2 m s-1) (ISO, 2007). When air 
    # movement is present, or when the body moves, this will affect the vapour 
    # resistance (typically lowering it), in which case it is referred to as the 
    # resultant or dynamic basic water vapour resistance (ISO, 2007).
    # rcvo can be calculated for input using the equation:

    # rcvo = Icl.*0.1555.*0.18.*18400; %ISO (2007) conversion from clo to s/m.

    # Where Icl is found in tables for specific clothing ensembles from ISO (2007), 
    # and 18400 is a conversion factor from Re,cl (vapour resistance, in m2kPaW-1) 
    # to rcvo, using Lv = 2.5×106 J kg-1, rho = 1.16 kg m-3, and Pa = 98kPa. 
    # See Kenny et al. (2009b) for further reference. 

    vr = COMFA_vr(vw, va)

    rcv = rcvo * (-0.80 * (1 - np.exp(-vr / 1.095))+1)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Alternative method for standing individual when clothing isn't known: -
    # UTCI regression equation (Havenith 2012 IJB). 
    
    # Calculates the clothing resistance for input for the COMFA energy balance 
    # equation with inputs of Tair (degrees C) to calculate the Icl (clo)based 
    # on air temperature from the UTCI equation. See Vanos, 2012b. 
    # For use on standing individuals only. Can be slowly walking or standing.

    # Written Oct 2009 by Jenni Vanos
    # Because we are assuming little to no activity speed (va) then rc = rco.

    # Icl = 1.372 - (0.01866*Tair) - (0.0004849*(Tair.^2)) - (0.000009333*(Tair.^3)); %Icl here is in clo.
    # Icl equation is from UTCI, and more info can be found in Vanos et al.
    # (2012b). 
    # rcvo = Icl.*0.1555.*0.18.*18400; %ISO (2007) conversion from clo to s/m.

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    return rcv

def COMFA_rav(vw, va):    
    # Calculates the resistance of the boundary air layer to vapour transport,
    # with input of relative wind velocity (m/s)

    ra = COMFA_ra(vw, va)
    # (s/m)
    rav = 0.92 * ra

    # Adapated April 2009 by Jenni Vanos to incorporate vr (m/s), which includes va and vw (m/s).
    # Units of rav are in s/m. 

    return rav

def COMFA_qsk(Mact, Ta, RH, vw, va, rco, age, kid):
    # Calculates the saturation specific humidity (kg/kg) at skin temperature for input into the
    # COMFA energy balance equation with inputs Mact (metabolic activity), Tair
    # (degrees C), RH (%), rco, va and Vr. 
    # Assuming atmospheric pressure (Pa) is 1013 mb (101.3kPa)

    Tsk = COMFA_Tsk(Mact, Ta, RH, vw, va, rco, age, kid)
    
    esk = 0.61365 * np.exp((17.502 * Tsk) / (240.97 + Tsk))

    qsk = 0.622 * (esk / (101.3 - esk))

    return qsk

def COMFA_qa(Ta, RH):
    # Calculates saturation specific humidity (kg/kg) at air temperatre for input
    # into the COMFA energy balance equation, with inputs Tair (air
    # temperature, degrees C), RH (relative humidity, %).
    # assuming atmospheric pressure is 1013 mb (101.3kPa)

    e = COMFA_e(Ta, RH)
    
    qa = 0.622 * e / (101.3 - e)
    
    # svp = Tetens_svp(Tair)
    #  
    # qa = svp.*(RH./100)
    
    # Written July, 2006 by Natasha Kenny, modified to reflect the changes to
    # this equation. See usage also in Kenny et al. (2009b) 
    # Re-written April 2009 by Jenni Vanos. 

    return qa

def COMFA_Ei(Mact, Ta, RH, vw, va, rco, rcvo, age, kid):
    # Calculates insensible evaporative heat losses in the COMFA energy balance
    # equation with inputs Tair (air temp, degrees C), RH (relative humidity %)
    # vw (wind velocity (m/s), va (acivity velocity), Mact (Metabolic activity, W/m2), 
    # rco(clothing resistance, s/m), rcvo (clothing vapour resistance s/m)

    # Density of air in kg/m3
    rho = 1.16
    # Latent heat of vapourization (J/kg)
    Lv = 2.442e+6
    # Resistance of skin to vapour transfer (s/m) - Campbell 1977
    rskv = 7.7e+3
    # clothing vapour resistance (s/m)
    rcv = COMFA_rcv(vw, va, rcvo)
    # Resistance of boundary air layer to vapour transport
    rav = COMFA_rav(vw, va)
    # update
    qsk = COMFA_qsk(Mact, Ta, RH, vw, va, rco, age, kid)
    # update
    qa = COMFA_qa(Ta, RH) 

    Ei = rho * Lv * ((qsk - qa) / (rcv + rav + rskv))
    
    # note Ei is very unimportant.
    
    # Written July, 2006 by Natasha Kenny, adapted April 2009 by Jenni
    # Vanos to incorporate changes with vr and rcv from changes made by
    # Kenny et al. (2009). 

    return Ei

def COMFA_Etot(Mact, Ta, RH, vw, va, rco, rcvo, age, kid):
    # Calculates the total evaporative heat loss for input into the COMFA
    # equation with inputs (Tair air temp in C, RH Relative Humidity %, vw, 
    # wind velocity (m/s), Mact, metabolic activity W/m2, rco (clothing resistance, s/m), va activity velocity (m/s), rcvo (clothing vapour resistance, s/m) not to be confused with
    # EVAP which is the final input in the equation

    Ei = COMFA_Ei(Mact, Ta, RH, vw, va, rco, rcvo, age, kid)

    Es = COMFA_Es(Mact, Ta, RH, age, kid)

    Etot = Ei + Es
    
    # Written July, 2006 by Natasha Kenny. Adapted April 2009 by Jenni Vanos to
    # include vr, and subtract P based on Kenny et al. (2009). 

    return Etot

def COMFA_Em(Mact, Ta, RH, vw, va, rco, rcvo, age, kid):
    # Calculates evaporative maximum, Em (W/m2) (where Em = the maximum possible evapoaration 
    # possible given the environment and clothing, and assumes a skin wettedness value (w) of 1)
    # for input into COMFA equation.
     
    # Inputs vr (relative wind velocity, m/s) (vr = vw when va=0), Mact (metabolic
    # activity), Tair (air temperature, degrees C), RH (relative humdity, %)
    
    # Density of air kg m-3
    rho = 1.16
    # latent heat of vaporization J kg-1
    Lv = 2.442e+6
    
    rav = COMFA_rav(vw, va)
    
    rcv = COMFA_rcv(vw, va, rcvo)
    
    qsk = COMFA_qsk(Mact, Ta, RH, vw, va, rco, age, kid)
    
    qa = COMFA_qa(Ta, RH)
       
    Em = rho * Lv * ((qsk - qa) / (rcv + rav))

    # Written by Natasha Kenny, June 2006 for input into the COMFA equation.
    # Adapted by Jenni Vanos April 2009 for addition of changes by Kenny et al.
    # (2009). (-P and + vr, va, vw)

    return Em

def COMFA_EVAP(Mact, Ta, RH, vw, va, rco, rcvo, age, kid):

    # Calculates the total evaporative heat losses for input into the COMFA
    # energy balance model with inputs Ei (insensible heat loss through the skin
    # and Es (sensible heat losses through the skin) W/m2

    Etot = COMFA_Etot(Mact, Ta, RH, vw, va, rco, rcvo, age, kid)

    Em = COMFA_Em(Mact, Ta, RH, vw, va, rco, rcvo, age, kid)
 
    if Etot <= Em:
        EVAP = Etot
    else:
        EVAP = Em

    # EVAP = np.zeros((Em.shape[0]))
    
    # ind = Etot <= Em
 
    # EVAP[ind] = Etot[ind]

    # The lower of E or Em is used in the energy budget equation.
    
    # Initially written by Kenny, modified by Jenni Vanos, April 2009. 

    return EVAP

def COMFA_Ts(Mact, Ta, RH, vw, va, rco, age, kid):
    # Calculates the average surface temperature of a clothed person (degrees C) with inputs Tair (air
    # temperature, degrees C), rco (clothing resistance,  s/m), vw (air velocity, s/m), Mact (Metabolic Activity, W/m2)

    rc = COMFA_rc(va, rco)
    # Resistance of boundary layer (s/m)
    ra = COMFA_ra(vw, va)
    # skin temperature (degC)
    Tsk = COMFA_Tsk(Mact, Ta, RH, vw, va, rco, age, kid)
    # Average surface temperature of a person
    Ts = ((Tsk - Ta) / (rc + ra)) * ra + Ta

    # Adapted from Natasha on April 2009 to include vr, va and take out P. 
    # Modified April, 2007 by Natasha Kenny to revise the heat flow from the skin surface to the
    # air rather than from the core temperature to the air.

    return Ts

def COMFA_TREMIT(Mact, Ta, RH, vw, va, rco, age, kid):
    # Calculates the terrestrial radiation losses from the body for input into the COMFA
    # energy balance equation with inputs of Tair (air temp, degrees C), rco (clothing
    # insulation (s/m)), wind velocity (m/s), activity velocity (m/s)), and Mact,
    # (metabolic activity (w/m2)).
    
    # surface temperature (C)
    Ts = COMFA_Ts(Mact, Ta, RH, vw, va, rco, age, kid)

    # 0.95 is the emissivity of human (Campbell & Norman 1998)
    # 0.78 is the radiative surface area of a human standing (Kerslake, 1972)
    # 0.7 for sitting/on bicycle
    TREMIT = 0.78 * (0.95 * 5.67e-8 * (Ts + 273.15)**4) 

    # Adapted April 2009 by J. Vanos to incorporate new vr geomtric equation (vr =
    # COMFA_vr(va, vw)), and take out P, which is no longer needed in the COMFA model
    # (Kenny et al., 2009a,b). 

    # modified April, 2007 by Natasha Kenny to incorportate emissivity of 0.95 and radiative
    # surface of a human 0.78 (this is for when walking or running). 

    return TREMIT

def COMFA_BUDGET(Mact, Ta, RH, vw, va, rco, rcvo, weight, height, age, kid):
    # Calculates the four energy fluxes MET, CONV, EVAP, and TREMIT in the COMFA energy budget, with inputs 
    # Tair (air temperature, degrees C), Metabolic activity (W/m2), Relative Humidty (%),
    # Wind Velocity (m/s), activity velocity (va, m s-1), static clothing resistance (rco, sm-1), static clothing vapor reistance (rcvo, s m-1). 
    # Note can also make Rabs an input. See below.

    ### This COFMA model is general for people at or near comfort, and is not
    ### tested yet on extreme conditions. 

    # %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% 
    # %%%% CLOTHING - generally can estimate Icl value (clo) and convert them to
    # %%%% resistances as follows: 

    # Rco: rco = Icl * 186.6. Or convert from conductivity based on 1 Clo = 0.1555 (m 2K W-1).		

    # NOTE: Icl values can be found in clo or conductivity values from ISO 9920 2007. 

    # Rcvo: Can convert from a conductivity value of Icl, where Icl is found in tables for specific clothing ensembles from ISO (2007). 
    # 1) convert from Icl to Re,cl (m2 kPa W-1):  Re,cl = Icl (m 2K W-1)*0.18 (constant from ISO 9920 pg 12 based on 1 or 2-layer clothing ensembles). 
    # 2) Convert Re,cl to rcvo:  rcvo = Re,cl*18,400, where 18400 is a conversion factor from Re,cl (vapour resistance, in m2kPaW-1) 		
    # to rcvo, using Lv = 2.5*106 J kg-1, rho = 1.16 kg m-3, and Pa = 98kPa. 		
    # %%%%%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% 
    
    # W m-2 %%%%random value ... Here you need to link to another dataset that calcs Rabs or bring o link to more doce. 
    # Rabs = 400
    
    MET = COMFA_MET(Mact, Ta, RH)

    CONV = COMFA_CONV(Mact, Ta, RH, vw, va, rco, weight, height, age, kid)

    EVAP = COMFA_EVAP(Mact, Ta, RH, vw, va, rco, rcvo, age, kid)
    
    TREMIT = COMFA_TREMIT(Mact, Ta, RH, vw, va, rco, age, kid)
    
    # BUDGET = MET + Rabs - CONV - EVAP - TREMIT

    # all are in W m-2

    # return BUDGET

    return MET, CONV, EVAP, TREMIT