
def rpgen(strike, dip, rake, gamma, sigma, TKO, AZM):

    '''
    RPGEN Calculate radiation pattern using shear-tensile source model.
    rpgen(strike,dip,rake,gamma,sigma, TKO, AZM) calculates P-wave, S-wave,
    SH-wave and SV-wave radiation pattern using shear-tensile source model
    presented in [see references 1, 2, 3 for details]. All input angles 
    (strike, dip, rake of the fault, tensile angle gamma, takeoff angle 
    TKO and azimuth from the source to the observation point AZM) should 
    be in degrees. The takeoff angle is measure from bottom. The function 
    returns matrices of the same size as input TKO and AZM matrices. 

    Input parameters:

      strike, dip, rake: fault plane parameters (degrees).
      gamma:  tensile angle in degrees (0 degrees for pure shear faulting, 
              90 degrees for pure tensile opening).
      sigma:  Poisson's ratio.
      TKO:    matrix of takeoff angles for which to calculate the correspo-
              nding radiation pattern coefficients (degrees, the takeoff 
              angles are measured from bottom).
      AZM:    matrix of corresponding azimuths (in degrees) for which the 
              radiation pattern coefficients should be calculated.

    Output parameters:
       
      Gp, Gs, Gsh, Gsv - P-wave, S-wave, SH-wave, and SV-wave radiation 
      pattern coefficients calculated for corresponding takeoff angles 
      and azimuths specified in TKO and AZM matrices.

    References:

      [1] Kwiatek, G. and Y. Ben-Zion (2013). Assessment of P and S wave 
          energy radiated from very small shear-tensile seismic events in 
          a deep South African mine. J. Geophys. Res. 118, 3630-3641, 
          doi: 10.1002/jgrb.50274
      [2] Ou, G.-B., 2008, Seismological Studies for Tensile Faults. 
          Terrestrial, Atmospheric and Oceanic Sciences 19, 463.
      [3] Vavryèuk, V., 2001. Inversion for parameters of tensile 
          earthquakes.” J. Geophys. Res. 106 (B8): 1633916355. 
          doi: 10.1029/2001JB000372.

    Copyright 2012-2013 Grzegorz Kwiatek.
    $Revision: 1.3 $  $Date: 2013/09/15 $
    '''
    
    #Import required modules
    import numpy as np

    #Alter inputs
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    AZM, TKO = np.meshgrid(np.radians(AZM),np.radians(TKO))
    
    #P radiation pattern
    Gp = np.cos(TKO)*(np.cos(TKO)*(np.sin(gamma)*(2*np.cos(dip)**2 - (2*sigma)/(2*sigma - 1)) + np.sin(2*dip)*np.cos(gamma)*np.sin(rake)) - np.cos(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) + np.sin(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma))) + np.sin(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma)) + np.cos(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) + np.sin(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) - np.sin(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) + np.cos(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)));

    #S radiation pattern
    Gs = ((np.sin(AZM)*np.sin(TKO)*(np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma)) + np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(TKO)*(np.sin(TKO)*(np.sin(gamma)*(2*np.cos(dip)**2 - (2*sigma)/(2*sigma - 1)) + np.sin(2*dip)*np.cos(gamma)*np.sin(rake)) + np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) - np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma))) + np.cos(AZM)*np.sin(TKO)*(np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) + np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.sin(dip)**2*np.sin(strike)**2))))**2 + (np.cos(TKO)*(np.cos(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma)) + np.sin(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike))) - np.sin(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.cos(AZM)*(np.cos(gamma)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) + np.cos(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(gamma)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)) + np.cos(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma))))**2)**(1/2);

    #SH radiation pattern
    Gsh = np.cos(TKO)*(np.cos(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma)) + np.sin(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike))) - np.sin(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.cos(AZM)*(np.cos(gamma)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) + np.cos(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(gamma)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)) + np.cos(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)))

    #SV radiation pattern
    Gsv = np.sin(AZM)*np.sin(TKO)*(np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma)) + np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(TKO)*(np.sin(TKO)*(np.sin(gamma)*(2*np.cos(dip)**2 - (2*sigma)/(2*sigma - 1)) + np.sin(2*dip)*np.cos(gamma)*np.sin(rake)) + np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) - np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma))) + np.cos(AZM)*np.sin(TKO)*(np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) + np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)));

    return Gp, Gs, Gsh, Gsv

