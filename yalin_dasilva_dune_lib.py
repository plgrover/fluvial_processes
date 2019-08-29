#!/usr/bin/env python

import numpy as np
import math
import sed_trans


_nu = 1.0e-6
_g  = 9.81
_gamma_s = 16186.5
_kappa = 0.4

_rho_water = 1000.
_rho_particle = 2650.0


_maxiter = 100
_tol = 1.e-4

def reconstruct_flow_from_dunes2(D50, Lambda, Delta, rho_particule=_rho_particle):
    ubar_low = None
    S_low = None
    ubar_high = None
    S_high = None

    # Check if the current steepness is greater than the maximum steepness
    delta_d = Delta / Lambda
    h = get_water_depth_from_length(Lambda)
    cf = dunelib.get_cf(h, D50)
    c =  dunelib.get_total_chezy(h,D50, Lambda, Delta)
    Z = h / D50
    delta_d_max = get_delta_d_max(h/D50)
    slope = get_slope_at_max_steepness(h_test, D50, rho_particule)

    print('The flow depth: {0} m and Z: {1}'.format(h, Z))
    print('The cf = {0} and c {1} '.format(cf,c))
    print('The maximum dune steepness: {0} and the observed steepness: {1}'.format(delta_d_max, delta_d))

    if delta_d < delta_d_max:

        _tol = 1.e-8

        slope1 = slope * 0.9
        slope0 = slope * 2.05

        k = 0
        while k <= _maxiter and abs((slope1 - slope0)) >= _tol:
            print(slope1, slope0)
            h0 = get_water_depth_from_length(Lambda, D50, slope0)
            h1 = get_water_depth_from_length(Lambda, D50, slope1)

            fx0 = get_steepness(h0, slope0, D50, rho_particule)[0] - delta_d
            fx1 = get_steepness(h1, slope1, D50, rho_particule)[0] - delta_d

            slope = slope0 - fx0 * (slope0 - slope1) / (fx0 - fx1)

            slope1 = slope0
            slope0 = slope

            k = k + 1

        if k > _maxiter:
            print('Error: exceeded {0} iterations'.format(k))
    else:
        print('The steepness is greater than then maximum dune steepness')






def reconstruct_flow_from_dunes(D50, Lambda, Delta, rho_particule=_rho_particle):
    
    # Check if the current steepness is greater than the maximum steepness
    delta_d = Delta/Lambda
    h_test = Lambda / (6.)
    delta_d_max = get_delta_d_max(h_test / D50)
    
    slope = get_slope_at_max_steepness(h_test, D50, rho_particule)
    
    print('delta_d: {0}   delta_d_max: {1}'.format(delta_d,delta_d_max))
    
    print('Slope at max steepness: {0}'.format(slope))
    
    
    
    if delta_d < delta_d_max:
                
        _tol = 1.e-8
        
        slope1 = slope * 0.9
        slope0 = slope * 2.05
        
        k = 0
        while k <= _maxiter and abs( (slope1 - slope0) ) >= _tol:
            
            print(slope1,slope0) 
            h0 = Lambda/6. #get_water_depth_from_length(Lambda, D50, slope0)
            h1 = Lambda/6. #get_water_depth_from_length(Lambda, D50, slope1)
    
            fx0 = get_steepness(h0, slope0, D50, rho_particule)[0] - delta_d
            fx1 = get_steepness(h1, slope1, D50, rho_particule)[0] - delta_d

    
            slope = slope0 - fx0 * ( slope0 - slope1 ) / ( fx0 - fx1 )
            
            slope1 = slope0
            slope0 = slope
    
            k = k + 1
    
        if k > _maxiter:
            print('Error: exceeded {0} iterations'.format(k))
    else:
        print('The steepness is greater than then maximum dune steepness')
        
    h = get_water_depth_from_length(Lambda, D50, slope)
    return slope,h
    
    
def get_total_chezy(h, D50, Lambda, Delta):
    cf =   (1./_kappa)  * math.log(0.368 * h / (D50*2.)) + 8.5
    
    
    c =( (1./cf**2) + 0.5 * ( (Delta/Lambda)**2 ) * (Lambda/h) )**(-0.5)
    
    return c
    
    
def get_c_delta(h,Lambda, Delta):
    
    retval = 0.5 * ( (Delta/Lambda)**2 ) * (Lambda/h)
    
    return 
    
    
def get_cf(h,D50, Bs=8.5):
    cf =   (1./_kappa)  * math.log(0.368 * h / (D50*2.)) + Bs
    
    return cf
    
    

def get_steepness(h, slope, D50, rho_particule=_rho_particle):
    '''
    Eq. 2.16
    '''
    
    Z = h/D50    
    X = sed_trans.get_X(h,slope,D50)
    
    phi_d = 1. - math.exp( -1. * (X/10.)**2  )
    
    m_gamma = get_m_gamma(Z)
    nu_star_d = get_nu_star_d(Z)
    nu_star = sed_trans.get_nu_star(slope, h, D50, rho_particule)   
    zeta_d = (nu_star-1)/(nu_star_d - 1) # Eq 2.17
    delta_d_max = get_delta_d_max(Z)

    delta_d = delta_d_max * (zeta_d * math.exp(1.- zeta_d))**m_gamma
    delta_d = delta_d * phi_d
    
    return delta_d, delta_d_max
    
    
def get_delta_d_max(Z):
    '''
    Eq 2.14
    '''
    
    retval = 0.00047 * Z**1.2 * math.exp( -0.17 * Z**0.47 )    
    retval += 0.04 * (1. - math.exp(-0.002 * Z))
    
    return retval

def get_m_gamma(Z):
    '''
    Eq 2.18
    '''    
    retval = 1. + 0.6 * math.exp(-0.1 * (5. - math.log10(Z))**3.6  )
    
    return retval

def get_nu_star_d(Z):
    '''
    Eq. 2015
    '''
    retval = 35.*(1. - math.exp(-0.074*Z**0.4)) - 5.0
    return retval
 
def get_slope_at_max_steepness(h, D50, rho_particule=_rho_particle):  
    Ycr = sed_trans.get_Ycr(D50,rho_particule)
    gamma_s = sed_trans.get_gamma_s(rho_particule)
    tau_cr = gamma_s*D50*Ycr
    
    Z = h/D50
    nu_star_d = get_nu_star_d(Z)
    
    slope_max_gamma = nu_star_d*tau_cr / (_rho_water*_g*h)
    
    return slope_max_gamma

def get_water_depth_from_length(Lambda=10.0, D50=0.0005, Slope=0.01):
    return Lambda/6.

    '''
    Uses equation 2.11 from Yalin and da Silva
    Z = h/D50
    '''
    h = 0.
    h0 = Lambda / (5.)
    h1 = Lambda / (7.)


    k = 0

    while k <= _maxiter and abs( h1 - h0 ) >= _tol:

        fx0 = _eq_2_11(D50,h0,Slope,Lambda)
        fx1 = _eq_2_11(D50,h1,Slope,Lambda)

        h = h0 - fx0 * ( h0 - h1 ) / ( fx0 - fx1 )
        h1 = h0
        h0 = h
        k = k + 1

    if k > _maxiter:
        print('Error: exceeded {0} iterations'.format(k))

    return h


def get_dune_length(D50,h,slope):
    Z = h/D50
    X =  sed_trans.get_X(h,slope,D50)

    mDelta = 0.055*math.sqrt(Z) + 0.04*X

    retval = 6.0*Z*( 1. + 0.01 * math.exp(-1. * mDelta ) * ((Z-40.)*(Z-400.)) /Z  )

    return retval * D50

def _eq_2_11(D50,h,slope,Lambda):
    Z = h/D50
    X =sed_trans.get_X(h,slope,D50)

    mDelta = 0.055*math.sqrt(Z) + 0.04*X

    retval = 6.0*Z*( 1. + 0.01 * math.exp(-1. * mDelta ) * ((Z-40.)*(Z-400.)) /Z  ) - Lambda/D50

    return retval



if __name__ == "__main__":
    Lambda=42.977
    D50=0.0005
    slope=1.39e-04
    h = get_water_depth_from_length(Lambda, D50, slope)
    
    Z = h/D50

    print('Depth: {0}, Z: {1}, L/D: {2}'.format(h,Z,Lambda/D50) )
    print('Expected Dune length: {0}'.format(get_Dune_Length(D50,h,slope)))


    steep,max_steep = get_steepness(h, slope, D50)
    print('The steepness of dune: {0} and max {1}'.format(steep,max_steep))
    print('The height of dune: {0}'.format(steep*Lambda))

    print('------------- 0.5 mm -------------------')
    slope,h = reconstruct_flow_from_dunes(D50, Lambda, 1.6 )
    print('Expect: {0}, got {1} - height: {2}'.format(0.000139254159320774,slope,h))
    
    print('------------- 1.0 mm -------------------')
    slope,h = reconstruct_flow_from_dunes(0.001, Lambda, 1.6 )
    print('Expect: {0}, got {1} - height: {2}'.format(0.000314646615653066,slope,h))
    
    
    print('------------- 1.5 mm -------------------')
    slope,h = reconstruct_flow_from_dunes(0.0015, Lambda, 1.6 )
    print('Expect: {0}, got {1} - height: {2}'.format(0.000529601593646557,slope,h))
    
    
    print('------------- Final -------------------')
    slope,h = reconstruct_flow_from_dunes(0.0005, 36.4, 1.15 )
    
    
    
    
    