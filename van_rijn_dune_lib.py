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

_tol = 1.e-4
_maxiter = 100


def reconstruct_flow_from_dunes(D50, Lambda, Delta, rho_particule=_rho_particle):
    
    # Check if the current steepness is greater than the maximum steepness
    gamma_d = Delta/Lambda
    h = get_water_depth_from_length(Lambda)
    gamma_d_max = get_gamma_d_max(h,D50)    
    
    slope = get_slope_at_max_steepness(h, D50, rho_particule)
    
    print('gamma_d: {0}   gamma_d_max: {1}'.format(gamma_d,gamma_d_max))
    
    print('Slope at max steepness: {0}'.format(slope))
    
    
    
    if gamma_d < gamma_d_max:   
                
        _tol = 1.e-8
        
        slope1 = slope * 0.9
        slope0 = slope * 2.05
        
        k = 0
        while k <= _maxiter and abs( (slope1 - slope0) ) >= _tol:
            
            print(slope1,slope0) 
            h0 = get_water_depth_from_length(Lambda)
            h1 = get_water_depth_from_length(Lambda)
    
            fx0 = get_steepness(h0, slope0, D50, rho_particule) - gamma_d
            fx1 = get_steepness(h1, slope1, D50, rho_particule) - gamma_d

    
            slope = slope0 - fx0 * ( slope0 - slope1 ) / ( fx0 - fx1 )
            
            slope1 = slope0
            slope0 = slope
    
            k = k + 1
    
        if k > _maxiter:
            print('Error: exceeded {0} iterations'.format(k))
    else:
        print('The steepness is greater than then maximum dune steepness')
        
    h = get_water_depth_from_length(Lambda)
    return slope,h


def get_slope_at_max_steepness(h, D50, rho_particule=_rho_particle):  
    Ycr = sed_trans.get_Ycr(D50,rho_particule)
    gamma_s = sed_trans.get_gamma_s(rho_particule)
    u_cr_sqr =  Ycr*gamma_s*D50/_rho_water 
    
    T = 4.81247       
    ustar_max_gamma_sqr =  T * u_cr_sqr + u_cr_sqr 
    slope_max_gamma = ustar_max_gamma_sqr / (_g*h)
    
    return slope_max_gamma


def get_gamma_d_max(h,D50):
    
    Z = h / D50
    # The maximum steepness occurs at this value of T
    # Performed in Wolfram Alpha
    # d/dx((1-exp(-0.5 x)) (25-x)) = e^(-0.5 x) (13.5-0.5 x)-1    
    T = 4.81247
    gamma_d = 0.015 * (1./Z)**0.5 
    gamma_d *= (1. - math.exp(-0.5 * T)) * (25. - T)
    
    return gamma_d
    

def get_water_depth_from_length(Lambda):
    return Lambda / 7.3

def get_steepness(h, slope, D50, rho_particule):
    
    T = get_T(h, slope, D50, rho_particule)
    Z = h / D50
    
    gamma_d = 0.015 * (1./Z)**0.5 
    gamma_d *= (1. - math.exp(-0.5 * T)) * (25. - T)
    
    return gamma_d

def get_T(h, slope, D50, rho_particule=_rho_particle):
    u_star = math.sqrt(_g*slope*h)
    gamma_s= get_gamma_s(rho_particule)
    u_star_cr = sed_trans.get_Ycr(D50,rho_particule) * gamma_s * D50
    
    T = ( (u_star**2) - (u_star_cr**2) )/(u_star_cr**2)
    
    return T


if __name__ == "__main__":
    Lambda=42.977
    D50=0.0005
    slope=1.39e-04


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
    
    