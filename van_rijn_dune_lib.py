#!/usr/bin/env python

import numpy as np
import math
import sed_trans
import dunelib
from scipy.optimize import fsolve

_nu = 1.0e-6
_g  = 9.81
_gamma_s = 16186.5
_kappa = 0.4

_rho_water = 1000.
_rho_particle = 2650.0

# This is T where the max steepness occurs
_T_max = 4.81247

_maxiter = 100
_tol = 1.e-4

def reconstruct_flow_from_dunes(D50, Lambda, Delta, rho_particule=_rho_particle):
    ubar_low = None
    S_low = None
    ubar_high = None
    S_high = None

    # Check if the current steepness is greater than the maximum steepness
    delta_d = Delta / Lambda
    h = get_water_depth_from_length(Lambda)
    cf = dunelib.get_cf(h, D50)
    c =  dunelib.get_total_chezy(h,D50, Lambda, Delta)
    uprime_star_cr = get_uprime_star_cr(D50, rho_particule, cf)
    Z = h / D50

    delta_d_max = get_delta_d_max(h, D50)

    print('The flow depth: {0} m and Z: {1}'.format(h, Z))
    print('The cf = {0} and c {1} '.format(cf,c))
    print('The maximum dune steepness: {0} and the observed steepness: {1}'.format(delta_d_max, delta_d))

    if delta_d < delta_d_max:
        # Solve for T
        Tlow, Thigh = solve_for_T(delta_d, Z)

        # Solve for ubars
        ubar_low = get_uprime_star(Tlow, uprime_star_cr) * cf
        S_low = ((ubar_low/c)**2)/(_g*h)
        ubar_high = get_uprime_star(Thigh, uprime_star_cr) * cf
        S_high = ((ubar_high / c) ** 2) / (_g * h)


    else:
        print('The steepness of the dune is greater than the maximum steepness')
        print('Calculating the flow velocity and slope using the value derived from Tmax')
        ubar_low = get_uprime_star(_T_max, uprime_star_cr) * cf
        S_low = ((ubar_low / c) ** 2) / (_g * h)

    return h, ubar_low, S_low, ubar_high, S_high

def solve_for_T(delta_d_target, Z):

    func = lambda T: delta_d_target - (0.015 * (1./Z)**0.3)*(1. - math.exp(-0.5 * T)) * (25. - T)

    T_initial_guess_low = 0.5
    T_initial_guess_high = 11.0
    T_low = fsolve(func, T_initial_guess_low)
    T_high = fsolve(func, T_initial_guess_high)

    return T_low, T_high


def get_uprime_star(T, uprime_star_cr):
    return math.sqrt(T*(uprime_star_cr**2.) + (uprime_star_cr**2.))


def get_uprime_star_cr(D50,rho_particule, cf):
    gamma_s = sed_trans.get_gamma_s(rho_particule)
    u_star_cr = sed_trans.get_Ycr(D50, rho_particule) * gamma_s * D50

    return  u_star_cr/cf

def get_water_depth_from_length(Lambda):
    return Lambda / 7.3



def get_delta_d_max(h,D50):

    Z = h / D50
    # The maximum steepness occurs at this value of T
    # Performed in Wolfram Alpha
    # d/dx((1-exp(-0.5 x)) (25-x)) = e^(-0.5 x) (13.5-0.5 x)-1

    return get_delta_d(_T_max, Z)


def get_delta_d(T, Z):

    delta_d = (0.015 * (1./Z)**0.3)*(1. - math.exp(-0.5 * T)) * (25. - T)

    return delta_d




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
    
    