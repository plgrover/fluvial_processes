#!/usr/bin/env python

import numpy as np
import math

_nu = 1.0e-6
_g  = 9.81
_gamma_s = 16186.5
_kappa = 0.41
_Bs = 8.5

_rho_water = 1000.
_rho_particle = 2650.0


def get_X(h,slope,D50):
    X =  math.sqrt(_g*h*slope)*D50/_nu
    return X 

def get_gamma_s(rho_particule=_rho_particle):
    ''' 
    '''
    return _g*(rho_particule-_rho_water)    

    
def get_Y(h, slope, D50, rho_particule=_rho_particle):
    '''
    '''
    tau =  _rho_water * _g * slope *  h
    gamma_s = get_gamma_s(rho_particule)
    retval = tau/(gamma_s * D50)
    
    return retval
    
def get_Ycr(D50,rho_particule):
    '''
    Equation 1.34
    '''
    xi = get_Xi(D50,rho_particule=2650.)
    ycr = 0.13*(xi**(-0.392))*math.exp(-0.015*xi**2.) + 0.045*(1-math.exp(-0.068*xi))
    return ycr

def get_Xi(D50,rho_particule):
    '''
    Equation 1.31
    '''
    gammaS = get_gamma_s(rho_particule)
    xi = ((gammaS * D50**3.)/(_rho_water * _nu**2.))**(1./3.)
    return xi



def get_nu_star(slope, h, D50, rho_particule=_rho_particle):
    '''
    '''
    Y = get_Y(slope, h, D50, rho_particule)
    Ycr = get_Ycr(D50, rho_particule)
    nu_star = Y / Ycr
    
    return nu_star
