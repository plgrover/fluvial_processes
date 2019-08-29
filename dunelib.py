
import math

_nu = 1.0e-6
_g  = 9.81
_gamma_s = 16186.5
_kappa = 0.4
_rho_water = 1000.
_rho_particle = 2650.0


def get_total_chezy(h, D50, Lambda, Delta):
    cf = get_cf(h, D50, 8.5)

    c = ((1. / cf ** 2) + 0.5 * ((Delta / Lambda) ** 2) * (Lambda / h)) ** (-0.5)

    return c


def get_c_delta(h, Lambda, Delta):
    retval = 0.5 * ((Delta / Lambda) ** 2) * (Lambda / h)

    return retval


def get_cf(h, D50, Bs=8.5):
    cf = (1. / _kappa) * math.log(0.368 * h / (D50 * 2.)) + Bs

    return cf