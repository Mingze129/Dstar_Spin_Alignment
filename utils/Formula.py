import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def BetheBlochAlephNP(p, m, params, masses):
    bg = p/m
    beta = bg/np.sqrt(1.+ bg*bg)
    aa   = beta**params[3]
    bb   = bg**(-params[4])
    bb   = np.log(params[2]+bb)
    charge_factor = 1. # params[6]         # params[5] = mChargeFactor, params[6] = mMIP
    final = (params[1]-aa-bb)*params[0]*charge_factor/aa
    return final

def gausscale(x, mean, sigma, scale):
    return (scale/0.4)*np.exp(-0.5*((x-mean)/sigma)**2)

def linear(x, a, b):
    return a*x + b

def gausgauslin(x, mu1, sigma1, scale1, mu2, sigma2, scale2, a, b):
    return gausscale(x, mu1, sigma1, scale1) + gausscale(x, mu2, sigma2, scale2) + (a*x + b)

def gauslin(x, mu1, sigma1, scale1, a, b):
    return gausscale(x, mu1, sigma1, scale1) + (a*x + b)

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def constant_func(x, pars):
    """
    """
    return pars[0]
