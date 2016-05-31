#!/usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt

def get_h_r(l_e, h_m, phi):
    #Kosinussatz
    return sqrt(l_e **2 + h_m**2 - 2*l_e*h_m * cos(phi))
    
def get_alpha(l_e, phi, h_r):
    # Sinussatz
    return arcsin(l_e * sin(phi)/h_r)
    
def get_beta(alpha, psi):
    # Analysis
    return pi - alpha - psi

def get_h(psi, h_b, beta):
    # Sinussatz
    return sin(psi) * h_b / sin(beta)

def get_gamma(lr,hr):
    return arcsin(lr/hr)

def get_third(beta,gamma):
    return pi-beta-gamma

def get_eps(beta):
    return pi - beta

def get_b_bottom(h, beta, gamma):
    # Sinussatz
    delta=get_third(beta,gamma)
    return (sin(gamma)/sin(delta))*h

def get_b_top(h,beta,gamma):
    eps = get_eps(beta)
    zeta = get_third(gamma,eps)
    return (sin(gamma)/sin(zeta))*h
    
def get_b(h,beta,gamma):
    return get_b_top(h,beta,gamma) + get_b_bottom(h,beta,gamma)

def vorwaerts(
    delta_phi = 0,
    x=[1, #l_e
    0.5, #l_r
    2, #h_m
    radians(180), #phi_0
    pi/2, #psi
    4] #h_b
	):
    
    phi = delta_phi + x[3]
    h_r = get_h_r(x[0], x[2], phi)
    alpha = get_alpha(x[0], phi, h_r)
    beta = get_beta(alpha, x[4])
    h = get_h(x[4], x[5], beta)
    gamma = get_gamma(x[1], h_r)
    b_bottom = get_b_bottom(h,beta,gamma)
    b_top = get_b_top(h,beta,gamma)
    b = get_b(h, beta, gamma)
    return b
