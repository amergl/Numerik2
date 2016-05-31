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
    
def get_b2(l_r, l_e, h_m, h_b, phi, psi):
    return (sin(arcsin(l_r/sqrt(l_e **2 + h_m**2 - 2*l_e*h_m * cos(phi))))/sin(pi-arcsin(l_r/sqrt(l_e **2 + h_m**2 - 2*l_e*h_m * cos(phi)))-(pi-(pi - arcsin(l_e * sin(phi)/sqrt(l_e **2 + h_m**2 - 2*l_e*h_m * cos(phi))) - psi))) + sin(arcsin(l_r/sqrt(l_e **2 + h_m**2 - 2*l_e*h_m * cos(phi))))/sin(pi-(pi - arcsin(l_e * sin(phi)/sqrt(l_e **2 + h_m**2 - 2*l_e*h_m * cos(phi))) - psi)-arcsin(l_r/sqrt(l_e **2 + h_m**2 - 2*l_e*h_m * cos(phi))))) * (sin(psi) * h_b / sin((pi - arcsin(l_e * sin(phi)/sqrt(l_e **2 + h_m**2 - 2*l_e*h_m * cos(phi))) - psi)))
    
    
    
def make_plot(h_b, h_m, l_e, l_r, phi, h_r, h, b, b_top, b_bottom, psi, alpha, beta, gamma):
    plt.clf()
    plt.close() 
    plt.axis((-0.5,h_b+0.5,-1.5*b,1.5*b))
    # untere Bodenlinie
    pointpsix, pointpsiy = plot_line(0,0,radians(0),h_b, name="hb", color='k-')
    plt.plot(h_b,0,'kx')
    plt.text(h_b,-0.3,'h_b')
    # Wand
    topx,topy = plot_line(pointpsix,pointpsiy,pi-psi,1.3*b, plot=False)
    plot_line(topx,topy,2*pi-psi,2.6*b, color='k-')
    # Abstand bis Drehteller vom Ursprung
    pointphix,pointphiy = plot_line(0,0,radians(0),h_m, name="hm", color='k-')
    plt.plot(h_m,0,'kx')
    plt.text(h_m,-0.3,'h_m')
    #plot_line(h_m,0,phi,l_e)
    # Winkelhalbierende
    pointbx, pointby = plot_line(0,0,alpha,h, name="h", color='g-')
    #Winkelhalbierende bis Kreismittelpunkt
    pointrx, pointry = plot_line(0,0,alpha,h_r, name="hr", color = 'b-')
    # b an Wand
    btopx,btopy = plot_line(pointbx,pointby,pi-psi,b_top, plot=False)
    plot_line(btopx,btopy,2*pi-psi,b, name="b", color='r-')
    bbottomx,bbottomy = plot_line(pointbx,pointby,2*pi-psi,b_bottom, plot=False)
    # unterer Scheitel
    plt.plot((0,bbottomx),(0,bbottomy),'c-.')
    plt.plot(bbottomx,bbottomy, 'kx')
    # oberer Scheitel
    plt.plot((0,btopx),(0,btopy),'c-.')
    plt.plot(btopx,btopy, 'kx')
    # Drehteller
    circlele=plt.Circle((pointphix,pointphiy),l_e,color='k',linestyle='dashed',fill=False)
    # Zylinder
    circler=plt.Circle((pointrx,pointry),l_r,color='k',fill=False)
    # Abstand Mitte Drehteller, Mitte Zylinder
    plt.plot((h_m, pointrx), (0, pointry), 'y-', label='l_e')
    plt.plot(pointrx,pointry, 'kx')
    plt.text(pointrx,pointry-0.3,'M')
    fig = plt.gcf()
    ax = plt.gca()
    fig.gca().add_artist(circlele)
    fig.gca().add_artist(circler)
    plt.axis('equal')
    plt.legend()
    plt.show()

def plot_line(pointx, pointy, angle, distance, plot=True, name=None, color='-'):
    x = cos(angle)*distance+pointx
    y = sin(angle)*distance+pointy
    if plot:
        plt.plot((pointx,x),(pointy,y),color, label=name)
    return x,y

def main(
    l_e = 1.,
    l_r = 0.5,
    h_m = 2.,
    phi_0 = radians(180),
    delta_phi = radians(0),
    psi = pi/2,
    h_b = 4.):
    
    phi = delta_phi + phi_0
    h_r = get_h_r(l_e, h_m, phi)
    alpha = get_alpha(l_e, phi, h_r)
    beta = get_beta(alpha, psi)
    h = get_h(psi, h_b, beta)
    gamma = get_gamma(l_r, h_r)
    b_bottom = get_b_bottom(h,beta,gamma)
    b_top = get_b_top(h,beta,gamma)
    b = get_b(h, beta, gamma)
    print "b:  ", b
    b2 = get_b2(l_r, l_e, h_m, h_b, phi, psi)
    print "b2: ", b2
    make_plot(h_b, h_m, l_e, l_r, phi, h_r, h, b, b_top, b_bottom, psi, alpha, beta, gamma)
    print abs(b2-b)
    print "-----------------"
    
    

if __name__ == "__main__":
    main(l_e=2., l_r=1., h_m=8., phi_0=radians(135), delta_phi=radians(0), psi=radians(75), h_b=20.)
    main(l_e=2., l_r=1., h_m=8., phi_0=radians(135), delta_phi=radians(45), psi=radians(75), h_b=20.) 
    main(l_e=2., l_r=1., h_m=8., phi_0=radians(135), delta_phi=radians(90), psi=radians(75), h_b=20.)
    main(l_e=2., l_r=1., h_m=8., phi_0=radians(135), delta_phi=radians(135), psi=radians(75), h_b=20.)
    main(l_e=2., l_r=1., h_m=8., phi_0=radians(135), delta_phi=radians(180), psi=radians(75), h_b=20.)
    main(l_e=2., l_r=1., h_m=8., phi_0=radians(135), delta_phi=radians(225), psi=radians(75), h_b=20.)
    main(l_e=2., l_r=1., h_m=8., phi_0=radians(135), delta_phi=radians(270), psi=radians(75), h_b=20.)
    main(l_e=2., l_r=1., h_m=8., phi_0=radians(135), delta_phi=radians(315), psi=radians(75), h_b=20.)
