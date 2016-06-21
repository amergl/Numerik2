#!/usr/bin/env python

from jacobi import jacobi
from vorwaerts import vorwaerts
from scipy.optimize import minimize, fmin
from Vorwaertsproblem import Vorwaertsproblem
import numpy as np
from time import time
from sys import argv


def newton(F,x0,eps=1e-5,alpha=1.05e-4):
    it=int(1e2)
    x=np.copy(x0)
    Fx=F(x)

    breakit=int(argv[2]) if len(argv) > 2  else 5
    residuum2=0
    residuum=np.linalg.norm(Fx,np.inf)
    diff=abs(residuum-residuum2)
    while residuum > eps and it > 0 and diff > eps and breakit > 0:
        dFx = jacobi(F, x)
        
        # Tikhonov regularisation
        A = dFx
        R = np.dot(np.linalg.inv(np.dot(A.T, A) + alpha * np.eye(A.shape[1])), A.T)
        xalpha = np.dot(R, Fx)
        x -= xalpha
        Fx = F(x)
        
        residuum2=np.linalg.norm(Fx)
        diff=abs(residuum-residuum2)
	if residuum2 > residuum:
            x+=xalpha
            Fx=F(x)
            alpha/=2
            alpha=min(1,abs(alpha))
            breakit-=1
	else:
            residuum=residuum2
            breakit=int(argv[2]) if len(argv) > 2  else 5
        it-=1

    return x, residuum

def rueckwaerts(n, z, l_r):
    inp = np.array([[np.radians(float(i)/n * 360), vorwaerts(np.radians(float(i)/n * 360), x=z, l_r=l_r)] for i in range(n+1)], dtype=np.float64)
    F = lambda x: np.array([abs(vorwaerts(i[0], x, l_r)-i[1]) for i in inp], dtype=np.float64)
    x_0 = np.array([10, 30, np.radians(90), np.radians(90), 90], dtype=np.float64)
    
    #if compare:
        #print("Method \t\tResiduum \t\t\tx = [l_e, h_m, phi_0, psi, h_b] \t\t\t\t\tDuration")
        #print("".ljust(155, "-"))

    for method in [ "CG", "BFGS", "COBYLA",]:
        t0=time()
        res = minimize(lambda x: np.linalg.norm(F(x), np.inf), x_0, method=method)
        duration=time()-t0
        if compare:
            print("%-20s \tResiduum = %e, \tx = %s, \tduration=%4fs"%("Minimize (%s):"%method, np.linalg.norm(F(res.x), np.inf),res.x,duration))
    t0=time()
    res = fmin(lambda x: np.linalg.norm(F(x), np.inf), x_0, disp=0)
    duration=time()-t0
    if compare:
        print("%-20s \tResiduum = %e, \tx = %s, \tduration=%4fs"%("Fmin:",np.linalg.norm(F(res), np.inf), res,duration))
    t0=time()
    x, residuum = newton(F, x_0, alpha=5e-4)
    duration=time()-t0
    if compare:
        print "%-20s \tResiduum = %e, \tx = %s, \tduration=%4fs"%("NumII-Team:",residuum,x,duration)
    return x, residuum
    
def print_x(title, x):    
    print(title)
    print("-------------------------------")
    print("l_e: \t{l_e}".format(l_e=x[0]))
    print("h_m: \t{h_m}".format(h_m=x[1]))
    print("phi_0: \trad={phi_0}, deg={deg}".format(phi_0=x[2], deg=np.degrees(x[2])))
    print("psi: \trad={psi}, deg={deg}".format(psi=x[3], deg=np.degrees(x[3])))
    print("h_b: \t{h_b}".format(h_b=x[4]))
    print("-------------------------------\n")

if __name__ == "__main__":
    compare=False
    if len(argv) > 1 and argv[1] == "compare":
        compare=True
    n = 30
    if not compare:
        print("Measurements: \t{n}".format(n=n))    
    delta_phi = 360/n
    if not compare:
        print("Delta Phi: \t\t{delta} degrees\n".format(delta=delta_phi))
    l_r = 0.3
    l_e = 11
    h_m = 27.32    
    phi_0 = 92.111
    psi = 90.7
    h_b = 89.8
    x_exact = [l_e, h_m , np.radians(phi_0), np.radians(psi), h_b]
    
    tstart = time()
    x_calc, residuum = rueckwaerts(n, x_exact, l_r)
    duration = time() - tstart
    if not compare:
        print_x("\nExact x:", x_exact)
        print_x("Calculated x:", x_calc)
        print("Residuum: \t{residuum}".format(residuum=residuum))
        print("Duration: \t{duration} seconds".format(duration=duration))
    
    #vorwaertsproblem = Vorwaertsproblem(l_e, l_r, h_m, h_b, phi_0, delta_phi, psi, show_legend=False)
    #vorwaertsproblem.plot()
