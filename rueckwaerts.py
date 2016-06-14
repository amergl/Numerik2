#!/usr/bin/env python

from jacobi import jacobi
from vorwaerts import vorwaerts
from scipy.optimize import minimize, fmin
from Vorwaertsproblem import Vorwaertsproblem
import numpy as np
import time


def newton(F, x0, eps=1e-5, alpha=1.05e-4):
    it = int(1e2)
    x = np.copy(x0)
    Fx = F(x)

    residuum = np.linalg.norm(Fx, np.inf)
    while residuum > eps and it > 0:
        dFx = jacobi(F, x)
        
        # Tikhonov regularisation
        Fx = F(x)
        A = dFx
        R = np.dot(np.linalg.inv(np.dot(A.T, A) + alpha * np.eye(A.shape[1])), A.T)
        xalpha = np.dot(R, Fx)
        x -= xalpha

        residuum2 = np.linalg.norm(Fx)
	if residuum2 - residuum > 0:
            alpha /= 2
        if abs(residuum2/residuum) > 0.9:
	    alpha /= 2
        elif residuum > residuum2:
	    alpha *= 4
#         alpha = min(5e-5,residuum2/residuum)
	residuum = residuum2
        # svd = linalg.svd(A,compute_uv=0)
#	print residuum #,svd
        it -= 1

    return x, residuum

def rueckwaerts(n, z, l_r):
    inp = np.array([[np.radians(float(i)/n * 360), vorwaerts(np.radians(float(i)/n * 360), x=z, l_r=l_r)] for i in range(n+1)], dtype=np.float64)
    F = lambda x: np.array([abs(vorwaerts(i[0], x, l_r)-i[1]) for i in inp], dtype=np.float64)
    x_0 = np.array([10, 30, np.radians(90), np.radians(90), 90], dtype=np.float64)

    for method in [ "CG", "BFGS", "COBYLA",]:
        res = minimize(lambda x: np.linalg.norm(F(x), np.inf), x_0, method=method)
        print("Minimize ({method}): \tResiduum = {residuum}, \tx = {x}".format(method=method, residuum=np.linalg.norm(F(res.x), np.inf), x=res.x))
    res = fmin(lambda x: np.linalg.norm(F(x), np.inf), x_0, disp=0)
    print("Fmin: \t\tResiduum = {residuum}, \tx = {x}".format(residuum=np.linalg.norm(F(res), np.inf), x=res))
    x, residuum = newton(F, x_0, alpha=5e-4)
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
    n = 30
    print("Measurements: \t{n}".format(n=n))    
    delta_phi = 360/n
    print("Delta Phi: \t\t{delta} degrees\n".format(delta=delta_phi))
    l_r = 0.3
    l_e = 11
    h_m = 27.32    
    phi_0 = 92.111
    psi = 90.7
    h_b = 89.8
    x_exact = [l_e, h_m , np.radians(phi_0), np.radians(psi), h_b]
    
    tstart = time.time()
    x_calc, residuum = rueckwaerts(n, x_exact, l_r)
    duration = time.time() - tstart
    print_x("\nExact x:", x_exact)
    print_x("Calculated x:", x_calc)
    print("Residuum: \t{residuum}".format(residuum=residuum))
    print("Duration: \t{duration} seconds".format(duration=duration))
    
    vorwaertsproblem = Vorwaertsproblem(l_e, l_r, h_m, h_b, phi_0, delta_phi, psi, show_legend=False)
    vorwaertsproblem.plot()
