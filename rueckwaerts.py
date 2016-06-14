#!/usr/bin/env python

from jacobi import *
from vorwaerts import *

from scipy.optimize import minimize,fmin

from time import *

def newton(F,x0,eps=1e-5,alpha=1.05e-4):
    it=int(1e2)
    x=copy(x0)
    Fx=F(x)

    residuum=linalg.norm(Fx,inf)
    while residuum > eps and it > 0:
        dFx=jacobi(F,x)
        
        #tikhonov regularisation
        Fx=F(x)
        A=dFx
        R=dot(linalg.inv(dot(A.T,A) + alpha*eye(A.shape[1])),A.T)
        xalpha=dot(R,Fx)
        x-=xalpha

        residuum2=linalg.norm(Fx)
	if residuum2-residuum > 0:
            alpha/=2
        if abs(residuum2/residuum) > 0.9:
	    alpha/=2
        elif residuum > residuum2:
	    alpha*=4
#        alpha=min(5e-5,residuum2/residuum)
	residuum=residuum2
        #svd=linalg.svd(A,compute_uv=0)
#	print residuum#,svd
        it-=1

    return x, residuum

def rueckwaerts(n,z,l_r):
    inp=array([[radians(float(i)/n * 360),vorwaerts(radians(float(i)/n * 360),x=z, l_r=l_r)] for i in range(n+1)],dtype=float64)
    F=lambda x: array([abs(vorwaerts(i[0],x,l_r)-i[1]) for i in inp],dtype=float64)
    x0=array([10,30,radians(90),radians(90),90],dtype=float64)

    for method in [ "CG","BFGS","COBYLA",]:
        res=minimize(lambda x:linalg.norm(F(x),inf),x0,method=method)
        print "Minimize(%s): Residuum = %f, x = %s"%(method,linalg.norm(F(res.x),inf),res.x)
    res=fmin(lambda x:linalg.norm(F(x),inf),x0,disp=0)
    print "Fmin: Residuum = %f, x = %s"%(linalg.norm(F(res),inf),res)
    x,residuum=newton(F,x0, alpha=5e-4)
    print x

    return residuum
    
def print_x(x):    
#    print("Param \tValue")
    print("-------------------------------")
    print("l_e: \t{l_e}".format(l_e=x[0]))
    print("h_m: \t{h_m}".format(h_m=x[1]))
    print("phi_0: \trad={phi_0}, deg={deg}".format(phi_0=x[2], deg=degrees(x[2])))
    print("psi: \trad={psi}, deg={deg}".format(psi=x[3], deg=degrees(x[3])))
    print("h_b: \t{h_b}".format(h_b=x[4]))
    print("")

if __name__ == "__main__":
    n=30
    #x=[1,2,radians(90),radians(87),4]
    #l_r=0.5
    l_r = 0.3
    l_e=11
    h_m=27.32
    phi0=radians(92.111)
    psi=radians(90.7)
    h_b=89.8
    x=[l_e,h_m,phi0,psi,h_b]
    print_x(x)
       
    residuum=rueckwaerts(n,x,l_r)
    print "X actually is",x
    print "Residuum is",residuum
    
