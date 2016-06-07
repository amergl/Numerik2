#!/usr/bin/env python

from numpy import *

from jacobi import *
from vorwaerts import *

from time import *

def newton(F,x0,eps=1e-5,alpha=1.05e-4):
    it=int(1e5)
    x=copy(x0)
    Fx=F(x)

    residuum=linalg.norm(Fx,inf)
    while residuum > eps and it > 0:
        dFx=jacobi(F,x)

        #tikhonov regularisation

        A=dFx
        R=dot(linalg.inv(dot(A.T,A) + alpha*eye(A.shape[1])),A.T)
        x-=dot(R,Fx)

        Fx=F(x)
        it-=1
        residuum=linalg.norm(Fx,inf)
    if it == 0:
        residuum=-1
    return x, residuum

def rueckwaerts(inp,exact):
    F=lambda x: array([vorwaerts_explicit(i[0],x)-i[1] for i in inp],dtype=float64)
    x0=array([1,0.5,2.,radians(88),radians(88),2.],dtype=float64)
    x,residuum=newton(F,x0)
    print x
    return residuum

if __name__ == "__main__":
    n=10
    x=[1,0.5,2,radians(90),radians(87),4]
    inp=array([[radians(float(i)/n * 360),vorwaerts(radians(float(i)/n * 360),x)] for i in range(n+1)],dtype=float64)
    residuum=rueckwaerts(inp,x)
    print "X actually is",x
    print "Residuum is",residuum
    
