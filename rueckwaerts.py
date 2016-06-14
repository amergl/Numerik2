#!/usr/bin/env python

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
        Fx=F(x)
        A=dFx
        R=dot(linalg.inv(dot(A.T,A) + alpha*eye(A.shape[1])),A.T)
        xalpha=dot(R,Fx)
        x-=xalpha
        
        residuum2=linalg.norm(Fx-dot(dFx,xalpha))
#	if residuum2-residuum > 0:
#            alpha/=10
#        if abs(residuum2/residuum) > 0.9:
#	    alpha/=5
#        elif residuum/residuum2 < 1e2:
#	    alpha*=5
#        alpha=min(5e-5,residuum2/residuum)
	residuum=residuum2
        #svd=linalg.svd(A,compute_uv=0)
	print residuum#,svd
        it-=1

    if it == 0:
        residuum=-1
    return x, residuum

def rueckwaerts(n,x,l_r):
    inp=array([[radians(float(i)/n * 360),vorwaerts(radians(float(i)/n * 360),x=x, l_r=l_r)] for i in range(n+1)],dtype=float64)
    F=lambda x: array([vorwaerts(i[0],x,l_r)-i[1] for i in inp],dtype=float64)
    x0=array([30,30,radians(90),radians(90),90],dtype=float64)
    #x0=array([1.7,6.,radians(120),radians(88),15],dtype=float64)
    
    x,residuum=newton(F,x0, alpha=5e-4)
    print x
    return residuum

if __name__ == "__main__":
    n=30
    #x=[1,2,radians(90),radians(87),4]
    #l_r=0.5
    l_r = 0.3
    x=[29,29,radians(90.5),radians(90.5),91]
    residuum=rueckwaerts(n,x,l_r)
    print "X actually is",x
    print "Residuum is",residuum
    
