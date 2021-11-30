# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:13:50 2016

@author: camporea
"""

# example of Poisson solvers in 1D
# solve laplacian(f)= b
# conjugate-gradient method

import numpy as np
from scipy.sparse import spdiags
import scipy.sparse
from scipy.sparse import linalg
import matplotlib.pyplot as plt

plt.close('all')




tol=1e-09 

# method of manufacturing solution 
# solution is both periodic and with Dirichlet BC 

    
f, axarr = plt.subplots(1, sharex=False)

err_cg = np.zeros(200) # red-black Gauss-seidel -2nd order
err_rb_gs = np.zeros(200) # red-black Gauss-seidel -2nd order


N=250
u_cg = np.zeros(N)
u_rb_gs = np.zeros(N)

x=np.linspace(-1,1,N)
dx = 2.0/(N-1)
an_sol = np.sin(2.0 * np.pi * x) * x 
b = -4*np.pi * np.pi * x * np.sin(2*np.pi*x) + 4*np.pi*np.cos(2*np.pi*x)
# create tridiagonal matrix
Poisson = spdiags(np.array(np.ones((1,N-2))),np.array([-1]),N-2,N-2) + spdiags(np.array(np.ones((1,N-2))),np.array([1]),N-2,N-2) + spdiags(np.array(-2* np.ones((1,N-2))),np.array([0]),N-2,N-2)

for maxiter in np.arange(8,200):

 temp = scipy.sparse.linalg.cg(Poisson, b[1:N-1]*dx**2,0*u_cg[1:N-1],tol,maxiter)
 u_cg[1:N-1]=temp[0]
 
 for i in np.arange(1,N-1,2):    
      u_rb_gs[i] = (-dx**2*b[i] + u_rb_gs[i+1] + u_rb_gs[i-1])/2.0

 for i in np.arange(2,N-1,2):    
      u_rb_gs[i] = (-dx**2*b[i] + u_rb_gs[i+1] + u_rb_gs[i-1])/2.0

 err_cg[maxiter]= (np.sum(np.abs(u_cg-an_sol)))/N
 err_rb_gs[maxiter]= (np.sum(np.abs(u_rb_gs-an_sol)))/N


 axarr.plot(maxiter,err_cg[maxiter],'k.')
 axarr.plot(maxiter,err_rb_gs[maxiter],'r.')
 axarr.legend(['CG','Red-black Gauss-Seidel'])
 
 axarr.set_yscale('log')
    
 plt.pause(0.001)

