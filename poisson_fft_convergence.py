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
from scipy.fftpack import dst
from scipy.fftpack import idst


plt.close('all')




tol=1e-09 
maxiter=100     

# method of manufacturing solution 
# solution is both periodic and with Dirichlet BC 

    

f, axarr = plt.subplots(1, sharex=False)

err_cg = np.zeros(1000) # red-black Gauss-seidel -2nd order
err_fft = np.zeros(1000) # red-black Gauss-seidel -2nd order

z=np.linspace(8,100)
plt.plot(z,1/(z**2),'k')

for grid_it in np.arange(9,50,2): # odd numbers
 N=grid_it
# k =np.concatenate((np.arange(0,N/2+1),np.arange(N/2,0,-1)))
 k=np.arange(1,N-1)
 u_cg = np.zeros(N)
 u_fft = np.zeros(N)
 
 x=np.linspace(-1,1,N)
 dx = 2.0/(N-1)
 an_sol = np.sin(2.0 * np.pi * x) * x 
 b = -4*np.pi * np.pi * x * np.sin(2*np.pi*x) + 4*np.pi*np.cos(2*np.pi*x)
 b_fft =dst(b[1:N-1],type=1)
 b_fft = b_fft/(-k**2 *np.pi*np.pi)
# b_fft[0] = 0;
 u_fft[1:N-1] = idst(b_fft,type=1)*2/(N-1)
# create tridiagonal matrix
 Poisson = spdiags(np.array(np.ones((1,N-2))),np.array([-1]),N-2,N-2) + spdiags(np.array(np.ones((1,N-2))),np.array([1]),N-2,N-2) + spdiags(np.array(-2* np.ones((1,N-2))),np.array([0]),N-2,N-2)

 temp = scipy.sparse.linalg.cg(Poisson, b[1:N-1]*dx**2,u_cg[1:N-1],tol,maxiter)
 u_cg[1:N-1]=temp[0]
 err_cg[N] = (np.sum(np.abs(u_cg-an_sol)))/N
 err_fft[N]= (np.sum(np.abs(u_fft-an_sol)))/N

 axarr.plot(N,err_cg[N],'k.')
 axarr.plot(N,err_fft[N],'r.')
 
 axarr.set_yscale('log')
 axarr.set_xscale('log')    
    
 plt.pause(0.005)

