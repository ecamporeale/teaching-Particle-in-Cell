# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:13:50 2016

@author: camporea
"""

# example of Poisson solvers in 1D
# solve laplacian(f)= b

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')


N = 17 # number of grid points -- nodes based

# method of manufacturing solution 
# solution is both periodic and with Dirichlet BC 


# iterative red-black gauss-seidel

#Convergence to 2nd order is achieved only for large enough number of iterations
n_it = 5000 # number of iterations. 

f, axarr = plt.subplots(1, sharex=False)

err_rb_gs = np.zeros(1000) # red-black Gauss-seidel -2nd order

z=np.linspace(8,100)
plt.plot(z,1/(z**2),'k')

for grid_it in np.arange(8,100):
 N=grid_it
 u_rb_gs = np.zeros(N) # red-black Gauss-Seidel

 x=np.linspace(-1,1,N)
 dx = 2.0/(N-1)
 an_sol = np.sin(2.0 * np.pi * x) * x 
 b = -4*np.pi * np.pi * x * np.sin(2*np.pi*x) + 4*np.pi*np.cos(2*np.pi*x)

 for it in np.arange(0,n_it):

    for i in np.arange(1,N-1,2):    
      u_rb_gs[i] = (-dx**2*b[i] + u_rb_gs[i+1] + u_rb_gs[i-1])/2.0
 
    for i in np.arange(2,N-1,2):    
      u_rb_gs[i] = (-dx**2*b[i] + u_rb_gs[i+1] + u_rb_gs[i-1])/2.0

 
 err_rb_gs[N]= (np.sum(np.abs(u_rb_gs-an_sol)))/N


 axarr.plot(N,err_rb_gs[N],'k.')
 
 axarr.set_yscale('log')
 axarr.set_xscale('log')    
    
 plt.pause(0.005)

