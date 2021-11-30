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


N = 33 # number of grid points -- nodes based

# method of manufacturing solution 
# solution is both periodic and with Dirichlet BC 

x=np.linspace(-1,1,N)
dx = 2.0/(N-1)
an_sol = np.sin(2.0 * np.pi * x) * x 
b = -4*np.pi * np.pi * x * np.sin(2*np.pi*x) + 4*np.pi*np.cos(2*np.pi*x)

# iterative jacobi method

n_it = 200 # number of iterations
u = np.zeros(N) # initial guess
u_new = np.zeros(N) 
u_new_w = u_new
u_jacobi = np.zeros(N) 
u_gs = np.zeros(N) # Gauss-Seidel
u_rb_gs = np.zeros(N) # red-black Gauss-Seidel
u_sor = np.zeros(N) # successive over-relaxation of rd_gs
omega=1.3

f, axarr = plt.subplots(3, sharex=False)

err_jacobi = np.zeros(n_it)
err_gs = np.zeros(n_it) # Gauss-seidel
err_rb_gs = np.zeros(n_it) # red-black Gauss-seidel
err_sor = np.zeros(n_it) # SOR

for it in np.arange(0,n_it):
    for i in np.arange(1,N-1):    
      u_new[i] = (-dx**2*b[i] + u_jacobi[i+1] + u_jacobi[i-1])/2.0
      u_gs[i] = (-dx**2*b[i] + u_gs[i+1] + u_gs[i-1])/2.0

    for i in np.arange(1,N-1,2):    
      u_rb_gs[i] = (-dx**2*b[i] + u_rb_gs[i+1] + u_rb_gs[i-1])/2.0

    for i in np.arange(2,N-1,2):    
      u_rb_gs[i] = (-dx**2*b[i] + u_rb_gs[i+1] + u_rb_gs[i-1])/2.0

    for i in np.arange(1,N-1):    
     u_sor[i] = omega * u_rb_gs[i] +(1-omega)* u_sor[i]
     
    for i in np.arange(1,N-1):    
      u_jacobi[i] = u_new[i]  

    axarr[0].cla()
    axarr[1].cla()
      
    axarr[0].plot(x,u_jacobi)
    axarr[0].plot(x,u_gs)
    axarr[0].plot(x,u_rb_gs)
    axarr[0].plot(x,u_sor)
    
    axarr[0].plot(x,an_sol)

    axarr[0].legend(['Jacobi','Gauss-seidel','Red-black GS','SOR','Solution'])
    axarr[0].set_title(it)
 
    axarr[1].plot(x[1:N-1],np.abs(u_jacobi[1:N-1]-an_sol[1:N-1]))
    axarr[1].plot(x[1:N-1],np.abs(u_gs[1:N-1]-an_sol[1:N-1]))
    axarr[1].plot(x[1:N-1],np.abs(u_rb_gs[1:N-1]-an_sol[1:N-1]))
    axarr[1].plot(x[1:N-1],np.abs(u_sor[1:N-1]-an_sol[1:N-1]))
    axarr[1].set_ylim([1e-4,0.5])
    axarr[1].set_xlim([-1+dx, 1-dx])

    axarr[1].set_yscale('log')

    axarr[1].legend(['Jacobi','Gauss-seidel','Red-black GS','SOR'])


    err_jacobi[it]= np.linalg.norm(u_jacobi-an_sol,ord=2)
    err_gs[it]= np.linalg.norm(u_gs-an_sol,ord=2)
    err_rb_gs[it]= np.linalg.norm(u_rb_gs-an_sol,ord=2)
    err_sor[it]= np.linalg.norm(u_sor-an_sol,ord=2)
    
    axarr[2].plot(it,err_jacobi[it],'k.')
    axarr[2].plot(it,err_gs[it],'r.')
    axarr[2].plot(it,err_rb_gs[it],'g.')
    axarr[2].plot(it,err_sor[it],'b.')
    
    axarr[2].set_yscale('log')
    axarr[2].legend(['Jacobi','Gauss-seidel','Red-black GS','SOR'])

    plt.pause(0.005)

