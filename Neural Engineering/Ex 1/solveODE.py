# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 10:24:06 2022

@author: pwerginz
"""

### ------------------------------------------------
# Implementation of different ways of solving the ODE y' = -20*y
# The analytical solution of the ODE is y(t) = exp(-20*t)
# The computation is done analytically or by a fixed step size forward or backward Euler method


# Import necessary packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


### Main function to solve model
def solveODE():
    ### ------------- Step 1: PARAMETERS -------------
    y0 = 1 # Initial condition
    tStop = 1 # Total duration of simulation
    tDt = 0.025 # Time step
    
    
    ### ------------- Step 2: COMPUTE ADDITIONAL PARAMETERS -------------
    # Compute time step
    timeSteps = int(tStop/tDt)+1
    timeStep = np.linspace(0,tStop,int(tStop/tDt)+1)
    
    
    ### ------------- Step 3: ANALYTICAL -------------
    tDtAn = 0.001 # Fine time step for reference solution
    timeStepAn = np.linspace(0,tStop,int(tStop/tDtAn)+1)
    yAnVec = np.exp(-20*timeStepAn)
    
    
    ### ------------- Step 4: FORWARD EULER ----------
    yFEVec = np.zeros(timeSteps) # Allocate memory
    yFEVec[0] = y0 # Set initial condition
    # Loop over time
    for t in range(timeSteps-1):
        yFEVec[t+1] = yFEVec[t] + (-20*yFEVec[t])*tDt
    
    
    ### ------------- Step 5: BACKWARD EULER ---------
    yBEVec = np.zeros(timeSteps) # Allocate memory
    yBEVec[0] = y0 # Set initial condition
    
    # Implicit BE to explicit form
    # y1=y0+(-20*y1)*h
    # y1+20*y1*h=y0
    # y1*(1+20h)=y0
    # y1=y0/(1+20h)
    
    # Loop over time
    for t in range(timeSteps-1):
        yBEVec[t+1] = yBEVec[t]/(1+20*tDt)
        
    def odefun(t, y): return -20 * y
    sol = integrate.solve_ivp(odefun, [0, 1], [y0], method='RK45')
    
    ### --------- Step 6: VISUALIZE RESULTS ---------
    # y over time
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.plot(timeStepAn,yAnVec,label='Analytical')
    ax.plot(timeStep,yFEVec,label='FE')
    ax.plot(timeStep,yBEVec,label='BE')
    ax.plot(sol.t,sol.y[0],label='RK45')
    plt.legend()
    plt.xlim(0,tStop)
    plt.ylim(min(-0.1,min(yFEVec)),max(1.1,max(yFEVec)))
    plt.xlabel('Time')
    plt.ylabel('Y')
    plt.show()

#
if __name__ == '__main__':
    solveODE()
    