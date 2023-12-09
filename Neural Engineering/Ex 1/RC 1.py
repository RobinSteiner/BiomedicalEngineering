# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 13:43:11 2022

@author: paulwerginz
"""

# ------------------------------------------------
# RC-circuit (passive) single-compartment (local) model
# The computation is done by a fixed step size forward or backward Euler method
# One stimulation mode: Current injection into the compartment ('iClamp')
# ------------------------------------------------
#
#       _||_
#       \  /
#   ---- \/ ----
#   |          |
#   ------------
#
#
# Import necessary packages
import numpy as np
import matplotlib.pyplot as plt


### Main function to solve model
def RC():
    ### ------------- Step 1: PARAMETERS -------------
    # Solver
    solver = 'BE' # Type of solver, 'FE'==Forward Euler or 'BE'==Backward Euler
    
    # Stimulus parameters
    I = 1 # Stimulus amplitude, in nA
    
    # Temporal parameters
    tStop = 200 # Total duration of simulation, in ms
    tDel = 50 # Delay until stimulus starts, in ms
    tDur = 100 # Duration of stimulus ON, in ms
    tDt = 0.025 # Time step, in ms
    
    # Membrane parameters
    R = 30 # Resistance, in MOhm
    C = 0.5 # Capacitance, in nF


    ### ------------- Step 2: COMPUTE ADDITIONAL PARAMETERS -------------
    # Compute time step
    timeSteps = int(tStop/tDt)+1
    timeStep = np.linspace(0,tStop,timeSteps)

    # Allocate memory for v
    vVec = np.zeros(timeSteps)
    
    
    ### --------- Step 3: SOLVE ODE & POSTPROCESSING ---------  
    for t in range(0,timeSteps-1):
        # Stimulus current
        IStim = 0
        if t>=int(tDel/tDt) and t<int((tDel+tDur)/tDt):
            IStim = I # in nA
        
        # Compute change of v
        if solver=='FE':
            vVec[t+1] = vVec[t]+((-vVec[t]/R+IStim)/C)*tDt # in (nA/nF)*ms==mV
        elif solver=='BE':
            vVec[t+1] = (vVec[t]+IStim*(tDt/C))/(1+tDt/(R*C)) # in ((nA*ms)/nF)/(ms/(MOhm*nF))==mV
        
    
    ### --------- Step 4: VISUALIZE RESULTS ---------
    # Membrane voltage over time
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.plot(timeStep,vVec,'b',label='Response')
    plt.xlim(0,tStop)
    plt.ylim(min(vVec)-55,max(vVec)+55)
    ax.plot([0,tDel],[min(vVec)-50,min(vVec)-50],'r',label='Input')
    ax.plot([tDel,tDel],[min(vVec)-50,min(vVec)-30],'r')
    ax.plot([tDel,tDel+tDur],[min(vVec)-30,min(vVec)-30],'r')
    ax.plot([tDel+tDur,tDel+tDur],[min(vVec)-30,min(vVec)-50],'r')
    ax.plot([tDel+tDur,tStop],[min(vVec)-50,min(vVec)-50],'r')
    ax.legend()
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane voltage (mV)')
    plt.show()

#
if __name__ == '__main__':
    RC()