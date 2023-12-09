# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 10:24:06 2022

@author: pwerginz
"""

# ------------------------------------------------
# Hodgkin & Huxley single-compartment (local) model
# Run: timeStep,vVec,mVec,hVec,nVec = HH_single()
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
import time as _time
import sys
import matplotlib.pyplot as plt


### ------------- DEFINE FUNCTIONS FOR ALPHAS AND BETAS -------------
def alphaM(v,kTemp):
    return kTemp*(0.1*(v+40))/(1-np.exp(-(v+40)/10))
def betaM(v,kTemp):
    return kTemp*4*np.exp(-(v+65)/18)
def alphaH(v,kTemp):
    return kTemp*0.07*np.exp(-(v+65)/20)
def betaH(v,kTemp):
    return kTemp*1/(1+np.exp(-(v+35)/10))
def alphaN(v,kTemp):
    return kTemp*(0.01*(v+55))/(1-np.exp(-(v+55)/10))
def betaN(v,kTemp):
    return kTemp*0.125*np.exp(-(v+65)/80)


### Main function to solve model
def HH_single():
    ### ------------- Step 1: PARAMETERS -------------
    # Solver
    solver = 'FE' # Type of solver, 'FE'==Forward Euler or 'BE'==Backward Euler

    # Stimulus parameters
    I = 20 # Stimulus amplitude, in uA/cm2

    # Temporal parameters
    tStop = 10 # Total duration of simulation, in ms
    tDel = 1 # Time when stimulus starts, in ms
    tDur = 0.5 # Time of stimulus ON, in ms
    tDt = 0.1 # Time step, in ms

    # Biophysics
    c = 1 # Membrane specific capacitance, in uF/cm2, original=1

    # Temperature
    temp = 6.3 # Model temperature, in Celsius, original=6.3

    # HH parameters
    vInit = -65 # Membrane voltage to initialize the simulation, in mV
    gNa = 120 # Sodium channel maximum conductivity, in mS/cm2, original=120
    gK = 36 # Potassium channel maximum conductivity, in mS/cm2, original=26
    gL = 0.3 # Leak channel maximum conductivity, in mS/cm2, original=0.3
    eNa = 50 # Sodium reversal/equilibrium potential, in mV, original=50
    eK = -77 # Potassium reversal/equilibrium potential, in mV, original=-77
    eL = -54.3 # Leak reversal/equilibrium potential, in mV, original=-54.3

    # Plot flag
    plotV = 1 # Membrane voltage over time


    ### ------------- Step 2: COMPUTE ADDITIONAL PARAMETERS -------------
    # Compute time step
    timeSteps = int(tStop/tDt)+1
    timeStep = np.linspace(0,tStop,timeSteps)

    # Temperature adjustment
    kTemp = 3**((temp-6.3)/10)

    # Other constants
    vAdd = 0.001

    # Compute initial values
    v0 = vInit
    m0 = alphaM(v0,kTemp)/(alphaM(v0,kTemp)+betaM(v0,kTemp))
    h0 = alphaH(v0,kTemp)/(alphaH(v0,kTemp)+betaH(v0,kTemp))
    n0 = alphaN(v0,kTemp)/(alphaN(v0,kTemp)+betaN(v0,kTemp))

    # Allocate memory for v, m, h and n
    vVec = np.zeros(timeSteps)
    mVec = np.zeros(timeSteps)
    hVec = np.zeros(timeSteps)
    nVec = np.zeros(timeSteps)

    # Set initial values
    vVec[0] = v0
    mVec[0] = m0
    hVec[0] = h0
    nVec[0] = n0


    ### --------- Step 3: SOLVE ODE & POSTPROCESSING ---------
    Tic = _time.time()
    for t in range(0,timeSteps-1):

        # States at current time step
        vT = vVec[t]
        mT = mVec[t]
        hT = hVec[t]
        nT = nVec[t]

        # Stimulus current
        iStim = 0
        if t>=int(tDel/tDt) and t<int((tDel+tDur)/tDt):
            iStim = I # in uA/cm2

        # Ionic currents
        # Sodium
        iNa = gNa*mT**3*hT*(vT-eNa) # in (mS/cm2)*mV==uA/cm2
        # Potassium
        iK = gK*nT**4*(vT-eK) # in (mS/cm2)*mV==uA/cm2
        # Leak
        iL = gL*(vT-eL) # in (mS/cm2)*mV==uA/cm2
        # Sum
        iIon = iNa+iK+iL # in uA/cm2

        # Update v, m, h and n
        if solver=='FE':
            # Compute change of v
            vVec[t+1] = vT+(-iIon+iStim)*(tDt/c) # in mV

            # Update gating variables with new v
            mVec[t+1] = mT+(alphaM(vT,kTemp)*(1-mT)-betaM(vT,kTemp)*mT)*tDt
            hVec[t+1] = hT+(alphaH(vT,kTemp)*(1-hT)-betaH(vT,kTemp)*hT)*tDt
            nVec[t+1] = nT+(alphaN(vT,kTemp)*(1-nT)-betaN(vT,kTemp)*nT)*tDt
        elif solver=='BE':
            # Additonal ionic contribution needed for BE
            # Sodium
            iNaAux = gNa*mT**3*hT*(vT+vAdd-eNa) # in (mS/cm2)*mV==uA/cm2
            # Potassium
            iKAux = gK*nT**4*(vT+vAdd-eK) # in (mS/cm2)*mV==uA/cm2
            # Leak
            iLAux = gL*(vT+vAdd-eL) # in (mS/cm2)*mV==uA/cm2
            # Sum
            rhsdidv = (iNaAux-iNa+iKAux-iK+iLAux-iL)/vAdd # in (uA/cm2)/mV

            # Compute membrane voltage at next time step
            vVec[t+1] = (vT+(-iIon+rhsdidv*vT+iStim)*(tDt/c))/(1+rhsdidv*(tDt/c)) # in mV

            # Update gating variables with new v
            mVec[t+1] = (mT+tDt*alphaM(vVec[t+1],kTemp))/(1+tDt*(alphaM(vVec[t+1],kTemp)+betaM(vVec[t+1],kTemp)))
            hVec[t+1] = (hT+tDt*alphaH(vVec[t+1],kTemp))/(1+tDt*(alphaH(vVec[t+1],kTemp)+betaH(vVec[t+1],kTemp)))
            nVec[t+1] = (nT+tDt*alphaN(vVec[t+1],kTemp))/(1+tDt*(alphaN(vVec[t+1],kTemp)+betaN(vVec[t+1],kTemp)))
        else:
            print('--- Unknown solver type!')
            sys.exit()

    Toc = _time.time()
    print('--- Solving time was %.3f seconds' % (Toc - Tic))


    ### --------- Step 4: VISUALIZE RESULTS ---------
    # Membrane voltage over time
    if plotV:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid()
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.plot(timeStep,vVec)
        ca = plt.gca()
        ca.add_patch(plt.Rectangle((tDel,-100),tDur,10,facecolor='r'))
        plt.xlim(0,tStop)
        plt.ylim(-100, 60)
        plt.xlabel('Time (ms)')
        plt.ylabel('Membrane voltage (mV)')
        plt.show()


    ### --------- Step 5: RETURN STATE VARIABLES AND CURRENT DENSITIES ---------
    return timeStep, vVec, mVec, hVec, nVec


#
if __name__ == '__main__':
    HH_single()
