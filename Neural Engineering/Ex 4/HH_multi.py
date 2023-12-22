# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 10:43:38 2022

@author: pwerginz
"""

# ------------------------------------------------
# Hodgkin & Huxley multi-compartment (stick) model
# Run: timeStep,vMat,mMat,hMat,nMat,iNaMat,iKMat,iLMat = HH_multi()
# The computation is done by a fixed step size forward or backward Euler method
# Two stimulation modes: (1) Current injection into a single compartment ('iClamp')
#                        (2) Extracellular stimulation ('extracellular')
# ------------------------------------------------
#             
#                       _
#                  (2) / \     ||   
#                      \_/    _||_ (1)
#                             \  /                    
#   -------------------------- \/ --------------------------
#   |    1     |    2     |    3     |    4     |    5     |
#   --------------------------------------------------------
#
#
# Import necessary packages
import time as _time
import matplotlib.pyplot as plt
import numpy as np
import sys
import csv
# import scipy.sparse
# import scipy.sparse.linalg
    

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
def HH_multi():
    ### ------------- Step 1: PARAMETERS -------------
    # Solver
    solver = 'BE' # Type of solver, 'FE'==Forward Euler or 'BE'==Backward Euler
    
    # Stimulus parameters
    mode = 'extracellular' # Stimulus mode, 'iClamp' or 'extracellular'
    elecType = 'FEM' # Electrode type, 'point' or 'disk' or 'FEM', for 'extracellular' only
    dDisk = 50 # Electrode diameter, in um, for 'disk' only
    rhoExt = 300 # Extracellular resistivity, in Ohm*cm, for 'point and 'disk' only
    FEMfile = 'fiber-vertical.csv' # File with FEM Solution, for 'FEM' only
    I = 100 # Stimulus amplitude, in pA ('iClamp') or uA ('extracellular')
    cellX = 0 # X-shift of cell relative to the electrode (at 0/0)), in um
    cellY = 25 # Y-shift of cell relative to the electrode (at 0/0)), in um

    # Temporal parameters
    tStop = 10 # Total duration of simulation, in ms
    tDel = 1 # Time when stimulus starts, in ms
    tDur = 0.5 # Time of stimulus ON, in ms
    tDt = 0.0125 # Time step, in ms
    
    # Stick parameters, total length = (nComp-1) x lComp
    nComp = 201 # Number of compartments
    lComp = 10 # Compartment length, in um
    rComp = 1 # Compartment radius, in um

    # Biophysics
    c = 1 # Membrane specific capacitance, in uF/cm2
    rhoA = 100 # Axial/Intracellular conductivity, in Ohm*cm, original=100
    
    # Temperature
    temp = 6.3 # Model temperature, in Celsius, original=6.3
    
    # Hodgkin & Huxley parameters
    vInit = -65 # Membrane voltage to initialize the simulation, in mV
    gNa = 120 # Sodium channel maximum conductivity, in mS/cm2, original=120
    gK = 36 # Potassium channel maximum conductivity, in mS/cm2, original=36
    gL = 0.3 # Leak channel maximum conductivity, in mS/cm2, original=0.3
    eNa = 50 # Sodium reversal/equilibrium potential, in mV, original=50
    eK = -77 # Potassium reversal/equilibrium potential, in mV, original=-77
    eL = -54.3 # Leak reversal/equilibrium potential, in mV, original=-54.3

    # Plot flags
    plotV = 1 # Membrane voltage over time
    
    
    ### ------------- Step 2: COMPUTE ADDITIONAL PARAMETERS -------------
    # Central compartment of stick
    centerComp = int((nComp-1)/2+1)
    
    # Axial resistance of each compartment, in kOhm
    R = ((2*rhoA*lComp*1e-4)/(2*(rComp*1e-4)**2*np.pi))*1e-3
    # Compartment surface, in cm2
    A = (2*rComp*np.pi*lComp)*1e-8
    # Compartment membrane capacitance, in uF
    C = c*A
    
    # Time step
    timeSteps = int(tStop/tDt)+1
    timeStep = np.linspace(0,tStop,timeSteps)
    
    # Temperature adjustment
    kTemp = 3**((temp-6.3)/10)
    
    # Other constants
    vAdd = 0.001

    # Set up tridiagonal matrix for fast computation of iAxial (BE)
    matInvDiag = np.concatenate(([1+(tDt/C)*(1/R)],np.ones(nComp-2)+(tDt/C)*(2/R),[1+(tDt/C)*(1/R)])) # in (ms/uF)*(1/kOhm)==1/V
    matInvOffDiag = np.ones(nComp-1)*-(tDt/C)*(1/R)
    matInv = np.diag(matInvDiag,0) + np.diag(matInvOffDiag,-1) + np.diag(matInvOffDiag,1)
    
    # Set up tridiagonal matrix for fast iAxial (FE) and iStim (extracellular) computation
    matAxDiag = np.concatenate(([-1/R],np.ones(nComp-2)*-2/R,[-1/R])) # in 1/kOhm
    matAxOffDiag = np.ones(nComp-1)*1/R
    matAxial = np.diag(matAxDiag,0) + np.diag(matAxOffDiag,-1) + np.diag(matAxOffDiag,1)
    
    # Stimulus currents depending on stimulus type
    if mode=='iClamp': # Simple current conversion
        iStim = np.zeros(nComp)
        iStim[centerComp-1] = 1e-6*I/A # in 1e-6*pA/cm2==uA/cm2
    elif mode=='extracellular':
        # Compute potentials at compartment centers, electrode is located at 0/0
        x = np.linspace(-(centerComp-1)*lComp,(centerComp-1)*lComp,nComp)+cellX
        y = np.ones(nComp)*cellY
        # Compute extracellular potentials depending on electrode type
        if elecType=='point':
            # Euklidean distance for each compartment center
            compDist = 1e-4*np.sqrt(x**2+y**2) # in cm
            # Analytical potentials (point source) for given distance
            potentials = 1e-3*(rhoExt*I)/(4*np.pi*compDist) # in 1e-3*(Ohm*cm*uA)/cm==mV
        elif elecType=='disk':
            # Axial and radial distance for each compartment center
            r = 1e-4*x # in cm
            z = 1e-4*y # in cm
            # Analytical potentials (disk electrode) for given distance
            rDisk = 1e-4*dDisk/2 # in cm
            potentials = 1e-3*(2.0*rhoExt*I)/(4.0*np.pi*rDisk) * np.arcsin((2*rDisk)/(np.sqrt((r-rDisk)**2+z**2)+np.sqrt((r+rDisk)**2+z**2))) # in 1e-3*(Ohm*cm*uA)/cm==mV
        elif elecType=='FEM':
            # Numerically computed potentials (arbitrary shape)
            file = open(FEMfile) # Load exported solution
            csv_reader = csv.reader(file,delimiter=';')
            next(csv_reader) # Skip header
            xx = []
            Ve = []
            for row in csv_reader:
                xx.append(float(row[10])*1e6) # FEM x-coordinates, conversion to um
                Ve.append(float(row[9])*1e3) # FEM potentials, conversion to mV
            xx = np.array(xx)
            Ve = np.array(Ve)
            potentials = I*np.interp(x,xx,Ve) # Interpolate to compartment centers
        # Matrix x vector for stimulus current
        iStim = np.matmul(matAxial,potentials)/A # in ((1/kOhm)*mV)/cm2==uA/cm2
    else:
        print('--- Unknown stimulus mode!')
        sys.exit()
    
    # Compute initial values
    v0 = vInit
    m0 = alphaM(v0,kTemp)/(alphaM(v0,kTemp)+betaM(v0,kTemp))
    h0 = alphaH(v0,kTemp)/(alphaH(v0,kTemp)+betaH(v0,kTemp))
    n0 = alphaN(v0,kTemp)/(alphaN(v0,kTemp)+betaN(v0,kTemp))
    
    # Allocate memory for v, m, h and n
    vMat = np.zeros((nComp,timeSteps))
    mMat = np.zeros((nComp,timeSteps))
    hMat = np.zeros((nComp,timeSteps))
    nMat = np.zeros((nComp,timeSteps))
    
    # Set initial values
    vMat[:,0] = v0
    mMat[:,0] = m0
    hMat[:,0] = h0
    nMat[:,0] = n0
    

    ### --------- Step 3: SOLVE ODE & POSTPROCESSING ---------
    Tic = _time.time()
    for t in range(0,timeSteps-1):
        
        # States at current time step
        vVecT = vMat[:,t]
        mVecT = mMat[:,t]
        hVecT = hMat[:,t]
        nVecT = nMat[:,t]
    
        # Stimulus current
        iStimVec = np.zeros(nComp)
        if t>=int(tDel/tDt) and t<int((tDel+tDur)/tDt):
            iStimVec = iStim # in uA/cm2
            
        # Ionic currents
        # Sodium
        iNaVec = gNa*mVecT**3*hVecT*(vVecT-eNa) # in (mS/cm2)*mV==uA/cm2
        # Potassium
        iKVec = gK*nVecT**4*(vVecT-eK) # in (mS/cm2)*mV==uA/cm2
        # Leak
        iLVec = gL*(vVecT-eL) # in (mS/cm2)*mV==uA/cm2
        # Sum
        iIonVec = iNaVec+iKVec+iLVec # in uA/cm2    
        
        # Update gating variables with new v
        if solver=='FE':
            # Axial currents
            iAxialVec = np.matmul(matAxial,vVecT)/A # in uA/cm2

            # Compute change of v
            vMat[:,t+1] = vVecT+(-iIonVec+iAxialVec+iStimVec)*(tDt/c) # in mV

            # Update gating variables with new v
            mMat[:,t+1] = mVecT+(alphaM(vVecT,kTemp)*(1-mVecT)-betaM(vVecT,kTemp)*mVecT)*tDt
            hMat[:,t+1] = hVecT+(alphaH(vVecT,kTemp)*(1-hVecT)-betaH(vVecT,kTemp)*hVecT)*tDt
            nMat[:,t+1] = nVecT+(alphaN(vVecT,kTemp)*(1-nVecT)-betaN(vVecT,kTemp)*nVecT)*tDt
        elif solver=='BE':
            # Additonal ionic contribution needed for BE
            # Sodium
            iNaAuxVec = gNa*mVecT**3*hVecT*(vVecT+vAdd-eNa) # in (mS/cm2)*mV==uA/cm2
            # Potassium
            iKAuxVec = gK*nVecT**4*(vVecT+vAdd-eK) # in (mS/cm2)*mV==uA/cm2
            # Leak
            iLAuxVec = gL*(vVecT+vAdd-eL) # in (mS/cm2)*mV==uA/cm2
            # Sum
            rhsdidvVec = (iNaAuxVec-iNaVec+iKAuxVec-iKVec+iLAuxVec-iLVec)/vAdd # in (uA/cm2)/mV

            # Compute change of v
            # Right hand side of matrix equation to be solved
            RHS = vVecT+(-iIonVec+rhsdidvVec*vVecT+iStimVec)*(tDt/c)
            # Add ionic current contribution to left hand side
            np.fill_diagonal(matInv,matInvDiag+rhsdidvVec*(tDt/c))
            
            # Solve matrix equation
            vMat[:,t+1] = np.linalg.solve(matInv,RHS) # in mV
            # vMat[:,t+1] = scipy.sparse.linalg.spsolve(scipy.sparse.csr_matrix(matInv),RHS) # in mV

            # Restore inverse Matrix (not very elegant solution...)
            np.fill_diagonal(matInv,matInvDiag-rhsdidvVec*(tDt/c))
            
            # Update gating variables with new v
            mMat[:,t+1] = (mVecT+tDt*alphaM(vMat[:,t+1],kTemp))/(1+tDt*(alphaM(vMat[:,t+1],kTemp)+betaM(vMat[:,t+1],kTemp)))
            hMat[:,t+1] = (hVecT+tDt*alphaH(vMat[:,t+1],kTemp))/(1+tDt*(alphaH(vMat[:,t+1],kTemp)+betaH(vMat[:,t+1],kTemp)))
            nMat[:,t+1] = (nVecT+tDt*alphaN(vMat[:,t+1],kTemp))/(1+tDt*(alphaN(vMat[:,t+1],kTemp)+betaN(vMat[:,t+1],kTemp)))
        else:
            print('--- Unknown solver type!')
            sys.exit()
            
    Toc = _time.time()
    print('--- Solving time was %.3f seconds' % (Toc - Tic))
    
    # Compute ionic current densities by using the state variables v, m, h and n
    iNaMat = gNa*mMat**3*hMat*(vMat-eNa) # in (mS/cm2)*mV==uA/cm2
    iKMat = gK*nMat**4*(vMat-eK) # in (mS/cm2)*mV==uA/cm2
    iLMat = gL*(vMat-eL) # in (mS/cm2)*mV==uA/cm2
    
    
    ### --------- Step 4: VISUALIZE RESULTS ---------
    if plotV:
        # Membrane voltage over time
        fig = plt.figure()

        # Membrane voltage over time (spatial)
        ax2 = fig.add_subplot()
        ax2.spines["right"].set_visible(False)
        ax2.spines["top"].set_visible(False)
        offset = -vInit;
        for i in range(nComp):
            if i==centerComp-1:
                ax2.plot(timeStep,vMat[i,:]+offset,'r');
            else:
                ax2.plot(timeStep,vMat[i,:]+offset,'b',lw=0.25);
            offset = offset-lComp;
        ca = plt.gca()
        ca.add_patch(plt.Rectangle((tDel,offset-150),tDur,50,facecolor='r'))
        plt.xlim(0,tStop)
        plt.ylim([offset-150,150])
        plt.xlabel('Time (ms)')
        plt.ylabel('Location along stick (um)')
        plt.show()
 

    ### --------- Step 5: RETURN STATE VARIABLES AND CURRENT DENSITIES --------- 
    return timeStep, vMat, mMat, hMat, nMat, iNaMat, iKMat, iLMat


#
if __name__ == '__main__':
    HH_multi()
