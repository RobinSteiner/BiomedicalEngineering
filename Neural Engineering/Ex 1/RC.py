# Import necessary packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

### Main function to solve model
def RC():
    ### ------------- Step 1: PARAMETERS -------------
    # Solver
    solvers = [#'BE', 'FE',
        'RK45',
        #'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA'
    ]  # Type of solvers: 'BE'==Backward Euler, 'FE'==Forward Euler, 'RK45'==Built-in solver
    showTimeSteps = True

    # Stimulus parameters
    I = 1  # Stimulus amplitude, in nA

    # Temporal parameters
    tStop = 200  # Total duration of simulation, in ms
    tDel = 50  # Delay until stimulus starts, in ms
    tDur = 100  # Duration of stimulus ON, in ms
    tDt = 5  # Time step, in ms

    # Membrane parameters
    R = 30  # Resistance, in MOhm
    C = 0.5  # Capacitance, in nF

    ### ------------- Step 2: COMPUTE ADDITIONAL PARAMETERS -------------
    # Compute time step
    timeSteps = int(tStop / tDt) + 1
    timeStep = np.linspace(0, tStop, timeSteps)

    # Allocate memory for v for each solver
    vVec_solvers = {solver: np.zeros(timeSteps) for solver in solvers}
    timeSteps_solvers = {solver: np.zeros(timeSteps) for solver in solvers}

    ### --------- Step 3: SOLVE ODE & POSTPROCESSING ---------
    for solver in solvers:
        if solver != 'FE' and solver != 'BE':
            # Use scipy's solve_ivp to solve the ODE system for the built-in solver
            sol = solve_ivp(
                lambda t, v: ode_system(t, v, I, R, C, tDel, tDur, tDt, solver),
                [0, tStop],
                [0],
                rtol=5e-5,
                atol=5e-7,
                method=solver
            )
            # Extract the solution
            vVec_solvers[solver] = sol.y[0]
            timeSteps_solvers[solver] = sol.t
        else:
            # Solve the ODE system for Forward Euler and Backward Euler
            vVec = np.zeros(timeSteps)
            for t in range(0, timeSteps - 1):
                IStim = 0
                if t >= int(tDel / tDt) and t < int((tDel + tDur) / tDt):
                    IStim = I  # in nA

                if solver == 'FE':
                    vVec[t + 1] = vVec[t] + ((-vVec[t] / R + IStim) / C) * tDt
                elif solver == 'BE':
                    vVec[t + 1] = (vVec[t] + IStim * (tDt / C)) / (1 + tDt / (R * C))

            vVec_solvers[solver] = vVec
            timeSteps_solvers[solver] = timeStep

    ### --------- Step 4: VISUALIZE RESULTS ---------
    # Membrane voltage over time for each solver
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    for solver in solvers:
        plot, = ax.plot(timeSteps_solvers[solver], vVec_solvers[solver], label=solver+' Nr. Timesteps='+str(timeSteps_solvers[solver].size))
        if showTimeSteps:
            for t_step in timeSteps_solvers[solver]:
                ax.axvline(x=t_step, color=plot.get_color(), linestyle='--', linewidth=0.8)
    
    ax.plot([0, tDel], [min(vVec_solvers[solver]) - 50, min(vVec_solvers[solver]) - 50], 'r', label='Input')
    ax.plot([tDel, tDel], [min(vVec_solvers[solver]) - 50, min(vVec_solvers[solver]) - 30], 'r')
    ax.plot([tDel, tDel + tDur], [min(vVec_solvers[solver]) - 30, min(vVec_solvers[solver]) - 30], 'r')
    ax.plot([tDel + tDur, tDel + tDur], [min(vVec_solvers[solver]) - 30, min(vVec_solvers[solver]) - 50], 'r')
    ax.plot([tDel + tDur, tStop], [min(vVec_solvers[solver]) - 50, min(vVec_solvers[solver]) - 50], 'r')
    ax.legend()
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane voltage (mV)')
    plt.show()


def ode_system(t, v, I, R, C, tDel, tDur, tDt, solver):
    # Stimulus current
    IStim = 0
    if t >= tDel and t < tDel + tDur:
        IStim = I  # in nA

    # Compute change of v
    if solver == 'FE':
        return (-v / R + IStim) / C
    elif solver == 'BE':
        return (v + IStim * tDt / C) / (1 + tDt / (R * C))
    else:
        return (-v / R + IStim) / C


# Run the model
if __name__ == '__main__':
    RC()
    