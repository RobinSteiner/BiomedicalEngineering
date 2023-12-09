# Import necessary packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from sympy import symbols, Function, dsolve, Eq

### Main function to solve model
def solveODE():
    ### ------------- Step 1: PARAMETERS -------------
    y0 = 9  # Initial condition
    tStop = 10  # Total duration of simulation
    tDt = 0.025  # Time step
    
    ### ------------- Step 2: COMPUTE ADDITIONAL PARAMETERS -------------
    # Compute time step
    timeSteps = int(tStop/tDt)+1
    timeStep = np.linspace(0, tStop, int(tStop/tDt)+1)
    
    ### ------------- Step 3: ANALYTICAL -------------
    t_sym = symbols('t')
    y_sym = Function('y')
    ode = Eq(y_sym(t_sym).diff(t_sym), -3*y_sym(t_sym) + 9*t_sym)
    analytical_solution = dsolve(ode, y_sym(t_sym), ics={y_sym(0): y0})
    
    # Convert the analytical solution to a function for numerical evaluation
    y_analytical = analytical_solution.rhs
    y_analytical_func = lambda t: y_analytical.subs(t_sym, t)
    
    ### ------------- Step 4: FORWARD EULER ----------
    yFEVec = np.zeros(timeSteps)  # Allocate memory
    yFEVec[0] = y0  # Set initial condition
    # Loop over time
    for t in range(timeSteps-1):
        yFEVec[t+1] = yFEVec[t] + (-3*yFEVec[t] + 9*timeStep[t])*tDt
    
    ### ------------- Step 5: BACKWARD EULER ---------
    yBEVec = np.zeros(timeSteps)  # Allocate memory
    yBEVec[0] = y0  # Set initial condition
    # Loop over time
    for t in range(timeSteps-1):
        yBEVec[t+1] = (yBEVec[t] + 9*timeStep[t+1]*tDt) / (1 + 3*tDt)

    ### ------------- Step 6: RK45 SOLVER ---------
    def odefun(t, y): return -3 * y + 9 * t
    sol = integrate.solve_ivp(odefun, [0, tStop], [y0], method='RK45')

    ### ------------- Step 7: VISUALIZE RESULTS ---------
    # y over time
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.plot(timeStep, yFEVec, label='Forward Euler')
    ax.plot(timeStep, yBEVec, label='Backward Euler')
    ax.plot(sol.t, sol.y[0], label='RK45')
    ax.plot(timeStep, [y_analytical_func(t) for t in timeStep], label='Analytical', linestyle='--')
    plt.legend()
    plt.xlim(0, tStop)
    plt.ylim(min(-0.1, min(yFEVec)), max(10, max(yFEVec)))
    plt.xlabel('Time')
    plt.ylabel('Y')
    plt.show()

# Run the solver
if __name__ == '__main__':
    solveODE()
    