import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from HH_single import HH_single

# Logistic function
def logistic_function(x, L, k, x0):
    return L / (1 + np.exp(-k * (x - x0)))

# Parameters
amplitudes = np.arange(0, 30, 1)  # Stimulus amplitudes from 0 to 30 μA/cm² in 1 μA/cm² steps
num_simulations = 50
tDur = 0.5

# Function to run simulation and compute spiking probability
def run_simulation_and_compute_probability(amplitude):
    spiking_count = 0

    for _ in range(num_simulations):
        _, v, _, _, _ = HH_single(I=amplitude, tDur=tDur)
        if np.max(v) > 0:
            spiking_count += 1

    spiking_probability = spiking_count / num_simulations * 100
    return spiking_probability

# Compute spiking probabilities
spiking_probabilities = [run_simulation_and_compute_probability(amplitude) for amplitude in amplitudes]

# Fit a logistic curve to the spiking probabilities
params, _ = curve_fit(logistic_function, amplitudes, spiking_probabilities, p0=[100, 1, 15])

# Plot results
plt.figure()
plt.scatter(amplitudes, spiking_probabilities, label='Simulated Data', color='blue')
plt.plot(amplitudes, logistic_function(amplitudes, *params), '--', label='Logistic Curve Fit', color='red')
plt.xlabel('Stimulus Amplitude (μA/cm²)')
plt.ylabel('Spiking Probability (%)')
plt.title('Spiking Probability vs. Stimulus Amplitude')
plt.legend()
plt.grid(True)
plt.show()

print(f'Logistic Curve Fit Parameters: L = {params[0]}, k = {params[1]}, x0 = {params[2]}')
