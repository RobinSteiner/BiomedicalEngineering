import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from HH_single import HH_single

def binary_search_threshold(time_step, stim_duration, pulse_amplitude):
    low, high = 0, 1000  # Initial search range for threshold
    threshold = None

    while high - low > 1e-6:
        current_amplitude = (low + high) / 2
        _, v, _, _, _ = HH_single(I=current_amplitude, tDur=stim_duration)
        max_voltage = np.max(v)

        if max_voltage > 0:
            high = current_amplitude
            threshold = current_amplitude
        else:
            low = current_amplitude

    return threshold

# Parameters
pulse_durations = np.logspace(-2, 2, 20)  # Pulse durations from 0.01 to 100 ms in 20 logarithmic steps
stim_duration_long_pulse = 100  # Duration of the long pulse for rheobase calculation

# Compute SD curve
thresholds = []

for duration in pulse_durations:
    threshold = binary_search_threshold(0.025, duration, 15)
    thresholds.append(threshold)

# Plot SD curve
plt.figure()
plt.semilogx(pulse_durations, thresholds, marker='o', linestyle='-', color='b')
plt.xlabel('Pulse Duration (ms)')
plt.ylabel('Threshold Current (uA/cm2)')
plt.title('Strength-Duration Curve')
plt.grid(True)
plt.show()

# Extract rheobase and chronaxie
rheobase = binary_search_threshold(0.025, stim_duration_long_pulse, 15)
double_rheobase = 2 * rheobase
index_double_rheobase = np.argmin(np.abs(np.array(thresholds) - double_rheobase))
chronaxie = pulse_durations[index_double_rheobase]

print(f'Rheobase: {rheobase} uA/cm2')
print(f'Chronaxie: {chronaxie} ms')
