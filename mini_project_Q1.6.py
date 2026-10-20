import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

### Geometry and printing
R_c = 0.05      # m
N = 100         # number of segments
t_print = 20.0  # s
h = 0.0005      # filament square side [m]
T_dep = 200     # temperature of deposited filament [°C]
T_b = 60.0      # temperature of the heated bed [°C]
T_e = 25.0      # ambient temperature [°C]

### Filament properties
rho = 1240.0  # kg/m^3
c_p = 1800.0  # J/(kg K)
k = 0.13      # W/(m K)
h_b = 100.0   # W/(m^2 K)
h_c = 10.0    # W/(m^2 K)

T_eq = (T_b, T_e, h_b, h_c, k, R_c)  #equilibrium temperatures [°C]

dq = 2 * np.pi * R_c / N


def rhs(t, T, closed=False):
    # Conduction between neighboring segments
    if closed:
        lap = (np.roll(T, -1) - 2*T + np.roll(T, 1)) / dq**2
    else:
        lap = np.zeros(len(T))
        if len(T) > 1:
            lap[0] = (T[1] - T[0]) / dq**2
            lap[-1] = (T[-2] - T[-1]) / dq**2
            lap[1:-1] = (T[2:] - 2*T[1:-1] + T[:-2]) / dq**2
    return (k / (rho * c_p))*lap + (h_b / (rho * c_p * h))*(T_b - T) - (3 * h_c / (rho * c_p * h))*(T - T_e)


print_times = []
print_diff = []
# Printing: add one hot segment every t_print/N seconds
T = np.array([T_dep])
for i in range(1, N):
    sol = solve_ivp(rhs, (0, t_print/N), T)
    T = np.append(sol.y[:, -1], T_dep)
    # Temperature distributions at 25%, 50%, and 75% of the printing process for plotting
    if i == N//4 - 1:
        T_quarter = T.copy()
    if i == N//2 - 1:
        T_half = T.copy()
    if i == 3*N//4 - 1:
        T_three_quarters = T.copy()
    print_times.append((i + 1) * t_print / N)
    print_diff.append(np.max(T) - np.min(T))

# Finish cooling the last deposited segment up to t_print
T = solve_ivp(rhs, (0, t_print/N), T).y[:, -1]
T_full = T.copy()

# After printing: close the circle and continue until nearly uniform
cool_times = np.linspace(0, 80, 800)
sol = solve_ivp(lambda t, x: rhs(t, x, closed=True),
                (0, 200), T, t_eval=cool_times)

max_min_diff = np.max(sol.y, axis=0) - np.min(sol.y, axis=0)

reached = np.where(max_min_diff < 1.0)[0]
j = reached[0]
print(f'Printing completed in {t_print:.1f} s')
print(f'Temperature nearly uniform after {t_print + cool_times[j]:.1f} s total')


# Plot temperature distribution during and after printing
q = (np.arange(N) + 0.5) * dq
plt.figure(1)
plt.plot(q[:len(T_quarter)], T_quarter, label='25% printed')
plt.plot(q[:len(T_half)], T_half, label='50% printed')
plt.plot(q[:len(T_three_quarters)], T_three_quarters, label='75% printed')
plt.plot(q, T_full, label='100% printed')
plt.plot(q, sol.y[:, j], label='Nearly uniform')
plt.xlabel('Arc length q [m]')
plt.ylabel('Temperature [°C]')
plt.legend()
plt.grid()
plt.show(block=False)


# Plot of max-min temperature difference after printing
time_total = np.concatenate((print_times, t_print + sol.t))
diff_total = np.concatenate((print_diff, max_min_diff))

plt.figure(2)
plt.plot(time_total, diff_total, label='Temperature difference')

plt.axhline(1.0, color='r', linestyle='--',
            label='Uniformity threshold = 1.0 °C')

plt.axvline(t_print, color='k', linestyle='--',
            label='Printing completed')

plt.xlabel('Time [s]')
plt.ylabel('Max - Min Temperature [°C]')
plt.legend()
plt.grid()
plt.show(block=False)

plt.show()
