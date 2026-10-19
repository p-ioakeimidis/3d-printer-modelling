import numpy as np
import matplotlib.pyplot as plt
import control as ct


##### System parameters #####
# Mechanical parameters
m = 0.5       #kg
J_m = 5e-6    #kgm^2
r_m = 0.006   #m
r_b = 0.006   #m  

## Viscous friction coefficients (per bearing)
mu_r = 1e-6  #N m s/rad
mu_lx = 0.5  #N s/m
mu_ly = 0.5  #N s/m

# DC motor parameters
R_m = 2.0    #Ω
L_m = 0.002  #mH
K_m = 0.05   #N m/A

# Equivalent inertia and friction coefficients
J_eq = J_m + m*r_m**2/2  #kgm^2

b_x = 4*mu_r/r_b**2 + 2*mu_lx
b_y = 2*mu_r/r_b**2 + 2*mu_ly



##### Decoupled state-space descriptions #####
A_x = np.array([
    [0, 1, 0],
    [0, -r_m**2*b_x/(2*J_eq), K_m*r_m/(2*J_eq)],
    [0, -2*K_m/(L_m*r_m), -R_m/L_m]
])

A_y = np.array([
    [0, 1, 0],
    [0, -r_m**2*b_y/(2*J_eq), K_m*r_m/(2*J_eq)],
    [0, -2*K_m/(L_m*r_m), -R_m/L_m]
])

B = np.array([[0], [0], [1/L_m]])

C = np.eye(3)  # All 3 states are returned
D = np.zeros((3, 1))

sys_x = ct.ss(A_x, B, C, D)
sys_y = ct.ss(A_y, B, C, D)


##### Elliptical path #####
# Semiaxes
a = 0.08  #m
b = 0.05  #m

# Center of ellipse
x_c = 0.10  #m
y_c = 0.10  #m

# Period and angular frequency
T = 5              #sec
omega = 2*np.pi/T  #rad/s

t = np.linspace(0, 5, 5000)


# Positions
x_ellipse = x_c + a*np.cos(omega*t)
y_ellipse = y_c + b*np.sin(omega*t)

# Time derivatives
xd_ellipse = -a*omega*np.sin(omega*t)
xdd_ellipse = -a*omega**2*np.cos(omega*t)
xddd_ellipse = a*omega**3*np.sin(omega*t)

yd_ellipse = b*omega*np.cos(omega*t)
ydd_ellipse = -b*omega**2*np.sin(omega*t)
yddd_ellipse = -b*omega**3*np.cos(omega*t)


# Coefficients to simplify notation
C3 = 2*L_m*J_eq/(K_m*r_m)

C2_x = L_m*r_m*b_x/K_m + 2*R_m*J_eq/(K_m*r_m)
C2_y = L_m*r_m*b_y/K_m + 2*R_m*J_eq/(K_m*r_m)

C1_x = R_m*r_m*b_x/K_m + 2*K_m/r_m
C1_y = R_m*r_m*b_y/K_m + 2*K_m/r_m


# Decoupled voltages
u_x = C3*xddd_ellipse + C2_x*xdd_ellipse + C1_x*xd_ellipse
u_y = C3*yddd_ellipse + C2_y*ydd_ellipse + C1_y*yd_ellipse

# Left and right motor voltages
u_l = (u_x + u_y)/2
u_r = (u_x - u_y)/2

# Currents at t = 0
ix_0 = (2*J_eq/(K_m*r_m))*xdd_ellipse[0] + (r_m*b_x/K_m)*xd_ellipse[0]

iy_0 = (2*J_eq/(K_m*r_m))*ydd_ellipse[0] + (r_m*b_y/K_m)*yd_ellipse[0]

# Initial states of the decoupled systems
X0_x = np.array([x_ellipse[0], xd_ellipse[0], ix_0])
X0_y = np.array([y_ellipse[0], yd_ellipse[0], iy_0])


# Simulation of the system
response_x = ct.forced_response(sys_x, T=t, U=u_x, initial_state=X0_x)

response_y = ct.forced_response(sys_y, T=t, U=u_y, initial_state=X0_y)


# Positions x, y
x_sim = response_x.states[0]
y_sim = response_y.states[0]

# Plot of elliptical path
plt.figure(1)
plt.plot(x_ellipse, y_ellipse, 'k--', label='Reference ellipse')
plt.plot(x_sim, y_sim, 'b', label='Simulated ellipse')

plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Elliptical path')

plt.axis('equal')
plt.grid(True)
#plt.legend()
plt.show(block=False)

# Plot of input voltages for ellipse
plt.figure(2)

plt.plot(t, u_l, label='Left motor voltage')
plt.plot(t, u_r, label='Right motor voltage')

plt.xlabel('Time [s]')
plt.ylabel('Voltage [V]')
#plt.title('Motor inputs for elliptical trajectory')

plt.grid(True)
plt.legend()
plt.show(block=False)




##### Lissajous curve #####
A_l = 0.08  #m
B_l = 0.1   #m
alpha_l = 3
beta_l = 2 
delta = np.pi/2
t = np.linspace(0,10,10000)

# Positions
x_l = x_c + A_l*np.sin(alpha_l*t + delta)
y_l = y_c+0.05+ B_l*np.sin(beta_l*t)

# Time derivatives
x_ld = -A_l*alpha_l*np.cos(alpha_l*t+delta)
x_ldd = -A_l*alpha_l**2*np.sin(alpha_l*t+delta)
x_lddd = -A_l*alpha_l**3*np.cos(alpha_l*t+delta)

y_ld = B_l*beta_l*np.cos(beta_l*t)
y_ldd = -B_l*beta_l**2*np.sin(beta_l*t)
y_lddd = -B_l*beta_l**3*np.cos(beta_l*t)

# Decoupled voltages
u_x_lissajous = C3*x_lddd + C2_x*x_ldd + C1_x*x_ld
u_y_lissajous = C3*y_lddd + C2_y*y_ldd + C1_y*y_ld

# Left and right motor voltages
u_l_lissajous = (u_x_lissajous + u_y_lissajous)/2
u_r_lissajous = (u_x_lissajous - u_y_lissajous)/2

# Currents at t = 0
ix_0_lissajous = (2*J_eq/(K_m*r_m))*x_ldd[0] + (r_m*b_x/K_m)*x_ld[0]

iy_0_lissajous = (2*J_eq/(K_m*r_m))*y_ldd[0] + (r_m*b_y/K_m)*y_ld[0]

# Initial states of the decoupled systems
X0_x_lissajous = np.array([x_l[0], x_ld[0], ix_0_lissajous])
X0_y_lissajous = np.array([y_l[0], y_ld[0], iy_0_lissajous])


# Simulation of the system
response_x_lissajous = ct.forced_response(sys_x, T=t, U=u_x_lissajous, initial_state=X0_x_lissajous)

response_y_lissajous = ct.forced_response(sys_y, T=t, U=u_y_lissajous, initial_state=X0_y_lissajous)


# Positions x, y
x_sim_lissajous = response_x_lissajous.states[0]
y_sim_lissajous = response_y_lissajous.states[0]

# Plot of elliptical path
plt.figure(3)
#plt.plot(x_ellipse, y_ellipse, 'k--', label='Reference ellipse')
plt.plot(x_sim_lissajous, y_sim_lissajous, 'b', label='Simulated Lissajous curve')

plt.xlabel('x [m]')
plt.ylabel('y [m]')
#plt.title('Elliptical path')

plt.axis('equal')
plt.grid(True)
#plt.legend()
plt.show(block=False)

# Plot of input voltages for ellipse
plt.figure(4)

plt.plot(t, u_l_lissajous, label='Left motor voltage')
plt.plot(t, u_r_lissajous, label='Right motor voltage')

plt.xlabel('Time [s]')
plt.ylabel('Voltage [V]')
#plt.title('Motor inputs for elliptical trajectory')

plt.grid(True)
plt.legend()
plt.show(block=False)

plt.show()