import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass   

@dataclass
class PrinterParams:
    """Parameters for the 3D Printer Hybrid System."""

    m: float = 0.5   #kg    
    J_m: float = 5e-6  #kgm^2 b
    r_m: float = 0.006   #m
    r_b: float = 0.006 #m 

    ## Viscous friction coefficients (per bearing)
    mu_r:float = 1e-6  #N m s/rad
    mu_lx :float = 0.5  #N s/m
    mu_ly :float = 0.5  #N s/m

    # DC motor parameters
    R_m :float= 2.0    #Ω
    L_m :float = 0.002  #mH
    K_m :float= 0.05   #N m/A
    
    # Thermal Parameters (Nozzle)
    L_n: float = 0.01    # H
    C_n: float = 0.01    # F
    C_tn: float = 50.0   # Thermal capacity of the nozzle (Joules/Celsius)
    R_tn: float = 1.1    # Thermal resistance to environment (Celsius/Watt)
    R_n: float = 1.0     # Electrical resistance of the nozzle heater (Ohms)
    T_e: float = 25.0    #  Environment temperature (Celsius)
    T_print: float = 235.0  #operating/high-temperature threshold (degree)
    T_low: float = 30.0     # cooling threshold (degree)
    u_h: float = 14.2    # Heater voltage applied when ON (Volts)

    # motor input voltage for the simulation
    u_print: float = 0.3         

###physical modelling###
def f(t, x, mode, p,side):
    theta_l = x[0] 
    theta_r = x[1]
    dx_1 = x[2] # Motor velocities
    dx_2 = x[3] # Motor velocities
    i_l = x[4]# Motor currents
    i_r= x[5]   # Motor currents
    i_L = x[6] # Nozzle L 
    u_C = x[7] # Nozzle C
    T_n =  x[8] # Nozzle thermal state

    if mode == "q0":
        if side == 0:       # Right
            u_l = p.u_print
            u_r = p.u_print

        elif side == 1:     # Up
            u_l = p.u_print
            u_r = -p.u_print

        elif side == 2:     # Left
            u_l = -p.u_print
            u_r = -p.u_print

        else :             # Down
            u_l = -p.u_print
            u_r = p.u_print

    u_n = 0.0 if mode in ("q4","q5") else p.u_h # when u_n is 0 , the nozzle input voltage is 0
    

    # Dynamic mechanical parameters
    p.J_eq = p.J_m + p.m * p.r_m**2 / 2.0
    p.b_x = 4.0 * p.mu_r / p.r_b**2 + 2.0 * p.mu_lx
    p.b_y = 2.0 * p.mu_r / p.r_b**2 + 2.0 * p.mu_ly

    # Motor Mechanical Dynamics
    dx_3 = -(p.r_m**2 / (p.J_eq * 4.0)) * ((p.b_x + p.b_y) * dx_1 + (p.b_x - p.b_y) * dx_2) + (p.K_m * i_l) / p.J_eq
    dx_4 = -(p.r_m**2 / (p.J_eq * 4.0)) * ((p.b_x - p.b_y) * dx_1 + (p.b_x + p.b_y) * dx_2) + (p.K_m * i_r) / p.J_eq

    # Motor Electrical Dynamics
    dx_5 = (u_l - p.R_m * i_l - p.K_m * dx_1) / p.L_m
    dx_6 = (u_r - p.R_m * i_r - p.K_m * dx_2) / p.L_m

    # Q1.5 Nozzle Electrical & Thermal Dynamics
    dx_7 = (u_n - u_C) / p.L_n
    dx_8 = (1.0 / p.C_n) * (i_L - u_C / p.R_n)
    dx_9 = (1.0 / p.C_tn) * ((u_C**2 / p.R_n) - (T_n - p.T_e) / p.R_tn)

    # Cartesian Velocities & Odometer Integration (dx_10)
    v_x = (p.r_m / 2.0) * (dx_1 + dx_2)
    v_y = (p.r_m / 2.0) * (dx_1 - dx_2)
    
    if mode in ("q0","q3") :
        dx_10 = np.sqrt(v_x**2 + v_y**2)
    else:
        dx_10 = 0.0

    return np.array([dx_1, dx_2, dx_3, dx_4, dx_5, dx_6, dx_7, dx_8, dx_9, dx_10])

def get_xy(x, p):
    pos_x = (p.r_m / 2.0) * (x[0] + x[1])
    pos_y = (p.r_m / 2.0) * (x[0] - x[1])
    return pos_x, pos_y

def g_corner_0(t, x, mode, p,side):

    pos_x, pos_y = get_xy(x, p)

    if side == 0:
        return pos_x - (x_start + length) 
    # tool head position in x axis will equal to x_start + length = 0.05 + 10 = 0.15, so it will move to (0.15,0) starting from x_start, because for side= 0 we give the same input voltage in motors. The same logic applies to rest.
g_corner_0.terminal = True
g_corner_0.direction = 1
def g_corner_1(t, x, mode, p,side):
    pos_x, pos_y = get_xy(x, p)
    if side == 1:
        return pos_y - (y_start + length)
g_corner_1.terminal = True
g_corner_1.direction = 1
def g_corner_2(t, x, mode, p,side):
    pos_x, pos_y = get_xy(x, p)
    if side == 2:
        return x_start - pos_x
g_corner_2.terminal = True
g_corner_2.direction = 1
def g_corner_3(t, x, mode, p,side):
    pos_x, pos_y = get_xy(x, p)
    if side == 3:
        return y_start - pos_y
g_corner_3.terminal = True
g_corner_3.direction = 1

p = PrinterParams()
t_curr, t_end = 0.0, 100.0
length = 0.10 # the length of the each side of the square
x_start = 0.05 # we start at (x,y) = 0.05 m
y_start = 0.05
theta_l0 = (x_start + y_start )/ p.r_m
theta_r0 = (x_start - y_start)/ p.r_m
x_curr = np.array([theta_l0, theta_r0, 0.0, 0.0,0.0, 0.0, p.u_h / p.R_n, p.u_h, p.T_print, 0.0], dtype=float) # 
current_side = 0 # we first start by printing the bottom side by moving only on the x - axis
current_fsm = "q0" # we printing 

ts_log, ys_log, state_log = [], [], []

print("Starting Hybrid Simulation Loop...")
while  current_side <4 and t_curr < t_end:
    if current_side ==0:
        event = g_corner_0
    elif current_side ==1:
        event = g_corner_1
    elif current_side ==2 :
        event = g_corner_2
    else :
        event = g_corner_3
    sol = solve_ivp(
        f,
        (t_curr, t_end),
        x_curr,
        args=(current_fsm, p, current_side),
        events=event
    )

    ts_log.append(sol.t)
    ys_log.append(sol.y)
    state_log.extend([current_side] * len(sol.t))
    t_curr = sol.t[-1]
    x_curr = sol.y[:, -1].copy()
    if sol.status == 1:
        x_curr[2:6] = 0.0
        current_side = current_side + 1

t_total = np.concatenate(ts_log)
y_total = np.concatenate(ys_log, axis=1)
x_pos, y_pos = get_xy(y_total, p)

plt.figure()
plt.plot(x_pos, y_pos)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Square path")
plt.axis("equal")
plt.grid()
plt.show()

