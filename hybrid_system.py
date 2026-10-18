import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass   

@dataclass
class PrinterParams:
    """Parameters for the 3D Printer Hybrid System."""
    # Mechanical & Motor Parameters (CoreXY)
    r_m: float = 0.01   # Motor pulley radius (m)
    J_eq: float = 0.1   # Equivalent inertia
    b_x: float = 0.05   # Friction in x-direction
    b_y: float = 0.05   # Friction in y-direction
    K_m: float = 0.1    # Motor constant
    L_m: float = 0.01   # Motor inductance (H)
    R_m: float = 1.0    # Motor resistance (Ohms)
    
    # Thermal Parameters (Nozzle)
    C_tn: float = 50.0  # Thermal capacity of the nozzle (Joules/Celsius)
    R_tn: float = 1.1   # Thermal resistance to environment (Celsius/Watt)
    R_n: float  = 1.0   # Electrical resistance of the nozzle heater (Ohms)
    T_e: float  = 20.0  # Environment temperature (Celsius)
    u_h: float  = 12.0  # Heater voltage applied when ON (Volts)

###physical modelling###
def f(x, u, q, p):
    theta_l = x[0]
    theta_r = x[1]
    theta_l_dot = x[2]
    theta_r_dot = x[3]

    i_l = x[4]
    i_r = x[5]
    T_n = x[6]

    u_l = u[0]
    u_r = u[1]


    dx_1 = theta_l_dot
    dx_2 = theta_r_dot

    dx_3 = -(p.r_m**2/(p.J_eq*4))*((p.b_x+p.b_y)*dx_1+(p.b_x-p.b_y)*dx_2)+(p.K_m*i_l)/p.J_eq
    dx_4 = -(p.r_m**2/(p.J_eq*4))*((p.b_x-p.b_y)*dx_1+(p.b_x+p.b_y)*dx_2)+(p.K_m*i_r)/p.J_eq

    dx_5 = (u_l-p.R_m*i_l-p.K_m*dx_1)/p.L_m
    dx_6 = (u_r-p.R_m*i_r-p.K_m*dx_2)/p.L_m

    #cartesian velocities 
    v_x = (p.r_m / 2) * (x[2] + x[3])
    v_y = (p.r_m / 2) * (x[2] - x[3])

    if q == 1:  # Heater is on
        dT_n = p.u_h**2 /(p.R_n * p.C_tn)  - (T_n - p.T_e) / (p.C_tn * p.R_tn)
    else:  # Heater is off (u_h = 0)
        dT_n = -(T_n - p.T_e) / (p.C_tn * p.R_tn)   

    #Integrate distance only when printing or wiping (u_l and u_r are active)
    if q==1 and (abs(u[0])>0 or abs(u[1])>0):
        dx_7 = np.sqrt(v_x**2 + v_y**2)
    else :
        dx_7 = 0.0

    return np.array([dx_1, dx_2, dx_3, dx_4, dx_5, dx_6, dT_n,dx_7]) 
    #x[7] is the distance printed of the filament 
    #dT_n has different eguation depending the q 
    ###return dx/dt = f(x,u)

###################################
    ##guards##
###################################

#F
#solve_ivp integrates dx_7 (velocity) starting from 0.0 , we take x[7] position/distance and when the value is equal to MAX_FILAMENT, then this guard is triggered
#the solver updates the value of x[7] in every step , actually x[7] is a discrete value , while dx_7 is a continiues one (dx/ dt)
def g_filament_finished(t,x,u,q,p):
    MAX_FILAMENT = 10.0
    return x[7] - MAX_FILAMENT
g_filament_finished.terminal = True 


## trigger when the cartesian position to the origin is 0
#i_reaced_0_0
def g_reached_0_0(t,x,u,q,p):
    pos_x = (p.r_m/2)*(x[0]+x[1])
    pos_y= (p.r_m/2)*(x[0]-x[1])
    return np.sqrt(pos_x**2 + pos_y**2) - 0.001 
g_reached_0_0.terminal = True


#here we measure the distance traveled by the nozzle , so that it can do the wiped-process (to clean itsself from old filament)
#so again x[7] has the total traveled distance of the nozzle 
#when this distance reaches WIPED_THRESHOLD, then this guard is triggered
def g_filament_wiped(t,x,u,q,p):
    WIPED_THRESHOLD = 0.05
    return x[7] - WIPED_THRESHOLD
g_filament_wiped.terminal = True 
    
#with q == 0, which means that theres not voltage to the heater, so the dx_6 term is continiuously being decreasing (we have negative sign for q=0) from ex. 300oC to the TARGET_TEMPERATURE
def g_temp_low(t,x,u,q,p):
    TARGET_TEMPERATURE = 30.0
    return x[6] - TARGET_TEMPERATURE
g_temp_low.Terminal = True

#last_state_t is the exact timestamp the temp_low state finished
#t is the the master clock of the simulation, in solve_ivp
# when we enter in q5 we capture the current time (t_curr)
def g_loading_done(t,x,u,q,p):
    last_state_t = u[2]
    return t - (last_state_t - 5.0)
g_loading_done.terminal = True

#here for q=1 we have the heater on so we have positive sign and then the temperature will be increasing until the TARGET_TEMPERATURE§
def g_temp_high(t,x,u,q,p):
    TARGET_TEMPERATURE = 200.0
    return x[6] - TARGET_TEMPERATURE
g_temp_high.terminal = True


###################################
##jumps##
###################################





##initial_mode###
##simulate using  solve_ivp()####



# p = PrinterParams()  # Create an instance of the PrinterParams dataclass
# x0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, p.T_e])  # Initial state vector
# fsm_start_state = 'q0_printing'