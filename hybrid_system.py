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
    T_high: float = 235.0  #operating/high-temperature threshold (degree)
    T_low: float = 30.0     # cooling threshold (degree)
    u_h: float = 14.2    # nozzle input voltage voltage applied 

    #hybrid-system event thresholds 
    print_distance: float = 5.0       # the available length of filament that can be printed  (meters)
    wipe_distance: float = 0.10     # the assumed length of filament that it is need to be printed to execute the wipe process(meters)
    position_tol: float = 0.001       # the targeted position of the head (meters)
    loading_time: float = 10.0         # the assumed loading time of the filament

    # motor input voltages for the simulation
    u_print: float = 0.3         # cosine input for back-and-forth printing (volts)
    print_period: float = 12.0       # seconcs
    u_move: float = 0.2             #constant input for moving and wiping (volts)

###physical modelling###
def f(t, x, mode, p):
    theta_l = x[0] 
    theta_r = x[1]
    dx_1 = x[2] # Motor velocities
    dx_2 = x[3] # Motor velocities
    i_l = x[4]# Motor currents
    i_r= x[5]   # Motor currents
    i_L = x[6] # Nozzle L 
    u_C = x[7] # Nozzle C
    T_n =  x[8] # Nozzle thermal state

    if mode == "q0" :
        u_l  = u_r = p.u_print * np.cos(2.0 * np.pi * t / p.print_period)
    elif mode == "q2":
        u_l  = u_r = -p.u_move
    elif mode =="q3":
        u_l = u_r = p.u_move
    elif mode == "q7":
        u_l = u_r = direction * p.u_move 
    else:
        u_l = u_r = 0.0

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

    # Cartesian Velocities & accumulated path length (dx_10)
    v_x = (p.r_m / 2.0) * (dx_1 + dx_2)
    v_y = (p.r_m / 2.0) * (dx_1 - dx_2)
    
    if mode in ("q0","q3") :
        dx_10 = np.sqrt(v_x**2 + v_y**2)
    else:
        dx_10 = 0.0

    return np.array([dx_1, dx_2, dx_3, dx_4, dx_5, dx_6, dx_7, dx_8, dx_9, dx_10])

#helper function
#Convert motor angles to Cartesian toolhead position
def get_xy(x, p):
    pos_x = (p.r_m / 2.0) * (x[0] + x[1])
    pos_y = (p.r_m / 2.0) * (x[0] - x[1])
    return pos_x, pos_y


##############
##guards###
##############

#solve_ivp integrates dx_7 (velocity) starting from 0.0 , we take x[9] position/distance and when the value is equal to p.print_distance, then this guard is triggered
#the solver updates the value of x[9] in every step

def g_filament_finished(t, x, mode, p):
    return x[9] - p.print_distance# 
g_filament_finished.terminal = True
g_filament_finished.direction = 1

def g_colour_change (t, x, mode, p):
    return t - colour_change_time
g_colour_change.terminal = True
g_colour_change.direction = 1 


## trigger when the cartesian position to the origin is 0
def g_reached_0_0(t, x, mode, p):
    pos_x, pos_y = get_xy(x, p)
    return np.sqrt(pos_x**2 + pos_y**2) - p.position_tol
g_reached_0_0.terminal = True
g_reached_0_0.direction = -1


#here we measure the distance traveled by the nozzle , so that it can do the wiped-process (to clean itsself from old filament)
#so again x[9] has the total traveled distance of the nozzle 
#when this distance reaches p.wipe_distance, then this guard is triggered
def g_filament_wiped(t, x, mode, p):
    return x[9] - p.wipe_distance 
g_filament_wiped.terminal = True
g_filament_wiped.direction = 1


def g_temp_low(t, x, mode, p):
    return x[8] - p.T_low# Cooled down to 30°C
g_temp_low.terminal = True
g_temp_low.direction = -1


def g_loading_done(t, x, mode, p):
    return t - (entry_time + p.loading_time) # 5s loading timer
g_loading_done.terminal = True
g_loading_done.direction = 1

def g_temp_high(t, x, mode, p):
    return x[8] - p.T_high # Heated up and above  235°C
g_temp_high.terminal = True
g_temp_high.direction = 1


#here we used the saved location of the head the moment the printer stopped printing
#and using the eucleidian eguation we compute the distance between the current location and saved location
#when this distance reaches almost 0, then the head will be positioned it stopped printed in the first place.
def g_reached_saved(t, x, mode, p):
    pos_x, pos_y = get_xy(x, p)
    return np.sqrt((pos_x - saved_x)**2 + (pos_y - saved_y)**2) - p.position_tol
g_reached_saved.terminal = True
g_reached_saved.direction = -1


p = PrinterParams()
t_curr, t_end = 0.0, 1000.0
theta_l0 = 0.145 / p.r_m
theta_r0 = 0.145 / p.r_m
direction  = 1.0
x_curr = np.array([theta_l0, theta_r0, 0.0, 0.0,0.0, 0.0, p.u_h / p.R_n, p.u_h, p.T_high, 0.0], dtype=float) # initial values for the state space parameters 
current_fsm = "q0"
entry_time = 0.0
colour_change_time = 250.0
colour_change_pending = True
saved_x, saved_y = 0.0, 0.0
ts_log, ys_log, state_log = [], [], []

print("Starting Hybrid Simulation Loop...")

#simulate hybrid system
while t_curr < t_end:
    
    # Instantaneous Jump for Saving Location (q1)
    if current_fsm == "q1":
        saved_x, saved_y = get_xy(x_curr, p)
        print(f"saved position {saved_x:.2f},{saved_y:.2f}")
        current_fsm = "q2"
        continue

   
    if current_fsm == "q0": # we assume that the printer is printing 
        if colour_change_pending:
            event =[g_filament_finished,g_colour_change]
        else :  
            event = g_filament_finished

    elif current_fsm == "q2":
        event = g_reached_0_0

    elif current_fsm == "q3":
        event = g_filament_wiped
        
    elif current_fsm == "q4":
        event = g_temp_low
    elif current_fsm == "q5":
        entry_time = t_curr
        event = g_loading_done
    elif current_fsm == "q6":
        event = g_temp_high
    elif current_fsm == "q7":
        event = g_reached_saved

    # we call solver 
    sol = solve_ivp(
        f,
        (t_curr, t_end),
        x_curr,
        args=(current_fsm, p),
        events=event
    )

    ts_log.append(sol.t)
    ys_log.append(sol.y)
    
    state_mapping = {"q0":0, "q1":1, "q2":2, "q3":3, "q4":4, "q5":5, "q6":6, "q7":7}
    state_log.extend([state_mapping[current_fsm]] * len(sol.t))

    t_curr = sol.t[-1]
    x_curr = sol.y[:, -1]

    #  jumps logic 

    # here the logic is that from the simulation in the while loop , when an event is triggered starting from q0,  
    # then the solver  has reached its target (ex. g_reaced_0_0.terminal = True) and then the integration stops 
    #so the status of the solver at this exact moment is 1, 
    #and then we jump to the next event and so on, until we hit the break and then the simulation ends.
    # Jumps logic 
    if sol.status == 1:
        print(f"[{t_curr:.2f}s] Event triggered. Leaving {current_fsm}")
        
        if current_fsm == "q0":
            if colour_change_pending:
                if len(sol.t_events[1]) > 0 :
                    colour_change_pending = False 
                    print(f"New colour requested")
            current_fsm = "q1"
            
        elif current_fsm == "q2":
            x_curr[2:6] = 0.0     # R(q2, q3) we zero out motor velocities and currents
            x_curr[9] = 0.0       # we reset odometer distance
            current_fsm = "q3"
            
        elif current_fsm == "q3":
            x_curr[2:6] = 0.0     # we zero out motor velocities and currents
            current_fsm = "q4"
            
        elif current_fsm == "q4":
            current_fsm = "q5"
            
        elif current_fsm == "q5":
            x_curr[9] = 0.0       # R(q5, q6) we reset odometer distance for the new filament
            current_fsm = "q6"
            
        elif current_fsm == "q6":
            pos_x, pos_y = get_xy(x_curr, p)
            direction = 1.0 if saved_x >= pos_x else -1.0 # complemetary line of code  to do a 2 - cycle simulation
            current_fsm = "q7"

            
        elif current_fsm == "q7":
            x_curr[2:6] = 0.0     # R(q7, q0) we zero out motor velocities and currents upon arrival
            print(f"[{t_curr:.2f}s] Cycle Complete! Returned to ({saved_x:.3f}, {saved_y:.3f}).")
            break


### Process and plot results ###
t_total = np.concatenate(ts_log)
y_total = np.concatenate(ys_log, axis=1)
x_pos, y_pos = get_xy(y_total, p)

plt.figure(figsize=(10, 8))

plt.subplot(3, 1, 1)
plt.plot(t_total, x_pos, label="movement in x")
plt.plot(t_total, y_pos, label="movement in y")
plt.ylabel("Position (m)")
plt.legend()
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(t_total, y_total[8])
plt.ylabel("Temperature (C)")
plt.grid()

plt.subplot(3, 1, 3)
plt.step(t_total, state_log, where="post")
plt.yticks(range(8), [f"q{i}" for i in range(8)])
plt.ylabel("State")
plt.xlabel("Time (s)")
plt.grid()

plt.tight_layout()
plt.show()