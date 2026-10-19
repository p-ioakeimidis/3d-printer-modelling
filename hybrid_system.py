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
    T_e: float = 20.0    #  Environment temperature (Celsius)
    u_h: float = 15.0    # Heater voltage applied when ON (Volts)
    
###physical modelling###
def f(t, x, u, q, p):
    theta_l = x[0] 
    theta_r = x[1]
    dx_1 = x[2] # Motor velocities
    dx_2 = x[3] # Motor velocities
    i_l = x[4]# Motor currents
    i_r= x[5]   # Motor currents
    i_L = x[6]
    u_C = x[7]
    T_n =  x[8] # Nozzle RLC & thermal states

    u_l, u_r = u[0], u[1]
    u_n = p.u_h if q == 1 else 0.0 

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
    
    if q == 1 and (abs(u_l) > 0 or abs(u_r) > 0):
        dx_10 = np.sqrt(v_x**2 + v_y**2)
    else:
        dx_10 = 0.0

    return np.array([dx_1, dx_2, dx_3, dx_4, dx_5, dx_6, dx_7, dx_8, dx_9, dx_10])

# ==========================================
# 2. DISCRETE EVENT GUARDS
# ==========================================

#solve_ivp integrates dx_7 (velocity) starting from 0.0 , we take x[7] position/distance and when the value is equal to MAX_FILAMENT, then this guard is triggered
#the solver updates the value of x[7] in every step , actually x[7] is a discrete value , while dx_7 is a continiues one (dx/ dt)
def g_filament_finished(t, x, u, q, p):
    return x[9] - 5.0 # 
g_filament_finished.terminal = True; g_filament_finished.direction = 1

## trigger when the cartesian position to the origin is 0
def g_reached_0_0(t, x, u, q, p):
    pos_x = (p.r_m / 2.0) * (x[0] + x[1])
    pos_y = (p.r_m / 2.0) * (x[0] - x[1])
    return np.sqrt(pos_x**2 + pos_y**2) - 0.001
g_reached_0_0.terminal = True; g_reached_0_0.direction = -1


#here we measure the distance traveled by the nozzle , so that it can do the wiped-process (to clean itsself from old filament)
#so again x[9] has the total traveled distance of the nozzle 
#when this distance reaches WIPED_THRESHOLD, then this guard is triggered
def g_filament_wiped(t, x, u, q, p):
    return x[9] - 0.05 # Traveled 5 cm of scrub distance
g_filament_wiped.terminal = True; g_filament_wiped.direction = 1


def g_temp_low(t, x, u, q, p):
    return x[8] - 30.0 # Cooled down to 30°C
g_temp_low.terminal = True; g_temp_low.direction = -1


#last_state_t is the exact timestamp the temp_low state finished
#t is the the master clock of the simulation, in solve_ivp
# when we enter in q5 we capture the current time (t_curr)
def g_loading_done(t, x, u, q, p):
    entry_time = u[2]
    return t - (entry_time + 5.0) # 5s loading timer
g_loading_done.terminal = True; g_loading_done.direction = 1

def g_temp_high(t, x, u, q, p):
    return x[8] - 200.0 # Heated up and above  200°C
g_temp_high.terminal = True; g_temp_high.direction = 1


#here we used the saved location of the head the moment the printer stopped printing
#and using the eucleidian eguation we compute the distance between the current location and saved location
#when this distance reaches almost 0, then the head will be where it stopped printed in the first place.


def g_reached_saved(t, x, u, q, p):
    pos_x = (p.r_m / 2.0) * (x[0] + x[1])
    pos_y = (p.r_m / 2.0) * (x[0] - x[1])
    saved_x, saved_y = u[2], u[3]
    return np.sqrt((pos_x - saved_x)**2 + (pos_y - saved_y)**2) - 0.001
g_reached_saved.terminal = True; g_reached_saved.direction = -1


p = PrinterParams()
t_curr, t_end = 0.0, 300.0
x_curr = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, p.T_e, 0.0]) # initial values for the state space parameters 
current_fsm = "q0"

saved_x, saved_y = 0.0, 0.0
ts_log, ys_log, state_log = [], [], []

print("Starting Hybrid Simulation Loop...")

#simulate hybrid system
while t_curr < t_end:
    
    # Instantaneous Jump for Saving Location (q1)
    if current_fsm == "q1":
        saved_x = (p.r_m / 2.0) * (x_curr[0] + x_curr[1])
        saved_y = (p.r_m / 2.0) * (x_curr[0] - x_curr[1])
        print(f"[{t_curr:.2f}s] q1: Saved location ({saved_x:.3f}, {saved_y:.3f})")
        current_fsm = "q2"
        continue

   
    if current_fsm == "q0": # we assume that the printer is printing 
        u_active = (2.5, 2.5)  # moving diagonally
        q, active_events = 1, [g_filament_finished] # we only care for this guard , in order to move to the next state q1 
        
    elif current_fsm == "q2":
        u_active = (-2.5, -2.5)  #reverse voltage back to (0,0) (as it initialy started at 0,0 before we gave voltage)
        q, active_events = 1, [g_reached_0_0]
        
    elif current_fsm == "q3":
        u_active = (1.5, -1.5) # Wiggle movement at origin
        q, active_events = 1, [g_filament_wiped]
        
    elif current_fsm == "q4":
        u_active = (0.0, 0.0) # Motors off, heater off
        q, active_events = 0, [g_temp_low]
        
    elif current_fsm == "q5":
        u_active = (0.0, 0.0, t_curr) # Motors off, pass entry time
        q, active_events = 0, [g_loading_done]
        
    elif current_fsm == "q6":
        u_active = (0.0, 0.0) # Motors off, heater on
        q, active_events = 1, [g_temp_high]
        
    elif current_fsm == "q7":
        u_active = (2.5, 2.5, saved_x, saved_y) # constant positice voltage to go back to saved point
        q, active_events = 1, [g_reached_saved]

    # we call solver 
    sol = solve_ivp(f, (t_curr, t_end), x_curr, args=(u_active, q, p),
                    events=active_events, dense_output=True, max_step=0.05)

    ts_log.append(sol.t)
    ys_log.append(sol.y)
    
    state_mapping = {"q0":0, "q1":1, "q2":2, "q3":3, "q4":4, "q5":5, "q6":6, "q7":7}
    state_log.extend([state_mapping[current_fsm]] * len(sol.t))

    t_curr = sol.t[-1]
    x_curr = sol.y[:, -1]

    #  jumps logic 
    if sol.status == 1:
        print(f"[{t_curr:.2f}s] Event triggered. Leaving {current_fsm}")
        if current_fsm == "q0": 
            current_fsm = "q1"
        elif current_fsm == "q2":
            x_curr[2:6], x_curr[9] = 0.0, 0.0 # Reset motor speeds and odometer distance
            current_fsm = "q3"
        elif current_fsm == "q3":
            x_curr[0], x_curr[1] = 0.0, 0.0   # Align exact origin
            x_curr[2:6] = 0.0
            current_fsm = "q4"
        elif current_fsm == "q4": 
            current_fsm = "q5"
        elif current_fsm == "q5":
            x_curr[9] = 0.0 # Reset odometer distance for new filament s
            current_fsm = "q6"
        elif current_fsm == "q6": 
            current_fsm = "q7"
        elif current_fsm == "q7":
            x_curr[2:6] = 0.0
            print(f"[{t_curr:.2f}s] Cycle Complete! Returned to ({saved_x:.3f}, {saved_y:.3f}).")
            break

t_total = np.concatenate(ts_log)
y_total = np.concatenate(ys_log, axis=1)
state_total = np.array(state_log)

x_pos = (p.r_m / 2.0) * (y_total[0, :] + y_total[1, :])
y_pos = (p.r_m / 2.0) * (y_total[0, :] - y_total[1, :])
temp_tot = y_total[8, :]

fig = plt.figure(figsize=(14, 9))

# Plot 1: Cartesian Trajectory
ax1 = plt.subplot(2, 2, 1)
ax1.plot(x_pos, y_pos, label="Tool Head Path", color="blue", lw=2)
ax1.scatter([0], [0], color="red", s=60, zorder=5, label="Wipe Tower (0,0)")
ax1.scatter([saved_x], [saved_y], color="green", s=60, zorder=5, label=f"Saved Pos ({saved_x:.2f}, {saved_y:.2f})")
ax1.set_title("Cartesian Trajectory", fontsize=11, fontweight='bold')
ax1.set_xlabel("X Position (m)"); ax1.set_ylabel("Y Position (m)")
ax1.legend(); ax1.grid(True)

# Plot 2: Nozzle Temperature
ax2 = plt.subplot(2, 2, 2)
ax2.plot(t_total, temp_tot, color="orange", linewidth=2.5)
ax2.axhline(y=30, color='b', linestyle='--', alpha=0.5, label="heater off (30°C)")
ax2.axhline(y=200, color='r', linestyle='--', alpha=0.5, label="heater on (200°C)")
ax2.set_title("Nozzle Temperature", fontsize=11, fontweight='bold')
ax2.set_xlabel("Time (s)"); ax2.set_ylabel("Temperature (°C)")
ax2.legend(); ax2.grid(True)

# Plot 3: Discrete State Timeline
ax3 = plt.subplot(2, 1, 2)
ax3.step(t_total, state_total, where='post', color='purple', linewidth=2.5)
ax3.set_yticks(range(8))
ax3.set_yticklabels([f"q{i}" for i in range(8)])
ax3.set_title("Hybrid Automaton State Transitions vs Time", fontsize=11, fontweight='bold')
ax3.set_xlabel("Time (s)"); ax3.set_ylabel("Active Discrete Mode")
ax3.grid(True, linestyle=':', alpha=0.7)
plt.tight_layout()
plt.show()