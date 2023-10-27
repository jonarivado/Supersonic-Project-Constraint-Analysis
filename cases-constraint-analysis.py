# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ambiance as amb
import math

# Aircraft variables & estimates in SI units
W = 25  # weight in kg
S = 2  # wing area in m^2
b = 2  # wing span in m
c = 0.5  # mean aerodynamic chord in m
AR = b**2/S  # aspect ratio
e = 0.8  # oswald efficiency factor
V = 100  # cruise speed in m/s
V_stall = 10   # stall speed in m/s
rho = float(amb.Atmosphere(0).density)  # air density in kg/m^3
mu = float(amb.Atmosphere(0).dynamic_viscosity)  # dynamic viscosity in kg/m/s
k = 1/(math.pi*e*AR)  # induced drag constant, k1
k2 = 0  # coefficient in lift-drag polar # TODO can be set to zero?
CD0 = 0.005  # zero lift drag coefficient
CL = 1.2   # lift coefficient
CD = CD0 + k*CL**2  # drag coefficient
CDR = 0
g0 = float(amb.Atmosphere(0).grav_accel)  # gravitational acceleration in m/s^2
q = 0.5*rho*V**2  # dynamic pressure in N/m^2
ROC = 10  # rate of climb in m/s
TR = 100  # turn radius in m
n = math.sqrt(1+(V**2/(g0*TR))**2)  # load factor
h_climb = 1000  # TODO should be climb_rate (dh/dt)
# climb_rate = (V * (T-D)) / (W * g0)  # climbing distance in m # TODO needs to be defined?
Vc = 1  # # TODO climb rate, climb velocity
theta = np.arctan(Vc/V)  # climb angle

alpha = 1  # thrust lapse factor
beta = 1  # how much the weight of the aircraft is reduced by the fuel burn compared to MTOW

# Case 1: Constant speed cruise at constant altitude (P_s = 0)
# Minimum Weight to Area ratio WTO_S_1 (Wing loading)
WTO_S_1 = (q / beta) * (np.sqrt((CD0 + CDR)/k))
# Minimum Thrust to Weight ratio TSL_WTO_1 (Thrust loading)
TSL_WTO_1 = (beta / alpha) * (2*np.sqrt(k*(CD0 + CDR)) + k2)

# Case 2: Constant speed climb (kinetic part of P_s = 0)
# Wing loading
WTO_S_2 = (q / beta) * (np.sqrt((CD0 + CDR)/k))
# Thrust loading
TSL_WTO_2 = (beta / alpha) * (2*np.sqrt(k*(CD0 + CDR)) + k2 + (1 / V) * h_climb)

# Case 3: Constant altitude, speed turn (P_s = 0)
# Wing loading
WTO_S_3 = (q / beta * n) * (np.sqrt((CD0 + CDR)/k))
# Thrust loading
TSL_WTO_3 = (beta * n / alpha) * (2*np.sqrt(k*(CD0 + CDR)) + k2 + (1 / V) * h_climb)

# Case 4: Service Ceiling (max. altitude performance)
# same as for constant speed climb

# Case 5: Constant altitude, constant acceleration

# Case 6: Climb angle
# Wing loading
WTO_S_6 = (q / beta) * (np.sqrt((CD0 + CDR)/k))
# Thrust loading
TSL_WTO_6 = (beta / alpha) * (2*np.sqrt(k*(CD0 + CDR)) + k2 + np.sin(theta))