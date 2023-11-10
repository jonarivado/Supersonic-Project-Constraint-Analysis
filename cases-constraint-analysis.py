import numpy as np
import matplotlib.pyplot as plt
import ambiance as amb
import math
import pandas as pd

# Aircraft variables & estimates in SI units
W = 10  # weight in kg
S = 0.46  # wing area in m^2
b = 0.91  # wing span in m
AR = 1.787  # aspect ratio
e = 0.8  # oswald efficiency factor
V = 330*0.85  # cruise speed in m/s
V_stall = 80/3.6  # stall speed in m/s
V_takeoff = 1.2 * V_stall  # take-off speed in m/s
rho = float(amb.Atmosphere(0).density)  # air density in kg/m^3
mu = float(amb.Atmosphere(0).dynamic_viscosity)  # dynamic viscosity in kg/m/s
k = 0.0043  # induced drag constant, k1
k2 = -0.022  # coefficient in lift-drag polar TODO can be set to zero?
CD0 = 0.033  # zero lift drag coefficient
CL = 0.76  # lift coefficient CLMax (Land @ MLW)
CD = CD0 + k * CL ** 2  # drag coefficient
CDR = 0 # Additional drag factors? TODO can be set to zero?
g0 = float(amb.Atmosphere(0).grav_accel)  # gravitational acceleration in m/s^2
q = 0.5 * rho * V ** 2  # dynamic pressure in N/m^2
ROC = 20  # rate of climb in m/s
TR = 1000  # turn radius in m
n = math.sqrt(1 + (V ** 2 / (g0 * TR)) ** 2)  # load factor
dv_dt = 10  # acceleration dv_dt = v_final - v_initial / delta_t TODO correct implementation / calculation

alpha = 1  # thrust lapse factor
beta = 1  # how much the weight of the aircraft is reduced by the fuel burn compared to MTOW

# W/S = 4330 N/m^2
# T/W = 0.343

# if ground distance must be calculated, use following
"""
T0 = 1
a = 1

A = g0 * ((T0 / W) - mu)  # T0: Total static thrust at V = 0
B = (g0 / W) * (0.5 * rho * S * (CD - mu * CL) + a)  # a: Constant
# take-off ground run distance
S = (1 / 2 * B) * np.log(A / (A - B * V_takeoff ** 2))
"""

# Stalling W/S STALL
WTO_S_stall = 0.5 * rho * V_stall ** 2 * CL

# Case 1: Straight and level flight, CRUISE
# simplified master equation
def TSL_WTO_CRUISE(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S))


# Case 2: CLIMB
# simplified master equation
def TSL_WTO_CLIMB(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + (1 / V) * ROC)

# Case: TAKE-OFF
def TSL_WTO_TO(WTO_S):
    return (beta / alpha) * ((CD + CDR - mu * CL) * (beta / q) * (WTO_S ** -1) + mu + (1 / g0) * dv_dt)

# Case: TURN
# simplified master equation
def TSL_WTO_TURN(WTO_S):
    return (beta / alpha) * (
                k * (n ** 2) * (beta / q) * WTO_S + k2 * n + ((CD0 + CDR) * q) / (beta * WTO_S))

WTO_S = np.linspace(0.1, 7000,2000)
# plot

#plot a point
# plt.plot(4330, 0.343, 'o', color='black', label='B737-100')

plt.plot(WTO_S, TSL_WTO_TO(WTO_S=WTO_S), ls='--', color='m', label='Take-off condition')
plt.axvline(x=WTO_S_stall, ls='--', label='Stalling condition')
plt.plot(WTO_S, TSL_WTO_CRUISE(WTO_S=WTO_S), ls='--', color='b', label='Cruise condition')
plt.plot(WTO_S, TSL_WTO_CLIMB(WTO_S=WTO_S), ls='--', color='r', label='Climb condition')
plt.plot(WTO_S, TSL_WTO_TURN(WTO_S=WTO_S), ls='--', color='g', label='Turn condition')
plt.xlabel('Wing loading [N/m^2]')
plt.ylabel('Thrust loading [-]')
plt.xlim(0, 500)
plt.ylim(0, 10)
plt.legend()
plt.show()

# T/W = 1.25, W/S = 310 N/m^2
# --> w must be 15.5kg and Engine must provide 200N


# other cases and optimum points
"""
# Case: Cruise
# Optimum points
# Minimum Weight to Area ratio WTO_S_1 (Wing loading)
WTO_S_1 = (q / beta) * np.sqrt((CD0 + CDR) / k)
# Minimum Thrust to Weight ratio TSL_WTO_1 (Thrust loading)
TSL_WTO_1 = (beta / alpha) * (2 * np.sqrt(k * (CD0 + CDR)) + k2)

# Case: Climb
# Optimum points
# Wing loading
WTO_S_2 = (q / beta) * np.sqrt((CD0 + CDR) / k)
# Thrust loading
TSL_WTO_2 = (beta / alpha) * (2 * np.sqrt(k * (CD0 + CDR)) + k2 + (1 / V) * dh_dt)

# Case 3: Turn radius
# Optimum points
# Wing loading
WTO_S_3 = (q / (beta * n)) * np.sqrt((CD0 + CDR) / k)
# Thrust loading
TSL_WTO_3 = (beta * n / alpha) * (2 * np.sqrt(k * (CD0 + CDR)) + k2)


# Case 4: Service Ceiling (max. altitude performance)
# same as for constant speed climb
# simplified master equation
def TSL_WTO_M4(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + (1 / V) * ROC)

# Optimum points
# Wing loading
WTO_S_4 = (q / beta) * np.sqrt((CD0 + CDR) / k)
# Thrust loading
TSL_WTO_4 = (beta / alpha) * (2 * np.sqrt(k * (CD0 + CDR)) + k2 + (1 / V) * ROC)


# Case 5: Constant altitude, horizontal acceleration
# simplified master equation
def TSL_WTO_M5(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + (1 / g0) * dv_dt)

# Optimum points
# Wing loading
WTO_S_5 = (q / (beta * n)) * np.sqrt((CD0 + CDR) / k)
# Thrust loading
TSL_WTO_5 = (beta * n / alpha) * (2 * np.sqrt(k * (CD0 + CDR)) + k2 + (1 / g0) * dv_dt)

# Case 6: Climb angle
# simplified master equation
def TSL_WTO_M6(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + np.sin(theta))


# Optimum points
# Wing loading
WTO_S_6 = (q / beta) * (np.sqrt((CD0 + CDR) / k))
# Thrust loading
TSL_WTO_6 = (beta / alpha) * (2 * np.sqrt(k * (CD0 + CDR)) + k2 + np.sin(theta))  # TODO if needed: theta = np.arctan(Vc / V)  # climb angle
"""