import numpy as np
import matplotlib.pyplot as plt
import ambiance as amb
import math
import pandas as pd

# Aircraft variables & estimates in SI units
W = 25  # weight in kg
S = 2  # wing area in m^2
b = 2  # wing span in m
c = 0.5  # mean aerodynamic chord in m
AR = b ** 2 / S  # aspect ratio
e = 0.8  # oswald efficiency factor
V = 100  # cruise speed in m/s
V_stall = 10  # stall speed in m/s
V_takeoff = 1.2 * V_stall  # take-off speed in m/s
rho = float(amb.Atmosphere(0).density)  # air density in kg/m^3
mu = float(amb.Atmosphere(0).dynamic_viscosity)  # dynamic viscosity in kg/m/s
k = 1 / (math.pi * e * AR)  # induced drag constant, k1
k2 = 0  # coefficient in lift-drag polar TODO can be set to zero?
CD0 = 0.005  # zero lift drag coefficient
CL = 1.2  # lift coefficient
CD = CD0 + k * CL ** 2  # drag coefficient
CDR = 0
g0 = float(amb.Atmosphere(0).grav_accel)  # gravitational acceleration in m/s^2
q = 0.5 * rho * V ** 2  # dynamic pressure in N/m^2
ROC = 10  # rate of climb in m/s
TR = 100  # turn radius in m
n = math.sqrt(1 + (V ** 2 / (g0 * TR)) ** 2)  # load factor
dv_dt = 10  # acceleration dv_dt = v_final - v_initial / delta_t TODO correct implementation / calculation

alpha = 1  # thrust lapse factor
beta = 1  # how much the weight of the aircraft is reduced by the fuel burn compared to MTOW

# ################################### Boeing 737-100 Validation ###################################
# Source: http://www.b737.org.uk/techspecsdetailed.htm, library: OpenAP, https://calculator.academy/aircraft-turn-radius-calculator/

# Aircraft variables & estimates in SI units
W = 44225  # weight in kg
S = 102  # wing area in m^2
b = 28.35  # wing span in m
AR = 8.83  # aspect ratio
e = 0.8  # oswald efficiency factor
V = 1000/3.6  # cruise speed in m/s
V_stall = 220/3.6  # stall speed in m/s
V_takeoff = 1.2 * V_stall  # take-off speed in m/s
rho = float(amb.Atmosphere(0).density)  # air density in kg/m^3
mu = float(amb.Atmosphere(0).dynamic_viscosity)  # dynamic viscosity in kg/m/s
k = 0.043  # induced drag constant, k1
k2 = 0  # coefficient in lift-drag polar TODO can be set to zero?
CD0 = 0.022  # zero lift drag coefficient
CL = 2.75  # lift coefficient CLMax (Land @ MLW)
CD = CD0 + k * CL ** 2  # drag coefficient
CDR = 0
g0 = float(amb.Atmosphere(0).grav_accel)  # gravitational acceleration in m/s^2
q = 0.5 * rho * V ** 2  # dynamic pressure in N/m^2
ROC = 12  # rate of climb in m/s
TR = 8000  # turn radius in m
n = math.sqrt(1 + (V ** 2 / (g0 * TR)) ** 2)  # load factor
dv_dt = 3  # acceleration dv_dt = v_final - v_initial / delta_t TODO correct implementation / calculation

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
print(WTO_S_stall)
# plot

#B737 val
#plot a point
plt.plot(4330, 0.343, 'o', color='black', label='B737-100')

plt.plot(WTO_S, TSL_WTO_TO(WTO_S=WTO_S), ls='--', color='m', label='Take-off condition')
plt.axvline(x=WTO_S_stall, ls='--', label='Stalling condition')
plt.plot(WTO_S, TSL_WTO_CRUISE(WTO_S=WTO_S), ls='--', color='b', label='Cruise condition')
plt.plot(WTO_S, TSL_WTO_CLIMB(WTO_S=WTO_S), ls='--', color='r', label='Climb condition')
plt.plot(WTO_S, TSL_WTO_TURN(WTO_S=WTO_S), ls='--', color='g', label='Turn condition')
plt.xlabel('Wing loading [N/m^2]')
plt.ylabel('Thrust loading [-]')
# plt.xlim(0, 2000)
plt.ylim(0, 1)
plt.legend()
plt.show()



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