import numpy as np
import matplotlib.pyplot as plt
import ambiance as amb
import math
import pandas as pd
import scipy.optimize as opt

# Aircraft variables & estimates in SI units
W = 10  # weight in kg
S = 0.5  # wing area in m^2
b = 1.5  # wing span in m
AR = 4.54  # aspect ratio
e = 0.8  # oswald efficiency factor
V = 330*0.85  # cruise speed in m/s
V_stall = 50/3.6  # stall speed in m/s
V_takeoff = 1.2 * V_stall  # take-off speed in m/s
rho = float(amb.Atmosphere(0).density)  # air density in kg/m^3
mu = float(amb.Atmosphere(0).dynamic_viscosity)  # dynamic viscosity in kg/m/s
k = 0.0946  # induced drag constant, k1
k2 = 0  # coefficient in lift-drag polar TODO can be set to zero?
CD0 = 0.00728  # zero lift drag coefficient
CL = 3  # lift coefficient CLMax (Land @ MLW)
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

safety_margin_TW = 0.1
safety_margin_WS = 10

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

# # Stalling W/S STALL
# WTO_S_stall = 0.5 * rho * V_stall ** 2 * CL

# # Case 1: Straight and level flight, CRUISE
# # simplified master equation
# def TSL_WTO_CRUISE(WTO_S):
#     return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S))


# # Case 2: CLIMB
# # simplified master equation
# def TSL_WTO_CLIMB(WTO_S):
#     return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + (1 / V) * ROC)

# # Case: TAKE-OFF
# def TSL_WTO_TO(WTO_S):
#     return (beta / alpha) * ((CD + CDR - mu * CL) * (beta / q) * (WTO_S ** -1) + mu + (1 / g0) * dv_dt)

# # Case: TURN
# # simplified master equation
# def TSL_WTO_TURN(WTO_S):
#     return (beta / alpha) * (
#                 k * (n ** 2) * (beta / q) * WTO_S + k2 * n + ((CD0 + CDR) * q) / (beta * WTO_S))

# #Optimization Algo

# safety_margin_TW = 0.1
# safety_margin_WS = 10

# # T/W = x[0] and W/S = x[1]
# constraints = [{'type': 'ineq', 'fun':  lambda x: x[0] - TSL_WTO_CRUISE(x[1])-safety_margin_TW},
#                 {'type': 'ineq', 'fun':  lambda x: x[0] - TSL_WTO_CLIMB(x[1])-safety_margin_TW},
#                 {'type': 'ineq', 'fun':  lambda x: x[0] - TSL_WTO_TO(x[1])-safety_margin_TW},
#                 {'type': 'ineq', 'fun':  lambda x: x[0] - TSL_WTO_TURN(x[1])-safety_margin_TW},
#                 {'type': 'eq', 'fun': lambda x: x[1] - WTO_S_stall + safety_margin_WS}]

# # initial guess
# x0 = [0.5, 100]
# # optimization
# res = opt.minimize(lambda x: x[0], x0, constraints=constraints,bounds=((0, 10), (0, 1000)))
# print("Opimal T/W: ", round(res.x[0],3))
# print("Optimal W/S: ", round(res.x[1],3))

# WTO_S = np.linspace(0.1, 7000,2000)
# # plot

# #plot the optimization result
# plt.plot(res.x[1], res.x[0], 'o', color='r', label='Optimum point')

# plt.plot(WTO_S, TSL_WTO_TO(WTO_S=WTO_S), ls='--', color='m', label='Take-off condition')
# plt.axvline(x=WTO_S_stall, ls='--', label='Stalling condition')
# plt.plot(WTO_S, TSL_WTO_CRUISE(WTO_S=WTO_S), ls='--', color='b', label='Cruise condition')
# plt.plot(WTO_S, TSL_WTO_CLIMB(WTO_S=WTO_S), ls='--', color='r', label='Climb condition')
# plt.plot(WTO_S, TSL_WTO_TURN(WTO_S=WTO_S), ls='--', color='g', label='Turn condition')
# plt.xlabel('Wing loading [N/m^2]')
# plt.ylabel('Thrust loading [-]')
# plt.xlim(0, 500)
# plt.ylim(0, 2)
# plt.legend()
# plt.title('Optima w/ margin: T/W: ' +str(round(res.x[1],3))+ ', W/S: '+str(round(res.x[0],3)))
# plt.show()
    
# other cases and optimum points
"""

# Case 4: Service Ceiling (max. altitude performance)
# same as for constant speed climb
# simplified master equation
def TSL_WTO_M4(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + (1 / V) * ROC)


# Case 5: Constant altitude, horizontal acceleration
# simplified master equation
def TSL_WTO_M5(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + (1 / g0) * dv_dt)


# Case 6: Climb angle
# simplified master equation
def TSL_WTO_M6(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + np.sin(theta))


"""