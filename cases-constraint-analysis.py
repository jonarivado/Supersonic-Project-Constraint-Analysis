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
rho = float(amb.Atmosphere(0).density)  # air density in kg/m^3
mu = float(amb.Atmosphere(0).dynamic_viscosity)  # dynamic viscosity in kg/m/s
k = 1 / (math.pi * e * AR)  # induced drag constant, k1
k2 = 0  # coefficient in lift-drag polar # TODO can be set to zero?
CD0 = 0.005  # zero lift drag coefficient
CL = 1.2  # lift coefficient
CD = CD0 + k * CL ** 2  # drag coefficient
CDR = 0
g0 = float(amb.Atmosphere(0).grav_accel)  # gravitational acceleration in m/s^2
q = 0.5 * rho * V ** 2  # dynamic pressure in N/m^2
ROC = 10  # rate of climb in m/s
TR = 100  # turn radius in m
n = math.sqrt(1 + (V ** 2 / (g0 * TR)) ** 2)  # load factor
h_climb = 1000  # TODO should be climb_rate (dh/dt)
# climb_rate = (V * (T-D)) / (W * g0)  # climbing distance in m # TODO needs to be defined?
Vc = 1  # # TODO climb rate, climb velocity
theta = np.arctan(Vc / V)  # climb angle

alpha = 1  # thrust lapse factor
beta = 1  # how much the weight of the aircraft is reduced by the fuel burn compared to MTOW

# Case 1: Constant speed cruise at constant altitude (P_s = 0)
# simplified master equation
def TSL_WTO_M1(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S))
# Optimum points
# Minimum Weight to Area ratio WTO_S_1 (Wing loading)
WTO_S_1 = (q / beta) * (np.sqrt((CD0 + CDR) / k))
# Minimum Thrust to Weight ratio TSL_WTO_1 (Thrust loading)
TSL_WTO_1 = (beta / alpha) * (2 * np.sqrt(k * (CD0 + CDR)) + k2)

# Case 2: Constant speed climb (kinetic part of P_s = 0)
# simplified master equation
def TSL_WTO_M2(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + (1 / V) * h_climb)
# Optimum points
# Wing loading
WTO_S_2 = (q / beta) * (np.sqrt((CD0 + CDR) / k))
# Thrust loading
TSL_WTO_2 = (beta / alpha) * (2 * np.sqrt(k * (CD0 + CDR)) + k2 + (1 / V) * h_climb)

# Case 3: Constant altitude, speed turn (P_s = 0)
# simplified master equation
def TSL_WTO_M3(WTO_S):
    return (beta / alpha) * (k * n**2 * (beta / q) * WTO_S + k2 * n + (CD0 + CDR) / ((beta / q) * WTO_S) + (1 / V) * h_climb)
# Optimum points
# Wing loading
WTO_S_3 = (q / beta * n) * (np.sqrt((CD0 + CDR) / k))
# Thrust loading
TSL_WTO_3 = (beta * n / alpha) * (2 * np.sqrt(k * (CD0 + CDR)) + k2 + (1 / V) * h_climb)

# Case 4: Service Ceiling (max. altitude performance)
# same as for constant speed climb
# simplified master equation
def TSL_WTO_M4(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + (1 / V) * h_climb)

# Case 5: Constant altitude, horizontal acceleration
# simplified master equation
# TODO implement

# Case 6: Climb angle
# simplified master equation
def TSL_WTO_M6(WTO_S):
    return (beta / alpha) * (k * (beta / q) * WTO_S + k2 + (CD0 + CDR) / ((beta / q) * WTO_S) + np.sin(theta))
# Optimum points
# Wing loading
WTO_S_6 = (q / beta) * (np.sqrt((CD0 + CDR) / k))
# Thrust loading
TSL_WTO_6 = (beta / alpha) * (2 * np.sqrt(k * (CD0 + CDR)) + k2 + np.sin(theta))

# write each T/W and W/S to a dataframe and print it
data = {'Case': ['Constant speed cruise at constant altitude', 'Constant speed climb', 'Constant altitude, speed turn',
                 'Service Ceiling', 'Constant altitude, constant acceleration', 'Climb angle'],
        'T/W': [TSL_WTO_1, TSL_WTO_2, TSL_WTO_3, TSL_WTO_2, TSL_WTO_2, TSL_WTO_6],
        'W/S': [WTO_S_1, WTO_S_2, WTO_S_3, WTO_S_2, WTO_S_2, WTO_S_6]}
df = pd.DataFrame(data)
print(df)

# plot of thrust and wing loading # TODO constant altitude, speed turn out of range, need to tune values
wing_loading = [WTO_S_1, WTO_S_2, WTO_S_3, WTO_S_2, WTO_S_2, WTO_S_6]
thrust_loading = [TSL_WTO_1, TSL_WTO_2, TSL_WTO_3, TSL_WTO_2, TSL_WTO_2, TSL_WTO_6]

# plt.axis((950, 1000, -0.5, 15))
WTO_S = np.arange(100, 1000, 10)
# plot with only simplified master equations
"""
plt.grid(True)
plt.xlabel('W/S - Wing Loading [N/m^2]')
plt.ylabel('T/W - Thrust Loading [-]')
plt.plot(WTO_S, TSL_WTO_M1(WTO_S=WTO_S), '-bo', WTO_S, TSL_WTO_M2(WTO_S=WTO_S), '-ro', WTO_S, TSL_WTO_M3(WTO_S=WTO_S), '-go', WTO_S, TSL_WTO_M4(WTO_S=WTO_S), '-mo', WTO_S, TSL_WTO_M6(WTO_S=WTO_S), '-ko')
plt.show()
"""
# subplots with simplified master equations and optimum points
fig, axs = plt.subplots(2)
fig.suptitle('Constraint Plots and optimal points')
axs[0].grid(True)
axs[1].grid(True)
axs[0].set_xlabel('W/S - Wing Loading [N/m^2]')
axs[0].set_ylabel('T/W - Thrust Loading [-]')
axs[0].plot(WTO_S, TSL_WTO_M1(WTO_S=WTO_S), '-bo', WTO_S, TSL_WTO_M2(WTO_S=WTO_S), '-ro', WTO_S, TSL_WTO_M3(WTO_S=WTO_S), '-go', WTO_S, TSL_WTO_M4(WTO_S=WTO_S), '-mo', WTO_S, TSL_WTO_M6(WTO_S=WTO_S), '-ko')
axs[1].plot(wing_loading, thrust_loading, 'c*')
plt.show()

# use for naming of points
"""case_index = 0
for t, i in zip(wing_loading, thrust_loading):
    plt.annotate(str(case_index), (t, i), xytext=(t+1, i+1))
    case_index += 1"""