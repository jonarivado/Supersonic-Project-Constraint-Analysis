from supersonic_analysis import ConstraintAnalysis
import ambiance as amb
import math

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


newCA = ConstraintAnalysis(W, S, b, AR, e, V, V_stall, V_takeoff, rho, mu, k, k2, CD0, CL, CD, CDR, g0, q, ROC, TR, n, dv_dt, alpha, beta,safety_margin_TW,safety_margin_WS)
newCA.optimize()
newCA.plot()