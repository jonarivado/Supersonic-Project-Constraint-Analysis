from supersonic_analysis import ConstraintAnalysis, DragPolar, MissionAnalysis
import ambiance as amb
import math

# Aircraft variables & estimates in SI units
W = 10  # weight in kg
WP = 1  # payload weight in kg
S = 0.463318  # wing area in m^2
b = 0.91  # wing span in m
AR = 1.78733  # aspect ratio
e = 0.8  # oswald efficiency factor
M = 0.85  # Mach number
a = 330  # speed of sound in m/s
V = a*M  # cruise speed in m/s
V_stall = 50/3.6  # stall speed in m/s
V_takeoff = 1.2 * V_stall  # take-off speed in m/s
rho = float(amb.Atmosphere(0).density)  # air density in kg/m^3
mu = float(amb.Atmosphere(0).dynamic_viscosity)  # dynamic viscosity in kg/m/s
k = 0.2129  # induced drag constant, k1
k2 = 0.0053  # coefficient in lift-drag polar
CD0 = 0.0059  # zero lift drag coefficient
CL = 1.1  # lift coefficient CLMax (Land @ MLW)
CD = CD0 + k * CL ** 2  # drag coefficient
CDR = 0 # Additional drag factors?
g0 = float(amb.Atmosphere(0).grav_accel)  # gravitational acceleration in m/s^2
q = 0.5 * rho * V ** 2  # dynamic pressure in N/m^2
ROC = 200  # rate of climb in m/s
TR = 250  # turn radius in m
n = math.sqrt(1 + (V ** 2 / (g0 * TR)) ** 2)  # load factor
dv_dt = 20  # acceleration dv_dt = v_final - v_initial / delta_t TODO correct implementation / calculation
delta_h = 400  # mission altitude in m
N = 1  # number of turns, counted in video
theta = 0.9863  # mission angle at mission altitude on standard day conditions
R = 600  # mission range in m

alpha = 1  # thrust lapse factor
beta = 1  # how much the weight of the aircraft is reduced by the fuel burn compared to MTOW

safety_margin_TW = 0.1
safety_margin_WS = 10

plot_max_x = 500
plot_max_y = 4

dragpolar = DragPolar(filename="testpolar.polar")
DPcoeff = dragpolar.calculate_coeff()
print(DPcoeff)

newCA = ConstraintAnalysis(W, WP, S, b, AR, e, V, V_stall, V_takeoff, rho, mu, k, k2, CD0, CL, CD, CDR, g0, q, ROC, TR, n, dv_dt, alpha, beta,safety_margin_TW,safety_margin_WS,plot_max_x,plot_max_y)
newCA.optimize()
print(newCA.load_factor())
# newCA.plot()

newMA = MissionAnalysis(WP, a, M, q, CL, CD, CDR, alpha, beta, newCA.optimize(), delta_h, n, theta, N, V, g0, R, mu, V_takeoff)
W_fuel = newMA.TOTAL_FUEL_WR()
print(W_fuel)
We_WTO, WTO_calc = newMA.TOTAL_EMPTY_WR(W_guess=23, W_fuel=W_fuel)
# W_guess=newCA.optimize()[1]*S/g0
newMA.THRUST(WTO=WTO_calc)
newMA.WING_AREA(WTO=WTO_calc)
