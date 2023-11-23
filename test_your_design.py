from supersonic_analysis import ConstraintAnalysis, DragPolar
import ambiance as amb
import math

#Drag Polar determination
dragpolar = DragPolar(filename="my modified design from anastase v2_DegenGeom_v3.polar")
DPcoeff = dragpolar.calculate_coeff()
print(DPcoeff)


# Aircraft variables & estimates in SI units
W = 10  # weight in kg
S = 0.7296  # wing area in m^2
b = 1.2  # wing span in m
AR = 1.96828  # aspect ratio
e = 0.8  # oswald efficiency factor
V = 330*0.85  # cruise speed in m/s
V_stall = 70/3.6  # stall speed in m/s
V_takeoff = 1.2 * V_stall  # take-off speed in m/s
rho = float(amb.Atmosphere(0).density)  # air density in kg/m^3
mu = float(amb.Atmosphere(0).dynamic_viscosity)  # dynamic viscosity in kg/m/s
k = DPcoeff['K1']  # induced drag constant, k1
k2 = DPcoeff['K2']  # coefficient in lift-drag polar TODO can be set to zero?
CD0 = DPcoeff['CD0']  # zero lift drag coefficient
CL = 1.1  # lift coefficient CLMax (Land @ MLW)
CD = CD0 + k * CL ** 2  # drag coefficient
CDR = 0 # Additional drag factors?
g0 = float(amb.Atmosphere(0).grav_accel)  # gravitational acceleration in m/s^2
q = 0.5 * rho * V ** 2  # dynamic pressure in N/m^2
ROC = 20  # rate of climb in m/s
TR = 300  # turn radius in m
n = math.sqrt(1 + (V ** 2 / (g0 * TR)) ** 2)  # load factor
dv_dt = 15  # acceleration dv_dt = v_final - v_initial / delta_t TODO correct implementation / calculation

alpha = 1  # thrust lapse factor
beta = 1  # how much the weight of the aircraft is reduced by the fuel burn compared to MTOW

safety_margin_TW = 0.1
safety_margin_WS = 10

plot_max_x = 500
plot_max_y = 4

newCA = ConstraintAnalysis(W=W, S=S, b=b, AR=AR, e=e, V=V, V_stall=V_stall, V_takeoff=V_takeoff, rho=rho, mu=mu, k=k, k2=k2, CD0=CD0, CL=CL, CD=CD, CDR=CDR, g0=g0, q=q, ROC=ROC, TR=TR, n=n, dv_dt=dv_dt, alpha=alpha, beta=beta,safety_margin_TW=safety_margin_TW,safety_margin_WS=safety_margin_WS,plot_max_x=plot_max_x,plot_max_y=plot_max_y)
newCA.optimize()
print(newCA.load_factor())
newCA.plot()
