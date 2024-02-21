from supersonic_analysis import ConstraintAnalysis, DragPolar, MissionAnalysis, MissionStep
import ambiance as amb
import math

####################################################################################################################
# Aircraft variables & estimates in SI units, PART WHICH CAN BE CHANGED
W = 10  # weight in kg
WP = 10  # payload weight in lb
S = 0.463318  # wing area in m^2
b = 0.91  # wing span in m
AR = 1.78733  # aspect ratio
e = 0.8  # oswald efficiency factor
M = 0.85  # Mach number
a = 330  # speed of sound in m/s
V_cruise = a*M  # cruise speed in m/s
V_stall = 50/3.6  # stall speed in m/s
V_takeoff = 1.2 * V_stall  # take-off speed in m/s
rho = float(amb.Atmosphere(0).density)  # air density in kg/m^3
mu = float(amb.Atmosphere(0).dynamic_viscosity)  # dynamic viscosity in kg/m/s
k = 0.2129  # induced drag constant, k1
k2 = 0.0053  # coefficient in lift-drag polar
CD0 = 0.0059  # zero lift drag coefficient
CL = 1.1  # lift coefficient CLMax (Land @ MLW)
CD = CD0 + k * CL ** 2  # drag coefficient
CDR = 0  # Additional drag factors?
g0 = float(amb.Atmosphere(0).grav_accel)  # gravitational acceleration in m/s^2
# q = 0.5 * rho * V_cruise ** 2  # dynamic pressure in N/m^2
ROC = 200  # rate of climb in m/s
TR = 250  # turn radius in m
n = math.sqrt(1 + (V_cruise ** 2 / (g0 * TR)) ** 2)  # load factor
dv_dt = 20  # acceleration dv_dt = v_final - v_initial / delta_t TODO correct implementation / calculation
# delta_h = 400  # mission altitude in m
# N = 1  # number of turns, counted in video
# theta = 0.9863  # mission angle at mission altitude on standard day conditions
# R = 1000  # mission range in m
####################################################################################################################

alpha = 1  # thrust lapse factor
beta = 1  # how much the weight of the aircraft is reduced by the fuel burn compared to MTOW

safety_margin_TW = 0.1
safety_margin_WS = 10

plot_max_x = 500
plot_max_y = 4

dragpolar = DragPolar(filename="testpolar.polar")
DPcoeff = dragpolar.calculate_coeff()
print(DPcoeff)

# Constraint analysis to obtain optimal thrust-to-weight ratio and optimal wing loading
newCA = ConstraintAnalysis(W, WP, S, b, AR, e, V_cruise, V_stall, V_takeoff, rho, mu, k, k2, CD0, CL, CD, CDR, g0, ROC, TR, n, dv_dt, alpha, beta, safety_margin_TW, safety_margin_WS, plot_max_x, plot_max_y)
# Thrust-to-weight ratio, Wing area-to-weight ratio
TWR, WSR = newCA.optimize()
print(newCA.load_factor())
# newCA.plot()

####################################################################################################################
# Mission steps definition, PART WHICH CAN BE CHANGED
mission = MissionAnalysis(WP, CD, CDR, CD0, mu, CL, TWR, WSR, g0, rho, a, S, e)
mission.TAKEOFF()
mission.CLIMB(M=0.3)
mission.TURN(M=0.7, theta=0.7, N=25, V_turn=120, TR=250)  # value for theta from Mattingly, depends on altitude
mission.CRUISE(M=0.85, theta=0.6, R=100, V_cruise=250)  # range R in km
mission.CLIMB(M=0.85)
mission.ACCELERATE(M=0.8)
mission.ACCELERATE(M=2.0)
mission.LANDING()
####################################################################################################################
# Mission analysis, details collected for printing, DO NOT CHANGE THIS PART
mission.analyze()
result = mission.create_dataframe()
selected_columns = result[['Mach number', 'Weight ratio', 'Total weight ratio']]
# print(selected_columns.columns)
selected_columns_indexed = selected_columns.reset_index()
print("-" * 80)
# print(selected_columns_indexed.columns)
selected_columns_indexed.rename(columns={'index': 'Missionstep'}, inplace=True)
# selected_columns_indexed.index.name = 'Step'
print(selected_columns_indexed)
print("-" * 80)
fuel_ratio = mission.TOTAL_FUEL_WR()
# ewf = mission.calculate_ewf(w_guess=55)
# print("Empty weight ratio: " + str(ewf))
takeoff_weight, final_ewf, final_ew, nr_iteration, fuel_weight = mission.TOTAL_WR(initial_w_guess=55, w_fuel=fuel_ratio)
print("Number of iterations: " + str(nr_iteration))
print("Final empty weight ratio: " + str(final_ewf))
print("Final empty weight: " + str(final_ew) + " lb")
print("Fuel ratio: " + str(fuel_ratio))
print("Final fuel weight: " + str(fuel_weight) + " lb")
print("-" * 80)
print("Final takeoff weight: " + str(takeoff_weight) + " lb")
thrust = mission.THRUST(TWR=TWR, WTO=mission.convert_lb_to_kg(input=takeoff_weight))
wing_area = mission.WING_AREA(WSR=WSR, WTO=mission.convert_lb_to_kg(input=takeoff_weight))
print("Thrust: " + str(thrust) + " N")
print("Wing area: " + str(wing_area) + " m^2")
print("-" * 80)


