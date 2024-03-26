from supersonic_analysis import ConstraintAnalysis, DragPolar, MissionAnalysis, MissionStep
import ambiance as amb
import math


####################################################################################################################
# Aircraft variables & estimates in SI units, PART WHICH CAN BE CHANGED
filename_dragpolar = "my modified design from anastase v2_DegenGeom_V6.polar" #the filename of the polar, in the form name.polar
# S = 1.6235442959655093  # wing area in m^2
S = 0.37 
# b = 0.91  # wing span in m
b = 0.32
M = 0.1  # Mach number
V_stall = 30/3.6  # stall speed in m/s 
CL_max = 1.1  # lift coefficient CLMax (Land @ MLW)
CDR = 0  # Additional drag factors?
ROC = 2  # rate of climb in m/s
TR = 250  # turn radius in m
dv_dt = 5  # acceleration dv_dt = v_final - v_initial / delta_t TODO correct implementation / calculation
e = 0.8  # oswald efficiency factor

WP = 2  # payload weight in kg



####################################################################################################################
# W = 10  # weight in kg
#AR = 1.78733  # aspect ratio
# q = 0.5 * rho * V_cruise ** 2  # dynamic pressure in N/m^2
# k1 = 0.2129  # induced drag constant, k1
# k2 = 0.0053  # coefficient in lift-drag polar
# CD0 = 0.0059  # zero lift drag coefficient
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

# Constraint analysis to obtain optimal thrust-to-weight ratio and optimal wing loading
newCA = ConstraintAnalysis(S, M, V_stall, CL_max, ROC, TR, dv_dt, filename_dragpolar, alpha, beta, CDR, safety_margin_TW, safety_margin_WS, plot_max_x, plot_max_y)
# Thrust-to-weight ratio, Wing area-to-weight ratio
TWR, WSR = newCA.optimize()
print(newCA.load_factor())
newCA.plot()

####################################################################################################################
# Mission steps definition, PART WHICH CAN BE CHANGED
mission = MissionAnalysis(WP, CDR, CL_max, TWR, WSR,  S, e, filename_dragpolar)
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