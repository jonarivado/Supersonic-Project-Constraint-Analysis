import numpy as np
import matplotlib.pyplot as plt
import ambiance as amb
import math
import pandas as pd
import scipy.optimize as opt

class ConstraintAnalysis:
    def __init__(self, W, WP, S, b, AR, e, V_cruise, V_stall, V_takeoff, rho, mu, k, k2, CD0, CL, CD, CDR, g0, ROC, TR, n, dv_dt, alpha, beta,safety_margin_TW=0,safety_margin_WS=0,plot_max_x=500,plot_max_y=2):
        self.W = W
        self.WP = WP
        self.S = S
        self.b = b
        self.AR = AR
        self.e = e
        self.V_cruise = V_cruise
        self.V_stall = V_stall
        self.V_takeoff = V_takeoff
        self.rho = rho
        self.mu = mu
        self.k = k
        self.k2 = k2
        self.CD0 = CD0
        self.CL = CL
        self.CD = CD
        self.CDR = CDR
        self.g0 = g0
        self.ROC = ROC
        self.TR = TR
        self.n = n
        self.dv_dt = dv_dt
        self.alpha = alpha
        self.beta = beta
        self.res = [0, 0]
        self.safety_margin_TW = safety_margin_TW
        self.safety_margin_WS = safety_margin_WS
        self.plot_max_x = plot_max_x
        self.plot_max_y = plot_max_y

    def TSL_WTO_CRUISE(self, WTO_S):
        q = 0.5 * self.rho * (self.V_cruise ** 2)
        return (self.beta / self.alpha) * (self.k * (self.beta / q) * WTO_S + self.k2 + (self.CD0 + self.CDR) / ((self.beta / q) * WTO_S))

    def TSL_WTO_CLIMB(self, WTO_S):
        q = 0.5 * self.rho * (self.V_cruise ** 2)  # climbing velocity the same as cruising velocity TODO ???
        return (self.beta / self.alpha) * (self.k * (self.beta / q) * WTO_S + self.k2 + (self.CD0 + self.CDR) / ((self.beta / q) * WTO_S) + (1 / self.V_cruise) * self.ROC)

    def TSL_WTO_TO(self, WTO_S):
        q = 0.5 * self.rho * (self.V_takeoff ** 2)
        return (self.beta / self.alpha) * ((self.CD + self.CDR - self.mu * self.CL) * (self.beta / q) * (WTO_S ** -1) + self.mu + (1 / self.g0) * self.dv_dt)

    def TSL_WTO_TURN(self, WTO_S):
        q = 0.5 * self.rho * (self.V_cruise ** 2)  # turn at cruise speed TODO ???
        return (self.beta / self.alpha) * (
                    self.k * (self.n ** 2) * (self.beta / q) * WTO_S + self.k2 * self.n + ((self.CD0 + self.CDR) * q) / (self.beta * WTO_S))
    
    def TSL_WTO_STALL(self):
        return 0.5 * self.rho * self.V_stall ** 2 * self.CL


    def optimize(self):
        constraints = [{'type': 'ineq', 'fun':  lambda x: x[0] - self.TSL_WTO_CRUISE(x[1])-self.safety_margin_TW},
                {'type': 'ineq', 'fun':  lambda x: x[0] - self.TSL_WTO_CLIMB(x[1])-self.safety_margin_TW},
                {'type': 'ineq', 'fun':  lambda x: x[0] - self.TSL_WTO_TO(x[1])-self.safety_margin_TW},
                {'type': 'ineq', 'fun':  lambda x: x[0] - self.TSL_WTO_TURN(x[1])-self.safety_margin_TW},
                {'type': 'eq', 'fun': lambda x: x[1] - self.TSL_WTO_STALL() + self.safety_margin_WS}]

        # initial guess
        x0 = [0.5, 100]
        # optimization
        self.res = opt.minimize(lambda x: x[0], x0, constraints=constraints, bounds=((0, 10), (0, 1000)))
        print("Optimal T/W: ", round(self.res.x[0], 3))
        print("Optimal W/S: ", round(self.res.x[1], 3))
        print("Required Thrust: ", round(self.res.x[0]*round(self.res.x[1]*self.S, 3), 3), " N")
        print("Optimal W: ", round(self.res.x[1]*self.S/self.g0, 3), " kg")

        return self.res.x[0], self.res.x[1]
    
    def plot(self):
        
        #assert that the optimize function has been run, otherwise throw an error
        assert hasattr(self, 'res') and hasattr(self.res, 'x'), "Optimization not performed. Please run optimize() first."

        WTO_S = np.linspace(0.1, 7000, 2000)
        # plot

        #plot the optimization result
        plt.plot(self.res.x[1], self.res.x[0], 'o', color='r', label='Optimum point')

        plt.plot(WTO_S, self.TSL_WTO_TO(WTO_S=WTO_S), ls='--', color='m', label='Take-off condition')
        plt.axvline(x=self.TSL_WTO_STALL(), ls='--', label='Stalling condition')
        plt.plot(WTO_S, self.TSL_WTO_CRUISE(WTO_S=WTO_S), ls='--', color='b', label='Cruise condition')
        plt.plot(WTO_S, self.TSL_WTO_CLIMB(WTO_S=WTO_S), ls='--', color='r', label='Climb condition')
        plt.plot(WTO_S, self.TSL_WTO_TURN(WTO_S=WTO_S), ls='--', color='g', label='Turn condition')
        plt.xlabel('Wing loading [N/m^2]')
        plt.ylabel('Thrust loading [-]')
        plt.xlim(0, self.plot_max_x)
        plt.ylim(0, self.plot_max_y)
        plt.legend()
        plt.title('Optima w/ margin: T/W: ' +str(round(self.res.x[0],3))+ ', W/S: '+str(round(self.res.x[1],3))+', W: '+str(round(self.res.x[1]*self.S/self.g0,3))+' kg')
        plt.show()

    def load_factor(self):
        return round(math.sqrt(1 + (self.V_cruise ** 2 / (self.g0 * self.TR)) ** 2),3)

class MissionStep:
    def __init__(self, step_type, details):
        self.step_type = step_type
        self.details = details

    def display_details(self):
        print(f"Mission part type: {self.step_type}")
        print("Details:")
        for key, value in self.details.items():
            print(f"{key}: {value}")
        print()

class MissionAnalysis:
    def __init__(self, WP, CD, CDR, CD0, mu, CL, TWR, WSR, g0, rho, a, S, e):
        self.mission_steps = []
        self.WP = WP  # payload weight
        self.CD = CD
        self.CDR = CDR
        self.CD0 = CD0
        self.mu = mu
        self.CL = CL
        self.TWR = TWR  # thrust to weight ratio
        self.WSR = WSR  # weight to wing area ratio
        self.g0 = g0
        self.rho = rho
        self.a = a
        self.S = S
        self.e = e
        self.total_WR = 1.0

    def add(self, step):
        self.mission_steps.append(step)

    def create_dataframe(self):
        data = {}
        for step in self.mission_steps:
            data[step.step_type] = step.details

        df = pd.DataFrame(data)
        pd.set_option('display.max_columns', None)  # Show all columns
        pd.set_option('display.max_rows', None)  # Show all rows
        return df.transpose()

    # Turbojet engine max power, eq. 3.55b, p. 71, Mattingly
    def TSFC(self, M, theta):
        return (1.5 + 0.23 * M) * np.sqrt(theta)

    # Equation 3.21, 3.22, p. 63, Mattingly
    def TAKEOFF(self):  # M=None, theta=None, V_takeoff=None, alpha=None, beta=None
        W_takeoff = 0.97  # in range 0.97 - 0.99
        self.total_WR = self.total_WR * W_takeoff
        step_details = {
            "Weight ratio": W_takeoff, "Total weight ratio": self.total_WR}  # "Theta": theta, "V_takeoff": V_takeoff, "Alpha": alpha, "Beta": beta
        step = MissionStep("Takeoff", step_details)
        self.mission_steps.append(step)
        # arguments = [M, theta, V_takeoff, alpha, beta]
        """if all(i is not None for i in arguments):
            C = self.TSFC(M=M, theta=theta) # [1/h]
            q = 0.5 * self.rho * V_takeoff**2
            xi = self.CD + self.CDR - self.mu * self.CL
            u = (xi * (q * beta) * ((self.WSR)**-1) + self.mu) * (beta / alpha) * self.TWR
            W_takeoff = np.exp(-C * np.sqrt(theta) / self.g0 * (V_takeoff / (1 - u)))
            print("W_takeoff: " + str(W_takeoff))"""
        # else:
            # raise ValueError("Error: Please provide corresponding values.")
        return W_takeoff

    # Equation 3.20, p 62, Mattingly
    def CLIMB(self, M=None):
        W_climb = None
        arguments = [M]  # theta, V_climb, delta_h, alpha, beta
        if all(i is not None for i in arguments):
            """theta = 0.9931
            C = self.TSFC(M=M, theta=theta)  # [1/h]
            u = ((self.CD + self.CDR) / self.CL) * (beta / alpha) * (self.TWR**-1)
            W_climb = np.exp((-C / self.a) * ((delta_h + (V_climb**2) / (2*self.g0)) / (1 - u)))
            print("W_climb: " + str(W_climb))"""
            if M > 0.2:
                W_climb = 1.0065 - 0.0325 * M  # eq. 6.9 from Raymer
                self.total_WR = self.total_WR * W_climb
                print("W_climb: " + str(W_climb))
                step_details = {"Mach number": M, "Weight ratio": W_climb, "Total weight ratio": self.total_WR}  # "Theta": theta, "V_climb": V_climb, "Mission altitude": delta_h, "Alpha": alpha, "Beta": beta
                step = MissionStep("Climb", step_details)
                self.mission_steps.append(step)
            else:
                raise ValueError("Error: Please provide larger Mach number")
        else:
            raise ValueError("Error: Please provide corresponding values.")
        return W_climb

    def ACCELERATE(self, M=None):
        W_accelerate = None
        arguments = [M]
        if all(i is not None for i in arguments):
            if M >= 1.0:  # supersonic
                W_accelerate = 0.991 - 0.007 * M - 0.01 * M**2
                self.total_WR = self.total_WR * W_accelerate
                print("W_accelerate (supersonic): " + str(W_accelerate))
                step_details = {"Mach number": M, "Weight ratio": W_accelerate, "Total weight ratio": self.total_WR}
                step = MissionStep("Accelerate at supersonic", step_details)
                self.mission_steps.append(step)
            else:  # subsonic
                W_accelerate = 1.0065 - 0.0325 * M
                self.total_WR = self.total_WR * W_accelerate
                print("W_accelerate (subsonic): " + str(W_accelerate))
                step_details = {"Mach number": M, "Weight ratio": W_accelerate, "Total weight ratio": self.total_WR}
                step = MissionStep("Accelerate at subsonic", step_details)
                self.mission_steps.append(step)
        else:
            raise ValueError("Error: Please provide corresponding values.")
        return W_accelerate


    # Equation 3.25, p. 64, Mattingly  CONVERT TO IMPERIAL UNITS!!!
    def TURN(self, M=None, theta=None, N=None, V_turn=None, TR=None):
        W_turn = None
        arguments = [M, theta, N, V_turn, TR]
        if all(i is not None for i in arguments):
            C = self.TSFC(M=M, theta=theta)  # [1/h]
            n = math.sqrt(1 + ((V_turn * 3.281) ** 2 / (self.g0 * 3.281) * TR * 0.54 * 6076) ** 2)
            W_turn = np.exp(-(C / 3600) * np.sqrt(theta) * ((self.CD + self.CDR) * n / self.CL) * (2*np.pi * N * V_turn) / (self.g0 * np.sqrt(n**2 - 1)))
            self.total_WR = self.total_WR * W_turn
            print("W_turn: " + str(W_turn))
            step_details = {"Mach number": M, "Weight ratio": W_turn, "Total weight ratio": self.total_WR}  # "Theta": theta, "Number of turns": N, "V_turn": V_turn, "Turn radius": TR
            step = MissionStep("Turn", step_details)
            self.mission_steps.append(step)
        else:
            raise ValueError("Error: Please provide corresponding values.")
        return W_turn

    # Equation 3.23, p.63, Mattingly
    """def FUEL_WR_CRUISE(self):
        C = 1.684  # [1/h]
        W_cruise = np.exp(- (C / (self.a * np.sqrt(self.theta)) * ((self.CD + self.CDR) / self.CL) * self.delta_s))
        return W_cruise
        pass"""

    # Brequet Range equation, Raymer  CONVERT TO IMPERIAL UNITS!!!
    def CRUISE(self, M=None, theta=None, R=None, V_cruise=None):
        W_cruise = None
        arguments = [M, theta, R, V_cruise]
        if all(i is not None for i in arguments):
            C = self.TSFC(M=M, theta=theta)
            q = 0.5 * self.rho * V_cruise ** 2
            # L_D_ratio = 4*(self.M+3)/self.M  # at supersonic 4(M+3)/M
            L_D_ratio = 1 / (((q * 0.02089 * self.CD0) / self.WSR) + (self.WSR * (1 / ((q * 0.02089) * np.pi * (self.S * 10.753) * self.e))))
            W_cruise = np.exp(- ((R * 0.54 * 6076) * (C / 3600)) / ((V_cruise * 3.281) * L_D_ratio))
            self.total_WR = self.total_WR * W_cruise
            print("W_cruise: " + str(W_cruise))
            step_details = {"Mach number": M, "Weight ratio": W_cruise, "Total weight ratio": self.total_WR}  # "Theta": theta, "Range": R, "V_cruise": V_cruise
            step = MissionStep("Cruise", step_details)
            self.mission_steps.append(step)
        else:
            raise ValueError("Error: Please provide corresponding values.")
        return W_cruise

    def LANDING(self):  # use fixed ratio
        W_landing = 0.993  # range 0.992 - 0.997
        self.total_WR = self.total_WR * W_landing
        step_details = {"Weight ratio": W_landing, "Total weight ratio": self.total_WR}
        step = MissionStep("Landing", step_details)
        self.mission_steps.append(step)
        return W_landing

    def analyze(self):
        print("")
        for step in self.mission_steps:
            step.display_details()


    def TOTAL_FUEL_WR(self):
        W_fuel = 1.06 * (1 - self.total_WR)
        # print('W_fuel: ' + str(W_fuel))
        return W_fuel

    # definition of empty weight ratio for UAV small
    def calculate_ewf(self, w_guess):
        ewf = (0.93 * (w_guess ** -0.06)) * 0.95
        # ewf = 1.495 * (w_guess ** -0.1)
        ew = ewf * w_guess
        print("Empty weight: " + str(ew))
        return ewf

    def TOTAL_WR(self, initial_w_guess, w_fuel):
        w_guess = initial_w_guess
        count = 1
        while True:
            print(" ")
            print("Iteration " + str(count))
            ewf = self.calculate_ewf(w_guess=w_guess)
            wf = w_fuel * w_guess
            print("Fuel weight: " + str(wf))

            WTO_calc = self.WP / (1 - w_fuel - ewf)

            if abs(WTO_calc - w_guess) < 0.001:
                break

            w_guess = WTO_calc
            count += 1
            print("Current takeoff weight: " + str(w_guess))
            print(" ")

        return w_guess

    """def WTO(self):  # TODO determine payload weight WP
        return self.WP / (1 - self.TOTAL_FUEL_WR() - self.TOTAL_EMPTY_WR(W_guess))"""

    def THRUST(self, WTO):
        thrust = self.res[0] * WTO
        print('Required thrust: ' + str(thrust))
        return thrust

    def WING_AREA(self, WTO):
        wing_area = self.res[1] / WTO
        print('Required wing area: ' + str(wing_area))
        return wing_area


class DragPolar:
    def __init__(self,filename):
        self.filename = filename
        self.df = pd.read_csv(self.filename,delim_whitespace=True)
        
        
    def dataframe(self):
        return self.df
    
    def calculate_coeff(self):
        """
        coeffcients = [K1,K2,CD0]
        """
        x = self.df['CL'].values.flatten()
        y = self.df['CDtot'].values.flatten()
        coefficients = np.polyfit(x, y, 2)
        K1 = coefficients[0]
        K2 = coefficients[1]
        CD0 = coefficients[2]
        return {'K1':K1,'K2':K2,'CD0':CD0}