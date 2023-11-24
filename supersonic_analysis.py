import numpy as np
import matplotlib.pyplot as plt
import ambiance as amb
import math
import pandas as pd
import scipy.optimize as opt

class ConstraintAnalysis:
    def __init__(self, W, WP, S, b, AR, e, V, V_stall, V_takeoff, rho, mu, k, k2, CD0, CL, CD, CDR, g0, q, ROC, TR, n, dv_dt, alpha, beta,safety_margin_TW=0,safety_margin_WS=0,plot_max_x=500,plot_max_y=2):
        self.W = W
        self.WP = WP
        self.S = S
        self.b = b
        self.AR = AR
        self.e = e
        self.V = V
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
        self.q = q
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
        return (self.beta / self.alpha) * (self.k * (self.beta / self.q) * WTO_S + self.k2 + (self.CD0 + self.CDR) / ((self.beta / self.q) * WTO_S))

    def TSL_WTO_CLIMB(self, WTO_S):
        return (self.beta / self.alpha) * (self.k * (self.beta / self.q) * WTO_S + self.k2 + (self.CD0 + self.CDR) / ((self.beta / self.q) * WTO_S) + (1 / self.V) * self.ROC)

    def TSL_WTO_TO(self, WTO_S):
        return (self.beta / self.alpha) * ((self.CD + self.CDR - self.mu * self.CL) * (self.beta / self.q) * (WTO_S ** -1) + self.mu + (1 / self.g0) * self.dv_dt)

    def TSL_WTO_TURN(self, WTO_S):
        return (self.beta / self.alpha) * (
                    self.k * (self.n ** 2) * (self.beta / self.q) * WTO_S + self.k2 * self.n + ((self.CD0 + self.CDR) * self.q) / (self.beta * WTO_S))
    
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
        return round(math.sqrt(1 + (self.V ** 2 / (self.g0 * self.TR)) ** 2),3)


class MissionAnalysis:
    def __init__(self, WP, a, M, q, CL, CD, CDR, alpha, beta, res, delta_h, n, theta, N, V, g0, R, mu, V_takeoff):
        self.WP = WP
        self.a = a
        self.M = M
        self.q = q
        self.CL = CL
        self.CD = CD
        self.CDR = CDR
        self.alpha = alpha
        self.beta = beta
        self.res = res
        self.delta_h = delta_h
        self.n = n
        self.theta = theta
        self.N = N
        self.V = V
        self.g0 = g0
        self.R = R
        self.mu = mu
        self.V_takeoff = V_takeoff

    # Turbojet engine max power, eq. 3.55b, p. 71, Mattingly
    def TSFC(self, M, theta):
        return (1.5 + 0.23 * M) * np.sqrt(theta)

    # Equation 3.20, p 62, Mattingly
    def FUEL_WR_CLIMB(self):
        """theta = 0.9931
        C = self.TSFC(M=self.M, theta=theta)  # [1/h]
        u = ((self.CD + self.CDR) / self.CL) * (self.beta / self.alpha) * (self.res[0]**-1)
        W_climb = np.exp((-C / self.a) * ((self.delta_h + (self.V**2) / (2*self.g0)) / (1 - u)))
        return W_climb"""
        return 0.920

    # Equation 3.25, p. 64, Mattingly
    def FUEL_WR_TURN(self):
        """C = self.TSFC(M=self.M, theta=self.theta)  # [1/h]
        W_turn = np.exp(-C * np.sqrt(self.theta) * ((self.CD + self.CDR) * self.n / self.CL) * (2*np.pi * self.N * self.V) / (self.g0 * np.sqrt(self.n**2 - 1)))
        return W_turn"""

    # Equation 3.23, p.63, Mattingly
    def FUEL_WR_CRUISE(self):
        """C = 1.684  # [1/h]
        W_cruise = np.exp(- (C / (self.a * np.sqrt(self.theta)) * ((self.CD + self.CDR) / self.CL) * self.delta_s))
        return W_cruise"""
        pass

    # Brequet Range equation, Raymer
    def FUEL_WR_CRUISE(self):
        C = self.TSFC(M=self.M, theta=self.theta)
        V = 280.5
        L_D_ratio = 4*(self.M+3)/self.M  # at supersonic 4(M+3)/M
        W_cruise = np.exp(- (self.R*C) / (self.V*L_D_ratio))
        return W_cruise

    # Equation 3.21, 3.22, p. 63, Mattingly
    def FUEL_WR_TAKEOFF(self):
        """theta = 1.000
        C = self.TSFC(M=0.05, theta=theta) # [1/h]
        xi = self.CD + self.CDR - self.mu * self.CL
        u = (xi * (self.q * self.beta) * ((self.res[1])**-1) + self.mu) * (self.beta / self.alpha) * self.res[0]
        W_takeoff = np.exp(-C * np.sqrt(theta) / self.g0 * (self.V_takeoff / (1 - u)))
        return W_takeoff"""
        # Historical value
        return 0.995

    def TOTAL_FUEL_WR(self):
        w_x = self.FUEL_WR_CLIMB() * self.FUEL_WR_CRUISE() * self.FUEL_WR_TAKEOFF()
        print('Climb: ' + str(self.FUEL_WR_CLIMB()))
        print('Turn: ' + str(self.FUEL_WR_TURN()))
        print('Cruise: ' + str(self.FUEL_WR_CRUISE()))
        print('Takeoff: ' + str(self.FUEL_WR_TAKEOFF()))
        print('W_x: ' + str(w_x))
        W_fuel = 1.06 * (1 - w_x)
        print('W_fuel: ' + str(W_fuel))
        return W_fuel

    # definition of empty weight ratio for UAV small
    def calculate_ewf(self, w_guess):
        ewf = 0.97 * w_guess ** -0.06
        return ewf

    def TOTAL_WR(self, initial_w_guess, w_fuel):
        w_guess = initial_w_guess
        while True:
            ewf = self.calculate_ewf(w_guess=w_guess)
            WTO_calc = self.WP / (1 - w_fuel - ewf)

            if abs(WTO_calc - w_guess) < 0.001:
                break

            w_guess = WTO_calc
            print(w_guess)

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