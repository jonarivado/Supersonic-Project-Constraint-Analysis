import numpy as np
import matplotlib.pyplot as plt
import ambiance as amb
import math
import pandas as pd
import scipy.optimize as opt

class ConstraintAnalysis:
    def __init__(self,W, S, b, AR, e, V, V_stall, V_takeoff, rho, mu, k, k2, CD0, CL, CD, CDR, g0, q, ROC, TR, n, dv_dt, alpha, beta,safety_margin_TW=0,safety_margin_WS=0):
        self.W = W
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
        self.res = [0,0]
        self.safety_margin_TW = safety_margin_TW
        self.safety_margin_WS = safety_margin_WS

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
        self.res = opt.minimize(lambda x: x[0], x0, constraints=constraints,bounds=((0, 10), (0, 1000)))
        print("Opimal T/W: ", round(self.res.x[0],3))
        print("Optimal W/S: ", round(self.res.x[1],3))
        print("Required Thrust: ", round(self.res.x[0]*round(self.res.x[1]*self.S,3),3), " N")
        print("Optimal W: ", round(self.res.x[1]*self.S/self.g0,3), " kg")

        return self.res.x[0], self.res.x[1]
    
    def plot(self):
        
        #assert that the optimize function has been run, otherwise throw an error
        assert hasattr(self, 'res') and hasattr(self.res, 'x'), "Optimization not performed. Please run optimize() first."

        WTO_S = np.linspace(0.1, 7000,2000)
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
        plt.xlim(0, 500)
        plt.ylim(0, 2)
        plt.legend()
        plt.title('Optima w/ margin: T/W: ' +str(round(self.res.x[1],3))+ ', W/S: '+str(round(self.res.x[0],3)))
        plt.show()
