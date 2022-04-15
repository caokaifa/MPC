import numpy as np
import pdb
import datetime
from Utilities import Curvature, getAngle


class Simulator():
    def __init__(self, map, lap = 0, flagLMPC = 0):
        self.map = map
        self.laps = lap
        self.flagLMPC = flagLMPC

    def Sim(self, ClosedLoopData, Controller, LMPCprediction=0):
        x = ClosedLoopData.x
        x_glob = ClosedLoopData.x_glob
        u = ClosedLoopData.u

        SimulationTime = 0
        for i in range(0, int(ClosedLoopData.Points)):
            Controller.solve(x[i, :])

            u[i, :] = Controller.uPred[0,:]

            if LMPCprediction != 0:
                Controller.LapTime = i
                LMPCprediction.PredictedStates[:,:,i, Controller.it]   = Controller.xPred
                LMPCprediction.PredictedInputs[:, :, i, Controller.it] = Controller.uPred
                LMPCprediction.SSused[:, :, i, Controller.it]          = Controller.SS_PointSelectedTot
                LMPCprediction.Qfunused[:, i, Controller.it]           = Controller.Qfun_SelectedTot
                x[i + 1, :], x_glob[i + 1, :] = dyModel(x[i, :], x_glob[i, :], u[i, :], np, ClosedLoopData.dt, self.map.PointAndTangent)
            else:
                x[i + 1, :], x_glob[i + 1, :] = dyModelori(x[i, :], x_glob[i, :], u[i, :], np, ClosedLoopData.dt, self.map.PointAndTangent)

            SimulationTime = i + 1

            if i <= 5:
                print("Linearization time: %.4fs Solver time: %.4fs" % (Controller.linearizationTime.total_seconds(), Controller.solverTime.total_seconds()))
                print("Time: ", i * ClosedLoopData.dt, "Current State and Input: ", x[i, :], u[i, :])

            if Controller.feasible == 0:
                print("Unfeasible at time ", i*ClosedLoopData.dt)
                print("Cur State: ", x[i, :], "Iteration ", Controller.it)
                break

            if self.flagLMPC == 1:
                Controller.addPoint(x[i, :], u[i, :])

            if (self.laps == 1) and (int(np.floor(x[i+1, 4] / (self.map.TrackLength))))>0:
                print("Simulation terminated: Lap completed")
                break

        ClosedLoopData.SimTime = SimulationTime
        print("Number of laps completed: ", int(np.floor(x[-1, 4] / (self.map.TrackLength))))


class PID:
    def __init__(self, vt):
        self.vt = vt
        self.uPred = np.zeros([1,2])

        startTimer = datetime.datetime.now()
        endTimer = datetime.datetime.now()
        deltaTimer = endTimer - startTimer
        self.solverTime = deltaTimer
        self.linearizationTime = deltaTimer
        self.feasible = 1

    def solve(self, x0):
        vt = self.vt
        self.uPred[0, 0] = - 0.6 * x0[5] - 0.9 * x0[3] + np.max([-0.9, np.min([np.random.randn() * 0.25, 0.9])])
        self.uPred[0, 1] = 1.5 * (vt - x0[0]) + np.max([-0.2, np.min([np.random.randn() * 0.10, 0.2])])


def dyModel(x, x_glob, u, np, dt, PointAndTangent):
    m = 1.98 # 1.98
    lf = 0.125
    lr = 0.125
    Iz = 0.024
    Df = 0.8 * m * 9.81 / 2.0
    Cf = 1.25 # 1.25
    Bf = 0.6 # 1.0
    Dr = 0.8 * m * 9.81 / 2.0
    Cr = 1.25
    Br = 1.0 # 1.0

    deltaT = 0.001
    x_next = np.zeros(x.shape[0])
    cur_x_next = np.zeros(x.shape[0])

    delta = u[0]
    a = u[1]

    psi = x_glob[3]
    X = x_glob[4]
    Y = x_glob[5]

    vx = x[0]
    vy = x[1]
    wz = x[2]
    epsi = x[3]
    s = x[4]
    ey = x[5]

    i = 0
    while (i+1) * deltaT <= dt:
        alpha_f = delta - np.arctan2(vy + lf * wz, vx)
        alpha_r = - np.arctan2(vy - lf * wz, vx)

        Fyf = Df * np.sin(Cf * np.arctan(Bf * alpha_f))
        Fyr = Dr * np.sin(Cr * np.arctan(Br * alpha_r))

        x_next[0] = vx + deltaT * (a - 1 / m * Fyf * np.sin(delta) + wz*vy)
        x_next[1] = vy + deltaT * (1 / m * (Fyf * np.cos(delta) + Fyr) - wz * vx)
        x_next[2] = wz + deltaT * (1 / Iz *(lf * Fyf * np.cos(delta) - lr * Fyr))
        x_next[3] = psi + deltaT * (wz)
        x_next[4] = X + deltaT * (vx * np.cos(psi) - vy * np.sin(psi))
        x_next[5] = Y + deltaT * (vx * np.sin(psi) + vy * np.cos(psi))

        cur = Curvature(s, PointAndTangent)
        cur_x_next[0] = x_next[0]
        cur_x_next[1] = x_next[1]
        cur_x_next[2] = x_next[2]
        cur_x_next[3] = epsi + deltaT * (wz - (vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey) * cur)
        cur_x_next[4] = s + deltaT * ((vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey))
        cur_x_next[5] = ey + deltaT * (vx * np.sin(epsi) + vy * np.cos(epsi))

        psi = x_next[3]
        X = x_next[4]
        Y = x_next[5]

        vx = cur_x_next[0]
        vy = cur_x_next[1]
        wz = cur_x_next[2]
        epsi = cur_x_next[3]
        s = cur_x_next[4]
        ey = cur_x_next[5]

        if s < 0:
            print("Start Point: ", x, " Input: ", u)
            print("x_next: ", x_next)
        i = i+1

    noise_vx = np.max([-0.05, np.min([np.random.randn() * 0.01, 0.05])])
    noise_vy = np.max([-0.1, np.min([np.random.randn() * 0.01, 0.1])])
    noise_wz = np.max([-0.05, np.min([np.random.randn() * 0.005, 0.05])])

    cur_x_next[0] = cur_x_next[0] + 0.1*noise_vx
    cur_x_next[1] = cur_x_next[1] + 0.1*noise_vy
    cur_x_next[2] = cur_x_next[2] + 0.1*noise_wz

    return cur_x_next, x_next


def dyModelori(x, x_glob, u, np, dt, PointAndTangent):
    m = 1.98 # 1.98
    lf = 0.125
    lr = 0.125
    Iz = 0.024
    Df = 0.8 * m * 9.81 / 2.0
    Cf = 1.25 # 1.25
    Bf = 1.0 # 1.0
    Dr = 0.8 * m * 9.81 / 2.0
    Cr = 1.25
    Br = 1.0 # 1.0

    deltaT = 0.001
    x_next = np.zeros(x.shape[0])
    cur_x_next = np.zeros(x.shape[0])

    delta = u[0]
    a = u[1]

    psi = x_glob[3]
    X = x_glob[4]
    Y = x_glob[5]

    vx = x[0]
    vy = x[1]
    wz = x[2]
    epsi = x[3]
    s = x[4]
    ey = x[5]

    i = 0
    while (i+1) * deltaT <= dt:
        alpha_f = delta - np.arctan2(vy + lf * wz, vx)
        alpha_r = - np.arctan2(vy - lf * wz, vx)

        Fyf = Df * np.sin(Cf * np.arctan(Bf * alpha_f))
        Fyr = Dr * np.sin(Cr * np.arctan(Br * alpha_r))

        x_next[0] = vx + deltaT * (a - 1 / m * Fyf * np.sin(delta) + wz*vy)
        x_next[1] = vy + deltaT * (1 / m * (Fyf * np.cos(delta) + Fyr) - wz * vx)
        x_next[2] = wz + deltaT * (1 / Iz *(lf * Fyf * np.cos(delta) - lr * Fyr))
        x_next[3] = psi + deltaT * (wz)
        x_next[4] = X + deltaT * (vx * np.cos(psi) - vy * np.sin(psi))
        x_next[5] = Y + deltaT * (vx * np.sin(psi) + vy * np.cos(psi))

        cur = Curvature(s, PointAndTangent)
        cur_x_next[0] = x_next[0]
        cur_x_next[1] = x_next[1]
        cur_x_next[2] = x_next[2]
        cur_x_next[3] = epsi + deltaT * (wz - (vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey) * cur)
        cur_x_next[4] = s + deltaT * ((vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey))
        cur_x_next[5] = ey + deltaT * (vx * np.sin(epsi) + vy * np.cos(epsi))

        psi = x_next[3]
        X = x_next[4]
        Y = x_next[5]

        vx = cur_x_next[0]
        vy = cur_x_next[1]
        wz = cur_x_next[2]
        epsi = cur_x_next[3]
        s = cur_x_next[4]
        ey = cur_x_next[5]

        if s < 0:
            print("Start Point: ", x, " Input: ", u)
            print("x_next: ", x_next)
        i = i+1

    noise_vx = np.max([-0.05, np.min([np.random.randn() * 0.01, 0.05])])
    noise_vy = np.max([-0.1, np.min([np.random.randn() * 0.01, 0.1])])
    noise_wz = np.max([-0.05, np.min([np.random.randn() * 0.005, 0.05])])

    cur_x_next[0] = cur_x_next[0] + 0.1*noise_vx
    cur_x_next[1] = cur_x_next[1] + 0.1*noise_vy
    cur_x_next[2] = cur_x_next[2] + 0.1*noise_wz

    return cur_x_next, x_next


def wrap(angle):
    if angle < -np.pi:
        w_angle = 2 * np.pi + angle
    elif angle > np.pi:
        w_angle = angle - 2 * np.pi
    else:
        w_angle = angle

    return w_angle