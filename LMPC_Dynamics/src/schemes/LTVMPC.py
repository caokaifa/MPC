from scipy import linalg
import numpy as np
from cvxopt.solvers import qp
from cvxopt import spmatrix, matrix, solvers
from Utilities import Curvature
import datetime
from numpy import linalg as la

solvers.options['show_progress'] = False


class LTV_MPC:

    def __init__(self, Q, R, N, vt, n, d, x, u, dt, map, Fx, bx, Fu, bu):
        self.A = []
        self.B = []
        self.C = []
        self.N = N
        self.n = n
        self.d = d
        self.vt = vt
        self.Q = Q
        self.R = R
        self.Datax = x
        self.Datau = u
        self.LinPoints = x[0:N+1, :]
        self.dt = dt
        self.map = map
        self.Fx = Fx
        self.bx = bx
        self.Fu = Fu
        self.bu = bu

        self.M, self.q = self.buildCost()
        self.F, self.b = self.buildIneqConst()
        self.G = []
        self.E = []

    def solve(self, x0):
        startTimer = datetime.datetime.now()
        self.A, self.B, self. C = EstimateABC(self)
        self.G, self.E, self.L = self.buildEqConst()
        endTimer = datetime.datetime.now()
        deltaTimer = endTimer - startTimer
        self.linearizationTime = deltaTimer

        M = self.M
        q = self.q
        G = self.G
        L = self.L
        E = self.E
        F = self.F
        b = self.b
        n = self.n
        N = self.N
        d = self.d

        startTimer = datetime.datetime.now()
        sol = qp(M, matrix(q), F, matrix(b), G, E * matrix(x0) + L)
        endTimer = datetime.datetime.now()
        deltaTimer = endTimer - startTimer
        self.solverTime = deltaTimer
        if sol['status'] == 'optimal':
            self.feasible = 1
        else:
            self.feasible = 0

        self.xPred = np.squeeze(np.transpose(np.reshape((np.squeeze(sol['x'])[np.arange(n * (N + 1))]), (N + 1, n))))
        self.uPred = np.squeeze(np.transpose(np.reshape((np.squeeze(sol['x'])[n * (N + 1) + np.arange(d * N)]), (N, d))))
        self.LinPoints = np.concatenate((self.xPred.T[1:, :], np.array([self.xPred.T[-1, :]])), axis=0)
        self.xPred = self.xPred.T
        self.uPred = self.uPred.T

    def buildIneqConst(self):
        rep_a = [self.Fx] * (self.N)
        Mat = linalg.block_diag(*rep_a)
        NoTerminalConstr = np.zeros((np.shape(Mat)[0], self.n))
        Fxtot = np.hstack((Mat, NoTerminalConstr))
        bxtot = np.tile(np.squeeze(self.bx), self.N)

        rep_b = [self.Fu] * (self.N)
        Futot = linalg.block_diag(*rep_b)
        butot = np.tile(np.squeeze(self.bu), self.N)

        rFxtot, cFxtot = np.shape(Fxtot)
        rFutot, cFutot = np.shape(Futot)
        Dummy1 = np.hstack((Fxtot, np.zeros((rFxtot, cFutot))))
        Dummy2 = np.hstack((np.zeros((rFutot, cFxtot)), Futot))
        F = np.vstack((Dummy1, Dummy2))
        b = np.hstack((bxtot, butot))
        F_sparse = spmatrix(F[np.nonzero(F)], np.nonzero(F)[0].astype(int), np.nonzero(F)[1].astype(int), F.shape)

        return F_sparse, b

    def buildCost(self):
        b = [self.Q] * (self.N)
        Mx = linalg.block_diag(*b)

        c = [self.R] * (self.N)
        Mu = linalg.block_diag(*c)

        M0 = linalg.block_diag(Mx, self.Q, Mu)
        xtrack = np.array([self.vt, 0, 0, 0, 0, 0])
        q = - 2 * np.dot(np.append(np.tile(xtrack, self.N + 1), np.zeros(self.R.shape[0] * self.N)), M0)
        M = 2 * M0
        M_sparse = spmatrix(M[np.nonzero(M)], np.nonzero(M)[0].astype(int), np.nonzero(M)[1].astype(int), M.shape)

        return M_sparse, q

    def buildEqConst(self):
        Gx = np.eye(self.n * (self.N + 1))
        Gu = np.zeros((self.n * (self.N + 1), self.d * (self.N)))

        E = np.zeros((self.n * (self.N + 1), self.n))
        E[np.arange(self.n)] = np.eye(self.n)

        L = np.zeros((self.n * (self.N + 1), 1))

        for i in range(0, self.N):
            ind1 = self.n + i * self.n + np.arange(self.n)
            ind2x = i * self.n + np.arange(self.n)
            ind2u = i * self.d + np.arange(self.d)

            Gx[np.ix_(ind1, ind2x)] = -self.A[i]
            Gu[np.ix_(ind1, ind2u)] = -self.B[i]
            L[ind1, :] = self.C[i]

        G = np.hstack((Gx, Gu))

        G_sparse = spmatrix(G[np.nonzero(G)], np.nonzero(G)[0].astype(int), np.nonzero(G)[1].astype(int), G.shape)
        E_sparse = spmatrix(E[np.nonzero(E)], np.nonzero(E)[0].astype(int), np.nonzero(E)[1].astype(int), E.shape)
        L_sparse = spmatrix(L[np.nonzero(L)], np.nonzero(L)[0].astype(int), np.nonzero(L)[1].astype(int), L.shape)

        return G_sparse, E_sparse, L_sparse


def EstimateABC(MPC):
    LinPoints = MPC.LinPoints
    N = MPC.N
    n = MPC.n
    d = MPC.d
    x = MPC.Datax
    u = MPC.Datau
    PointAndTangent = MPC.map.PointAndTangent
    dt = MPC.dt

    Atv = []; Btv = []; Ctv = []
    for i in range(0, N + 1):
        MaxNumPoint = 500
        x0 = LinPoints[i, :]

        Ai = np.zeros((n, n))
        Bi = np.zeros((n, d))
        Ci = np.zeros((n, 1))

        h = 2
        stateFeatures = [0, 1, 2]
        inputFeatures = [1]
        lamb = 0.0
        yIndex = 0
        scaling = np.array([[1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0]])

        Ai[yIndex, stateFeatures], Bi[yIndex, inputFeatures], Ci[yIndex], _ = LocLinReg(h, x, u, x0, yIndex,
                                                                                        stateFeatures,
                                                                                        inputFeatures, scaling, qp,
                                                                                        matrix, lamb, MaxNumPoint)
        h = 2
        stateFeatures = [0, 1, 2]
        inputFeatures = [0]
        lamb = 0.0
        yIndex = 1

        scaling = np.eye(len(stateFeatures))

        Ai[yIndex, stateFeatures], Bi[yIndex, inputFeatures], Ci[yIndex], _ = LocLinReg(h, x, u, x0, yIndex, stateFeatures,
                                                                             inputFeatures, scaling, qp, matrix, lamb, MaxNumPoint)

        h = 2
        stateFeatures = [0, 1, 2]
        inputFeatures = [0]
        lamb = 0.0
        yIndex = 2
        scaling = np.eye(len(stateFeatures))

        Ai[yIndex, stateFeatures], Bi[yIndex, inputFeatures], Ci[yIndex], _ = LocLinReg(h, x, u, x0, yIndex, stateFeatures,
                                                                             inputFeatures, scaling, qp, matrix, lamb, MaxNumPoint)

        vx = x0[0]
        vy = x0[1]
        wz = x0[2]
        epsi = x0[3]
        s = x0[4]
        ey = x0[5]

        if s < 0:
            print("s is negative, here the state: \n", LinPoints)

        cur = Curvature(s, PointAndTangent)
        den = 1 - cur *ey
        depsi_vx = -dt * np.cos(epsi) / den * cur
        depsi_vy = dt * np.sin(epsi) / den * cur
        depsi_wz = dt
        depsi_epsi = 1 - dt * (-vx * np.sin(epsi) - vy * np.cos(epsi)) / den * cur
        depsi_s = 0
        depsi_ey = dt * (vx * np.cos(epsi) - vy * np.sin(epsi)) / (den**2) * cur * (-cur)

        Ai[3, :] = [depsi_vx, depsi_vy, depsi_wz, depsi_epsi, depsi_s, depsi_ey]
        Ci[3] = epsi + dt * (wz - (vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey) * cur) - np.dot(Ai[3, :], x0)

        ds_vx = dt * (np.cos(epsi) / den)
        ds_vy = -dt * (np.sin(epsi) / den)
        ds_wz = 0
        ds_epsi = dt * (-vx * np.sin(epsi) - vy * np.cos(epsi)) / den
        ds_s = 1
        ds_ey = -dt * (vx * np.cos(epsi) - vy * np.sin(epsi)) / (den * 2) * (-cur)

        Ai[4, :] = [ds_vx, ds_vy, ds_wz, ds_epsi, ds_s, ds_ey]
        Ci[4] = s + dt * ((vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey)) - np.dot(Ai[4, :], x0)

        dey_vx = dt * np.sin(epsi)
        dey_vy = dt * np.cos(epsi)
        dey_wz = 0
        dey_epsi = dt * (vx * np.cos(epsi) - vy * np.sin(epsi))
        dey_s = 0
        dey_ey = 1

        Ai[5, :] = [dey_vx, dey_vy, dey_wz, dey_epsi, dey_s, dey_ey]
        Ci[5] = ey + dt * (vx * np.sin(epsi) + vy * np.cos(epsi)) - np.dot(Ai[5, :], x0)

        Atv.append(Ai)
        Btv.append(Bi)
        Ctv.append(Ci)

    return Atv, Btv, Ctv


def LocLinReg(h, x, u, x0, yIndex, stateFeatures, inputFeatures, scaling, qp, matrix, lamb, MaxNumPoint):
    oneVec = np.ones( (x.shape[0]-1, 1) )
    x0Vec = (np.dot( np.array([x0[stateFeatures]]).T, oneVec.T )).T
    diff  = np.dot(( x[0:-1, stateFeatures] - x0Vec ), scaling)

    norm = la.norm(diff, 1, axis=1)
    indexTot =  np.squeeze(np.where(norm < h))
    if indexTot.shape[0] >= MaxNumPoint:
        index = np.argsort(norm)[0:MaxNumPoint]
    else:
        index = indexTot

    K = (1 - ( norm[index] / h )**2) * 3/4
    X0 = np.hstack( ( x[np.ix_(index, stateFeatures)], u[np.ix_(index, inputFeatures)]))
    M = np.hstack((X0, np.ones((X0.shape[0],1))))
    y = x[np.ix_(index+1, [yIndex])]
    b = matrix(-np.dot(np.dot(M.T, np.diag(K)), y))

    Q0 = np.dot(np.dot(M.T, np.diag(K)), M)
    Q = matrix(Q0 + lamb * np.eye(Q0.shape[0]))

    res_cons = qp(Q, b)
    Result = np.squeeze(np.array(res_cons['x']))
    A = Result[0:len(stateFeatures)]
    B = Result[len(stateFeatures):(len(stateFeatures)+len(inputFeatures))]
    C = Result[-1]

    return A, B, C, index