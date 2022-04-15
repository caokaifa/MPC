from scipy import linalg
import numpy as np
from cvxopt.solvers import qp
from cvxopt import spmatrix, matrix, solvers
import datetime, pdb

solvers.options['show_progress'] = False


class LTI_MPC:

    def __init__(self, A, B, Q, R, N, vt, Fx, bx, Fu, bu):
        self.A = A
        self.B = B
        self.n = A.shape[0]
        self.d = B.shape[1]
        self.N = N
        self.Q = Q
        self.R = R
        self.vt = vt
        self.Fx = Fx
        self.bx = bx
        self.Fu = Fu
        self.bu = bu

        startTimer = datetime.datetime.now()
        endTimer = datetime.datetime.now()
        deltaTimer = endTimer - startTimer
        self.solverTime = deltaTimer
        self.linearizationTime = deltaTimer

        self.M, self.q = self.buildCost()
        self.F, self.b = self.buildIneqConst()
        self.G, self.E = self.buildEqConst()

    def solve(self, x0):
        M = self.M
        F = self.F
        G = self.G
        E = self.E
        q = self.q
        b = self.b
        n = self.n
        d = self.d
        N = self.N

        startTimer = datetime.datetime.now()
        sol = qp(M, matrix(q), F, matrix(b), G, E * matrix(x0))
        endTimer = datetime.datetime.now()
        deltaTimer = endTimer - startTimer
        self.solverTime = deltaTimer

        if sol['status'] == 'optimal':
            self.feasible = 1
        else:
            self.feasible = 0

        self.xPred = np.squeeze(np.transpose(np.reshape((np.squeeze(sol['x'])[np.arange(n * (N + 1))]), (N + 1, n)))).T
        self.uPred = np.squeeze(np.transpose(np.reshape((np.squeeze(sol['x'])[n * (N + 1) + np.arange(d * N)]), (N, d)))).T

    def buildEqConst(self):
        Gx = np.eye(self.n * (self.N + 1))
        Gu = np.zeros((self.n * (self.N + 1), self.d * (self.N)))
        E = np.zeros((self.n * (self.N + 1), self.n))
        E[np.arange(self.n)] = np.eye(self.n)

        for i in range(0, self.N):
            ind1 = self.n + i * self.n + np.arange(self.n)
            ind2x = i * self.n + np.arange(self.n)
            ind2u = i * self.d + np.arange(self.d)

            Gx[np.ix_(ind1, ind2x)] = -self.A
            Gu[np.ix_(ind1, ind2u)] = -self.B

        G = np.hstack((Gx, Gu))
        G_sparse = spmatrix(G[np.nonzero(G)], np.nonzero(G)[0].astype(int), np.nonzero(G)[1].astype(int), G.shape)
        E_sparse = spmatrix(E[np.nonzero(E)], np.nonzero(E)[0].astype(int), np.nonzero(E)[1].astype(int), E.shape)

        return G_sparse, E_sparse

    def buildIneqConst(self):
        rep_a = [self.Fx] * (self.N)
        Mat = linalg.block_diag(*rep_a)
        NoTerminalConstr = np.zeros((np.shape(Mat)[0], self.n))  # No need to constraint also the terminal point
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
        M = 2 * M0  # Need to multiply by two because CVX considers 1/2 in front of quadratic cost
        M_sparse = spmatrix(M[np.nonzero(M)], np.nonzero(M)[0].astype(int), np.nonzero(M)[1].astype(int), M.shape)
        return M_sparse, q



