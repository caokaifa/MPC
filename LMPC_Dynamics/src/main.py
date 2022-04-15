import sys
sys.path.append('/home/robert/LMPC论文/LMPC_Dynamics/src/schemes')
from SysModel import Simulator, PID
from Classes import ClosedLoopData, LMPCprediction
from LTVMPC import LTV_MPC
from LTIMPC import LTI_MPC
from Track import Map, unityTestChangeOfCoordinates
from LMPC import LMPCplus
from Utilities import Regression
from plot import plotTrajectory, plotClosedLoopLMPC, animation_xy, animation_states, saveGif_xyResults, Save_statesAnimation, plotCost, plotC, plotQC, plotQt, plotQcost
import numpy as np
import matplotlib.pyplot as plt
import pdb
import pickle


def main():

    # Parameter initialization
    dt = 1.0 / 10.0   # Controller discretization time
    Time = 100        # Simulation time for PID
    TimeMPC = 100     # Time for LTI-MPC
    TimeMPC_tv = 100  # Time for LTV-MPC
    TimeLMPC = 400    # Time for LMPC
    vt = 0.8          # Reference velocity for path controllers
    v0 = 0.5          # Initial velocity at lap 0
    N = 12            # Horizon
    dim_state = 6     # State dimension
    dim_input = 2     # Input dimension

    Q = np.diag([1.0, 1.0, 1, 1, 0.0, 100.0])             # vx, vy, wz, epsi, s, ey
    R = np.diag([1.0, 10.0])                              # delta, a
    Q_lmpc = np.diag([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) * 0  # vx, vy, wz, epsi, s, ey
    R_lmpc = np.diag([1.0, 1.0]) * 0                      # delta, a
    Qf = np.array([0, 10]) * 1
    QterminalSlack = np.diag([10, 1, 1, 1, 10, 1]) * 20
    dR_LMPC = np.array([1.0, 10.0]) * 10

    inputConstr = np.array([[0.5, 0.5],
                            [10.0, 10.0]])

    LMPC_Solver = "CVX"   # Can pick CVX for cvxopt or OSQP. For OSQP uncomment line 14 in LMPC.py
    numSS_it = 4          # Number of trajectories used at each iteration to build the safe set
    numSS_Points = 40     # Number of points to select from each trajectory to build the safe set

    Laps = 46 + numSS_it  # Total LMPC laps (50 laps)

    map = Map(0.4)                                             # Initialize the map
    model = Simulator(map)                                     # Initialize the MPC model
    LMPCmodel = Simulator(map, 1, 1)                           # Initialize the LMPC model

    # State constraints for LTI-MPC and LTV-MPC
    Fx_MPC = np.array([[1., 0., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 0., 1.],
                   [0., 0., 0., 0., 0., -1.]])

    bx_MPC = np.array([[10],
                   [2.],
                   [2.]])

    # Input constraints for LTI-MPC and LTV-MPC
    Fu_MPC = np.array([[1., 0.],
                       [-1., 0.],
                       [0., 1.],
                       [0., -1.]])

    bu_MPC = np.array([[inputConstr[0, 0]],
                       [inputConstr[0, 1]],
                       [inputConstr[1, 0]],
                       [inputConstr[1, 1]]])

    # State constraints for LMPC
    Fx = np.array([[0., 0., 0., 0., 0., 1.],
                   [0., 0., 0., 0., 0., -1.]])
    bx = np.array([[map.halfWidth],
                   [map.halfWidth]])

    # Input constraints for LMPC
    Fu = np.array([[1., 0.],
                   [-1., 0.],
                   [0., 1.],
                   [0., -1.]])
    bu = np.array([[inputConstr[0,0]],
                   [inputConstr[0,1]],
                   [inputConstr[1,0]],
                   [inputConstr[1,1]]])

    print("Starting PID")
    ClosedLoopDataPID = ClosedLoopData(dt, Time, v0)
    PIDController = PID(vt)
    model.Sim(ClosedLoopDataPID, PIDController)

    file_data = open(sys.path[0]+'\data\ClosedLoopDataPID.obj', 'wb')
    pickle.dump(ClosedLoopDataPID, file_data)
    file_data.close()
    print("===== PID terminated")

    print("Starting LTI-MPC")
    lamb = 0.0000001
    A, B, Error = Regression(ClosedLoopDataPID.x, ClosedLoopDataPID.u, lamb)
    ClosedLoopDataLTI_MPC = ClosedLoopData(dt, TimeMPC, v0)
    LTIMPC = LTI_MPC(A, B, Q, R, N, vt, Fx_MPC, bx_MPC, Fu_MPC, bu_MPC)
    model.Sim(ClosedLoopDataLTI_MPC, LTIMPC)

    file_data = open(sys.path[0] + '\data\ClosedLoopDataLTI_MPC.obj', 'wb')
    pickle.dump(ClosedLoopDataLTI_MPC, file_data)
    file_data.close()
    print("===== LTI-MPC terminated")

    print("Starting LTV-MPC")
    ClosedLoopDataLTV_MPC = ClosedLoopData(dt, TimeMPC_tv, v0)
    LTVMPC = LTV_MPC(Q, R, N, vt, dim_state, dim_input, ClosedLoopDataPID.x, ClosedLoopDataPID.u, dt, map, Fx_MPC, bx_MPC, Fu_MPC, bu_MPC)
    model.Sim(ClosedLoopDataLTV_MPC, LTVMPC)

    file_data = open(sys.path[0] + 'data\ClosedLoopDataLTV_MPC.obj', 'wb')
    pickle.dump(ClosedLoopDataLTV_MPC, file_data)
    file_data.close()
    print("===== LTV-MPC terminated")

    print("Starting LMPC")
    ClosedLoopLMPC =  ClosedLoopData(dt, TimeLMPC, v0)
    LMPCOpenLoopData = LMPCprediction(N, dim_state, dim_input, TimeLMPC, numSS_Points, Laps)

    LMPC = LMPCplus(numSS_Points, numSS_it, N, QterminalSlack, Qf, Q_lmpc, R_lmpc, dR_LMPC, dt, map, Laps, TimeLMPC, LMPC_Solver, Fx, bx, Fu, bu)
    LMPC.addTrajectory(ClosedLoopDataPID)
    LMPC.addTrajectory(ClosedLoopDataLTV_MPC)
    LMPC.addTrajectory(ClosedLoopDataPID)
    LMPC.addTrajectory(ClosedLoopDataLTI_MPC)

    x0 = np.zeros((1, dim_state))
    x0_glob = np.zeros((1, dim_state))
    x0[0, :] = ClosedLoopLMPC.x[0, :]
    x0_glob[0, :] = ClosedLoopLMPC.x_glob[0, :]

    for it in range(numSS_it, Laps):
        ClosedLoopLMPC.updateInitialConditions(x0, x0_glob)
        LMPCmodel.Sim(ClosedLoopLMPC, LMPC, LMPCOpenLoopData)
        LMPC.addTrajectory(ClosedLoopLMPC)

        if LMPC.feasible == 0:
            break
        else:
            # Reset Initial Conditions
            x0[0, :] = ClosedLoopLMPC.x[ClosedLoopLMPC.SimTime, :] - np.array([0, 0, 0, 0, map.TrackLength, 0])
            x0_glob[0, :] = ClosedLoopLMPC.x_glob[ClosedLoopLMPC.SimTime, :]

    file_data = open(sys.path[0] + '\data\LMPController.obj', 'wb')
    pickle.dump(ClosedLoopLMPC, file_data)
    pickle.dump(LMPC, file_data)
    pickle.dump(LMPCOpenLoopData, file_data)
    file_data.close()
    print("===== LMPC terminated")

    laptimes = np.zeros((50, 2))

    # Laptime Plot
    for i in range(0, LMPC.it):
        print("Lap time at iteration ", i, " is ", LMPC.Qfun[0, i] * dt, "s")
        laptimes[i, 0] = LMPC.Qfun[0, i] * dt
        laptimes[i, 1] = i
    plt.figure(3)
    plt.plot(laptimes[:, 1], laptimes[:, 0], '-o')
    plt.ylabel('Lap Time (sec)')
    plt.xlabel('Lap Number')

    print("===== Start Plotting")

    plotTrajectory(map, ClosedLoopDataPID.x, ClosedLoopDataPID.x_glob, ClosedLoopDataPID.u)

    plotTrajectory(map, ClosedLoopDataLTI_MPC.x, ClosedLoopDataLTI_MPC.x_glob, ClosedLoopDataLTI_MPC.u)

    plotTrajectory(map, ClosedLoopDataLTV_MPC.x, ClosedLoopDataLTV_MPC.x_glob, ClosedLoopDataLTV_MPC.u)

    plotCost(LMPC.Qfun, int(TimeLMPC / dt) + 1)

    plotC(LMPC.Qfun_SelectedTot, numSS_it)

    plotQC(LMPC.Qcost, numSS_Points)

    plotQt(LMPC.qq)

    plotQcost(LMPC.costSolved)

    plotClosedLoopLMPC(LMPC, map)

    animation_xy(map, LMPCOpenLoopData, LMPC, Laps - 2)

    animation_states(map, LMPCOpenLoopData, LMPC, 10)

    unityTestChangeOfCoordinates(map, ClosedLoopDataPID)
    unityTestChangeOfCoordinates(map, ClosedLoopDataLTI_MPC)
    unityTestChangeOfCoordinates(map, ClosedLoopLMPC)

    saveGif_xyResults(map, LMPCOpenLoopData, LMPC, Laps-1)
    Save_statesAnimation(map, LMPCOpenLoopData, LMPC, 5)

    plt.show()



if __name__ == "__main__":
    main()