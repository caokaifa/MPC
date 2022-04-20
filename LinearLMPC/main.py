import numpy as np
from FTOCP import FTOCP
from LMPC import LMPC
import pdb
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import copy
import pickle

def main():
	# Define system dynamics and cost
	A = np.array([[1,1],[0,1]])
	B = np.array([[0],[1]])
	Q = np.diag([1.0, 1.0]) #np.eye(2)
	R = 1.0#np.array([[1]])

	print("Computing first feasible trajectory")
	
	# Initial Condition
	x0 = [-15.0, 0.0];

	# Initialize FTOCP object
	N_feas = 10
	ftocp_for_mpc  = FTOCP(N_feas, A, B, 0.01*Q, R)

	# ====================================================================================
	# Run simulation to compute feasible solution
	# ====================================================================================
	xcl_feasible = [x0]
	ucl_feasible =[]
	xt           = x0
	time         = 0

	# time Loop (Perform the task until close to the origin)
	while np.dot(xt, xt) > 10**(-15):
		xt = xcl_feasible[time] # Read measurements

		ftocp_for_mpc.solve(xt, verbose = False) # Solve FTOCP

		# Read input and apply it to the system
		ut = ftocp_for_mpc.uPred[:,0][0]
		ucl_feasible.append(ut)
		xcl_feasible.append(ftocp_for_mpc.model(xcl_feasible[time], ut))
		time += 1

	xcl_N = np.array(xcl_feasible).T
	plt.plot(xcl_N[0,:], xcl_N[1,:], 'sr', label='sampled Safe Set')
	# plt.show()


	print(np.round(np.array(xcl_feasible).T, decimals=2))
	print(np.round(np.array(ucl_feasible).T, decimals=2))
	# ====================================================================================

	# ====================================================================================
	# Run LMPC
	# ====================================================================================

	# Initialize LMPC object
	N_LMPC = 3 # horizon length
	ftocp = FTOCP(N_LMPC, A, B, Q, R) # ftocp solved by LMPC
	lmpc = LMPC(ftocp, CVX=True) # Initialize the LMPC (decide if you wanna use the CVX hull)
	lmpc.addTrajectory(xcl_feasible, ucl_feasible) # Add feasible trajectory to the safe set
	
	totalIterations = 20 # Number of iterations to perform
	plt.figure()


	# run simulation
	# iteration loop
	print("Starting LMPC")
	for it in range(0,totalIterations):
		# Set initial condition at each iteration
		xcl = [x0] 
		ucl =[]
		time = 0
		# time Loop (Perform the task until close to the origin)
		while np.dot(xcl[time], xcl[time]) > 10**(-10):
			
			# Read measurement
			xt = xcl[time] 

			# Solve FTOCP
			lmpc.solve(xt, verbose = False) 
			# Read optimal input
			ut = lmpc.uPred[:,0][0]

			# Apply optimal input to the system
			ucl.append(ut)
			xcl.append(lmpc.ftocp.model(xcl[time], ut))
			time += 1

		# Add trajectory to update the safe set and value function
		lmpc.addTrajectory(xcl, ucl)

	# =====================================================================================

	# plt.figure()
	xcl_T = np.array(xcl).T
	plt.plot(xcl_T[0,:], xcl_T[1,:], '-ob', label='sampled Safe Set')

	

	# ====================================================================================
	# Compute optimal solution by solving a FTOCP with long horizon
	# ====================================================================================
	N = 100 # Set a very long horizon to fake infinite time optimal control problem
	ftocp_opt = FTOCP(N, A, B, Q, R)
	ftocp_opt.solve(xcl[0])
	xOpt = ftocp_opt.xPred
	uOpt = ftocp_opt.uPred
	costOpt = lmpc.computeCost(xOpt.T.tolist(), uOpt.T.tolist())
	print("Optimal cost is: ", costOpt[0])
	# Store optimal solution in the lmpc object
	lmpc.optCost = costOpt[0]
	lmpc.xOpt    = xOpt

	# Save the lmpc object
	filename = 'lmpc_object.pkl'
	filehandler = open(filename, 'wb')
	pickle.dump(lmpc, filehandler)
	xcl = np.array(lmpc.SS[0]).T
	plt.plot(xOpt[0,:], xOpt[1,:], '--g', label='Initial feasible trajectory')
	plt.show()


if __name__== "__main__":
  main()
