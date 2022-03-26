import osqp
import numpy as np
import scipy as sp
from scipy import sparse

# Generate problem data
sp.random.seed(1)
m = 30
n = 20
Ad = sparse.random(m, n, density=0.7, format='csc')
b = np.random.randn(m)

# OSQP data
P = sparse.block_diag([sparse.csc_matrix((n, n)), sparse.eye(m)], format='csc')
q = np.zeros(n+m)
A = sparse.vstack([
        sparse.hstack([Ad, -sparse.eye(m)]),
        sparse.hstack([sparse.eye(n), sparse.csc_matrix((n, m))])], format='csc')
l = np.hstack([b, np.zeros(n)])
u = np.hstack([b, np.ones(n)])

# Create an OSQP object
prob = osqp.OSQP()

# Setup workspace
prob.setup(P, q, A, l, u)

# Solve problem
res = prob.solve()