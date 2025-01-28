import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.integrate import solve_ivp

# Problem parameters
L = 1.0  # Length of the domain
nx = 50  # Number of spatial points
dx = L / (nx - 1)  # Spatial resolution
D = 1.0  # Diffusion coefficient

# Create spatial grid
x = np.linspace(0, L, nx)

# Initial condition
T0 = np.sin(np.pi * x)

# Boundary conditions
T_left = 0
T_right = 0

# Time parameters for ODE solver
t_start = 0
t_end = 1
dt = 0.01


# 1. Linear Matrix Solver for Steady State
def linear_solver_steady_state():
    # Create matrix A for the discretized Laplacian
    diagonals = [
        -2 * np.ones(nx),
        np.ones(nx - 1),
        np.ones(nx - 1),
    ]
    A = diags(diagonals, [0, -1, 1], format="csr") / (dx ** 2)

    # Enforce boundary conditions by modifying A and b
    A = A.tolil()
    A[0, :] = 0
    A[0, 0] = 1
    A[-1, :] = 0
    A[-1, -1] = 1
    A = A.tocsr()

    b = np.zeros(nx)
    b[0] = T_left
    b[-1] = T_right

    # Solve A * T = b
    T_steady = spsolve(A, b)
    return T_steady


# 2. ODE Solver for Transient Behavior
def transient_solver(T0):
    # Define the RHS of the ODE system
    def rhs(t, T):
        T_new = np.zeros_like(T)
        T_new[1:-1] = D * (T[:-2] - 2 * T[1:-1] + T[2:]) / (dx ** 2)
        return T_new

    # Solve the ODE
    sol = solve_ivp(rhs, (t_start, t_end), T0, t_eval=np.arange(t_start, t_end, dt))
    return sol.t, sol.y


# Run solvers
T_steady = linear_solver_steady_state()
time, T_transient = transient_solver(T0)

# Plot results
plt.figure(figsize=(12, 6))

# Steady-state solution
plt.plot(x, T_steady, label="Steady-State (Linear Solver)", color="red", linewidth=2)

# Transient solution over time
#for i in range(0, len(time), len(time) // 5):
 #   plt.plot(x, T_transient[:, i], label=f"Transient t={time[i]:.2f}")

plt.xlabel("x")
plt.ylabel("Temperature")
plt.title("Comparison of Transient and Steady-State Solutions")
plt.legend()
plt.grid()
plt.show()