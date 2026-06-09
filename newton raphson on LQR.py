# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 19:53:48 2026
@author: Dr_lecteur
"""

# ── Libraries ────────────────────────────────────────────────────────────────
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.linalg import solve_continuous_lyapunov

# ── Physical Parameters ──────────────────────────────────────────────────────
m_p = 0.009   # pole mass  [kg]
M_c = 0.020   # cart mass  [kg]
g   = 9.81    # gravity    [m/s²]
l   = 0.25    # pole length [m]

# ── State-Space Model ────────────────────────────────────────────────────────
def inverted_pendulum_ss(M, m, l, g):
    """
    Returns (A, B) — linearized inverted pendulum
    around the upright equilibrium.
    State: [theta, theta_dot, z, z_dot]
    """
    A = np.array([
        [0,  1,               0,  0],
        [0,  0,          -m*g/M,  0],
        [0,  0,               0,  1],
        [0,  0,  (M+m)*g/(l*M),  0]
    ])

    B = np.array([
        [0       ],
        [1/M     ],
        [0       ],
        [-1/(M*l)]
    ])

    return A, B

A, B = inverted_pendulum_ss(M_c, m_p, l, g)
print("A =\n", A)
print("B =\n", B)

# ── LQR Weights ──────────────────────────────────────────────────────────────
Q = np.diag([10.0, 1.0, 10.0, 1.0])
R = np.array([[0.01]])

# ── Riccati Residual ─────────────────────────────────────────────────────────
def riccati_residual(P, A, B, Q, R):
    """
    Computes R(P) = A'P + PA - PBR⁻¹B'P + Q
    Should be zero at the solution.
    """
    R_inv = np.linalg.inv(R)
    return A.T @ P + P @ A - P @ B @ R_inv @ B.T @ P + Q

# ── Newton Step via Lyapunov Solve ───────────────────────────────────────────
def newton_step(P, A, B, Q, R):
    """
    One Newton iteration:
      1. Find closed-loop matrix A_k = A - B R⁻¹ B' P
      2. Solve Lyapunov equation for ΔP:
             A_k' ΔP + ΔP A_k = -F(P)
      3. Give a result as return P + ΔP
    """
    R_inv  = np.linalg.inv(R)
    A_cl   = A - B @ R_inv @ B.T @ P                    # closed-loop A
    F_P    = riccati_residual(P, A, B, Q, R)             # residual
    dP     = solve_continuous_lyapunov(A_cl.T, -F_P)     # Lyapunov solve
    return P + dP
print("THe iteration will be : /")
# ── Newton Solver for the ARE ────────────────────────────────────────────────
def newton_are(A, B, Q, R, tol=1e-8, max_iter=100):
    """
    Solves R(P): A'P + PA - PBR⁻¹B'P + Q = 0
    using Newton's Method.
    """
    P = Q.copy()   # initial guess

    for k in range(max_iter):
        P_new = newton_step(P, A, B, Q, R)
        P_new = (P_new + P_new.T) / 2   # enforce symmetry

        err = np.linalg.norm(P_new - P, 'fro')
        print(f"  Iteration {k+1:3d} | residual = {err:.2e}")

        if err < tol:
            print(f"\n✓ Converged in {k+1} iterations.")
            return P_new

        P = P_new

    print("⚠ Did not converge within max iterations.")
    return P

# ── Solve & Compute Gain ─────────────────────────────────────────────────────
P = newton_are(A, B, Q, R)
K = np.linalg.inv(R) @ B.T @ P

print("\nP =\n", P)
print("\nK =\n", K)
