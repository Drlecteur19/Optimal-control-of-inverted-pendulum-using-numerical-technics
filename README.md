# Optimal-control-of-inverted-pendulum-using-numerical-technics
# Optimal Control of the Inverted Pendulum — Newton's Method for the Algebraic Riccati Equation

This project implements an **LQR controller** for the classic inverted pendulum system. Instead of using a built-in `lqr()` solver, the optimal gain matrix is computed by solving the **Algebraic Riccati Equation (ARE)** numerically via **Newton's Method** mater we will make more.

---

## Motivation

Standard control libraries (e.g., `scipy.linalg.solve_continuous_are`) solve the ARE internally using Schur decomposition. This project exposes the underlying mathematics step by step, making it a useful reference for students learning optimal control theory.

---

## System Model

The inverted pendulum is linearized around the unstable upright equilibrium and described by the standard state-space form:

$$\dot{x} = Ax + Bu$$

The state vector is:

$$x = \begin{bmatrix} \theta \\ \dot{\theta} \\ z \\ \dot{z} \end{bmatrix}$$

where $\theta$ is the pole angle, and $z$ is the cart position.

---

## LQR Problem

The objective is to minimize the infinite-horizon quadratic cost:

$$J = \int_0^{\infty} \left( x^T Q x + u^T R u \right) dt$$

The optimal control law is:

$$u^* = -Kx, \quad K = R^{-1} B^T P$$

where $P$ is the **symmetric positive definite** solution of the ARE.

---

## Algebraic Riccati Equation

$$A^T P + PA - PBR^{-1}B^T P + Q = 0$$

This is a **nonlinear matrix equation** in $P$. Newton's Method linearizes it at each iteration to find the solution efficiently.

---

## Newton's Method — Key Equations

At each iteration $k$, given the current estimate $P_k$, define the closed-loop matrix:

$$A_k = A - BR^{-1}B^T P_k$$

Then solve the **Lyapunov equation** for the update $\Delta P_k$:

$$A_k^T \Delta P_k + \Delta P_k A_k = -(A^T P_k + P_k A - P_k B R^{-1} B^T P_k + Q)$$

Update the solution:

$$P_{k+1} = P_k + \Delta P_k$$

Repeat until convergence:

$$\| P_{k+1} - P_k \|_F < \varepsilon$$

---

## Convergence & Numerical Notes

- Convergence is **quadratic** (the error squares at each step) once close to the solution.
- $P$ is re-symmetrized at every iteration to prevent floating-point drift: $P \leftarrow \frac{P + P^T}{2}$
- A stabilizing initial guess $P_0$ is required (e.g., $P_0 = Q$).
- Tolerance used: $\varepsilon = 10^{-8}$

---

## Project Structure

```
.
├── pendulum_model.py      # System matrices A, B
├── newton_riccati.py      # Newton's Method solver for the ARE
├── lqr_controller.py      # Computes gain K and simulates the system
├── main.py                # Entry point
└── README.md
```

---

## Requirements

```bash
pip install numpy scipy matplotlib
```

---

## Usage

```bash
python main.py
```

---

## References

- Lewis, F. L., Vrabie, D., & Syrmos, V. L. — *Optimal Control*, Wiley, 2012
- Khalil, H. K. — *Nonlinear Systems*, Prentice Hall, 2002
- Anderson & Moore — *Optimal Control: Linear Quadratic Methods*, Dover, 2007

