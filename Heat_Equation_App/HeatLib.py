import numpy as np
import scipy.linalg as lng
import numexpr as ne
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def Animation(dt, t, x, u):

    if dt < 0.001:
        t_for_animation = t[::10]
        u_for_animation = u[::10]
    else:
        t_for_animation = t
        u_for_animation = u

    fig, ax = plt.subplots()

    def animate(i):
        ax.clear()
        ax.plot(x, u_for_animation[i], color="blue")

        ymin, ymax = u.min(), u.max()

        if ymin == ymax:
            ymin -= 0.5
            ymax += 0.5

        ax.set_ylim([ymin, ymax])
        ax.set_xlim([x.min(), x.max()])
        ax.set_xlabel("x", fontsize=12)
        ax.set_title(f"Temperature — t = {t_for_animation[i]:.4f}",fontsize=14)
        ax.grid(True, linewidth=0.5, color="lightgray")

    ani = FuncAnimation(fig=fig,func=animate,frames=len(t_for_animation),interval=20,repeat=False)
    plt.show()

    return ani

    
# --- Central-difference spatial discretization using NumPy slicing. -------------

def CDS(r,n,u): return r*u[n,2:]+(1-2*r)*u[n,1:-1]+r*u[n,:-2]

# --------------------------------------------------------------------------------

# --- One FTCS Boundary formula --------------------------------------------------

def FTCS_boundary(u_boundary, u_adjacent, r, dx, a, k, f):

    return 2*r*u_adjacent + (1 - 2*r - 2*r*dx*(a/k))*u_boundary + 2*r*dx*(f/k)
    
# ---------------------------------------------------------------------------------

# --- Thomas Algorithm ------------------------------------------------------------

def thomas_solver(Aw, Ap, Ae, Q):
    
    n = len(Ap)
    Ap_star = Ap.copy()
    Q_star = Q.copy()

    for i in range(1, n):

        factor = Aw[i-1] / Ap_star[i-1]
        Ap_star[i] = Ap_star[i] - factor * Ae[i-1]
        Q_star[i] = Q_star[i] - (factor * Q_star[i - 1])

    solution = np.zeros(n)
    solution[-1] = Q_star[-1]/Ap_star[-1]

    for i in range(n - 2, -1, -1):
        solution[i] = (Q_star[i] - Ae[i] * solution[i + 1])/Ap_star[i]
        
    return solution

# --- Evaluate a user-defined expression -------------------------------

def evaluate_expression(expression, variable_name, variable):

    values = ne.evaluate(str(expression),local_dict={variable_name: variable,"pi": np.pi})

    # Handles constant expressions as well as arrays
    
    return np.broadcast_to(values, np.shape(variable)).astype(float).copy()

# --- Build boundary forcing arrays -----------------------------------

def boundary_forcing(t, initial_value, new_expression):

    # New prescribed boundary values for n >= 1
    f = evaluate_expression(new_expression, "t", t)

    # User-supplied initial boundary value at n = 0
    f[0] = float(initial_value)

    return f

# --- Forward Time - Central Difference Scheme ------------------------

def FTCS(TempFormula, a0, aN, f0_initial, fN_initial, f0, fN, k0, kN, dx, dt, L, T, a=1):

    # --- Input validation ---------------------------------------------

    if dx <= 0 or dt <= 0 or L <= 0 or T <= 0 or a <= 0: 
        raise ValueError("dx, dt, L, T and diffusivity must be positive.")

    if a0 == 0 and k0 == 0: raise ValueError("Invalid left boundary parameters.")

    if aN == 0 and kN == 0: raise ValueError("Invalid right boundary parameters.")

    # --- Grid ---------------------------------------------------------

    x = np.arange(0, L + dx/2, dx)
    t = np.arange(0, T + dt/2, dt)

    if len(x) < 3: raise ValueError("At least three spatial nodes are required.")

    if len(t) < 2: raise ValueError("At least two time levels are required.")

    if not np.isclose(x[-1], L): raise ValueError("L must be an integer multiple of dx.")

    if not np.isclose(t[-1], T): raise ValueError("T must be an integer multiple of dt.")

    r = a*dt / dx**2

    u = np.zeros((len(t), len(x)))

    # --- Initial temperature profile ----------------------------------

    u[0, :] = evaluate_expression(TempFormula, "x", x)

    # IMPORTANT:
    # The initial temperature profile is never overwritten.

    # --- Boundary forcing --------------------------------------------

    f0 = boundary_forcing(t, f0_initial, f0)
    fN = boundary_forcing(t, fN_initial, fN)

    # --- FTCS stability control ---------------------------------------

    r_max = 0.5

    if k0 != 0:

        if a0/k0 < 0:
            raise ValueError("Invalid left Robin coefficient ratio.")

        r_max = min(r_max, 1 / (2*(1 + dx*a0/k0)))

    if kN != 0:

        if aN/kN < 0:
            raise ValueError("Invalid right Robin coefficient ratio.")

        r_max = min(r_max,1 / (2*(1 + dx*aN/kN)))

    if r > r_max:
        raise ValueError(f"Unstable FTCS scheme: r = {r:.6f}, "f"maximum allowed r = {r_max:.6f}")

    # --- Time integration ---------------------------------------------

    for n in range(len(t)-1):

        # Interior nodes
        u[n+1, 1:-1] = CDS(r, n, u)

        # --- Left boundary --------------------------------------------

        if k0 == 0:
            # Dirichlet: new temperature at n+1
            u[n+1, 0] = f0[n+1] / a0

        else:
            # Neumann / Robin: forcing at n
            u[n+1, 0] = FTCS_boundary(u[n, 0],u[n, 1], r, dx, a0, k0, f0[n])

        # --- Right boundary -------------------------------------------

        if kN == 0:
            # Dirichlet: new temperature at n+1
            u[n+1, -1] = fN[n+1] / aN

        else:
            # Neumann / Robin: forcing at n
            u[n+1, -1] = FTCS_boundary(u[n, -1], u[n, -2], r, dx, aN, kN, fN[n])

    # --- Animation ----------------------------------------------------

    plt.close()
    ani = Animation(dt, t, x, u)

    return x, t, u

# --- Crank-Nicolson Scheme --------------------------------------------

def CrankNicolson(TempFormula, a0, aN, f0_initial, fN_initial, f0, fN, k0, kN, dx, dt, L, T, a=1):

    # --- Input validation ---------------------------------------------

    if dx <= 0 or dt <= 0 or L <= 0 or T <= 0 or a <= 0:
        raise ValueError("dx, dt, L, T and diffusivity must be positive.")

    if a0 == 0 and k0 == 0: raise ValueError("Invalid left boundary parameters.")

    if aN == 0 and kN == 0: raise ValueError("Invalid right boundary parameters.")

    # --- Grid ---------------------------------------------------------

    x = np.arange(0, L + dx/2, dx)
    t = np.arange(0, T + dt/2, dt)

    if len(x) < 3: raise ValueError("At least three spatial nodes are required.")

    if len(t) < 2: raise ValueError("At least two time levels are required.")

    if not np.isclose(x[-1], L): raise ValueError("L must be an integer multiple of dx.")

    if not np.isclose(t[-1], T): raise ValueError("T must be an integer multiple of dt.")

    r = a*dt / dx**2

    N = len(x)

    u = np.zeros((len(t), N))

    # --- Initial temperature profile ----------------------------------

    u[0, :] = evaluate_expression(TempFormula, "x", x)

    # The initial profile is preserved.

    # --- Boundary forcing arrays -------------------------------------

    f0 = boundary_forcing(t, f0_initial, f0)
    fN = boundary_forcing(t, fN_initial, fN)

    # --- Tridiagonal matrix A ----------------------------------------

    lower = -r * np.ones(N - 1)
    main = (2 + 2*r) * np.ones(N)
    upper = -r * np.ones(N - 1)

    # --- Left boundary coefficients ----------------------------------

    if k0 == 0:

        # Dirichlet: u[0] = f0/a0
        main[0] = 1.0
        upper[0] = 0.0

    else:

        # Neumann or Robin
        main[0] = 2 + 2*r + 2*r*dx*(a0/k0)
        upper[0] = -2*r

    # --- Right boundary coefficients ---------------------------------

    if kN == 0:

        # Dirichlet: u[-1] = fN/aN
        main[-1] = 1.0
        lower[-1] = 0.0

    else:

        # Neumann or Robin
        main[-1] = 2 + 2*r + 2*r*dx*(aN/kN)
        lower[-1] = -2*r

    # --- Time integration ---------------------------------------------

    for n in range(len(t) - 1):

        Q = np.zeros(N)

        # --- Interior nodes ------------------------------------------

        Q[1:-1] = (r*u[n, :-2] + (2 - 2*r)*u[n, 1:-1] + r*u[n, 2:])

        # --- Left boundary -------------------------------------------

        if k0 == 0:

            # Dirichlet at n+1
            Q[0] = f0[n+1] / a0

        else:

            # Neumann / Robin with trapezoidal forcing
            Q[0] = ((2 - 2*r - 2*r*dx*(a0/k0))*u[n, 0] + 2*r*u[n, 1] + 2*r*dx*(f0[n] + f0[n+1])/k0)

        # --- Right boundary ------------------------------------------

        if kN == 0:

            # Dirichlet at n+1
            Q[-1] = fN[n+1] / aN

        else:

            # Neumann / Robin with trapezoidal forcing
            Q[-1] = ((2 - 2*r - 2*r*dx*(aN/kN))*u[n, -1] + 2*r*u[n, -2] + 2*r*dx*(fN[n] + fN[n+1])/kN)

        # --- Solve tridiagonal system --------------------------------

        u[n+1, :] = thomas_solver( lower, main, upper, Q )

    # --- Animation ----------------------------------------------------

    plt.close()
    ani = Animation(dt, t, x, u)

    return x, t, u


if __name__ == "__main__":

    Initial_Temp = input("Enter initial temperature: ")

    f0_initial = float(input("Initial left boundary value (n=0): "))
    fN_initial = float(input("Initial right boundary value (n=0): "))

    f0 = input("New left boundary expression (n>=1): ")
    fN = input("New right boundary expression (n>=1): ")

    CrankNicolson(Initial_Temp,0,0,f0_initial,fN_initial,f0,fN,1, 1,0.1,0.001,1,1)

