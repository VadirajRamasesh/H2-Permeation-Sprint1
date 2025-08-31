# H2 permeation – Sprint 1 (1D diffusion baseline)
# simple FTCS implementation with Sieverts' law boundary conditions

import numpy as np
import matplotlib.pyplot as plt
import math

# ---- material parameters for 316/316L stainless steel (from Sandia refs)
R = 8.314462618              # gas constant [J/mol/K]
PHI0, E_PHI = 2.81e-4, 62.27e3
S0,   E_S   = 4.88e2,  8.65e3

def phi_T(T):   return PHI0 * math.exp(-E_PHI/(R*T))
def S_T(T):     return S0   * math.exp(-E_S/(R*T))
def D_T(T):     return phi_T(T) / S_T(T)
def sieverts(T, p):          # surface concentration from Sieverts' law
    return S_T(T) * math.sqrt(max(p, 0.0))

# ---- case setup
L   = 1e-3     # wall thickness [m]
Nx  = 201      # grid points
T   = 298      # temperature [K]
pL  = 0.1      # left H2 pressure [MPa] ~1 bar
pR  = 0.0      # right H2 pressure [MPa] vacuum
t_end = 3600   # total time [s]
snaps = 6      # how many profiles to save/plot

# ---- mesh + material
x  = np.linspace(0, L, Nx)
dx = x[1]-x[0]
S  = S_T(T)
D  = D_T(T)
cL = sieverts(T, pL)
cR = sieverts(T, pR)

# FTCS stability: r <= 0.5
dt = 0.4 * dx*dx / D
r  = D*dt/(dx*dx)

# initial concentration profile (linear)
C = np.linspace(cL, cR, Nx)

nsteps   = int(np.ceil(t_end/dt))
save_ids = np.unique(np.round(np.linspace(0, nsteps, snaps)).astype(int))
profiles = [(0.0, C.copy())]
times    = np.linspace(0, nsteps*dt, nsteps+1)
fluxL    = np.zeros(nsteps+1)

def ftcs(C, r, cL, cR):
    Cn = C
    C2 = Cn.copy()
    C2[1:-1] = Cn[1:-1] + r*(Cn[2:] - 2*Cn[1:-1] + Cn[:-2])
    C2[0], C2[-1] = cL, cR
    return C2

# ---- time loop
for k in range(1, nsteps+1):
    C = ftcs(C, r, cL, cR)
    fluxL[k] = -D*(C[1]-C[0])/dx
    if k in save_ids:
        profiles.append((times[k], C.copy()))

# analytic steady state (linear, constant D)
C_ss  = cL + (cR - cL)*(x/L)
l2err = float(np.sqrt(np.mean((C - C_ss)**2)))

# ---- plots
plt.figure()
for t, Ck in profiles:
    plt.plot(x, Ck, label=f"t={t:.0f}s")
plt.xlabel("x [m]"); plt.ylabel("C [mol/m^3]")
plt.legend(); plt.tight_layout()
plt.savefig("profiles.png", dpi=200)

plt.figure()
plt.plot(times, fluxL)
plt.xlabel("time [s]"); plt.ylabel("J_left [mol/m^2/s]")
plt.tight_layout()
plt.savefig("flux_time.png", dpi=200)

# ---- summary printout
print("\n=== Sprint 1 (1D diffusion) ===")
print(f"T={T} K, L={L:.1e} m, Nx={Nx}, dx={dx:.2e}")
print(f"S={S:.2e}, D={D:.2e}")
print(f"cL={cL:.2e}, cR={cR:.2e}")
print(f"dt={dt:.2e} s, r={r:.2f}  (<=0.5 ok)")
print(f"L2 error vs linear SS = {l2err:.2e}")
print("obs: FTCS stable, profile → linear, flux → plateau")


