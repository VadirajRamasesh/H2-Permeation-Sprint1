# H2 permeation — Sprint 1 (1D, FTCS + Sieverts)

import numpy as np
import matplotlib.pyplot as plt
import math

# 316/316L (units: phi mol m^-1 s^-1 Pa^-0.5, S mol m^-3 Pa^-0.5)
R = 8.314
phi0, e_phi = 2.81e-4, 6.227e4
s0,  e_s    = 4.88e2,  8.65e3

def phi(T): return phi0 * math.exp(-e_phi/(R*T))
def S(T):   return s0   * math.exp(-e_s/(R*T))
def D(T):   return phi(T) / S(T)
def sieverts(T, p_pa):  # p in Pa
    return S(T) * math.sqrt(p_pa if p_pa > 0.0 else 0.0)

# case
L   = 1e-3
nx  = 201
T   = 298
pL_MPa, pR_MPa = 0.1, 0.0   # MPa
t_end, snaps = 3600, 6

# pressures to Pa
pL = pL_MPa * 1e6
pR = pR_MPa * 1e6

# mesh + material
x  = np.linspace(0.0, L, nx)
dx = x[1] - x[0]
Dv = D(T)
cL = sieverts(T, pL)
cR = sieverts(T, pR)

# timestep (stable + enough steps)
dt_stable = 0.4 * dx*dx / Dv
nsteps    = max(400, int(np.ceil(t_end / dt_stable)))
dt        = t_end / nsteps
r         = Dv * dt / (dx*dx)   # <= 0.5

# IC (non-steady)
C = np.zeros(nx)

# store a few profiles
save_ids = np.unique(np.round(np.linspace(0, nsteps, snaps)).astype(int))
times = np.linspace(0.0, t_end, nsteps+1)
profiles = [(0.0, C.copy())]
fluxL = np.zeros(nsteps+1)

# FTCS step
for k in range(1, nsteps+1):
    Cn = C.copy()
    C[1:-1] = Cn[1:-1] + r*(Cn[2:] - 2*Cn[1:-1] + Cn[:-2])
    C[0], C[-1] = cL, cR
    fluxL[k] = -Dv * (C[1] - C[0]) / dx
    if k in save_ids:
        profiles.append((times[k], C.copy()))

# steady state (constant D -> linear)
C_ss  = cL + (cR - cL)*(x/L)
l2err = float(np.sqrt(np.mean((C - C_ss)**2)))

# plots (save artifacts)
plt.figure()
for t, Ck in profiles:
    plt.plot(x, Ck, label=f"t={int(t)}s")
plt.xlabel("x [m]"); plt.ylabel("C [mol/m^3]"); plt.legend(); plt.tight_layout()
plt.xlim(0, 5e-6)  # zoom first 5 µm to make the boundary layer visible at 298 K
plt.savefig("profiles.png", dpi=200); plt.close()

plt.figure()
plt.plot(times, fluxL)
plt.xlabel("time [s]"); plt.ylabel("J_left [mol/m^2/s]"); plt.tight_layout()
plt.savefig("flux_time.png", dpi=200); plt.close()

print(f"r={r:.3f}, D={Dv:.3e} m^2/s, cL={cL:.3e}, cR={cR:.3e}, L2={l2err:.2e}")


