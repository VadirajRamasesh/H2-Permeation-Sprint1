# H2 permeation through 316L – quick and dirty FTCS test
# summer holidy project 2025 – Vadiraj @ Rosenheim
# yeah I know FTCS is ancient but it's transparent and I wanted to see the boundary layer myself
# will replace with Crank-Nicolson later when I have time

import numpy as np
import matplotlib.pyplot as plt
import math

R = 8.314
phi0, e_phi = 2.81e-4, 6.227e4    # Reiter data for 316L
s0,  e_s    = 4.88e2,  8.65e3

def phi(T): return phi0 * math.exp(-e_phi/(R*T))
def S(T):   return s0   * math.exp(-e_s/(R*T))
def D(T):   return phi(T) / S(T)

def sieverts(T, p_pa):
    return S(T) * math.sqrt(max(p_pa, 0.0))   # no negative pressure nonsense

L   = 1e-3
nx  = 201
x   = np.linspace(0, L, nx)
dx  = x[1] - x[0]

T   = 298                     # room temp – diffusion is dog slow
pL  = 0.1e6                   # 0.1 MPa left side
pR  = 0.0

cL = sieverts(T, pL)
cR = sieverts(T, pR)
Dv = D(T)

# I want to see the first 5 µm clearly – diffusion length is tiny at 298 K
# yes I calculated it by hand on paper first

dt = 0.4 * dx*dx / Dv          # keep r ≤ 0.4 just to be safe
nsteps = int(np.ceil(3600 / dt)) + 1
t_end = 3600

C = np.zeros(nx)

# save some snapshots
snap_times = np.linspace(0, 3600, 7)
save_ids = [int(t/dt) for t in snap_times]

profiles = []
fluxL = []

for k in range(1, nsteps):
    Cn = C.copy()
    C[1:-1] = Cn[1:-1] + (Dv*dt/dx**2) * (Cn[2:] - 2*Cn[1:-1] + Cn[:-2])
    C[0] = cL
    C[-1] = cR
    
    J_left = -Dv * (C[1] - C[0]) / dx
    fluxL.append(J_left)
    
    if k in save_ids:
        profiles.append((k*dt, C.copy()))

# steady state for error check
C_ss = cL + (cR - cL) * (x/L)
L2 = np.sqrt(np.mean((C - C_ss)**2))

# plots – zoom first 5 µm because at 298 K nothing happens in the bulk
plt.figure(figsize=(8,5))
for t, prof in profiles:
    plt.plot(x*1e6, prof, label=f't = {int(t)} s')
plt.xlim(0, 5)
plt.xlabel('x [µm]')
plt.ylabel('c [mol/m³]')
plt.legend()
plt.tight_layout()
plt.savefig('H2_permeation_316L_298K_first5um.png', dpi=300)
plt.close()

print(f'D = {Dv:.2e} m²/s, r = {Dv*dt/dx**2:.3f}, L2 error = {L2:.2e}')
# L2 should be basically zero – if not I screwed up the boundary conditions again


