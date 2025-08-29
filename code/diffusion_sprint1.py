#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sprint 1 — 1D Hydrogen Diffusion in 316/316L Stainless Steel
Model: ∂C/∂t = ∂/∂x ( D(T) ∂C/∂x ),  with Dirichlet BCs from Sieverts’ law C = S(T)√p_H2.
Outputs: CSV, PNG, JSON under results/YYYYMMDD-HHMMSS.
"""

from __future__ import annotations
import argparse, json, math, datetime
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# ---------- constants ----------
R = 8.314462618  # J/(mol·K)

# ---------- material model ----------
class Material:
    """
    Φ(T) = Φ0 exp(-EΦ/RT)   [mol·m^-1·s^-1·MPa^-0.5]
    S(T) = S0 exp(-ES/RT)   [mol·m^-3·MPa^-0.5]
    D(T) = Φ/S              [m^2/s]
    Parameterization consistent with:
      - San Marchi, Somerday, Robinson, Int. J. Hydrogen Energy 32 (2007) 100–116.
      - Somerday et al., Sandia SAND2012-7321 (2012).
    """
    def __init__(self, name: str, phi0: float, E_phi: float, s0: float, E_s: float):
        self.name  = name
        self.phi0  = phi0
        self.E_phi = E_phi
        self.s0    = s0
        self.E_s   = E_s

    def phi(self, T: float) -> float:
        return self.phi0 * math.exp(-self.E_phi / (R * T))

    def S(self, T: float) -> float:
        return self.s0 * math.exp(-self.E_s / (R * T))

    def D(self, T: float) -> float:
        return self.phi(T) / self.S(T)

# 316/316L (MPa convention for Sieverts units)
SS316 = Material(
    name  = "316/316L Stainless Steel",
    phi0  = 2.81e-4,   # mol·m^-1·s^-1·MPa^-0.5
    E_phi = 62.27e3,   # J/mol
    s0    = 4.88e2,    # mol·m^-3·MPa^-0.5
    E_s   = 8.65e3     # J/mol
)

# ---------- numeric helpers ----------
def l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

# ---------- solvers ----------
def step_ftcs(C: np.ndarray, r: float, cL: float, cR: float) -> np.ndarray:
    Cn = C
    Cnew = Cn.copy()
    Cnew[1:-1] = Cn[1:-1] + r * (Cn[2:] - 2.0 * Cn[1:-1] + Cn[:-2])
    Cnew[0], Cnew[-1] = cL, cR
    return Cnew

def build_cn_matrices(Nx: int, r: float):
    # Interior nodes i = 1..Nx-2
    main = (1.0 + r) * np.ones(Nx - 2)
    off  = (-0.5 * r) * np.ones(Nx - 3)
    return main, off

def solve_tridiag(main: np.ndarray, off: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    # Thomas algorithm
    a = off.copy(); b = main.copy(); c = off.copy(); d = rhs.copy()
    n = b.size
    for i in range(1, n):
        w = a[i - 1] / b[i - 1]
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    x = np.empty(n)
    x[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    return x

def step_cn(C: np.ndarray, r: float, cL: float, cR: float, main: np.ndarray, off: np.ndarray) -> np.ndarray:
    Cn = C
    rhs = (1.0 - r) * Cn[1:-1] + 0.5 * r * (Cn[2:] + Cn[:-2])
    rhs[0]  += 0.5 * r * cL
    rhs[-1] += 0.5 * r * cR
    Cnew = Cn.copy()
    Cnew[1:-1] = solve_tridiag(main, off, rhs)
    Cnew[0], Cnew[-1] = cL, cR
    return Cnew

# ---------- simulation ----------
def simulate(
    L_m: float,
    Nx: int,
    T_K: float,
    p_left_MPa: float,
    p_right_MPa: float,
    t_end_s: float,
    dt_s: float,
    scheme: str,
    save_every: int,
    mat: Material
) -> dict:
    x = np.linspace(0.0, L_m, Nx)
    dx = x[1] - x[0]

    S_T = mat.S(T_K)      # mol·m^-3·MPa^-0.5
    D_T = mat.D(T_K)      # m^2/s
    cL  = S_T * math.sqrt(max(p_left_MPa,  0.0))  # mol·m^-3
    cR  = S_T * math.sqrt(max(p_right_MPa, 0.0))  # mol·m^-3

    C = np.linspace(cL, cR, Nx)  # linear IC

    # FTCS stability: r = D dt / dx^2 ≤ 0.5
    r = D_T * dt_s / (dx * dx)
    if scheme == "ftcs" and r > 0.5:
        dt_s = 0.5 * dx * dx / D_T
        r    = D_T * dt_s / (dx * dx)
        print(f"[stability] FTCS dt adjusted: dt={dt_s:.3e} s, r={r:.3f}")

    Nt   = int(math.ceil(t_end_s / dt_s))
    tvec = np.linspace(0.0, Nt * dt_s, Nt + 1)

    profiles = [C.copy()]
    times    = [0.0]
    flux     = np.zeros(Nt + 1)

    if scheme == "cn":
        main, off = build_cn_matrices(Nx, D_T * dt_s / (dx * dx))

    for n in range(1, Nt + 1):
        if scheme == "ftcs":
            C = step_ftcs(C, r, cL, cR)
        elif scheme == "cn":
            C = step_cn(C, D_T * dt_s / (dx * dx), cL, cR, main, off)
        else:
            raise ValueError("scheme ∈ {ftcs, cn}")

        # boundary flux at left face
        flux[n] = -D_T * (C[1] - C[0]) / dx

        if (n % save_every == 0) or (n == Nt):
            profiles.append(C.copy())
            times.append(tvec[n])

    # analytic steady state (constant D, Dirichlet)
    C_ss = cL + (cR - cL) * (x / L_m)
    err  = l2(C, C_ss)

    return dict(
        x=x, profiles=profiles, times=times, flux=flux, tvec=tvec,
        D=D_T, S=S_T, cL=cL, cR=cR,
        L=L_m, Nx=Nx, T=T_K,
        p_left_MPa=p_left_MPa, p_right_MPa=p_right_MPa,
        dt=dt_s, scheme=scheme, material=mat.name,
        C_final=C, C_ss=C_ss, l2_err=err
    )

# ---------- output ----------
def save_outputs(res: dict):
    outdir = ensure_dir(Path("results") / datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))

    # metadata
    meta = {
        "material": res["material"],
        "T_K": res["T"], "L_m": res["L"], "Nx": res["Nx"],
        "p_left_MPa": res["p_left_MPa"], "p_right_MPa": res["p_right_MPa"],
        "scheme": res["scheme"], "dt_s": res["dt"],
        "D_m2_s": res["D"], "S_mol_m3_MPa_m05": res["S"],
        "l2_final_vs_analytic": res["l2_err"]
    }
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2))

    # profiles.csv: x plus profiles at saved times
    hdr = "x_m," + ",".join([f"t={t:.6g}s" for t in res["times"]])
    arr = np.column_stack([res["x"]] + res["profiles"])
    np.savetxt(outdir / "profiles.csv", arr, delimiter=",", header=hdr, comments="")

    # flux_time.csv
    np.savetxt(
        outdir / "flux_time.csv",
        np.column_stack([res["tvec"], res["flux"]]),
        delimiter=",",
        header="time_s,boundary_flux_mol_m2_s",
        comments=""
    )

    # plots
    plt.figure()
    for k, Ck in enumerate(res["profiles"]):
        label = "t=0 s" if k == 0 else f"t={res['times'][k]:.3g} s"
        plt.plot(res["x"], Ck, label=label)
    plt.xlabel("x (m)"); plt.ylabel("C (mol m$^{-3}$)")
    plt.title("Hydrogen concentration profiles")
    plt.legend(); plt.tight_layout()
    plt.savefig(outdir / "profiles.png", dpi=200)

    plt.figure()
    plt.plot(res["tvec"], res["flux"])
    plt.xlabel("time (s)"); plt.ylabel("J_left (mol m$^{-2}$ s$^{-1}$)")
    plt.title("Boundary flux vs time")
    plt.tight_layout()
    plt.savefig(outdir / "flux_time.png", dpi=200)

    plt.figure()
    plt.plot(res["x"], res["C_final"], label="Final numeric")
    plt.plot(res["x"], res["C_ss"], "--", label="Analytic steady state")
    plt.xlabel("x (m)"); plt.ylabel("C (mol m$^{-3}$)")
    plt.title(f"Final vs analytic steady state (L2={res['l2_err']:.3e})")
    plt.legend(); plt.tight_layout()
    plt.savefig(outdir / "steady_state_check.png", dpi=200)

    print(f"[out] {outdir}")

# ---------- cli ----------
def build_parser():
    p = argparse.ArgumentParser(description="1D H diffusion (FTCS/CN) with Sieverts BCs")
    p.add_argument("--L", type=float, default=1e-3, help="Wall thickness [m]")
    p.add_argument("--Nx", type=int, default=201, help="Grid points")
    p.add_argument("--T", type=float, default=298.15, help="Temperature [K]")
    p.add_argument("--p_left", type=float, default=0.1, help="Left H2 pressure [MPa]")
    p.add_argument("--p_right", type=float, default=0.0, help="Right H2 pressure [MPa]")
    p.add_argument("--t_end", type=float, default=3600.0, help="Total simulated time [s]")
    p.add_argument("--dt", type=float, default=0.5, help="Time step [s] (guarded for FTCS)")
    p.add_argument("--scheme", choices=["ftcs", "cn"], default="cn", help="Numerical scheme")
    p.add_argument("--save_every", type=int, default=200, help="Save every N steps")
    return p

def main(argv=None):
    args = build_parser().parse_args(argv)
    res = simulate(
        L_m=args.L, Nx=args.Nx, T_K=args.T,
        p_left_MPa=args.p_left, p_right_MPa=args.p_right,
        t_end_s=args.t_end, dt_s=args.dt,
        scheme=args.scheme, save_every=args.save_every,
        mat=SS316
    )
    print(f"[run] {res['material']} | T={res['T']} K | D={res['D']:.3e} m^2/s | {res['scheme']} | L2={res['l2_err']:.3e}")
    save_outputs(res)

if __name__ == "__main__":
    main()

