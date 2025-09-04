# H₂ Permeation — Sprint 1

First sprint in a three-step project to build a hydrogen permeation simulator for stainless steel.  
This sprint implements a **1-D diffusion solver** (Fick’s law) with **Sieverts’ law** boundary conditions.

---

## About
- Material: 316/316L stainless steel  
- Physics: Fick’s 2nd law, constant-D; surface concentrations from **Sieverts’ law** (pressure in **Pa**)  
- Numerics: explicit **FTCS** scheme, with stability enforced (`r = D Δt / Δx² ≤ 0.5`)  
- Initial condition: **zero concentration** to show the transient; steady state is linear for constant D

**Two code versions (for transparency):**
- **Patched (validated):** [`code/diffusion_sprint1(patched).py`](code/diffusion_sprint1%28patched%29.py)  
  – fixes MPa→Pa, non-steady IC, and timestep control  
- **Original draft:** [`code/Sprint 1.py`](code/Sprint%201.py)  
  – first attempt kept for iteration history

---

## Results
See the folder: [`results/`](results)

- **Concentration profile:** [`results/profile.png`](results/profile.png)  
  Transient profiles evolving toward the expected **linear** steady-state gradient.
- **Flux vs time:** [`results/flux.png`](results/flux.png)  
  Flux rises from ~0 and approaches a plateau as steady diffusion sets in.

> Note: at **298 K**, hydrogen diffusivity in 316/316L is very small (≈10⁻¹⁶ m²/s), so in 1 h only a thin region near the entry face evolves—this is expected and reflected in the plots.

---

## How to run

### 1) Create a virtual environment (Python 3.12 recommended)
```bash
# Windows (PowerShell)
py -3.12 -m venv .venv
.\.venv\Scripts\activate

# macOS / Linux
python3 -m venv .venv
source .venv/bin/activate
