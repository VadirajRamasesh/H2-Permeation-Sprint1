# H2-Permeation-Sprint1

### Project Goal
Hydrogen always escapes. This project builds a **1D permeation simulator** to model hydrogen diffusion through stainless steel from first principles.

---

### What’s Inside
- **Derivations (`/docs`)** — Fick’s Laws + Sieverts’ Law, step-by-step
- **Code (`/code`)** — Python finite-difference (FTCS) implementation
- **Results (`/results`)** — Plots of concentration vs depth/time

---

### Physics & Numerics
- Fick’s 2nd Law: ∂C/∂t = D ∂²C/∂x²  
- Sieverts’ Law: C ∝ √p_H₂  
- Scheme: **FTCS** (transient diffusion), steady-state check

---

### Sprint 1 Scope
- 1D flat wall, thickness = **1 mm**  
- Boundary conditions: high H₂ pressure on one side, **vacuum** on the other  
- Output: concentration profile vs x and t

---

### Next
- Sprint 2: Multi-layer barrier coatings  
- Sprint 3: Fusion-relevant H–T systems & realistic geometries

---

**Author:** Vadiraj Br — M.Sc. Hydrogen Technology,

