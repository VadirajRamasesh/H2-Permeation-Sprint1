# H2-Permeation Sprint 1

This is my first sprint in a 3 step sprint process in building a hydrogen permeation simulator which would mimic the real world conditions sytems face.  
The idea is to start simple: **1-D hydrogen diffusion through a flat stainless steel wall**, using Fick’s law and Sieverts’ law as boundary conditions.  

---

### What I did
- Implemented an explicit finite-difference solver (FTCS).
- Applied Sieverts’ law at the boundaries (high pressure → vacuum).
- Checked stability (r ≤ 0.5) and compared with the analytic steady state.

---

### Results
The model runs and produces:
- Concentration profile across the wall.
- Flux vs time at the boundary.

 output:

<img width="697" height="207" alt="image" src="https://github.com/user-attachments/assets/230f84dc-0a92-4c90-9e5a-4dd14ba2eebc" />
<img width="553" height="412" alt="image" src="https://github.com/user-attachments/assets/ab54f49a-0030-48f5-824f-a54b70275c13" />
<img width="545" height="400" alt="image" src="https://github.com/user-attachments/assets/76fd2525-9fe1-4717-b9ff-aba8edb33232" />





---

### What’s next
- Sprint 2: add multi-layer barrier coatings.  
- Sprint 3: extend to fusion-relevant H–T systems and realistic geometries.  

---

*Author: Vadiraja ramasesh — M.Sc. Hydrogen Technology*

