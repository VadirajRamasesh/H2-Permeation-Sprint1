# H2-Permeation Sprint 1

This is my first sprint in a 3 step sprint process in building a hydrogen permeation simulator which would mimic the real world conditions sytems face.  
About this Project. In this sprint, I built 1D Python solver for diffusion using Fick’s law with Sieverts’ 
law as boundary conditions. This allowed me to simulate how hydrogen enters and travels through stainless steel over time. 

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
## Key Outcomes 

Produced concentration profiles and flux trends across the wall.
- Validated stability of the numerical solver against analytical steady state.
- - Gained practical experience with time-dependent monitoring and interpretation of diffusion processes.

### What’s next

- Sprint 2: Next, I plan to add temperature and pressure effects to make the model closer to real operating conditions.
- Sprint 3: Later, I’d like to try more complex geometries and connect diffusion with electrochemical effects.

---

*Author: Vadiraja ramasesh — M.Sc. Hydrogen Technology*

