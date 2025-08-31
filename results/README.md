This folder shows the main outputs from Sprint 1, where I modeled hydrogen diffusion through a stainless steel wall using Fick’s law and Sieverts’ law.  

The **concentration profile plot**  

![Concentration Profile](Concentration%20Profile.png)

 shows how hydrogen spreads through the wall. At the start (t=0s) the profile is flat, but after ~43,000s it settles into a straight line, which is the classic steady-state gradient.  

The **flux vs time plot** shows how the transport rate at the boundary changes over time. It starts at zero, rises steadily as hydrogen enters the wall, and then levels off toward equilibrium.  

I also included the **console log output**, which lists the parameters and checks used in this run. It shows the grid size, diffusion coefficient, boundary concentrations, and the stability condition (r=0.40 ≤ 0.5). The L² error compared to the analytical steady-state was ~2.5×10⁻¹⁶, confirming that the solver behaves as expected.  

Finally, there is a **raw_data.csv** file with the numerical values for concentration and flux. This makes it possible to reproduce the plots or use the data for further analysis.  

### Key outcomes  
- The solver produced stable results and matched the expected steady-state behavior.  
- The concentration profile evolves exactly as theory predicts.  
- Flux increases with time and trends toward a plateau.  
- These outputs give me a reliable baseline for Sprint 2, where I will add temperature and pressure effects.  
