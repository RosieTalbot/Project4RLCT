# Project4RLCT
Numerical And Analytical Bifurcation Analysis of Swarm Models - Project 4 2024/5


This is the github for my Project supervised by Andrew Krause and Denis Patterson.

It includes my code for the Kuramoto Model and the Vicsek Model.

Kuramoto Model 
---------
A model looking at synchronised oscillators, r is the order parameter.

- KuramotoModel = The Kuramoto Model; produces graphs : r-t, circle plot, histogram of final omegas, omega-t.

- KuramotoKVals = Runs the Kuramoto Model over a range of K values; produces graphs : r-t overlays, r-K.

- KuramotoConstW = Runs the Kuramoto Model; produces graphs : Moving circle plot, r-t. Uses PlotCircle and PlotR.

- KuramotoWVals = Runs the Kuramoto Model over a range of W values; produces graphs : r-t overlays, r-W.

K - K is the coupling strength
- 0 - fully unsychronised
- 5 - synchronised (or should be)

Omega - The distribution of natural frequencies
- randn + pi
- rand
- zeros + pi
- ones 

Initial Conditions for the 
- randn + pi, mod 2pi
- rand*2pi, mod 2pi
- rand, mod 2pi

Vicsek Model
---------
A model looking at the directions of particles in a swarm, W is the order parameter

- VicsekModel = The Vicsek Model; produces graphs : Moving plot, stationary overlay plot.

- VicsekNuVals = Runs the Vicsek Model over a range of Nu values; produces graphs : W-Nu, W-t overlays.

Nu - The intrinsic noise of each particle 
- 0 - No noise, fully synchronised
- 5 - Lots of noise, unsynchronised 
 
