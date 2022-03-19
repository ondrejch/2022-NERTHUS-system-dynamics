# 2021-NERTHUS-core
Scripts for our 2021 NERTHUS neutronic and depletion paper

~~~~~~~~~~~~~~~~~~~~~~~~
~~~~ File Directory ~~~~

hsat.m			Function needed for initial conditions of saturated region in parameters.m
hsh.m			Function needed for initial conditions of superheated region in parameters.m
kin_dyn_edit.txt	Depletion data needed for read_depletion.m
Nerthus_sim.slx		MATLAB Simulink simulation file
parameters.m		Collects all parameters for the simulation. Is run by transient file.
read_depletion.m	Needed to for depletion dependence option
transients.m 		Designer for the transients. Run this before running simulation
XSteam.m 		MATLAB program to calculate specific heat, density, and enthalpy of water

* All of these files need to be in the same workspace
NOTE: Add an empty file named "Results" to the workspace. Resulting .mat data files will be deposited there
