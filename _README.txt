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
~~~~~~~~~~~~~~~~~~~~~~~~
~~~~      TIPS      ~~~~

 - If the simulation crashes after a successful run without any changes inbetween it is recommended to close the simulation, delete the .slxc file and the slprj file, then restart Simulink.
 - To run the simulation faster go into the OTSG, disconnect the feedwater PID from the MUX block, and attach the feedwater constant block. The feedwater PID is necessary for load following. 
 - When running at very low (<100MW) and very high (>1000MW) power, the feedwater PID has been known to cause issues. It is recommended that the feedwater constant block be used instead of the PID, or it may be necessary to recalibrate the PID as needed.
~~~~~~~~~~~~~~~~~~~~~~~~
