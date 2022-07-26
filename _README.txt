~~~~~~~~~~~~~~~~~~~~~~~~
~~~~ File Directory ~~~~

flibe_PKPs.txt		FLiBe fuel salt depletion data needed for read_depletion.m
hsat.m			Function needed for initial conditions of saturated region in parameters.m
hsh.m			Function needed for initial conditions of superheated region in parameters.m
Nerthus_sim.slx		MATLAB Simulink simulation file
parameters.m		Collects all parameters for the simulation. It is run by transient file.
read_depletion.m	Needed to for depletion dependence option
thor_PKPs.txt		FNaBe LEU-Thorium fuel salt depletion data needed for read_depletion.m
transients.m 		Designer for the transients. Run this before running simulation
XSteam.m 		MATLAB program to calculate specific heat, density, and enthalpy of water

* All of these files need to be in the same workspace
~~~~~~~~~~~~~~~~~~~~~~~~
~~~~      TIPS      ~~~~

 - transients.m is designed to be the control panel
 - parameters.m is under the hood
 - Put a empty file named 'Results' in the directory. Data will be deposited there
 - If the simulation crashes after a successful run without any changes inbetween it is recommended to close the simulation, delete the .slxc file and the slprj file, then restart Simulink.
 - At very low to no xenon removal, use the saturation block. Controls are in parameters.m
~~~~~~~~~~~~~~~~~~~~~~~~

