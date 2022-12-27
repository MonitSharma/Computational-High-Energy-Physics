# eos_with_critical_point_qcd
Lattice QCD simulations and finding critical point in QCD transitions

This is an updated version of the BEST Collaboration 
program producing an EoS matching lattice QCD at muB=0, 
and containing a critical point in the 3D Ising universality 
class, in a parametrized form.
It allows for different choices of constraints on the 
shape of the critical line, which reduce the number of 
parameters, as well as the inclusion of no critical point,
corresponding to a Taylor expansion of lattice QCD result 
only.

Moreover, it allows for the choice of strangeness neutrality
conditions:
nS = 0 && nQ = 0.4 nB
or the case of only baryon number:
muS = muQ = 0

THE INPUT
The input required to run the program is:
- Input lattice data. They correspond to the 0th, 2nd 
	and 4th coefficients in the Taylor expansion of the 
	pressure. They are given both in the cases with and 
	without strangeness neutrality. They are contained 
	in the files:
	"Lattice_Data_5_821_dT1.dat"
	"Lattice_Data_5_821_dT1_SN.dat"
- Input HRG pressure data. This is a table containing 
	the values of the pressure, calculated making use of 
	a non-interacting Hadron Resonance Gas (HRG) model in 
	the computational range. The data are contained in the 
	files:
	"Press_HRG_MUB000601_T005300_dT1.dat".
	"Press_HRG_MUB000601_T005300_dT1_SN.dat".
- Input parameter file. This file contains all the information 
	needed to run the program. More follows.

THE BODY OF THE PROGRAM
NOTE: The computational grid of the program is T = 5 - 821 MeV, 
			muB = 0 - 601 MeV, with grid size of 1 MeV. The output of 
			the program is in the range T = 30 - 800 MeV, 
			muB = 0 - 450 MeV. However, it can be modified within this 
			range, and extended further to muB = 600 MeV.
The program goes through the following steps:<br />

0A.) 	Read options<br /><br />
0B.) 	Read parameter file<br /><br />
1.) 	Import lattice data, multiply by T^4 and store. <br />
NOTE: throughout the program, all quantities are dimensional. 
			At the end, thermodynamic quantities exported are 
			normalized with powers of the temperature.<br /><br />
2.) 	Given the choice of parameters, a coordinates file is created. 
			It contains the map between (T,muB) and (R,theta). <br /><br /> 
4.) 	Using the coordinates file, the Ising contribution to the 
			expansion coefficients is calculated.<br /><br />
5.) 	Using the coordinate file, the 3D Ising pressure and its derivatives 
			are calculated at all points in the grid.<br />
NOTE: the calculation makes use of coordinate points with negative 
			muB too, due to the symmetrization procedure described in 
			the paper.<br /><br />
5.b)  The correlation length is calculated over the grid.<br /><br />
6.) 	The Non-Ising contribution to the expansion coefficients is 
			calculated from the lattice coefficients and the Ising ones.<br /><br />
7.) 	The Non-Ising contribution to the pressure is calculated in 
			all points in the grid, as a Taylor expansion of the Non-Ising 
			coefficients.<br /><br />
8.) 	The Non-Ising pressure is filtered through a Gaussian filter. 
			This is done in order to reduce unphysical wiggles in the 
			thermodynamic functions at largu muB due to the truncation of 
			the Taylor series. The parameters of the filter do not need to 
			be entered manually.<br /><br />
9.) 	Numerical derivatives of the Non-Ising pressure are calculated.<br /><br />
10.)	Ising and Non-Ising contributions to the pressure ad its 
			derivatives are combined (summed).<br /><br />
11.)	HRG pressure is imported and multiplied by T^4. Numerical 
			derivatives of the pressure are calculated.<br /><br />
12.) 	HRG and Ising+Non-Ising contribution to the pressure and 
		 	derivatives are combined via a smooth merging. The merging 
			makes use of an hyperbolic tangent. The parameters of the 
			merging do not need to be input manually.<br /><br />
13.)	Once derivatives are defined, they are combined and normalized 
			to generate the thermodynamics. Pressure, entropy density, 
			baryon density, energy density, speed of sound and second baryon 
			number cumulant are exported in the specified range. They are all 
			normalized with powers of the temperature.
			In this new version, we export the equilibrium correlation length
			as well, which is defined up to a constant.<br /><br />


RUNNING THE CODE
###                                                            ###	
### This is one of the major edits made in this newer version. ###
###                                                            ###	
Simply run the code with:
./EoS [options] <parameter_file>


THE OPTIONS
###                                                            ###	
### This is one of the major edits made in this newer version. ###
###                                                            ###	
 -r<lowT>:<highT>    Change T-range within T = 30 - 800 MeV
 -R<lowMU>:<highMU>  Change mu-range within muB = 0 - 600 MeV
 -L                  Lattice mode: in this case the equation
                     of state is produced without a critical point,
                     and the files copied in directory "LATonly"
 -S                  Consider the strangeness neutral case. In this case
                     the directory created in "output" has "SN" appended
                     to its name
 -p                  Outputs all the data from intermediate steps too
 -v                  Verbose, prints messages regarding the stage of 
										 running reached
 -h                  Print help


THE PARAMETER FILE
###                                                            ###	
### This is one of the major edits made in this newer version. ###
###                                                            ###	

The parameter file is in the format:
Lattice_chi=/path/to/input/Lattice_Data_5_821_dT1.dat
Lattice_chi_SN=/path/to/input/Lattice_Data_5_821_dT1_SN.dat
HRG_press=/path/to/input/Press_HRG_MUB000601_T005300_dT1.dat
HRG_press_SN=/path/to/input/Press_HRG_MUB000601_T005300_dT1_SN.dat
Mode=<mode>
T0=<value>
kappa=<value>
TC=<value>
angle1=<value>
muBC=<value>
anglediff=<value>
ww=<value>
rho=<value>


Where:
Lattice_chi points to the file Lattice_Data_5_821_dT1.dat
Lattice_chi_SN points to the file Lattice_Data_5_821_dT1_SN.dat
HRG_press points to the file Press_HRG_MUB000601_T005300_dT1.dat
HRG_press_SN points to the file Press_HRG_MUB000601_T005300_dT1_SN.dat
Mode is FREE, PAR, or STR
T0 is chiral transition temperature at muB=0 
kappa is curvature of chiral transition line at muB=0
TC is critical temperature (for mode FREE)
angle1 is angle alpha1 
muBC is critical muB
anglediff is alpha2-alpha1 (in degrees)
ww is value of w
rho is value of rho

Not all entries are necessary, depending on the MODE. 
For example, TC and angle1 are not needed in Mode "PAR".

1.) "Mode": FREE
		This means the choice of the six parameters in the map is 
		completely free. In this case, the additional parameters 
		needed are: 
		TC, the temperature at the critical point; 
		muBC, the chemical potential at the critical point; 
		angle1, the angle between the h=0 axis and the T=0 axis (in degrees); 
		angle2, the angle between the r=0 axis and the T=0 axis (in degrees); 
		w, the global scaling parameter in the mapping; 
		rho, the relative scaling in the mapping. 
2.) "Mode": PAR
		This means the critical point will lie on a parabola. In 
		this case, the additional parameters needed are: 
		T0, the value of temperature at which the parabolic pseudo-critical 
				line crosses the T axis; 
		kappa, the curvature of the transition line at the T axis; 
		muBC, the chemical potential at the critical point: since 
				both TC and angle1 are determined by this choice, they do 
				not need to be input; 
		anglediff, the difference between angle2 and angle1 (in degrees); 
		w, rho as before. 
		NOTE: the values of T0 and kappa currently used are from http://inspirehep.net/record/1385115.
3.)	"Mode": STR
		This means the critical point will lie on a straight line. In 
		this case, the additional parameters needed are: 
		T0, the value of temperature at which the straight pseudo-critical 
				line crosses the T axis; 
		muBC, the chemical potential at the critical point; 
		angle1, the angle between the h=0 axis and the T=0 axis (in degrees); 
		angle2, the angle between the r=0 axis and the T=0 axis (in degrees); 
		w, the global scaling parameter in the mapping; 
		rho, the relative scaling in the mapping. 



NOTE ON THE PARAMETERS
This project is intended for use in the framework of the BES-II. 
The program has been tested for "realistic" values of the parameters.
For the location of the critical point, the code is intended to 
place it in the chemical potential range spanned by the BES-II 
program (muB < 500 MeV).
The values of the angles should not be multiples of 90 degrees. 
This would prevent the critical contribution from the Ising model 
to be reflected in the QCD thermodynamical quantities that are 
calculated.
The code was tested for values of w = 0.1 - 20, rho = 0.1 - 20. Larger 
values of w should not result in problems; however, the larger w is, 
the more the critical contribution is "washed away". 

NOTE ON THE FILE
The public version of the code does not contain the tables for the 
different thermodynamic quantities. In order to generate them, it is
sufficient to change the paths to your directories in the parameter 
file (input/param.example), set the parameters to the desired values, 
then run the code.

COMPILING
Simply running "make" in the main directory should not
result in problems.



WHEN USING THIS CODE
When using this code, the following works should be cited:
https://inspirehep.net/literature/1672952
https://inspirehep.net/literature/1851774

