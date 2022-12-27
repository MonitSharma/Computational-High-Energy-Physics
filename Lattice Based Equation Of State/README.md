## Equation of State with Baryon Number, Electric charge and Strangeness Chemical Potential
## Lattice Based QCD
Quantum Chromodynamics Equation of State using all three conserved quantities of Baryon number, Electric charge and
Strangeness Chemical Potential.

#### INPUT :
We need to provide the lattice simulated data of the coefficients of the taylor expansion, to calculate the value of
pressure , that can be used to calculate other thermodynamic values. 

The required input is a .dat file name "coefficients.dat" which contains the parameters for the parameterization 
of the taylor coefficients up to an order 4, in total there are 22 coefficients.

#### OUTPUT: :
After running the code, it automatically creates a folder by the name of "coefficients_check", where all the 22 coefficients
are printed as functions of the temperature in 30 < T < 800 MeV range, also we get the first and second order derivatives 
with respect to the temperature. Then a folder named "thermodynamics" is created which contain many files, namely :

  "EoS_Taylor_AllMu.dat" : In this file all the thermodynamic quantities are calculated on the whole
    4D computational grid. All these thermodynamic quantities are Temperature , Baryonic chemical potential, Electric
    charge chemical potential , Strangeness Chemical potential, pressure, entropy, density, Baryon density,
    Strangeness density, charge density, energy density and speed of sound, in that order.

  "AllTherm_No_QS_Taylor.dat" : In this file all the thermodynamic quantities are calculated on a 2D computational
    grid, rather than a 4D one in the previous case. The 2D being Temperature and Baryonic chemical potential. Here
    we take the assumption of Electric charge chemical potential and Strangeness Chemical Potential to be vanishing.
    It also prints the same thermodynamic quantities as above and in the same order.

  "AllTherm_StrNeutr_muBTConst_Taylor.dat" : In this file the thermodynamic quantities are calculated for the
    T = 30 to 800 MeV range, having a constant values for mu_B/T. These values are {0.5,1,1.5,2,2.5,3} and also considers
    the case of Strangeness neutraility. It gives us mu_B/T,  Temperature , Baryonic chemical potential, Electric
    charge chemical potential , Strangeness Chemical potential, pressure, entropy, density, Baryon density,
    Strangeness density, charge density, energy density and speed of sound, in that order.

  "AllTherm_NoQS-muBTConst_Taylor.dat" : In this file, the thermodynamic quantities are calculated in the same range 
    as above, and also among the constant values of mu_B/T , but here we consider the case of Charge chemical potential and Strangeness 
    chemical potential to be zero. It also gives us the same value as above, and in the same order.


#### COMPILING :
The code is written in C and it requires basic C libraries to run like the <stdio.h> , <stdlib.h> ,<math.h>
,<string.h> ,<time.h>. Also the libraries like stat.h and unistd.h header files were included to get the functions like
mkdir() and chdir() to make and change directory in between of compilation of code.

Since a make file is already there, the code can be simply run by:

```bash
 make
 ./EoS_BQS coefficients.dat
 ```


where the EoS_BQS is the make file and the coefficients.dat is the data file contain the lattice coefficients.


 This is an implementation of the paper "Lattice-based Equation of State at finite Baryon number, Electric charge
 and Strangeness chemical potential" by J.Noronha-Hostler, P.Parotto, C. Ratti and J.M Stafford

 Original code of the paper can be found from the paper : https://inspirehep.net/record/1720588
 This work should be cited whenever using the code.
 I changed and modified the code for my Thesis work. In the reference I added some papers and work I read for the work.
