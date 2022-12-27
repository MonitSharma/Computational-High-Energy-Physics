/*
    This is the main file where everything happens.
	This file produces a Taylor expanded EoS using Lattice QCD coefficients at muB=0.
*/ 

#define NRANSI

/* Import standard libraries, be sure to check they are all included in the path */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include "../include/nrD.h"
#include "../include/nrutilD.h"

/* Import additional library for many variables used in the program. */
#include "../include/Variables.h"
#include "../include/Functions.h"

/* Strangeness neutrality. */
#define NSN 2

/* Time variables for timing. */
clock_t start, end;
double cpu_time_used;

/* Functions whose zeroes we are seeking for strangeness neutrality. */
void funcvSN(int n,double x[],double f[]){
	f[1] = StrDensTaylor(Tval,muBval,x[1],x[2]);
	f[2] = ChDensTaylor(Tval,muBval,x[1],x[2]) - 0.4*BarDensTaylor(Tval,muBval,x[1],x[2]);
} 


/* The main body of the program. */
int main(int argc, char *argv[])
{
	char buff[FILENAME_MAX];

	/* Define time variables to measure time elapsed creating the tables.*/
	time_t start, stop;

	/* All vectors and matrices are initialized. */
	/* Vectors for coefficients. */
	// Order 0
	CHI000PAR=vector(1,21);
	// Order 2
	// Diagonal 
	CHI200PAR=vector(1,21); CHI020PAR=vector(1,21); CHI002PAR=vector(1,21); 
	// Mixed
	CHI110PAR=vector(1,21);	CHI101PAR=vector(1,21);	CHI011PAR=vector(1,21);	
	// Order 4
	// Diagonal
	CHI400PAR=vector(1,21);	CHI040PAR=vector(1,21);	CHI004PAR=vector(1,21);
	// Mixed 31
	CHI310PAR=vector(1,21);	CHI301PAR=vector(1,21);	CHI031PAR=vector(1,21);	
	// Mixed 13
	CHI130PAR=vector(1,21);	CHI103PAR=vector(1,21); CHI013PAR=vector(1,21);	
   // Mixed 22
   CHI220PAR=vector(1,21); CHI202PAR=vector(1,21); CHI022PAR=vector(1,21);	
   // Mied 112
   CHI211PAR=vector(1,21);	CHI121PAR=vector(1,21); CHI112PAR=vector(1,21);

	/* Matrix for coefficients: 22 coefficients, 21 parameters each. */
	parMatrix=matrix(1,22,1,21);
	
	/* For strangeness neutrality */
	xSN=vector(1,NSN);
 	fSN=vector(1,NSN);
	
	/* Assign the name of the main folder where the program lives and the files we wish to import are located.*/
	getcwd(buff,FILENAME_MAX);
	printf("Current working directory is: \n\n%s \n\n",buff);

   /* Start clock for importing lists */
  	start = clock();

	/* Parametrization parameters are read from the user-input file, and saved. */
  	FILE *ParametersIn = fopen(argv[1], "r");
  	if (ParametersIn == 0){
  		fprintf(stderr,"failed to open paremeters file\n");
  		exit(1);
  	}
  	for(i=1;fscanf(ParametersIn,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &A0, &A1, &A2, &A3, &A4, &A5, &A6, &A7, 
  	                                                                                &A8, &A9, &B0, &B1, &B2, &B3, &B4, &B5, &B6, &B7, &B8, &B9, &C0) !=EOF; i++){
       parMatrix[i][1] = A0;  	    parMatrix[i][2] = A1;  	    parMatrix[i][3] = A2;  	    parMatrix[i][4] = A3;
  	    parMatrix[i][5] = A4;  	    parMatrix[i][6] = A5;  	    parMatrix[i][7] = A6;  	    parMatrix[i][8] = A7;
  	    parMatrix[i][9] = A8;  	    parMatrix[i][10] = A9;
  	    parMatrix[i][11] = B0;  	    parMatrix[i][12] = B1;  	    parMatrix[i][13] = B2;  	    parMatrix[i][14] = B3;
  	    parMatrix[i][15] = B4;  	    parMatrix[i][16] = B5;  	    parMatrix[i][17] = B6;  	    parMatrix[i][18] = B7;
  	    parMatrix[i][19] = B8;  	    parMatrix[i][20] = B9;
  	    parMatrix[i][21] = C0;
  	    
  	    //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",A0,A1,A2,A3,A4,A5,A6,A7,B0,B1,B2,B3,B4,B5,B6,B7,C0);
  	}
  	fclose(ParametersIn);	
  	
 
	// Vectors for coefficients.   
	// All vectors are filled with parameters. 
	// Order 0
  	for(i=1;i<=21;i++) CHI000PAR[i] = parMatrix[1][i];
  	// Order 2
	// Diagonal 
  	for(i=1;i<=21;i++) CHI200PAR[i] = parMatrix[2][i];
  	for(i=1;i<=21;i++) CHI020PAR[i] = parMatrix[3][i];
	for(i=1;i<=21;i++) CHI002PAR[i] = parMatrix[4][i];
	// Mixed
	for(i=1;i<=21;i++) CHI110PAR[i] = parMatrix[5][i];
  	for(i=1;i<=21;i++) CHI101PAR[i] = parMatrix[6][i];
  	for(i=1;i<=21;i++) CHI011PAR[i] = parMatrix[7][i];
	// Order 4
   	// Diagonal	
   	for(i=1;i<=21;i++) CHI400PAR[i] = parMatrix[8][i];
  	for(i=1;i<=21;i++) CHI040PAR[i] = parMatrix[9][i];
  	for(i=1;i<=21;i++) CHI004PAR[i] = parMatrix[10][i];
  	// Mixed 31
  	for(i=1;i<=21;i++) CHI310PAR[i] = parMatrix[11][i];
  	for(i=1;i<=21;i++) CHI301PAR[i] = parMatrix[12][i];
  	for(i=1;i<=21;i++) CHI031PAR[i] = parMatrix[13][i];
  	// Mixed 13
  	for(i=1;i<=21;i++) CHI130PAR[i] = parMatrix[14][i];
  	for(i=1;i<=21;i++) CHI103PAR[i] = parMatrix[15][i];
  	for(i=1;i<=21;i++) CHI013PAR[i] = parMatrix[16][i];
  	// Mixed 22
  	for(i=1;i<=21;i++) CHI220PAR[i] = parMatrix[17][i];
  	for(i=1;i<=21;i++) CHI202PAR[i] = parMatrix[18][i];
  	for(i=1;i<=21;i++) CHI022PAR[i] = parMatrix[19][i];
  	// Mixed 211
  	for(i=1;i<=21;i++) CHI211PAR[i] = parMatrix[20][i];
  	for(i=1;i<=21;i++) CHI121PAR[i] = parMatrix[21][i];
  	for(i=1;i<=21;i++) CHI112PAR[i] = parMatrix[22][i];
  	
  	
  	
  	/* Create folder for coefficients checks. */
	mkdir("Coefficients_Check", S_IRWXU | S_IRWXG | S_IRWXO);
	chdir("Coefficients_Check");
  	
  	// Print values of all coefficients and derivatives thereof wrt T. (To check that everything is in order.)
  	FILE *CHIS = fopen("All_Chis.dat","w");
  	FILE *DCHISDT = fopen("All_DChisDT.dat","w");
  	FILE *D2CHISDT2 = fopen("All_D2ChisDT2.dat","w");
  	for(i=30;i<=800;i+=1){
  	    fprintf(CHIS,"%d    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf\n",
  	                i,CHI000(i),CHI200(i),CHI020(i),CHI002(i),CHI110(i),CHI101(i),CHI011(i),CHI400(i),CHI040(i),CHI004(i),CHI310(i),CHI301(i),CHI031(i),CHI130(i),CHI103(i),CHI013(i),
  	                CHI220(i),CHI202(i),CHI022(i),CHI211(i),CHI121(i),CHI112(i));
  	    fprintf(DCHISDT,"%d    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf\n",
  	                i,DCHI000DT(i),DCHI200DT(i),DCHI020DT(i),DCHI002DT(i),DCHI110DT(i),DCHI101DT(i),DCHI011DT(i),DCHI400DT(i),
  	                DCHI040DT(i),DCHI004DT(i),DCHI310DT(i),DCHI301DT(i),DCHI031DT(i),DCHI130DT(i),DCHI103DT(i),DCHI013DT(i),
  	                DCHI220DT(i),DCHI202DT(i),DCHI022DT(i),DCHI211DT(i),DCHI121DT(i),DCHI112DT(i));
  	    fprintf(D2CHISDT2,"%d    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf\n",
 	                i,D2CHI000DT2(i),D2CHI200DT2(i),D2CHI020DT2(i),D2CHI002DT2(i),D2CHI110DT2(i),D2CHI101DT2(i),D2CHI011DT2(i),D2CHI400DT2(i),
  	                D2CHI040DT2(i),D2CHI004DT2(i),D2CHI310DT2(i),D2CHI301DT2(i),D2CHI031DT2(i),D2CHI130DT2(i),D2CHI103DT2(i),D2CHI013DT2(i),
  	                D2CHI220DT2(i),D2CHI202DT2(i),D2CHI022DT2(i),D2CHI211DT2(i),D2CHI121DT2(i),D2CHI112DT2(i));            
  	}
  	fclose(CHIS);  	fclose(DCHISDT);  fclose(D2CHISDT2);

	chdir(buff);
  	
  	
  	/* Create folder for thermodynamic quantities. */
	mkdir("Thermodynamics", S_IRWXU | S_IRWXG | S_IRWXO);
	chdir("Thermodynamics");
  	
	/* (Unconstrained) thermodynamics for all T, muB, muS, muQ. */  	
  	FILE *All_Therm_Taylor = fopen("EoS_Taylor_AllMu.dat","w");
  	for(i=30;i<=800;i+=5){
  	    for(j=0;j<=450;j+=10){
      	    for(k=0;k<=450;k+=10){
   	      	    for(l=0;l<=450;l+=10){
                	Tval = i; muBval = j;  muQval = k; muSval = l;
                   
                	PressVal = PressTaylor(i,j,k,l);
   	            	EntrVal = EntrTaylor(i,j,k,l);
         	    	BarDensVal = BarDensTaylor(i,j,k,l);
	            	StrDensVal = StrDensTaylor(i,j,k,l);
	            	ChDensVal = ChDensTaylor(i,j,k,l);
	            	EnerDensVal = EntrVal - PressVal + muBval/Tval*BarDensVal + muQval/Tval*ChDensVal+ muSval/Tval*StrDensVal;
  	          	    SpSoundVal = SpSound(Tval,muBval,muQval,muSval);
	                
	                fprintf(All_Therm_Taylor,"%lf  %lf  %lf  %lf  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f\n", Tval, muBval, muQval, muSval, PressVal, EntrVal, 
	                        BarDensVal, StrDensVal, ChDensVal, EnerDensVal, SpSoundVal);                 
           	    }
            }	    
	    }
	}    
	fclose(All_Therm_Taylor);
	

	/* (Unconstrained) thermodynamics for muS = muQ = 0. */
  	FILE *AllTherm_No_QS = fopen("AllTherm_No_QS_Taylor.dat","w");
  	for(i=30;i<=800;i+=1){
  	    for(j=0;j<=450;j+=1){
  	        k = 0; l = 0;
  	        Tval = i; muBval = j;
  	        muQval = k; muSval = l;
  	        
  	        PressVal = PressTaylor(i,j,k,l);
   	     	EntrVal = EntrTaylor(i,j,k,l);
	        BarDensVal = BarDensTaylor(i,j,k,l);
	        StrDensVal = StrDensTaylor(i,j,k,l);
	        ChDensVal = ChDensTaylor(i,j,k,l);
	        EnerDensVal = EntrVal - PressVal + muBval/Tval*BarDensVal + muQval/Tval*ChDensVal+ muSval/Tval*StrDensVal;
  	        SpSoundVal = SpSound(Tval,muBval,muQval,muSval);
	        
	        fprintf(AllTherm_No_QS,"%d  %d  %d  %d  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f\n", i, j, k, l, PressVal, EntrVal, 
	                        BarDensVal, StrDensVal, ChDensVal, EnerDensVal, SpSoundVal);
	    }
	}    
	fclose(AllTherm_No_QS);
	

	/* Thermodynamics for Strangeness neutrality. */
  	FILE *AllTherm_StrNeutr = fopen("AllTherm_StrNeutr_Taylor.dat","w");
  	for(i=30;i<=800;i+=1){
  		for(j=0;j<=450;j+=1){
  			Tval = i; muBval = j;
	    	xSN[1] = 0.0; xSN[2] = 0.0;	
	    	newt(xSN,NSN,&check,funcvSN);
			funcvSN(NSN,xSN,fSN);
	    	muQval = xSN[1]; muSval = xSN[2];
	        
	    	PressVal = PressTaylor(Tval,muBval,muQval,muSval);
   	   		EntrVal = EntrTaylor(Tval,muBval,muQval,muSval);
	    	BarDensVal = BarDensTaylor(Tval,muBval,muQval,muSval);
	    	StrDensVal = StrDensTaylor(Tval,muBval,muQval,muSval);
	    	ChDensVal = ChDensTaylor(Tval,muBval,muQval,muSval);
	   		EnerDensVal = EntrVal - PressVal + muBval/Tval*BarDensVal + muQval/Tval*ChDensVal+ muSval/Tval*StrDensVal;	
	    	SpSoundVal = SpSound(Tval,muBval,muQval,muSval);
	   
	    	fprintf(AllTherm_StrNeutr,"%f  %f  %f  %f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f\n", Tval, muBval, muQval, muSval, PressVal, EntrVal, 
	                        BarDensVal, StrDensVal, ChDensVal, EnerDensVal, SpSoundVal);
	    }
	}    
	fclose(AllTherm_StrNeutr);
	

	/* Thermodynamics for muB/T = constant and Strangeness Neutrality. */
	FILE *AllTherm_StrNeutr_muBTConst = fopen("AllTherm_StrNeutr_muBTConst_Taylor.dat","w");
	for(Rat=0.5;Rat<=3;Rat+=0.5){  
	   	for(i=30;i<=800;i+=1){
	  		if(Rat*i < 451){    
  	    		Tval = i; muBval = Rat*Tval;
	        	xSN[1] = 0.0; xSN[2] = 0.0;	
	        	newt(xSN,NSN,&check,funcvSN);
		    	funcvSN(NSN,xSN,fSN);
	        	if (check) printf("Convergence problems.\n");		        
		    	muQval = xSN[1]; muSval = xSN[2];
	            
	        	PressVal = PressTaylor(Tval,muBval,muQval,muSval);
   	      		EntrVal = EntrTaylor(Tval,muBval,muQval,muSval);
	        	BarDensVal = BarDensTaylor(Tval,muBval,muQval,muSval);
	        	StrDensVal = StrDensTaylor(Tval,muBval,muQval,muSval);
	        	ChDensVal = ChDensTaylor(Tval,muBval,muQval,muSval);
	        	EnerDensVal = EntrVal - PressVal + muBval/Tval*BarDensVal + muQval/Tval*ChDensVal+ muSval/Tval*StrDensVal;
	        	SpSoundVal = SpSound(Tval,muBval,muQval,muSval);
	         
	        	fprintf(AllTherm_StrNeutr_muBTConst,"%f %f  %f  %f  %f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f\n", Rat, Tval, muBval, muQval, muSval, PressVal, EntrVal, 
	                        BarDensVal, StrDensVal, ChDensVal, EnerDensVal, SpSoundVal);
	        }
	    }
	}    
	fclose(AllTherm_StrNeutr_muBTConst);

	
	/* Thermodynamics for MuB/T = constant and muS = muQ = 0. */
  	FILE *AllTherm_NoQS_muBTConst = fopen("AllTherm_NoQS_muBTConst_Taylor.dat","w");
	for(Rat=0.5;Rat<=3;Rat+=0.5){  
	   	for(i=30;i<=800;i+=1){
    		if(Rat*i < 451){    
  	        	Tval = i; muBval = Rat*Tval;
	        	muQval = 0.0; muSval = 0.0;
	            
	        	PressVal = PressTaylor(Tval,muBval,muQval,muSval);
   	      		EntrVal = EntrTaylor(Tval,muBval,muQval,muSval);
	        	BarDensVal = BarDensTaylor(Tval,muBval,muQval,muSval);
	        	StrDensVal = StrDensTaylor(Tval,muBval,muQval,muSval);
	        	ChDensVal = ChDensTaylor(Tval,muBval,muQval,muSval);
	        	EnerDensVal = EntrVal - PressVal + muBval/Tval*BarDensVal + muQval/Tval*ChDensVal+ muSval/Tval*StrDensVal;
	        	SpSoundVal = SpSound(Tval,muBval,muQval,muSval);
	         
	        	fprintf(AllTherm_NoQS_muBTConst,"%lf %lf  %lf  %lf  %lf  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f\n", Rat, Tval, muBval, muQval, muSval, PressVal, EntrVal, 
	                        BarDensVal, StrDensVal, ChDensVal, EnerDensVal, SpSoundVal);
	        }
	    }
	}    
	fclose(AllTherm_NoQS_muBTConst);
	
	chdir(buff);

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Finished calculating in %lf seconds.\n", cpu_time_used);	
	
	return 0;
}


#undef NRANSI