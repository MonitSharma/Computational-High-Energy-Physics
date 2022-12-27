/* Copyright (c) 2018-2020, Paolo Parotto and Debora Mroczek, 
 * Department of Physics, University of Houston, Houston, TX 77204, US. */

/* This main produces an EoS matching Lattice QCD at muB=0, and containing 
 * a critical point in the 3D Ising universality class, in a parametrized form.
 * It allows for different choices of constraints on the shape of the critical 
 * line, which reduce the number of parameters. */ 

#define NRANSI

// --------- INCLUDE -------- // 
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <err.h>

#include "../include/nrD.h"
#include "../include/nrutilD.h"
#include "../include/Variables.h"
#include "../include/Functions.h"

// --------- DEFINE --------- //
#ifndef PI
#define PI 3.141592653589
#endif

#ifndef N
#define N 2
#endif 


/* default settings */
int verbose=0;
int print=0;
int strange=0;
int lattice=0;
char lat_data[FILENAME_MAX];
char lat_data_SN[FILENAME_MAX];
char HRG_press[FILENAME_MAX];
char HRG_press_SN[FILENAME_MAX];
char MODE[FILENAME_MAX];
int lowT=5, highT=821;
int lowMU=0, highMU=601;
int	lowT_out=30, highT_out=800;
int	lowMU_out=0, highMU_out=450;

void help() /*{{{*/
{
  puts("");
  puts(" Usage:");
  puts("        ./EoS [options] <param_file>");
  puts("");
  puts(" <param_file>");
  puts("        Lattice_chi:    path to file with lattice coefficients");
	puts("                        for muQ=muS=0");
  puts("        Lattice_chi_SN: path to file with lattice coefficients");
	puts("                        for strangeness neutrality");
  puts("        HRG_press:      path to fie with HRG pressure for muQ=muS=0");
  puts("        HRG_press_SN:   path to file with HRG pressure for");
	puts("                        strangeness neutrality");
  puts("        Mode:           whether CP is to be placed onto parabola,");
	puts("                        straight line, or freely.");
  puts("        T0:             chiral transition temperature at muB=0");
	puts("                        (only for mode PAR,STR");
  puts("        kappa:          curvature of chiral transition at muB=0");
	puts("                        (only for mode PAR)");
  puts("        TC:             critical temperature (only for mode FREE)");
  puts("        angle1:         angle alpha1 (only for modes FREE,STR)");
  puts("        muBC:           critical chemical potential");
  puts("        anglediff:      alpha2-alpha1 (in degrees)");
  puts("        ww:             parameter w");
  puts("        rho:            parameter rho");
  puts("");
  puts(" options:");
  puts("        -r<lowT>:<highT>    Change T-range (still within limits)");
  puts("        -R<lowMU>:<highMU>  Change mu-range (still within limits)");
  puts("        -L                  Lattice mode: produce EoS without CP");
  puts("        -S                  Consider the strangeness neutral case");
  puts("        -p                  Output data for intermediate steps too");
  puts("        -v                  Verbose, print more stuff");
  puts("        -h                  Print help");
  puts("");
} /*}}}*/

void allocate () /*{{{*/ 
{
	/* Vectors for the root finding */ 
	x=vector(1,N);	f=vector(1,N);
	// Matrix for Jacobian. 
	JJ=matrix(1,N,1,N);
	
	// Rank 3 tensor for coordinates.
	Coords=f3tensor(lowT,highT,0,2*highMU,1,N);
	// Vectors for Chi's at muB=0. 
	Chi0LatVec=vector(lowT,highT);			
	Chi2LatVec=vector(lowT,highT);			
	Chi4LatVec=vector(lowT,highT);			
	Chi0IsingVec=vector(lowT,highT);			
	Chi2IsingVec=vector(lowT,highT);			
	Chi4IsingVec=vector(lowT,highT);			
	Chi0NoIsingVec=vector(lowT,highT);		
	Chi2NoIsingVec=vector(lowT,highT);		
	Chi4NoIsingVec=vector(lowT,highT);		
	dChi0NoIsingdTVec=vector(lowT,highT);	
	dChi2NoIsingdTVec=vector(lowT,highT);	
	dChi4NoIsingdTVec=vector(lowT,highT);	
	d2Chi0NoIsingdT2Vec=vector(lowT,highT);	
	d2Chi2NoIsingdT2Vec=vector(lowT,highT);	
	d2Chi4NoIsingdT2Vec=vector(lowT,highT);	
	// Matrices for thermodynamic functions over the phase diagram. 
	// NonIsing
	PressNoIsingMat=matrix(lowT,highT,0,highMU);			
	dPressNoIsingdTMat=matrix(lowT,highT,0,highMU);				
	dPressNoIsingdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressNoIsingdT2Mat=matrix(lowT,highT,0,highMU);		
	d2PressNoIsingdmuB2Mat=matrix(lowT,highT,0,highMU);			
	d2PressNoIsingdTdmuBMat=matrix(lowT,highT,0,highMU);
	// NonIsing (after filtering)
	PressNoIsingFilterMat=matrix(lowT,highT,0,highMU);		
	dPressNoIsingFilterdTMat=matrix(lowT,highT,0,highMU);		
	dPressNoIsingFilterdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressNoIsingFilterdT2Mat=matrix(lowT,highT,0,highMU);	
	d2PressNoIsingFilterdmuB2Mat=matrix(lowT,highT,0,highMU);	
	d2PressNoIsingFilterdTdmuBMat=matrix(lowT,highT,0,highMU);
	// Ising
	PressIsingMat=matrix(lowT,highT,0,highMU);				
	dPressIsingdTMat=matrix(lowT,highT,0,highMU);				
	dPressIsingdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressIsingdT2Mat=matrix(lowT,highT,0,highMU);			
	d2PressIsingdmuB2Mat=matrix(lowT,highT,0,highMU);			
	d2PressIsingdTdmuBMat=matrix(lowT,highT,0,highMU);
	// Ising + NonIsing (after filtering)
	PressTotMat=matrix(lowT,highT,0,highMU);				
	dPressTotdTMat=matrix(lowT,highT,0,highMU);					
	dPressTotdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressTotdT2Mat=matrix(lowT,highT,0,highMU);			
	d2PressTotdmuB2Mat=matrix(lowT,highT,0,highMU);				
	d2PressTotdTdmuBMat=matrix(lowT,highT,0,highMU);
	// HRG
	PressHRGMat=matrix(lowT,highT,0,highMU);				
	dPressHRGdTMat=matrix(lowT,highT,0,highMU);					
	dPressHRGdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressHRGdT2Mat=matrix(lowT,highT,0,highMU);			
	d2PressHRGdmuB2Mat=matrix(lowT,highT,0,highMU);				
	d2PressHRGdTdmuBMat=matrix(lowT,highT,0,highMU);
	// Total + HRG
	PressTotHRGMat=matrix(lowT,highT,0,highMU);				
	dPressTotHRGdTMat=matrix(lowT,highT,0,highMU);				
	dPressTotHRGdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressTotHRGdT2Mat=matrix(lowT,highT,0,highMU);		
	d2PressTotHRGdmuB2Mat=matrix(lowT,highT,0,highMU);			
	d2PressTotHRGdTdmuBMat=matrix(lowT,highT,0,highMU);
	// Final (normalized)
	PressFinalMat=matrix(lowT,highT,0,highMU);				
	EntropyFinalMat=matrix(lowT,highT,0,highMU);				
	BarDensityFinalMat=matrix(lowT,highT,0,highMU);
	EnerDensityFinalMat=matrix(lowT,highT,0,highMU);		
	SpSoundFinalMat=matrix(lowT,highT,0,highMU);				
	Chi2FinalMat=matrix(lowT,highT,0,highMU);			
	/* For correlation length */
	CorrLengthEpsExpMat=matrix(lowT,highT,0,highMU);				
	CorrLengthAsymptMat=matrix(lowT,highT,0,highMU);				
	CorrLengthMergedMat=matrix(lowT,highT,0,highMU);				
	// If LAT only is chosen (not normalized)
	PressLATonlyMat=matrix(lowT,highT,0,highMU);			
	dPressLATonlydTMat=matrix(lowT,highT,0,highMU);				
	dPressLATonlydmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressLATonlydT2Mat=matrix(lowT,highT,0,highMU);		
	d2PressLATonlydmuB2Mat=matrix(lowT,highT,0,highMU);			
	d2PressLATonlydTdmuBMat=matrix(lowT,highT,0,highMU);
	PressLATonlyFilterMat=matrix(lowT,highT,0,highMU);
	// If LAT only is chosen (normalized)
	PressLATonlyNormMat=matrix(lowT,highT,0,highMU);		
	EntropyLATonlyNormMat=matrix(lowT,highT,0,highMU);			
	BarDensityLATonlyNormMat=matrix(lowT,highT,0,highMU);	
	EnerDensityLATonlyNormMat=matrix(lowT,highT,0,highMU);	
	SpSoundLATonlyNormMat=matrix(lowT,highT,0,highMU);			
	Chi2LATonlyNormMat=matrix(lowT,highT,0,highMU);
} /*}}}*/

void read_param(char * param_file) /*{{{*/
{
  /* Read from the parameter file the following:
   * Lattice_chi: lattice coefficients for muS=muQ=0
   * Lattice_chi_SN: lattice coefficients for strangeness neutrality
   * HRG_press: HRG pressure
	 * Mode: whether CP is to be placed onto parabola, straight line ...
	 * T0: chiral transition temp at muB=0 (only for mode PAR,STR)
	 * kappa: curvature of chiral transition at muB=0 (only for modes PAR)
	 * TC: critical temperature (only for mode FREE)
	 * angle1: angle alpha1 (only for modes FREE,STR)
	 * muBC: critical chemical potential
	 * anglediff: alpha2-alpha1
	 * ww: parameter w
	 * rho: parameter rho */
  fprintf(stderr,"Reading parameter file.\n");

  FILE * file=fopen(param_file,"r");
  char line[FILENAME_MAX];
  if(file==NULL) err(1,"Could not open %s\n", param_file);
  while(fgets(line, sizeof(line),file)){
    char *c=strchr(line,'#'); 
    if (c) *c=0; 

    if(sscanf(line,"Lattice_chi=%1024s ", lat_data)==1) continue;
    if(sscanf(line,"Lattice_chi_SN=%1024s ", lat_data_SN)==1) continue;
    if(sscanf(line,"HRG_press=%1024s ", HRG_press)==1) continue;
    if(sscanf(line,"HRG_press_SN=%1024s ", HRG_press_SN)==1) continue;
    if(sscanf(line,"Mode=%16s ", MODE)==1) continue;
    if(sscanf(line,"T0=%lf ", &T0)==1) continue;
    if(sscanf(line,"kappa=%lf ", &kappa)==1) continue;
    if(sscanf(line,"TC=%lf ", &TC_in)==1) continue;
    if(sscanf(line,"angle1=%lf ", &angle1_in)==1) continue;
    if(sscanf(line,"muBC=%lf ", &muBC)==1) continue;
    if(sscanf(line,"anglediff=%lf ", &anglediff)==1) continue;
    if(sscanf(line,"ww=%lf ", &ww)==1) continue;
    if(sscanf(line,"rho=%lf ", &rho)==1) continue;

    err(1,"Line: <%s> was not understood", line);
  }
  fclose(file);

	if(lattice) goto end_param;

	/* Import and store parameters */
 	fprintf(stderr,"MODE: %s\n", MODE);

	if(strcmp(MODE,"FREE") == 0){
		strcpy(MODESTR, "FREE");
		TC = TC_in;    
		angle1=angle1_in;
		angle2 = angle1+anglediff;    
		fprintf(stderr,"The parameters you entered are:\n"
									 "TC = %f\nmuBC = %f\nangle1 = %f\n"
									 "angle2 = %f\nw = %f\nrho = %f\n",
										TC,muBC,angle1,angle2,ww,rho);
	} else if(strcmp(MODE,"PAR") == 0){
		strcpy(MODESTR, "PAR");
		TC=T0+kappa/T0*muBC*muBC; 
		angle1 = 180/PI*fabs(atan(-2.0*kappa/T0*muBC)); 
		angle2 = angle1 + anglediff;
		fprintf(stderr,"PARABOLA\n");
		fprintf(stderr,"The parameters you entered are:\n"
									 "TC = %f\nmuBC = %f\nangle1 = %f\n"
									 "angle2 = %f\nw = %f\nrho = %f\n",
										TC,muBC,angle1,angle2,ww,rho);
	} else if(strcmp(MODE,"STR") == 0){
		strcpy(MODESTR, "STR");
		angle1 = angle1_in;
		TC=T0-atan(angle1*PI/180)*muBC;
		fprintf(stderr,"The parameters you entered are:\n"
									 "TC = %f\nmuBC = %f\nangle1 = %f\n"
									 "angle2 = %f\nw = %f\nrho = %f\n",
										TC,muBC,angle1,angle2,ww,rho);
	} else err(1,"Mode could not be set properly.\n");
	
	/* we determine dTC and dmuBC from TC, muBC, ww and rho. */
	dTC = ww*TC; 
	dmuBC = dTC*rho;
	
	end_param:
	fprintf(stderr,"Finished importing and setting parameters.\n");
} /*}}}*/

void read_lattice_data (char * lattice_data) /*{{{*/
{
	/* Chi's from Lattice are imported and stored in a vector */
	fprintf(stderr,"Importing lattice data.\n");
	FILE *LatIn = fopen(lattice_data,"r");
	if (LatIn == 0)	err(1, "Failed to open Lattice Data\n");
	fprintf(stderr,"File imported succesfully!\n");

 	/* Save lattice chis into vectors (make them dimension of 
	 * energy^4). */	
	for(i=lowT;fscanf(LatIn,"%lf %lf %lf %lf",&xIn1,&xIn2,&xIn3,&xIn4) !=EOF;i++){
  	Chi0LatVec[i] = xIn2*pow(i,4);
   	Chi2LatVec[i] = xIn3*pow(i,4);
   	Chi4LatVec[i] = xIn4*pow(i,4);
 	}
 	fclose(LatIn);
 	fprintf(stderr,"Successfully stored lattice data.\n");
} /*}}}*/

void do_lattice() /*{{{*/
{
	fprintf(stderr,"You chose to consider the lattice-only case," 
								 "without a Critical Point.\n");
  /* Creating directory and moving into it. */
	char LATdir[FILENAME_MAX];
	
	if (!strange) strcpy(LATdir,"LATonly"); 
	else strcpy(LATdir,"LATonly_SN");
	
	mkdir(LATdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	chdir(LATdir);

	/* Calculate Taylor expanded pressure. */
	for (i=lowT;i<=highT; i++) for(j=0; j<=highMU; j++){
		muBval = (double) j; 	
		Tval = (double) i;
  	PressLATonlyMat[i][j] = Chi0LatVec[i] 
														+ 1.0/2.0*Chi2LatVec[i]*pow(muBval/Tval,2) 
														+ 1.0/24.0*Chi4LatVec[i]*pow(muBval/Tval,4);
	}

	// Take all the needed derivatives of the pressure.
	// Wrt T
	deriv_matrix(PressLATonlyMat,dPressLATonlydTMat,
								1,lowT,highT,0,highMU,0);	
	deriv_matrix(dPressLATonlydTMat,d2PressLATonlydT2Mat,
								1,lowT,highT,0,highMU,0);
	// Wrt muB
	deriv_matrix(PressLATonlyMat,dPressLATonlydmuBMat,
								2,lowT,highT,0,highMU,0);
	deriv_matrix(dPressLATonlydmuBMat,d2PressLATonlydmuB2Mat,
								2,lowT,highT,0,highMU,1);
	// Wrt T and muB
	deriv_matrix(dPressLATonlydTMat,d2PressLATonlydTdmuBMat,
								2,lowT,highT,0,highMU,0);
	
	// That is it. Now we just need to combine these into the final observables (normalized).
	fprintf(stderr,"Calculating thermodynamics quantities. \n");

	// Create the files for export.
	FILE *FilePressLATonlyNorm = fopen("PressLATonlyNorm3D.dat", "w");		
	FILE *FileEntrLATonlyNorm = fopen("EntrLATonlyNorm3D.dat", "w");
	FILE *FileBarDensLATonlyNorm = fopen("BarDensLATonlyNorm3D.dat", "w");	
	FILE *FileEnerDensLATonlyNorm = fopen("EnerDensLATonlyNorm3D.dat", "w");
	FILE *FileSpSoundLATonlyNorm = fopen("SpsoundLATonlyNorm3D.dat", "w");	
	FILE *FileChi2LATonlyNorm = fopen("Chi2LATonlyNorm3D.dat", "w");
		
	// Export.
	for (i=lowT_out; i<=highT_out; i++) for(j=0;j<=highMU_out; j++){
		Tval = (double) i; muBval = (double) j;
		
		PressLATonlyNormMat[i][j] = PressLATonlyMat[i][j]/pow(Tval,4);
 		EntropyLATonlyNormMat[i][j] = dPressLATonlydTMat[i][j]/pow(Tval,3);
  	BarDensityLATonlyNormMat[i][j] = dPressLATonlydmuBMat[i][j]/pow(Tval,3);
  	EnerDensityLATonlyNormMat[i][j] = 
			(dPressLATonlydTMat[i][j]*Tval-PressLATonlyMat[i][j]
			 +muBval*dPressLATonlydmuBMat[i][j])/pow(Tval,4);
  	SpSoundLATonlyNormMat[i][j] = 
			(dPressLATonlydmuBMat[i][j]*dPressLATonlydmuBMat[i][j]
					*d2PressLATonlydT2Mat[i][j] 
			 - 2.0*dPressLATonlydTMat[i][j]*dPressLATonlydmuBMat[i][j]
					*d2PressLATonlydTdmuBMat[i][j] 
			 + dPressLATonlydTMat[i][j]*dPressLATonlydTMat[i][j]
					*d2PressLATonlydmuB2Mat[i][j])
			/(dPressLATonlydTMat[i][j]*Tval
			 +muBval*dPressLATonlydmuBMat[i][j])
			/(d2PressLATonlydT2Mat[i][j]*d2PressLATonlydmuB2Mat[i][j]
			 -d2PressLATonlydTdmuBMat[i][j]*d2PressLATonlydTdmuBMat[i][j]); 						 	Chi2LATonlyNormMat[i][j] = d2PressLATonlydmuB2Mat[i][j]/pow(Tval,2);

  	fprintf(FilePressLATonlyNorm,"%3.16f	%3.1f	%12.16f \n", 
			muBval, Tval, PressLATonlyNormMat[i][j]);
  	fprintf(FileEntrLATonlyNorm,"%3.16f	%3.1f	%12.16f \n", 
			muBval, Tval, EntropyLATonlyNormMat[i][j]);
  	fprintf(FileBarDensLATonlyNorm,"%3.16f	%3.1f	%12.16f \n", 
			muBval, Tval, BarDensityLATonlyNormMat[i][j]);
  	fprintf(FileEnerDensLATonlyNorm,"%3.16f	%3.1f	%12.16f \n", 
			muBval, Tval, EnerDensityLATonlyNormMat[i][j]);
 	  fprintf(FileSpSoundLATonlyNorm,"%3.16f	%3.1f	%12.16f \n", 
			muBval, Tval, SpSoundLATonlyNormMat[i][j]);
		fprintf(FileChi2LATonlyNorm,"%3.16f	%3.1f	%12.16f \n", 
			muBval, Tval, Chi2LATonlyNormMat[i][j]);
  	}
	// Close the files.
	fclose(FilePressLATonlyNorm);		
	fclose(FileEntrLATonlyNorm);		
	fclose(FileBarDensLATonlyNorm);
	fclose(FileEnerDensLATonlyNorm);		
	fclose(FileSpSoundLATonlyNorm);		
	fclose(FileChi2LATonlyNorm);	

	fprintf(stderr,"The procedure is completed. Exiting.\n");
} /*}}}*/

void write_filenames() /*{{{*/
{
	if (!strange)
		sprintf(nameFolder, "Files_%s_%d_%d_%d_%d_%d_%d",MODESTR,(int)TC,(int)muBC,
												(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	else 
		sprintf(nameFolder, "Files_%s_%d_%d_%d_%d_%d_%d_SN",MODESTR,(int)TC,(int)muBC,
												(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	
	sprintf(nameCoords, "Coords_%s_%d_%d_%d_%d_%d_%d.dat",MODESTR,(int)TC,(int)muBC,
												(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameChisIsing, "Chis_Ising_%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,(int)muBC,
												(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameChisNoIsing, "Chis_No_Ising_%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,
												(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedChisNoIsingdT, "dChis_No_Ising_dT_%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,
												(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2ChisNoIsingdT2, "d2Chis_No_Ising_dT2%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,
												(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressIsing3D, "Press_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,
												(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTIsing3D, "dPdT_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,
												(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBIsing3D, "dPdmuB_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,
												(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2Ising3D, "d2PdT2_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,
												(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2Ising3D, "d2PdmuB2_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,
												(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBIsing3D, "d2PdTdmuB_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressNoIsing3D, "Press_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTNoIsing3D, "dPdT_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBNoIsing3D, "dPdmuB_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2NoIsing3D, "d2PdT2_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2NoIsing3D, "d2PdmuB2_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBNoIsing3D, "d2PdTdmuB_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressIsingPlusNoIsing3D, "Press_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTIsingPlusNoIsing3D, "dPdT_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBIsingPlusNoIsing3D, "dPdmuB_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2IsingPlusNoIsing3D, "d2PdT2_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2IsingPlusNoIsing3D, "d2PdmuB2_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBIsingPlusNoIsing3D, "d2PdTdmuB_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressTotHRG3D, "Press_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTTotHRG3D, "dPdT_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBTotHRG3D, "dPdmuB_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2TotHRG3D, "d2PdT2_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2TotHRG3D, "d2PdmuB2_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBTotHRG3D, "d2PdTdmuB_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressFinal3D, "Press_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameEntrFinal3D, "Entr_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameBarDensFinal3D, "BarDens_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameEnerDensFinal3D, "EnerDens_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameSpsoundFinal3D, "SpSound_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameChi2Final3D, "Chi2_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameCorrLengthFinal3D, "CorrLength_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
} /*}}}*/

void calc_jacobian() /*{{{*/
{
	/* For calculations hereafter, angles are expressed in radians. */ 
	angle1 = angle1*PI/180;
	angle2 = angle2*PI/180;

	/* Jacobian matrix of the (r,h) -> (T,muB) map is generated. */  
	Jacobian(JJ,dTC,dmuBC,angle1,angle2);
	drdmuB = JJ[1][2];
	dhdmuB = JJ[2][2];
	drdT = JJ[1][1];
	dhdT = JJ[2][1];

	fprintf(stderr,"\nThe Jacobian calculation gave:\n"
								 "drdmuB = %12.16f\ndhdT = %12.16f\n"
								 "drdT = %12.16f\ndhmuB = %12.16f\n",
									drdmuB,dhdmuB,drdT,dhdT);
} /*}}}*/

void make_coords() /*{{{*/
{
	FILE *fp=fopen(nameCoords,"w");
	fprintf(stderr,"\nCreating the coordinates file...\n");
	// Time is taken at start and end to measure how long it took to generate the coordinates file.
	time(&start);
	for(Ti=lowT;Ti<=highT;Ti+=1.0) for(muBi=-highMU;muBi<=highMU;muBi+=1.0){
		h_i = hh(Ti,muBi);
		r_i  = rr(Ti, muBi);

    /* Condition for choosing init condition in root finding. */
   	x[1] = 1.0;
		if(h_i > 0) x[2] = 1.0;
		if(h_i <= 0) x[2] = -1.0;
						
		/* The actual root finding routine. */
		newt(x,N,&check,funcv);
		funcv(N,x,f);
		if (check) printf("Convergence problems.\n");
		/* Print to file the coordinates {T,muB,R,Theta} (in this order). */
		fprintf(fp,"%3.1f %3.1f %2.16f %2.16f\n",Ti,muBi,x[1],x[2]);
	} 
	time(&end);
	/* A control on whether the loops reached the end, and did not get stuck, otherwise file is eliminated. */	
	if(Ti == (highT+1) && muBi == (highMU+1))	fprintf(stderr,"\nThe file %s was created successfully in %d seconds\n", nameCoords, (int) difftime(end,start));
	else{
		remove(nameCoords);
		err(1,"The file %s was NOT created successfully,"
					" removing it.\n",nameCoords);
	} 								
	fclose(fp);

	// Coordinates are imported and stored in a Rank 3 tensor Coords. 
	fprintf(stderr,"Importing coordinates from file %s \n",nameCoords);
	fp=fopen(nameCoords,"r");
	for(j=0, x1int=lowT-1;fscanf(fp,"%lf %lf %lf %lf\n", &xIn1, &xIn2, &xIn3, &xIn4) != EOF; j++ ){
		x2int = j % (2*highMU+1);
		if (x2int == 0) x1int++;
 		Coords[x1int][x2int][1] = xIn3;
		Coords[x1int][x2int][2] = xIn4; 
	}
	fclose(fp);
  fprintf(stderr,"Imported coordinates successfully.\n");
} /*}}}*/

void make_Ising() /*{{{*/
{
	/* The coordinates corresponding to muB=0 are used 
	 * to generate the critical Chi's at different temperatures
	 * and exported. */
	fprintf(stderr,"Generating Ising at muB=0\n");

	FILE * FileChiIsing=NULL;
	FILE * FilePressIsing=NULL;				
	FILE * FiledPdTIsing=NULL;
	FILE * FiledPdmuBIsing=NULL;			
	FILE * Filed2PdT2Ising=NULL;
	FILE * Filed2PdmuB2Ising=NULL;		
	FILE * Filed2PdTdmuBIsing=NULL;
	if (print){
		FileChiIsing=fopen(nameChisIsing,"w");
		FilePressIsing=fopen(namePressIsing3D,"w");				
		FiledPdTIsing=fopen(namedPdTIsing3D,"w");
		FiledPdmuBIsing=fopen(namedPdmuBIsing3D,"w");			
		Filed2PdT2Ising=fopen(named2PdT2Ising3D,"w");
		Filed2PdmuB2Ising=fopen(named2PdmuB2Ising3D,"w");		
		Filed2PdTdmuBIsing=fopen(named2PdTdmuBIsing3D,"w");
	}

	for(i=lowT;i<=highT;i++){
		Tval = (double) i;
		Chi0IsingVec[i] = -G(Coords[i][highMU][1],Coords[i][highMU][2])*pow(TC,4);
		Chi2IsingVec[i] = -Tval*Tval*d2GdmuB2ConT(Coords[i][highMU][1],Coords[i][highMU][2])*pow(TC,4);
		Chi4IsingVec[i] = -Tval*Tval*Tval*Tval*d4GdmuB4ConT(Coords[i][highMU][1],Coords[i][highMU][2])*pow(TC,4);
	
		if (print){
			fprintf(FileChiIsing,"%3.1f %12.16f %12.16f %12.16f\n",(double) i,Chi0IsingVec[i],Chi2IsingVec[i],Chi4IsingVec[i]);
		}
	}
	fprintf(stderr,"Generation successful\n");

	/* The pressure is generated over the whole range of coordinates.
	 * NOTE: the pressure is symmetrized around muB=0 to ensure all odd order 
	 * derivatives vanish at muB=0. */
	fprintf(stderr,"Generating Ising in 3D\n");
	
	for(i=lowT;i<highT;i++) for(j=0;j<=highMU;j++){
		k = j + highMU;
		if(j==0){
			PressIsingMat[i][j] = - G(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
			dPressIsingdTMat[i][j] = - dGdTConmuB(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
			dPressIsingdmuBMat[i][j] = 0.0;
			d2PressIsingdT2Mat[i][j] = - d2GdT2ConmuB(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
			d2PressIsingdmuB2Mat[i][j] = - d2GdmuB2ConT(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
			d2PressIsingdTdmuBMat[i][j] = 0.0;
			CorrLengthEpsExpMat[i][j] = CalcCorrLengthEpsExp(Coords[i][k][1],Coords[i][k][2]);
			CorrLengthAsymptMat[i][j] = CalcCorrLengthAsympt(Coords[i][k][1],Coords[i][k][2]);
		} else{
	 		PressIsingMat[i][j] = - (G(Coords[i][k][1],Coords[i][k][2]) + G(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		dPressIsingdTMat[i][j] = - (dGdTConmuB(Coords[i][k][1],Coords[i][k][2]) + dGdTConmuB(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		dPressIsingdmuBMat[i][j] = - (dGdmuBConT(Coords[i][k][1],Coords[i][k][2]) - dGdmuBConT(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		d2PressIsingdT2Mat[i][j] = - (d2GdT2ConmuB(Coords[i][k][1],Coords[i][k][2]) + d2GdT2ConmuB(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		d2PressIsingdmuB2Mat[i][j] = - (d2GdmuB2ConT(Coords[i][k][1],Coords[i][k][2]) + d2GdmuB2ConT(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		d2PressIsingdTdmuBMat[i][j] = - (d2GdTdmuB(Coords[i][k][1],Coords[i][k][2]) - d2GdTdmuB(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
			CorrLengthEpsExpMat[i][j] = CalcCorrLengthEpsExp(Coords[i][k][1],Coords[i][k][2]);
			CorrLengthAsymptMat[i][j] = CalcCorrLengthAsympt(Coords[i][k][1],Coords[i][k][2]);
 		}
		if (print){	
 			fprintf(FilePressIsing,"%3.1f	%3.1f	%12.16f\n",
								(double) j,(double) i, PressIsingMat[i][j]);
 			fprintf(FiledPdTIsing,"%3.1f	%3.1f	%12.16f\n",
								(double) j,(double) i, dPressIsingdTMat[i][j]);
 			fprintf(FiledPdmuBIsing,"%3.1f	%3.1f	%12.16f\n",
								(double) j,(double) i, dPressIsingdmuBMat[i][j]);
 			fprintf(Filed2PdT2Ising,"%3.1f	%3.1f	%12.16f\n",
								(double) j,(double) i, d2PressIsingdT2Mat[i][j]);
			fprintf(Filed2PdmuB2Ising,"%3.1f	%3.1f	%12.16f\n",
								(double) j,(double) i, d2PressIsingdmuB2Mat[i][j]);
 			fprintf(Filed2PdTdmuBIsing,"%3.1f	%3.1f	%12.16f\n",
								(double) j,(double) i, d2PressIsingdTdmuBMat[i][j]);
		}
	}
	
	if (print){
		fclose(FileChiIsing);
		fclose(FilePressIsing); 
		fclose(FiledPdTIsing); 
		fclose(FiledPdmuBIsing); 
		fclose(Filed2PdT2Ising); 
		fclose(Filed2PdmuB2Ising);	
		fclose(Filed2PdTdmuBIsing);		
	}
} /*}}}*/

void merge_CorrLength() /*{{{*/
{
	Ttip = 145;				muBtip = 320;
	TMlow = 195.0; 		T0low = 58.0; 			klow = 0.00318462;
	T0high = 282.594; khigh = 0.00207197;
	TMmerg = 17.0; 		T0merg = 0.1; 			kmerg = 0.0075;

	for (i=lowT_out; i<=highT_out; i++) for(j=0;j<=highMU_out; j++){
		Tval = (double) i; muBval = (double) j;

		Tlow = TMlow - (TMlow - T0low)*exp(-klow*muBval);
		Thigh = T0high*exp(-khigh*muBval);
		Tmerg = T0merg + (TMmerg - T0merg)*exp(-kmerg*muBval);

		if(i<Ttip && j<=muBtip){   
			CorrLengthMergedMat[i][j] = CorrLengthAsymptMat[i][j]*0.5*(1.0 + tanh((Tval - Tlow)/Tmerg)) 
																+ CorrLengthEpsExpMat[i][j]*0.5*(1.0 - tanh((Tval - Tlow)/Tmerg));
		}
		else if(i>=Ttip && j<=muBtip){
			CorrLengthMergedMat[i][j] = CorrLengthAsymptMat[i][j]*0.5*(1.0 - tanh((Tval - Thigh)/Tmerg)) 
																+ CorrLengthEpsExpMat[i][j]*0.5*(1.0 + tanh((Tval - Thigh)/Tmerg));
		}
		else{
			 CorrLengthMergedMat[i][j] = CorrLengthEpsExpMat[i][j];
		}
	}
} /*}}}*/

void make_nonIsing() /*{{{*/
{
	/* Non-Ising Chi's are calculated, exported and stored. */ 
	fprintf(stderr,"Generating non-Ising at muB = 0\n");

	FILE * FileChiNoIsing;
	FILE * FilePressNoIsing;
	FILE * FilePressNoIsingFilter; 
	FILE * FiledPdTNoIsingFilter;
	FILE * FiledPdmuBNoIsingFilter; 
	FILE * Filed2PdT2NoIsingFilter; 
	FILE * Filed2PdmuB2NoIsingFilter; 
	FILE * Filed2PdTdmuBNoIsingFilter;
	if (print){ /* if printing intermediate steps is wanted */
		FileChiNoIsing=fopen(nameChisNoIsing,"w");
	 	FilePressNoIsing=fopen(namePressNoIsing3D,"w");			
		FilePressNoIsingFilter = fopen("NoIsingPressFilter.dat","w");
		FiledPdTNoIsingFilter = fopen("NoIsingdPdTFilter.dat","w");						
		FiledPdmuBNoIsingFilter = fopen("NoIsingdPdmuBFilter.dat","w");
		Filed2PdT2NoIsingFilter = fopen("NoIsingd2PdT2Filter.dat","w");				
		Filed2PdmuB2NoIsingFilter = fopen("NoIsingd2PdmuB2Filter.dat","w");
		Filed2PdTdmuBNoIsingFilter = fopen("NoIsingd2PdTdmuBFilter.dat","w");	
	}

	for (i=lowT; i<=highT; i++){
		Tval = (double) i;
		Chi0NoIsingVec[i] = Chi0LatVec[i] - Chi0IsingVec[i];
		Chi2NoIsingVec[i] = Chi2LatVec[i] - Chi2IsingVec[i];
		Chi4NoIsingVec[i] = Chi4LatVec[i] - Chi4IsingVec[i];
        
	  if(print)	fprintf(FileChiNoIsing,"%3.1f %12.16f %12.16f %12.16f\n", Tval, Chi0NoIsingVec[i],Chi2NoIsingVec[i],Chi4NoIsingVec[i]);
	}

	/* The non-Ising pressure in 3D is calculated from Taylor expansion, 
	 * stored and exported. */
 	fprintf(stderr,"Calculating non-Ising in 3D.\n");
	for (i=lowT; i<=highT; i++) for(j=0; j<=highMU; j++){
		muBval = (double) j; 
		Tval = (double) i;
    PressNoIsingMat[i][j] = Chi0NoIsingVec[i]
														+(1.0/2.0)*Chi2NoIsingVec[i]*pow(muBval/Tval,2)
														+(1.0/24.0)*Chi4NoIsingVec[i]*pow(muBval/Tval,4);

		if (print) fprintf(FilePressNoIsing,"%3.1f	%3.1f	%12.16f\n",
										muBval,Tval,PressNoIsingMat[i][j]);
	}

	/* Filtering of regular pressure here. Widths of the filter are also 
	 * set here. Result is stored and exported. */
	fprintf(stderr,"Filtering non-Ising Pressure in 3D\n");

	sigmax = 18.0;	sigmay = 1.0;

	for (i=lowT; i<=highT; i++)	for(j=0; j<=highMU; j++){
		muBval = (double) j; Tval = (double) i;

		kmin = i-sigmax; kmax = i+sigmax;
		lmin = j-sigmay; lmax = j+sigmay;	
 
		sum = 0.0; norm = 0.0;
		for(k=kmin;k<=kmax;k++) for(l=lmin;l<=lmax;l++){
			if(k<=lowT)		kint = lowT;
			else if(k>=highT)	kint = highT;	
			else		 	kint = k;
					
			if(l<=0)  		lint = -l;
			else if(l>=highMU)	lint = highMU;		
			else			lint = l;
									
			sum+=GausFunc(kint-i,lint-j,sigmax,sigmay)*PressNoIsingMat[kint][lint];
			norm+=GausFunc(kint-i,lint-j,sigmax,sigmay);
		}
		PressNoIsingFilterMat[i][j] = sum/norm;
		if (print) fprintf(FilePressNoIsingFilter, "%3.1f %3.1f %12.16f\n", 
													muBval,Tval,PressNoIsingFilterMat[i][j]);
	}

	/* Derivatives (numerical) of the filtered Non-Ising pressure 
	 * and all thermodynamic quantities of interest are calculated 
	 * and stored, over the phase diagram. */
	fprintf(stderr,"Calculating (numerical) derivatives"
								 " of the filtered Non-Ising Pressure\n");

	/* Take all the needed derivatives */
	/* Wrt T */
	deriv_matrix(PressNoIsingFilterMat,dPressNoIsingFilterdTMat,1,lowT,highT,0,highMU,0);	
	deriv_matrix(dPressNoIsingFilterdTMat,d2PressNoIsingFilterdT2Mat,1,lowT,highT,0,highMU,0);
	/* Wrt muB */
	deriv_matrix(PressNoIsingFilterMat,dPressNoIsingFilterdmuBMat,2,lowT,highT,0,highMU,0);
	deriv_matrix(dPressNoIsingFilterdmuBMat,d2PressNoIsingFilterdmuB2Mat,2,lowT,highT,0,highMU,1);
	/* Wrt T and muB */
	deriv_matrix(dPressNoIsingFilterdTMat,d2PressNoIsingFilterdTdmuBMat,2,lowT,highT,0,highMU,0);

	if(print){
		fclose(FileChiNoIsing);
		fclose(FilePressNoIsing);
		fclose(FilePressNoIsingFilter);
		fclose(FiledPdTNoIsingFilter);
		fclose(FiledPdmuBNoIsingFilter); 
		fclose(Filed2PdT2NoIsingFilter); 
		fclose(Filed2PdmuB2NoIsingFilter); 
		fclose(Filed2PdTdmuBNoIsingFilter);
	}
} /*}}}*/

void merge_Ising_nonIsing() /*{{{*/
{
	/* Calculate pressure and derivatives summing Ising and non-Ising 
	 * pressure, then normalize. Store and export. */
  fprintf(stderr,"Calculating other derivatives in 3D \n");
	FILE *FilePressTot;
	FILE *FiledPdTIsingPlusNoIsing;
	FILE *FiledPdmuBIsingPlusNoIsing;
	FILE *Filed2PdT2IsingPlusNoIsing;
	FILE *Filed2PdmuB2IsingPlusNoIsing;
	FILE *Filed2PdTdmuBIsingPlusNoIsing;
	if(print){
		FilePressTot = fopen(namePressIsingPlusNoIsing3D, "w");					
		FiledPdTIsingPlusNoIsing=fopen(namedPdTIsingPlusNoIsing3D,"w");
		FiledPdmuBIsingPlusNoIsing=fopen(namedPdmuBIsingPlusNoIsing3D,"w");		
		Filed2PdT2IsingPlusNoIsing=fopen(named2PdT2IsingPlusNoIsing3D,"w");
		Filed2PdmuB2IsingPlusNoIsing=fopen(named2PdmuB2IsingPlusNoIsing3D,"w");	
		Filed2PdTdmuBIsingPlusNoIsing=fopen(named2PdTdmuBIsingPlusNoIsing3D,"w");
  }

	for (i=lowT; i<=highT; i++) for(j=0;j<=highMU; j++){
		Tval = (double) i; muBval = (double) j;
		
		PressTotMat[i][j] = PressIsingMat[i][j] + PressNoIsingFilterMat[i][j];
   	dPressTotdTMat[i][j] = dPressIsingdTMat[i][j] + dPressNoIsingFilterdTMat[i][j];
   	dPressTotdmuBMat[i][j] = dPressIsingdmuBMat[i][j] + dPressNoIsingFilterdmuBMat[i][j];
   	d2PressTotdT2Mat[i][j] = d2PressIsingdT2Mat[i][j] + d2PressNoIsingFilterdT2Mat[i][j];
   	d2PressTotdmuB2Mat[i][j] = d2PressIsingdmuB2Mat[i][j] + d2PressNoIsingFilterdmuB2Mat[i][j];  		
   	d2PressTotdTdmuBMat[i][j] = d2PressIsingdTdmuBMat[i][j] + d2PressNoIsingFilterdTdmuBMat[i][j];   

		if(print){
   		fprintf(FilePressTot,"%3.16f	%3.1f	%12.16f \n", 
								muBval, Tval, PressTotMat[i][j]);
 			fprintf(FiledPdTIsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",
								muBval,Tval,dPressTotdTMat[i][j]);
 			fprintf(FiledPdmuBIsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",
								muBval,Tval,dPressTotdmuBMat[i][j]);
 			fprintf(Filed2PdT2IsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",
								muBval,Tval,d2PressTotdT2Mat[i][j]);
 			fprintf(Filed2PdmuB2IsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",		
								muBval,Tval,d2PressTotdmuB2Mat[i][j]);
 			fprintf(Filed2PdTdmuBIsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",
								muBval,Tval,d2PressTotdTdmuBMat[i][j]);
  	}
	}
	if(print){
		fclose(FilePressTot);					
		fclose(FiledPdTIsingPlusNoIsing); 		
		fclose(FiledPdmuBIsingPlusNoIsing); 
		fclose(Filed2PdT2IsingPlusNoIsing);	
		fclose(Filed2PdmuB2IsingPlusNoIsing);	
		fclose(Filed2PdTdmuBIsingPlusNoIsing);
	}	
} /*}}}*/

void merge_HRG() /*{{{*/
{
 	/* Now, the smooth joining with the HRG pressure starts. 
	 * Import HRG pressure and store. */
	fprintf(stderr,"Importing HRG Pressure \n");
	
	FILE * FilePressHRG;

	if (!strange) FilePressHRG=fopen(HRG_press,"r");
	else FilePressHRG=fopen(HRG_press_SN,"r"); 
	
	if (FilePressHRG==NULL) err(1, "Failed to open HRG Pressure.\n");

	for(i=0, x2int=0;fscanf(FilePressHRG,"%lf %lf %lf\n",&xIn1,&xIn2,&xIn3) !=EOF;i++){
		x1int = (i % 817) + 5;
		PressHRGMat[x1int][x2int] = xIn3*pow(x1int,4);
		if(x1int == 821) x2int++;
	}
	fclose(FilePressHRG); 

	/* Calculated derivatives of the HRG pressure */
	/* Wrt T */
	deriv_matrix(PressHRGMat,dPressHRGdTMat,1,lowT,highT,0,highMU,0);	
	deriv_matrix(dPressHRGdTMat,d2PressHRGdT2Mat,1,lowT,highT,0,highMU,0);
	/* Wrt muB */
	deriv_matrix(PressHRGMat,dPressHRGdmuBMat,2,lowT,highT,0,highMU,0);
	deriv_matrix(dPressHRGdmuBMat,d2PressHRGdmuB2Mat,2,lowT,highT,0,highMU,1);
	/* Wrt T and muB */
	deriv_matrix(dPressHRGdTMat,d2PressHRGdTdmuBMat,2,lowT,highT,0,highMU,0);
	
	/* Final pressure, joined with the HRG one is calculated, 
	 * stored and exported. The parameter deltaT for merging with 
	 * HRG is also set here. */
	fprintf(stderr,"Merging with HRG data.\n");

  /* Merging with HRG data. */
	FILE *FilePressTotHRG;
	FILE *FiledPdTTotHRG;
	FILE *FiledPdmuBTotHRG;
	FILE *Filed2PdT2TotHRG;
	FILE *Filed2PdmuB2TotHRG;
	FILE *Filed2PdTdmuBTotHRG;
	if (print){
		FilePressTotHRG = fopen(namePressTotHRG3D, "w");		
		FiledPdTTotHRG = fopen(namedPdTTotHRG3D,"w");
		FiledPdmuBTotHRG = fopen(namedPdmuBTotHRG3D,"w");		
		Filed2PdT2TotHRG = fopen(named2PdT2TotHRG3D,"w");
		Filed2PdmuB2TotHRG = fopen(named2PdmuB2TotHRG3D,"w");	
		Filed2PdTdmuBTotHRG = fopen(named2PdTdmuBTotHRG3D,"w");
	}


	// Set delta T
	deltaT = 17.0;

	/* The actual merging */
	for (i=lowT; i<=highT; i++) for(j=0;j<=highMU; j++){
		Tval = (double) i; muBval = (double) j;
  	Targum = (Tval - (T0 + kappa/T0*muBval*muBval - 23.0))/deltaT;

  	if (Targum>=4.5)				
			PressTotHRGMat[i][j] = PressTotMat[i][j];
		else if (Targum<=-4.5)	
			PressTotHRGMat[i][j] = PressHRGMat[i][j];
  	else 										
			PressTotHRGMat[i][j] = PressTotMat[i][j]*0.5*(1.0 + tanh(Targum)) 
                             	+ PressHRGMat[i][j]*0.5*(1.0 - tanh(Targum));
	
  	if(Targum >= 4.5)       
			dPressTotHRGdTMat[i][j] = dPressTotdTMat[i][j];
  	else if(Targum <= -4.5)      
			dPressTotHRGdTMat[i][j] = dPressHRGdTMat[i][j];
  	else                    
			dPressTotHRGdTMat[i][j] = dPressTotdTMat[i][j]*0.5*(1.0 + tanh(Targum))
                        + dPressHRGdTMat[i][j]*0.5*(1.0 - tanh(Targum))
                        + PressTotMat[i][j]*0.5/deltaT*pow(cosh(Targum),-2)
                        - PressHRGMat[i][j]*0.5/deltaT*pow(cosh(Targum),-2);       
  	if(Targum >= 4.5)       
			dPressTotHRGdmuBMat[i][j] = dPressTotdmuBMat[i][j];
  	else if(Targum <= -4.5)      
			dPressTotHRGdmuBMat[i][j] = dPressHRGdmuBMat[i][j];
  	else                    
			dPressTotHRGdmuBMat[i][j] = 
				dPressTotdmuBMat[i][j]*0.5*(1.0 + tanh(Targum))
        + dPressHRGdmuBMat[i][j]*0.5*(1.0 - tanh(Targum))
        - PressTotMat[i][j]*0.5*pow(cosh(Targum),-2)*(2*kappa/T0*muBval)/deltaT
        + PressHRGMat[i][j]*0.5*pow(cosh(Targum),-2)*(2*kappa/T0*muBval)/deltaT;

  	if(Targum >= 4.5)       
			d2PressTotHRGdT2Mat[i][j] = d2PressTotdT2Mat[i][j];
  	if(Targum <= -4.5)      
			d2PressTotHRGdT2Mat[i][j] = d2PressHRGdT2Mat[i][j];
  	else                    
			d2PressTotHRGdT2Mat[i][j] = 
				d2PressTotdT2Mat[i][j]*0.5*(1.0 + tanh(Targum)) 
        + d2PressHRGdT2Mat[i][j]*0.5*(1.0 - tanh(Targum))
        + dPressTotdTMat[i][j]/deltaT*pow(cosh(Targum),-2)
        - dPressHRGdTMat[i][j]/deltaT*pow(cosh(Targum),-2)
        - PressTotMat[i][j]*pow(cosh(Targum),-2)*tanh(Targum)/pow(deltaT,2) 
        + PressHRGMat[i][j]*pow(cosh(Targum),-2)*tanh(Targum)/pow(deltaT,2);   

  	if(Targum >= 4.5)       
			d2PressTotHRGdmuB2Mat[i][j] = d2PressTotdmuB2Mat[i][j];
  	if(Targum <= -4.5)      
			d2PressTotHRGdmuB2Mat[i][j] = d2PressHRGdmuB2Mat[i][j];
   	else                    
			d2PressTotHRGdmuB2Mat[i][j] = 
				d2PressTotdmuB2Mat[i][j]*0.5*(1.0 + tanh(Targum))
        + d2PressHRGdmuB2Mat[i][j]*0.5*(1.0 - tanh(Targum))
        - dPressTotdmuBMat[i][j]*pow(cosh(Targum),-2)*(2*kappa/T0*muBval)/deltaT
        + dPressHRGdmuBMat[i][j]*pow(cosh(Targum),-2)*(2*kappa/T0*muBval)/deltaT
        - PressTotMat[i][j]*0.5*pow(cosh(Targum),-2)
					*(2*tanh(Targum)*(2*kappa/T0*muBval)/deltaT*(2*kappa/T0*muBval)/deltaT
          	+ 2*kappa/T0/deltaT)
        + PressHRGMat[i][j]*0.5*pow(cosh(Targum),-2)
					*(2*tanh(Targum)*(2*kappa/T0*muBval)/deltaT*(2*kappa/T0*muBval)/deltaT
          	+ 2*kappa/T0/deltaT);

  	if(Targum >= 4.5)       
			d2PressTotHRGdTdmuBMat[i][j] = d2PressTotdTdmuBMat[i][j];
  	if(Targum <= -4.5)      
			d2PressTotHRGdTdmuBMat[i][j] = d2PressHRGdTdmuBMat[i][j];
  	else                    
			d2PressTotHRGdTdmuBMat[i][j] = 
				d2PressTotdTdmuBMat[i][j]*0.5*(1.0 + tanh(Targum))
        + d2PressHRGdTdmuBMat[i][j]*0.5*(1.0 - tanh(Targum))
        + dPressTotdmuBMat[i][j]*0.5/deltaT*pow(cosh(Targum),-2)
        - dPressHRGdmuBMat[i][j]*0.5/deltaT*pow(cosh(Targum),-2) 
        - dPressTotdTMat[i][j]*0.5/deltaT
					*(2*kappa/T0*muBval)*pow(cosh(Targum),-2)
        + dPressHRGdTMat[i][j]*0.5/deltaT
					*(2*kappa/T0*muBval)*pow(cosh(Targum),-2)
        + PressTotMat[i][j]*pow(cosh(Targum),-2)/pow(deltaT,2)
        	*tanh(Targum)*(2*kappa/T0*muBval)
        - PressHRGMat[i][j]*pow(cosh(Targum),-2)/pow(deltaT,2)
        	*tanh(Targum)*(2*kappa/T0*muBval);    

		/* Output is disabled now */
   	if (print){
			fprintf(FilePressTotHRG,"%3.16f	%3.1f	%12.16f \n", 
								muBval, Tval, PressTotHRGMat[i][j]);
  		fprintf(FiledPdTTotHRG,"%3.16f	%3.1f	%12.16f \n", 
								muBval, Tval, dPressTotHRGdTMat[i][j]);
   		fprintf(FiledPdmuBTotHRG,"%3.16f	%3.1f	%12.16f \n", 
								muBval, Tval, dPressTotHRGdmuBMat[i][j]);
   		fprintf(Filed2PdT2TotHRG,"%3.16f	%3.1f	%12.16f \n", 
								muBval, Tval, d2PressTotHRGdT2Mat[i][j]);
 			fprintf(Filed2PdmuB2TotHRG,"%3.16f	%3.1f	%12.16f \n", 
								muBval, Tval, d2PressTotHRGdmuB2Mat[i][j]);
  		fprintf(Filed2PdTdmuBTotHRG,"%3.16f	%3.1f	%12.16f \n", 
								muBval, Tval, d2PressTotHRGdTdmuBMat[i][j]);
		}
	}
	if (print){
		fclose(FilePressTotHRG);	
		fclose(FiledPdTTotHRG); 		
		fclose(FiledPdmuBTotHRG); 
		fclose(Filed2PdT2TotHRG);	
		fclose(Filed2PdmuB2TotHRG);		
		fclose(Filed2PdTdmuBTotHRG);
	}	

	fprintf(stderr,"Merged with HRG.\n");
} /*}}}*/

void out_thermodynamics() /*{{{*/
{
	/* Now all the derivatives have been merged with the HRG, we can finish with 
	 * the thermodynamics and export the quantities we want, in particular we 
	 * will normalize with suitable powers of the temperature to obtain 
	 * dimensionless quantities.*/
	fprintf(stderr,"Merging with HRG complete. \n"
								 "Calculating thermodynamics quantities.\n");
	/* Calculate thermodynamic quantities */
	FILE * FilePressFinal = fopen(namePressFinal3D, "w");			
	FILE * FileEntrFinal = fopen(nameEntrFinal3D, "w");
	FILE * FileBarDensFinal = fopen(nameBarDensFinal3D, "w");	
	FILE * FileEnerDensFinal = fopen(nameEnerDensFinal3D, "w");
	FILE * FileSpSoundFinal = fopen(nameSpsoundFinal3D, "w");	
	FILE * FileChi2Final = fopen(nameChi2Final3D, "w");
	FILE * FileCorrLengthFinal = fopen(nameCorrLengthFinal3D, "w");
	
	for (i=lowT_out; i<=highT_out; i++) for(j=0;j<=highMU_out; j++){
		Tval = (double) i; muBval = (double) j;

   	PressFinalMat[i][j] = PressTotHRGMat[i][j]/pow(Tval,4);
   	EntropyFinalMat[i][j] = dPressTotHRGdTMat[i][j]/pow(Tval,3);
   	BarDensityFinalMat[i][j] = dPressTotHRGdmuBMat[i][j]/pow(Tval,3);
   	EnerDensityFinalMat[i][j] = 
			(dPressTotHRGdTMat[i][j]*Tval - PressTotHRGMat[i][j] 
				+ muBval*dPressTotHRGdmuBMat[i][j])/pow(Tval,4);
   	SpSoundFinalMat[i][j] = 
			(dPressTotHRGdmuBMat[i][j]*dPressTotHRGdmuBMat[i][j]
				*d2PressTotHRGdT2Mat[i][j] 
			 - 2.0*dPressTotHRGdTMat[i][j]*dPressTotHRGdmuBMat[i][j]
				*d2PressTotHRGdTdmuBMat[i][j] 
			 + dPressTotHRGdTMat[i][j]*dPressTotHRGdTMat[i][j]
				*d2PressTotHRGdmuB2Mat[i][j])
   		/(dPressTotHRGdTMat[i][j]*Tval + muBval*dPressTotHRGdmuBMat[i][j])
	   	/(d2PressTotHRGdT2Mat[i][j]*d2PressTotHRGdmuB2Mat[i][j]
				-d2PressTotHRGdTdmuBMat[i][j]*d2PressTotHRGdTdmuBMat[i][j]);
    Chi2FinalMat[i][j] = d2PressTotHRGdmuB2Mat[i][j]/pow(Tval,2);

		fprintf(FilePressFinal,"%3.16f	%3.1f	%12.16f \n", 
							muBval, Tval, PressFinalMat[i][j]);
 		fprintf(FileEntrFinal,"%3.16f	%3.1f	%12.16f \n", 
							muBval, Tval, EntropyFinalMat[i][j]);
 		fprintf(FileBarDensFinal,"%3.16f	%3.1f	%12.16f \n", 
							muBval, Tval, BarDensityFinalMat[i][j]);
 		fprintf(FileEnerDensFinal,"%3.16f	%3.1f	%12.16f \n", 
							muBval, Tval, EnerDensityFinalMat[i][j]);
	 	fprintf(FileSpSoundFinal,"%3.16f	%3.1f	%12.16f \n", 
							muBval, Tval, SpSoundFinalMat[i][j]);
		fprintf(FileChi2Final,"%3.16f	%3.1f	%12.16f \n", 
							muBval, Tval, Chi2FinalMat[i][j]);
 		fprintf(FileCorrLengthFinal,"%3.16f	%3.1f	%12.16f \n", 
							muBval, Tval, CorrLengthMergedMat[i][j]);
 	}
	fclose(FilePressFinal);			
	fclose(FileEntrFinal);			
	fclose(FileBarDensFinal);
	fclose(FileEnerDensFinal);	
	fclose(FileSpSoundFinal);		
	fclose(FileChi2Final);			
	fclose(FileCorrLengthFinal);
} /*}}}*/


void free_alloc() /*{{{*/
{
	fprintf(stderr,"Free memory allocation.\n");
	/* Remove coordinate file */
	remove(nameCoords);
  /* Vectors for root finding */
	free_vector(x,1,N);
	free_vector(f,1,N);
  /* Jacobian */
	free_matrix(JJ,1,N,1,N);
  /* Tensor for coordinates*/
	free_f3tensor(Coords,lowT,highT,0,2*highMU,1,N);
  /* Vectors */
	free_vector(Chi0LatVec,lowT,highT);			
	free_vector(Chi2LatVec,lowT,highT);			
	free_vector(Chi4LatVec,lowT,highT);		
	free_vector(Chi0IsingVec,lowT,highT);		
	free_vector(Chi2IsingVec,lowT,highT);		
	free_vector(Chi4IsingVec,lowT,highT);	
	free_vector(Chi0NoIsingVec,lowT,highT);	
	free_vector(Chi2NoIsingVec,lowT,highT);	
	free_vector(Chi4NoIsingVec,lowT,highT);	
	/* For LAT-only mode */ 
	free_matrix(PressLATonlyMat,lowT,highT,0,highMU);					
	free_matrix(PressLATonlyFilterMat,lowT,highT,0,highMU);		
	free_matrix(dPressLATonlydTMat,lowT,highT,0,highMU);		
	free_matrix(dPressLATonlydmuBMat,lowT,highT,0,highMU);		
	free_matrix(d2PressLATonlydT2Mat,lowT,highT,0,highMU);		
	free_matrix(d2PressLATonlydmuB2Mat,lowT,highT,0,highMU);	
	free_matrix(d2PressLATonlydTdmuBMat,lowT,highT,0,highMU);	
  /* For LAT-only mode, normalized */
	free_matrix(PressLATonlyNormMat,lowT,highT,0,highMU);				
	free_matrix(EntropyLATonlyNormMat,lowT,highT,0,highMU);		
	free_matrix(BarDensityLATonlyNormMat,lowT,highT,0,highMU);
	free_matrix(EnerDensityLATonlyNormMat,lowT,highT,0,highMU);	
	free_matrix(SpSoundLATonlyNormMat,lowT,highT,0,highMU);		
	free_matrix(Chi2LATonlyNormMat,lowT,highT,0,highMU);
  /* Ising */ 
	free_matrix(PressNoIsingMat,lowT,highT,0,highMU);				
	free_matrix(dPressNoIsingdTMat,lowT,highT,0,highMU);			
	free_matrix(dPressNoIsingdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressNoIsingdT2Mat,lowT,highT,0,highMU);	
	free_matrix(d2PressNoIsingdmuB2Mat,lowT,highT,0,highMU);	
	free_matrix(d2PressNoIsingdTdmuBMat,lowT,highT,0,highMU);
  /* No-Ising filtered */
	free_matrix(PressNoIsingFilterMat,lowT,highT,0,highMU);				
	free_matrix(dPressNoIsingFilterdTMat,lowT,highT,0,highMU);		
	free_matrix(dPressNoIsingFilterdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressNoIsingFilterdT2Mat,lowT,highT,0,highMU);	
	free_matrix(d2PressNoIsingFilterdmuB2Mat,lowT,highT,0,highMU);	
	free_matrix(d2PressNoIsingFilterdTdmuBMat,lowT,highT,0,highMU);
  /* Ising */
	free_matrix(PressIsingMat,lowT,highT,0,highMU);			    
	free_matrix(dPressIsingdTMat,lowT,highT,0,highMU);			
	free_matrix(dPressIsingdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressIsingdT2Mat,lowT,highT,0,highMU);		
	free_matrix(d2PressIsingdmuB2Mat,lowT,highT,0,highMU);		
	free_matrix(d2PressIsingdTdmuBMat,lowT,highT,0,highMU);
  /* Total */
	free_matrix(PressTotMat,lowT,highT,0,highMU);			    
	free_matrix(dPressTotdTMat,lowT,highT,0,highMU);			
	free_matrix(dPressTotdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressTotdT2Mat,lowT,highT,0,highMU);		    
	free_matrix(d2PressTotdmuB2Mat,lowT,highT,0,highMU);		
	free_matrix(d2PressTotdTdmuBMat,lowT,highT,0,highMU);
  /* HRG */
	free_matrix(PressHRGMat,lowT,highT,0,highMU);			    
	free_matrix(dPressHRGdTMat,lowT,highT,0,highMU);			
	free_matrix(dPressHRGdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressHRGdT2Mat,lowT,highT,0,highMU);		    
	free_matrix(d2PressHRGdmuB2Mat,lowT,highT,0,highMU);		
	free_matrix(d2PressHRGdTdmuBMat,lowT,highT,0,highMU);
  /* Total+HRG */
	free_matrix(PressTotHRGMat,lowT,highT,0,highMU);		    
	free_matrix(dPressTotHRGdTMat,lowT,highT,0,highMU);			
	free_matrix(dPressTotHRGdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressTotHRGdT2Mat,lowT,highT,0,highMU);	    
	free_matrix(d2PressTotHRGdmuB2Mat,lowT,highT,0,highMU);		
	free_matrix(d2PressTotHRGdTdmuBMat,lowT,highT,0,highMU);
	/* Final output */
	free_matrix(PressFinalMat,lowT,highT,0,highMU);			    
	free_matrix(EntropyFinalMat,lowT,highT,0,highMU);			
	free_matrix(BarDensityFinalMat,lowT,highT,0,highMU);
	free_matrix(EnerDensityFinalMat,lowT,highT,0,highMU);	    
	free_matrix(SpSoundFinalMat,lowT,highT,0,highMU);			
	free_matrix(Chi2FinalMat,lowT,highT,0,highMU);	
	/* Correlation length */
	free_matrix(CorrLengthEpsExpMat,lowT,highT,0,highMU);
	free_matrix(CorrLengthAsymptMat,lowT,highT,0,highMU);
	free_matrix(CorrLengthMergedMat,lowT,highT,0,highMU);
} /*}}}*/



int get_options(int argc, char *argv[]) /*{{{*/
{	
	int o;
	while ((o=getopt(argc, argv, "r:R:vpLSh")) != -1) {
		switch (o) {
			case 'v': verbose++; break;
			case 'p': print++; break;
			case 'S': strange++; break;
			case 'L': lattice++; break;
			case 'r': sscanf(optarg,"%d%*[:,-]%d",&lowT_out,&highT_out); break;
			case 'R': sscanf(optarg,"%d%*[:,-]%d",&lowMU_out,&highMU_out); break;
			case 'h': help(); exit(EXIT_FAILURE);
			default:	fprintf(stderr,"Run help [-h] for help.\n");
								exit(EXIT_FAILURE);
		}
	}
	if (!verbose) if(freopen("/dev/null","w",stderr)==NULL){
		printf("Could not suppress stderr. Exiting. \n");
		exit(EXIT_FAILURE);
	}
	if (lowT_out > highT_out || lowMU_out > highMU_out){
		printf("Wrong ordering in required output range!. Exiting. \n");
		exit(EXIT_FAILURE);
	}
	if (lowT_out < 30 || highT_out > 800 || lowMU_out < 0 || highMU_out > 600){
		printf("Error in setting the output range. \nYou cannot extend the"
					" range beyond:\n"
					"               T = 30-800 MeV\n"
					"               mu = 0-600 MeV\n"
					"Exiting. \n");
		exit(EXIT_FAILURE);
	}

	return 0;
} /*}}}*/

// The main
int main(int argc, char *argv[])
{
	char buff[FILENAME_MAX];

	get_options(argc,argv);

	allocate(); /* allocations for vectors and matrices */


 	/* Assign the name of the main folder where the program lives and 
	 * the files we wish to import are located. */
	getcwd(buff,FILENAME_MAX);	/* store current working dir in buf */

	read_param(argv[optind]);		/* read parameter file */

	if (!strange) read_lattice_data(lat_data);	/* import lattice coefficients */
	else read_lattice_data(lat_data_SN); 
	// Move into folder for output
	chdir("output");

	if (lattice){
		do_lattice();	/* if in lattice mode, run it */
		return 0;
	} 
	write_filenames();			/* determine file names from parameter choice */

	// Create folder with name depending on parameters, and move to the folder.
	mkdir(nameFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	chdir(nameFolder);

	calc_jacobian();				/* calculate Jacobian of (r,h)->(T,mu) map */
	make_coords(); 					/* generate coordinates */
	make_Ising(); 					/* generate Ising fluctuations */
	merge_CorrLength();			/* Merge correlation legth expressions */
	make_nonIsing(); 				/* determine non-Ising coefficients */
	merge_Ising_nonIsing(); /* merge Ising and non-Ising */
	merge_HRG(); 						/* merge with HRG */
	out_thermodynamics();		/* output thermodynamic quantities */

	remove(nameCoords); 		/* remove coordinate file */
	free_alloc();						/* free allocations */

	return 0;
}
#undef NRANSI
