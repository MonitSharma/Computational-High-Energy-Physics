//  Variables.h

#ifndef Variables_h
#define Variables_h

#include <time.h>

/* Define constants appearing in the parametrization 
 * (R,Theta) -> (r,h), critical exponents, 
 * and combinations thereof (the c's coefficients) */
#define A -0.76201
#define B 0.00804
#define beta 0.326
#define delta 4.8
#define alpha 0.11
#define M0 0.605
#define h0 0.364
#define c0 0.0424369
#define c1 0.321329
#define c2 -1.20375
#define c3 -0.00126032
#define bede 1.5648
#define I -2.34390
#define epsilon 1
#define aa 1
#define nu 0.63

#define MEDIUM_INT 512
#define LONG_INT 1024

// Define some chars.
char MODESTR[MEDIUM_INT];

// Clock variables for keeping track ot time expense.
clock_t start, end;
clock_t start1, end1;
clock_t start2, end2;
clock_t start3, end3;
double cpu_time_used;


// Limits in the phase diagram
int lowT, highT;
int lowMU, highMU;
int lowT_out, highT_out;
int lowMU_out, highMU_out;


// Define strings for filenames; used for exporting files.
// Name for coordinate file.
char nameCoords[LONG_INT];
// Name for folder to be created, depending on the parameters.
char nameFolder[LONG_INT];
// Name for Chis at muB=0, and derivatives wrt T.
char nameChisIsing[LONG_INT]; 
char nameChisNoIsing[LONG_INT]; 
char namedChisNoIsingdT[LONG_INT]; 
char named2ChisNoIsingdT2[LONG_INT];
// Name for thermodynamic quantities - Ising.
char namePressIsing3D[LONG_INT]; 
char nameCorrLenIsing[LONG_INT]; 
char namedPdTIsing3D[LONG_INT]; 
char namedPdmuBIsing3D[LONG_INT]; 
char named2PdT2Ising3D[LONG_INT]; 
char named2PdmuB2Ising3D[LONG_INT]; 
char named2PdTdmuBIsing3D[LONG_INT];
// Name for thermodynamic quantities - Non Ising.
char namePressNoIsing3D[LONG_INT]; 
char namedPdTNoIsing3D[LONG_INT]; 
char namedPdmuBNoIsing3D[LONG_INT]; 
char named2PdT2NoIsing3D[LONG_INT]; 
char named2PdmuB2NoIsing3D[LONG_INT]; 
char named2PdTdmuBNoIsing3D[LONG_INT];
// Name for thermodynamic quantities - Total (Ising+Non Ising, not normalized).
char namePressIsingPlusNoIsing3D[LONG_INT]; 
char namedPdTIsingPlusNoIsing3D[LONG_INT]; 
char namedPdmuBIsingPlusNoIsing3D[LONG_INT];
char named2PdT2IsingPlusNoIsing3D[LONG_INT]; 
char named2PdmuB2IsingPlusNoIsing3D[LONG_INT]; 
char named2PdTdmuBIsingPlusNoIsing3D[LONG_INT];
// Name for thermodynamic quantities - Total + HRG.
char namePressTotHRG3D[LONG_INT]; 
char namedPdTTotHRG3D[LONG_INT]; 
char namedPdmuBTotHRG3D[LONG_INT]; 
char named2PdT2TotHRG3D[LONG_INT]; 
char named2PdmuB2TotHRG3D[LONG_INT]; 
char named2PdTdmuBTotHRG3D[LONG_INT];
// Name for thermodynamic quantities - Final (normalized).
char namePressFinal3D[LONG_INT]; 
char nameEnerDensFinal3D[LONG_INT]; 
char nameBarDensFinal3D[LONG_INT]; 
char nameEntrFinal3D[LONG_INT]; 
char nameSpsoundFinal3D[LONG_INT]; 
char nameChi2Final3D[LONG_INT]; 
char nameCorrLengthFinal3D[LONG_INT];

/* Dummy Variables for parameter input*/
double v1,v2,v3,v4,v5,v6;
char v0[LONG_INT], p0[LONG_INT];

/* Parameters of the Ising -> QCD map */
double dTC; double dmuBC; double TC; double muBC; double angle1; double angle2; double ww; double rho;

/* Parameters used in case the choice of parametrization is such that the CEP lies on a straight line or a parabola. */
double T0; double  kappa; double  anglediff; double  TC_in; double  angle1_in;

/* Parameters that take on the values obtained from the inverse Jacobian of the Ising -> QCD map. */
double drdmuB; double dhdmuB; double drdT; double dhdT;

/* Dummy variable for data i/o */
double xIn1; double xIn2; double xIn3; double xIn4; double xIn5;

/* Integers that will be used as indices. */
int i,j,k,l,x1int,x2int,kint,lint;

/* Vectors and integer used in the root finding routine. */
double *x,*f;
int check;

/* Variables used for the Gauss filter. */
int kmin,kmax,lmin,lmax;
double sigmax; double sigmay; double sum; double norm;

/* Variables for correlation length merging */
double Ttip, muBtip, Tlow, TMlow, T0low, klow, Thigh, T0high, khigh, Tmerg, TMmerg, T0merg, kmerg;

/* Define the dummy variables for the root finding "for" loop. */ 
double Ti; double  muBi; double  r_i; double  h_i;

/* Variables used to export/use T and muB as doubles. */
double Tval; double muBval; double cure; double mapsign;

/* Variables used in the merging with the HRG pressure. */
double Tprime; double Targum; double deltaT; double deltaTMax; double deltaTMin;

/* These variables are defined to be used inside the expression for the G and derivatives thereof, to speed up the process. */
double gg0Num; double  gg1Num; double  gg2Num; double  htil0Num; double  htil1Num; double  htil2Num;
double d2GdR2Num; double  d2GdTheta2Num; double  d2GdRdThetaNum; double  dGdRNum; double  dGdThetaNum;
double dMdRNum; double  dMdThetaNum;
double dGdrConhNum; double  dGdhConrNum;
double dThetadRConhNum; double  dRdThetaConhNum; double  dThetadRConrNum; double  dRdThetaConrNum;
double drdRConhNum; double  drdThetaConhNum; double  dhdRConrNum; double  dhdThetaConrNum; double  d2rdR2ConhNum; double  d2rdTheta2ConhNum; double  d2hdR2ConrNum; double  d2hdTheta2ConrNum;
double dRdrConhNum; double  dThetadrConhNum; double  dRdhConrNum; double  dThetadhConrNum; double  d2Rdr2ConhNum; double  d2Thetadr2ConhNum; double  d2Rdh2ConrNum; double  d2Thetadh2ConrNum;
double d2RdrdhNum; double  d2ThetadrdhNum;
double d2Gdr2ConhNum; double  d2Gdh2ConrNum; double  d2GdrdhNum;
double d3Gdr3ConhNum; double  d3Gdh3ConrNum; double  d3Gdr2dhNum; double  d3Gdrdh2Num;
double g0; double htilnum; double htilprimenum; double htilsecondnum; double htilthirdnum; double m1Th2; double ggnum; double ggprimenum; double ggsecondnum; double ggthirdnum; double 
		tau; double ggamma; double omega; double omegap; double ggammap; double zeta; double zetap; double Th2; double R1malpha; double R1momegap; 

/* These vectors are used to store the expansion coefficients of the pressure at different T's, 
		as well as derivatives of the non-Ising one wrt T. */
double *Chi0LatVec, *Chi2LatVec; double  *Chi4LatVec; double  *Chi6LatVec;
double *Chi0IsingVec; double  *Chi2IsingVec; double  *Chi4IsingVec; double  *Chi6IsingVec;
double *Chi0NoIsingVec; double  *Chi2NoIsingVec; double  *Chi4NoIsingVec; double  *Chi6NoIsingVec;
double *dChi0NoIsingdTVec; double  *dChi2NoIsingdTVec; double  *dChi4NoIsingdTVec; double  *dChi6NoIsingdTVec;
double *d2Chi0NoIsingdT2Vec; double  *d2Chi2NoIsingdT2Vec; double  *d2Chi4NoIsingdT2Vec; double  *d2Chi6NoIsingdT2Vec;

/* The variables are used to export the Chi's from the Ising model; double  and the 3D critical pressure. */
double Chi0Ising; double  Chi2Ising; double  Chi4Ising; double  Pij; double  dPdTij; double  dPdmuBij; double  d2PdT2ij; double  d2PdmuB2ij; double  d2PdTdmuBij;

/* Matrix used for the inverse Jacobian of the Ising -> QCD map. */
double **JJ;

/* Matrices to store/export thermodynamics. */ /*{{{*/
/* Non-Ising quantities */
double **PressNoIsingMat;
double **dPressNoIsingdTMat; 
double **d2PressNoIsingdT2Mat; 
double **dPressNoIsingdmuBMat; 
double **d2PressNoIsingdmuB2Mat;
double **d2PressNoIsingdTdmuBMat;
double **PressNoIsingFilterMat; 
/* Non-Ising quantities (after filtering) */
double **dPressNoIsingFilterdTMat; 
double **d2PressNoIsingFilterdT2Mat; 
double **dPressNoIsingFilterdmuBMat; 
double **d2PressNoIsingFilterdmuB2Mat;
double **d2PressNoIsingFilterdTdmuBMat; 
/* Ising quantities */
double **PressIsingMat;
double **dPressIsingdTMat; 
double **d2PressIsingdT2Mat;
double **dPressIsingdmuBMat;
double **d2PressIsingdmuB2Mat;
double **d2PressIsingdTdmuBMat; 
/* Total quantities */
double **PressTotMat; 
double **dPressTotdTMat; 
double **d2PressTotdT2Mat; 
double **dPressTotdmuBMat; 
double **d2PressTotdmuB2Mat; 
double **d2PressTotdTdmuBMat; 
/* HRG quantities */
double **PressHRGMat; 
double **dPressHRGdTMat; 
double **d2PressHRGdT2Mat; 
double **dPressHRGdmuBMat; 
double **d2PressHRGdmuB2Mat; 
double **d2PressHRGdTdmuBMat; 
/* Total + HRG quantities */
double **PressTotHRGMat; 
double **dPressTotHRGdTMat; 
double **d2PressTotHRGdT2Mat; 
double **dPressTotHRGdmuBMat; 
double **d2PressTotHRGdmuB2Mat; 
double **d2PressTotHRGdTdmuBMat; 
double **dPdTMat; 
double **dPdmuBMat; 
double **d2PdT2Mat; 
double **d2PdmuB2Mat; 
double **d2PdTdmuBMat; 
/* Final quantities */
double **PressFinalMat; 
double **EntropyFinalMat; 
double **BarDensityFinalMat; 
double **EnerDensityFinalMat; 
double **SpSoundFinalMat; 
double **Chi2FinalMat; 
/* Special case: correlation length */
double **CorrLengthEpsExpMat;
double **CorrLengthAsymptMat;
double **CorrLengthMergedMat;
/* Lattice-only quantities */
double **PressLATonlyMat; 
double **PressLATonlyFilterMat; 
double **dPressLATonlydTMat; 
double **dPressLATonlydmuBMat; 
double **d2PressLATonlydT2Mat; 
double **d2PressLATonlydmuB2Mat; 
double **d2PressLATonlydTdmuBMat; 
double **PressLATonlyNormMat; 
double **EntropyLATonlyNormMat; 
double **BarDensityLATonlyNormMat; 
double **EnerDensityLATonlyNormMat; 
double **SpSoundLATonlyNormMat; 
double **Chi2LATonlyNormMat; 
/*}}}*/

/* A Rank=3 tensor to store the coordinates. */
double ***Coords;

#endif /* Variables_h */
