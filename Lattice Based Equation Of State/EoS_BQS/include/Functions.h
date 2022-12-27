// Header file for the functions

#ifndef Functions_h
#define Functions_h

double **parMatrix;
double *par;

// ------------- THE FUNCTIONAL FORMS --------------- //
double coeff(double *par,double x); 
double coeffprime(double *par,double x);

// Theese are where the respective coefficients of the lattice data will be strored
// Here chi200 represents chi^b and chi310 represent chi^bq and similarly for the other points
// We have restricted ourselves to fourth order, hence for no value of i,j,k the following constraint
// gets violated i + j + k <= 4
double CHI000(double T);
double CHI200(double T);
double CHI020(double T);
double CHI002(double T);
double CHI110(double T);
double CHI101(double T);
double CHI011(double T);
double CHI400(double T);
double CHI040(double T);
double CHI004(double T);
double CHI310(double T);
double CHI301(double T);
double CHI031(double T);
double CHI130(double T);
double CHI103(double T);
double CHI013(double T);
double CHI220(double T);
double CHI202(double T);
double CHI022(double T);
double CHI211(double T);
double CHI121(double T);
double CHI112(double T);

// Here we store all the respective derivates with respect to temperature
// this is only for the first order derivatives
double DCHI000DT(double T);
double DCHI200DT(double T);
double DCHI020DT(double T);
double DCHI002DT(double T);
double DCHI110DT(double T);
double DCHI101DT(double T);
double DCHI011DT(double T);
double DCHI400DT(double T);
double DCHI040DT(double T);
double DCHI004DT(double T);
double DCHI310DT(double T);
double DCHI301DT(double T);
double DCHI031DT(double T);
double DCHI130DT(double T);
double DCHI103DT(double T);
double DCHI013DT(double T);
double DCHI220DT(double T);
double DCHI202DT(double T);
double DCHI022DT(double T);
double DCHI211DT(double T);
double DCHI121DT(double T);
double DCHI112DT(double T);

// Here we store the Double Derivatives wrt T

double D2CHI000DT2(double T);
double D2CHI200DT2(double T);
double D2CHI020DT2(double T);
double D2CHI002DT2(double T);
double D2CHI110DT2(double T);
double D2CHI101DT2(double T);
double D2CHI011DT2(double T);
double D2CHI400DT2(double T);
double D2CHI040DT2(double T);
double D2CHI004DT2(double T);
double D2CHI310DT2(double T);
double D2CHI301DT2(double T);
double D2CHI031DT2(double T);
double D2CHI130DT2(double T);
double D2CHI103DT2(double T);
double D2CHI013DT2(double T);
double D2CHI220DT2(double T);
double D2CHI202DT2(double T);
double D2CHI022DT2(double T);
double D2CHI211DT2(double T);
double D2CHI121DT2(double T);
double D2CHI112DT2(double T);


// Here are all the functions that will store all the 
// thermodynamical quantities derived further by the code
double PressTaylor(double T, double muB, double muQ, double muS);          // This one contains the pressure term
double EntrTaylor(double T, double muB, double muQ, double muS);           // the entropy term
double BarDensTaylor(double T, double muB, double muQ, double muS);        // the baryon density
double StrDensTaylor(double T, double muB, double muQ, double muS);        // the strangeness density
double ChDensTaylor(double T, double muB, double muQ, double muS);         // and so on 

double Chi2BTaylor(double T, double muB, double muQ, double muS);
double Chi2QTaylor(double T, double muB, double muQ, double muS);
double Chi2STaylor(double T, double muB, double muQ, double muS);
double Chi11BQTaylor(double T, double muB, double muQ, double muS);
double Chi11BSTaylor(double T, double muB, double muQ, double muS);
double Chi11QSTaylor(double T, double muB, double muQ, double muS);

double DBarDensDTTaylor(double T, double muB, double muQ, double muS);        // double derivative for various
double DChDensDTTaylor(double T, double muB, double muQ, double muS);         // quantities
double DStrDensDTTaylor(double T, double muB, double muQ, double muS);        


double DEntrDTTaylor(double T, double muB, double muQ, double muS);          // entropy double derivative

double SpSound(double T, double muB, double muQ, double muS);               // speed of sound 

#endif /* Functions_h */