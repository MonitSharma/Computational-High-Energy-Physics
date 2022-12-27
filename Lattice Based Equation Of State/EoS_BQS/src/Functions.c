#include <math.h>
#include <stdio.h>
#define NRANSI
#include "../include/Variables.h"


// ------------- THE FUNCTIONAL FORMS --------------- //
// For all CHIS except CHI200
double coeff(double *par,double x){
     return par[21] + (par[1] + par[2]*(154.0/x) + par[3]*(154.0/x)*(154.0/x) + par[4]*(154.0/x)*(154.0/x)*(154.0/x) + par[5]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[6]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[7]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[8]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[9]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[10]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x))
                  *1.0/(par[11] + par[12]*(154.0/x) + par[13]*(154.0/x)*(154.0/x) + par[14]*(154.0/x)*(154.0/x)*(154.0/x)
                + par[15]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[16]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[17]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[18]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)
                + par[19]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[20]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x));
}
double coeffprime(double *par,double x){
    return (- 1.0/x*par[2]*(154.0/x) - 2.0/x*par[3]*(154.0/x)*(154.0/x) - 3.0/x*par[4]*(154.0/x)*(154.0/x)*(154.0/x) - 4.0/x*par[5]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                - 5.0/x*par[6]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) - 6.0/x*par[7]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                - 7.0/x*par[8]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) - 8.0/x*par[9]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                - 9.0/x*par[10]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x))
                     *1.0/(par[11] + par[12]*(154.0/x) + par[13]*(154.0/x)*(154.0/x) + par[14]*(154.0/x)*(154.0/x)*(154.0/x)
                + par[15]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[16]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[17]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[18]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)
                + par[19]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[20]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x))
           - (par[1] + par[2]*(154.0/x) + par[3]*(154.0/x)*(154.0/x) + par[4]*(154.0/x)*(154.0/x)*(154.0/x) + par[5]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[6]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[7]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[8]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[9]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[10]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x))
                    *(- 1.0/x*par[12]*(154.0/x) - 2.0/x*par[13]*(154.0/x)*(154.0/x) - 3.0/x*par[14]*(154.0/x)*(154.0/x)*(154.0/x) - 4.0/x*par[15]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                - 5.0/x*par[16]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) - 6.0/x*par[17]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                - 7.0/x*par[18]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) - 8.0/x*par[19]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)
                - 9.0/x*par[20]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x))
                     *1.0/pow((par[11] + par[12]*(154.0/x) + par[13]*(154.0/x)*(154.0/x) + par[14]*(154.0/x)*(154.0/x)*(154.0/x)
                + par[15]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[16]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[17]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[18]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)
                + par[19]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
                + par[20]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)),2);
}
double coeffsecond(double *par, double x){
   return (-2*(-((154.0*par[2])/pow(x,2)) - (2*pow(154.0,2)*par[3])/pow(x,3) - (3*pow(154.0,3)*par[4])/pow(x,4) - (4*pow(154.0,4)*par[5])/pow(x,5) - (5*pow(154.0,5)*par[6])/pow(x,6) - 
        (6*pow(154.0,6)*par[7])/pow(x,7) - (7*pow(154.0,7)*par[8])/pow(x,8) - (8*pow(154.0,8)*par[9])/pow(x,9) - (9*pow(154.0,9)*par[10])/pow(x,10))*
      (-((154.0*par[12])/pow(x,2)) - (2*pow(154.0,2)*par[13])/pow(x,3) - (3*pow(154.0,3)*par[14])/pow(x,4) - (4*pow(154.0,4)*par[15])/pow(x,5) - (5*pow(154.0,5)*par[16])/pow(x,6) - 
        (6*pow(154.0,6)*par[17])/pow(x,7) - (7*pow(154.0,7)*par[18])/pow(x,8) - (8*pow(154.0,8)*par[19])/pow(x,9) - (9*pow(154.0,9)*par[20])/pow(x,10)))/
    pow(par[11] + (154.0*par[12])/x + (pow(154.0,2)*par[13])/pow(x,2) + (pow(154.0,3)*par[14])/pow(x,3) + (pow(154.0,4)*par[15])/pow(x,4) + (pow(154.0,5)*par[16])/pow(x,5) + 
      (pow(154.0,6)*par[17])/pow(x,6) + (pow(154.0,7)*par[18])/pow(x,7) + (pow(154.0,8)*par[19])/pow(x,8) + (pow(154.0,9)*par[20])/pow(x,9),2) + 
   ((2*154.0*par[2])/pow(x,3) + (6*pow(154.0,2)*par[3])/pow(x,4) + (12*pow(154.0,3)*par[4])/pow(x,5) + 
      (20*pow(154.0,4)*par[5])/pow(x,6) + (30*pow(154.0,5)*par[6])/pow(x,7) + (42*pow(154.0,6)*par[7])/pow(x,8) + (56*pow(154.0,7)*par[8])/pow(x,9) + (72*pow(154.0,8)*par[9])/pow(x,10) + 
      (90*pow(154.0,9)*par[10])/pow(x,11))/
    (par[11] + (154.0*par[12])/x + (pow(154.0,2)*par[13])/pow(x,2) + (pow(154.0,3)*par[14])/pow(x,3) + (pow(154.0,4)*par[15])/pow(x,4) + (pow(154.0,5)*par[16])/pow(x,5) + 
      (pow(154.0,6)*par[17])/pow(x,6) + (pow(154.0,7)*par[18])/pow(x,7) + (pow(154.0,8)*par[19])/pow(x,8) + (pow(154.0,9)*par[20])/pow(x,9)) + 
   (par[1] + (154.0*par[2])/x + (pow(154.0,2)*par[3])/pow(x,2) + (pow(154.0,3)*par[4])/pow(x,3) + (pow(154.0,4)*par[5])/pow(x,4) + (pow(154.0,5)*par[6])/pow(x,5) + 
      (pow(154.0,6)*par[7])/pow(x,6) + (pow(154.0,7)*par[8])/pow(x,7) + (pow(154.0,8)*par[9])/pow(x,8) + (pow(154.0,9)*par[10])/pow(x,9))*
    ((2*pow(-((154.0*par[12])/pow(x,2)) - (2*pow(154.0,2)*par[13])/pow(x,3) - (3*pow(154.0,3)*par[14])/pow(x,4) - (4*pow(154.0,4)*par[15])/pow(x,5) - (5*pow(154.0,5)*par[16])/pow(x,6) - 
           (6*pow(154.0,6)*par[17])/pow(x,7) - (7*pow(154.0,7)*par[18])/pow(x,8) - (8*pow(154.0,8)*par[19])/pow(x,9) - (9*pow(154.0,9)*par[20])/pow(x,10),2))/
       pow(par[11] + (154.0*par[12])/x + (pow(154.0,2)*par[13])/pow(x,2) + (pow(154.0,3)*par[14])/pow(x,3) + (pow(154.0,4)*par[15])/pow(x,4) + (pow(154.0,5)*par[16])/pow(x,5) + 
         (pow(154.0,6)*par[17])/pow(x,6) + (pow(154.0,7)*par[18])/pow(x,7) + (pow(154.0,8)*par[19])/pow(x,8) + (pow(154.0,9)*par[20])/pow(x,9),3) - 
      ((2*154.0*par[12])/pow(x,3) + (6*pow(154.0,2)*par[13])/pow(x,4) + (12*pow(154.0,3)*par[14])/pow(x,5) + (20*pow(154.0,4)*par[15])/pow(x,6) + (30*pow(154.0,5)*par[16])/pow(x,7) + 
         (42*pow(154.0,6)*par[17])/pow(x,8) + (56*pow(154.0,7)*par[18])/pow(x,9) + (72*pow(154.0,8)*par[19])/pow(x,10) + (90*pow(154.0,9)*par[20])/pow(x,11))/
       pow(par[11] + (154.0*par[12])/x + (pow(154.0,2)*par[13])/pow(x,2) + (pow(154.0,3)*par[14])/pow(x,3) + (pow(154.0,4)*par[15])/pow(x,4) + (pow(154.0,5)*par[16])/pow(x,5) + 
         (pow(154.0,6)*par[17])/pow(x,6) + (pow(154.0,7)*par[18])/pow(x,7) + (pow(154.0,8)*par[19])/pow(x,8) + (pow(154.0,9)*par[20])/pow(x,9),2));          
}


// For CHI200
double coeffMod(double *par, double x){
    return par[21] + exp(-par[1]*(200.0/x) - par[2]*(200.0/x)*(200.0/x))*par[3]*(1.0+tanh(par[4]*(x/200.0)+par[5]));
}
double coeffprimeMod(double *par, double x){
    return (1.0/x*par[1]*(200.0/x) + 2.0/x*par[2]*(200.0/x)*(200.0/x))*exp(-par[1]*(200.0/x) - par[2]*(200.0/x)*(200.0/x))*par[3]*(1.0+tanh(par[4]*(x/200.0)+par[5]))
            + exp(-par[1]*(200.0/x) - par[2]*(200.0/x)*(200.0/x))*par[3]*par[4]/200.0*pow(cosh(par[4]*(x/200.0)+par[5]),-2.0);    
}
double coeffsecondMod(double *par, double x){
    return 2.0*exp(- par[2]*(200.0/x)*(200.0/x) - par[1]*(200.0/x))*par[3]*par[4]*(2.0/x*par[2]*(200.0/x)*(200.0/x) + 1.0/x*par[1]*(200.0/x))
                  *pow(cosh(par[5] + par[4]*x/200.0),-2)/200.0 
            - 2.0*exp(- par[2]*(200.0/x)*(200.0/x) - par[1]*(200.0/x))*par[3]*par[4]*par[4]*sinh(par[5] + par[4]*x/200.0)*pow(cosh(par[5] + par[4]*x/200.0),-3)/(200.0*200.0) 
            + par[3]*exp(- par[2]*(200.0/x)*(200.0/x) - par[1]*(200.0/x))*(- 6.0/(x*x)*par[2]*(200.0/x)*(200.0/x) - 2.0/(x*x)*par[1]*(200.0/x) +pow(2.0/x*par[2]*(200.0/x)*(200.0/x) 
               + 1.0/x*par[1]*(200.0/x),2))*(1.0 + tanh(par[5] + par[4]*x/200.0));    
}


// ------------- THE COEFFICIENTS ------------------ //
double CHI000(double T){
    return coeff(CHI000PAR,T);
}
double CHI200(double T){
    return coeffMod(CHI200PAR,T);
}
double CHI020(double T){
    return coeff(CHI020PAR,T);
}
double CHI002(double T){
    return coeff(CHI002PAR,T);
}
double CHI110(double T){
    return coeff(CHI110PAR,T);
}
double CHI101(double T){
    return coeff(CHI101PAR,T);
}
double CHI011(double T){
    return coeff(CHI011PAR,T);
}
double CHI400(double T){
    return coeff(CHI400PAR,T);
}
double CHI040(double T){
    return coeff(CHI040PAR,T);
}
double CHI004(double T){
    return coeff(CHI004PAR,T);
}
double CHI310(double T){
    return coeff(CHI310PAR,T);
}
double CHI301(double T){
    return coeff(CHI301PAR,T);
}
double CHI031(double T){
    return coeff(CHI031PAR,T);
}
double CHI130(double T){
    return coeff(CHI130PAR,T);
}
double CHI103(double T){
    return coeff(CHI103PAR,T);
}
double CHI013(double T){
    return coeff(CHI013PAR,T);
}
double CHI220(double T){
    return coeff(CHI220PAR,T);
}
double CHI202(double T){
    return coeff(CHI202PAR,T);
}
double CHI022(double T){
    return coeff(CHI022PAR,T);
}
double CHI211(double T){
    return coeff(CHI211PAR,T);
}
double CHI121(double T){
    return coeff(CHI121PAR,T);
}
double CHI112(double T){
    return coeff(CHI112PAR,T);
}


// ---- Derivatives wrt T --------- //
double DCHI000DT(double T){
    return coeffprime(CHI000PAR,T);
}
double DCHI200DT(double T){
    return coeffprimeMod(CHI200PAR,T);
}
double DCHI020DT(double T){
    return coeffprime(CHI020PAR,T);
}
double DCHI002DT(double T){
    return coeffprime(CHI002PAR,T);
}
double DCHI110DT(double T){
    return coeffprime(CHI110PAR,T);
}
double DCHI101DT(double T){
    return coeffprime(CHI101PAR,T);
}
double DCHI011DT(double T){
    return coeffprime(CHI011PAR,T);
}
double DCHI400DT(double T){
    return coeffprime(CHI400PAR,T);
}
double DCHI040DT(double T){
    return coeffprime(CHI040PAR,T);
}
double DCHI004DT(double T){
    return coeffprime(CHI004PAR,T);
}
double DCHI310DT(double T){
    return coeffprime(CHI310PAR,T);
}
double DCHI301DT(double T){
    return coeffprime(CHI301PAR,T);
}
double DCHI031DT(double T){
    return coeffprime(CHI031PAR,T);
}
double DCHI130DT(double T){
    return coeffprime(CHI130PAR,T);
}
double DCHI103DT(double T){
    return coeffprime(CHI103PAR,T);
}
double DCHI013DT(double T){
    return coeffprime(CHI013PAR,T);
}
double DCHI220DT(double T){
    return coeffprime(CHI220PAR,T);
}
double DCHI202DT(double T){
    return coeffprime(CHI202PAR,T);
}
double DCHI022DT(double T){
    return coeffprime(CHI022PAR,T);
}
double DCHI211DT(double T){
    return coeffprime(CHI211PAR,T);
}
double DCHI121DT(double T){
    return coeffprime(CHI121PAR,T);
}
double DCHI112DT(double T){
    return coeffprime(CHI112PAR,T);
}


// ---- Double derivatives wrt T --------- //
double D2CHI000DT2(double T){
    return coeffsecond(CHI000PAR,T);
}
double D2CHI200DT2(double T){
    return coeffsecondMod(CHI200PAR,T);
}
double D2CHI020DT2(double T){
    return coeffsecond(CHI020PAR,T);
}
double D2CHI002DT2(double T){
    return coeffsecond(CHI002PAR,T);
}
double D2CHI110DT2(double T){
    return coeffsecond(CHI110PAR,T);
}
double D2CHI101DT2(double T){
    return coeffsecond(CHI101PAR,T);
}
double D2CHI011DT2(double T){
    return coeffsecond(CHI011PAR,T);
}
double D2CHI400DT2(double T){
    return coeffsecond(CHI400PAR,T);
}
double D2CHI040DT2(double T){
    return coeffsecond(CHI040PAR,T);
}
double D2CHI004DT2(double T){
    return coeffsecond(CHI004PAR,T);
}
double D2CHI310DT2(double T){
    return coeffsecond(CHI310PAR,T);
}
double D2CHI301DT2(double T){
    return coeffsecond(CHI301PAR,T);
}
double D2CHI031DT2(double T){
    return coeffsecond(CHI031PAR,T);
}
double D2CHI130DT2(double T){
    return coeffsecond(CHI130PAR,T);
}
double D2CHI103DT2(double T){
    return coeffsecond(CHI103PAR,T);
}
double D2CHI013DT2(double T){
    return coeffsecond(CHI013PAR,T);
}
double D2CHI220DT2(double T){
    return coeffsecond(CHI220PAR,T);
}
double D2CHI202DT2(double T){
    return coeffsecond(CHI202PAR,T);
}
double D2CHI022DT2(double T){
    return coeffsecond(CHI022PAR,T);
}
double D2CHI211DT2(double T){
    return coeffsecond(CHI211PAR,T);
}
double D2CHI121DT2(double T){
    return coeffsecond(CHI121PAR,T);
}
double D2CHI112DT2(double T){
    return coeffsecond(CHI112PAR,T);
}


// ------------ THERMODYNAMICS --------------------- //
double PressTaylor(double T, double muB, double muQ, double muS){
    return CHI000(T) + CHI110(T)*muB/T*muQ/T + CHI101(T)*muB/T*muS/T + CHI011(T)*muQ/T*muS/T + 1.0/2.0*CHI200(T)*muB/T*muB/T + 1.0/2.0*CHI020(T)*muQ/T*muQ/T + 1.0/2.0*CHI002(T)*muS/T*muS/T 
            + 1.0/2.0*CHI211(T)*muB/T*muB/T*muQ/T*muS/T + 1.0/2.0*CHI121(T)*muB/T*muQ/T*muQ/T*muS/T + 1.0/2.0*CHI112(T)*muB/T*muQ/T*muS/T*muS/T + 1.0/4.0*CHI220(T)*muB/T*muB/T*muQ/T*muQ/T 
            + 1.0/4.0*CHI202(T)*muB/T*muB/T*muS/T*muS/T + 1.0/4.0*CHI022(T)*muQ/T*muQ/T*muS/T*muS/T + 1.0/6.0*CHI310(T)*muB/T*muB/T*muB/T*muQ/T + 1.0/6.0*CHI130(T)*muB/T*muQ/T*muQ/T*muQ/T 
            + 1.0/6.0*CHI301(T)*muB/T*muB/T*muB/T*muS/T + 1.0/6.0*CHI103(T)*muB/T*muS/T*muS/T*muS/T + 1.0/6.0*CHI031(T)*muQ/T*muQ/T*muQ/T*muS/T + 1.0/6.0*CHI013(T)*muQ/T*muS/T*muS/T*muS/T 
            + 1.0/24.0*CHI400(T)*muB/T*muB/T*muB/T*muB/T + 1.0/24.0*CHI040(T)*muQ/T*muQ/T*muQ/T*muQ/T + 1.0/24.0*CHI004(T)*muS/T*muS/T*muS/T*muS/T;
}
double EntrTaylor(double T, double muB, double muQ, double muS){
    return  1.0/(24.0*T*T*T)*(96.0*T*T*T*CHI000(T) + 24.0*muS*muS*T*CHI002(T) + 48.0*muQ*muS*T*CHI011(T) + 24.0*muQ*muQ*T*CHI020(T) + 48.0*muB*muS*T*CHI101(T) + 48.0*muB*muQ*T*CHI110(T) 
            + 24.0*muB*muB*T*CHI200(T) + 24.0*T*T*T*T*DCHI000DT(T) + 12.0*muS*muS*T*T*DCHI002DT(T) + muS*muS*muS*muS*DCHI004DT(T) + 24.0*muQ*muS*T*T*DCHI011DT(T) 
            + 4.0*muQ*muS*muS*muS*DCHI013DT(T) + 12.0*muQ*muQ*T*T*DCHI020DT(T) + 6.0*muQ*muQ*muS*muS*DCHI022DT(T) + 4.0*muQ*muQ*muQ*muS*DCHI031DT(T) + muQ*muQ*muQ*muQ*DCHI040DT(T) 
            + 24.0*muB*muS*T*T*DCHI101DT(T) + 4.0*muB*muS*muS*muS*DCHI103DT(T) + 24.0*muB*muQ*T*T*DCHI110DT(T) + 12.0*muB*muQ*muS*muS*DCHI112DT(T) + 12.0*muB*muQ*muQ*muS*DCHI121DT(T) 
            + 4.0*muB*muQ*muQ*muQ*DCHI130DT(T) + 12.0*muB*muB*T*T*DCHI200DT(T) + 6.0*muB*muB*muS*muS*DCHI202DT(T) + 12.0*muB*muB*muQ*muS*DCHI211DT(T) + 6.0*muB*muB*muQ*muQ*DCHI220DT(T) 
            + 4.0*muB*muB*muB*muS*DCHI301DT(T) + 4.0*muB*muB*muB*muQ*DCHI310DT(T) + muB*muB*muB*muB*DCHI400DT(T));
}
double BarDensTaylor(double T, double muB, double muQ, double muS){
    return CHI110(T)*muQ/T + CHI101(T)*muS/T + CHI200(T)*muB/T + CHI211(T)*muB/T*muQ/T*muS/T + 1.0/2.0*CHI121(T)*muQ/T*muQ/T*muS/T + 1.0/2.0*CHI112(T)*muQ/T*muS/T*muS/T 
            + 1.0/2.0*CHI220(T)*muB/T*muQ/T*muQ/T + 1.0/2.0*CHI202(T)*muB/T*muS/T*muS/T + 1.0/2.0*CHI310(T)*muB/T*muB/T*muQ/T + 1.0/6.0*CHI130(T)*muQ/T*muQ/T*muQ/T 
            + 1.0/2.0*CHI301(T)*muB/T*muB/T*muS/T + 1.0/6.0*CHI103(T)*muS/T*muS/T*muS/T + 1.0/6.0*CHI400(T)*muB/T*muB/T*muB/T;
}
double StrDensTaylor(double T, double muB, double muQ, double muS){
    return CHI101(T)*muB/T + CHI011(T)*muQ/T + CHI002(T)*muS/T + 1.0/2.0*CHI211(T)*muB/T*muB/T*muQ/T + 1.0/2.0*CHI121(T)*muB/T*muQ/T*muQ/T + CHI112(T)*muB/T*muQ/T*muS/T 
            + 1.0/2.0*CHI202(T)*muB/T*muB/T*muS/T + 1.0/2.0*CHI022(T)*muQ/T*muQ/T*muS/T + 1.0/6.0*CHI301(T)*muB/T*muB/T*muB/T + 1.0/2.0*CHI103(T)*muB/T*muS/T*muS/T 
            + 1.0/6.0*CHI031(T)*muQ/T*muQ/T*muQ/T + 1.0/2.0*CHI013(T)*muQ/T*muS/T*muS/T + 1.0/6.0*CHI004(T)*muS/T*muS/T*muS/T;
}
double ChDensTaylor(double T, double muB, double muQ, double muS){
    return CHI110(T)*muB/T + CHI011(T)*muS/T + CHI020(T)*muQ/T + 1.0/2.0*CHI211(T)*muB/T*muB/T*muS/T + CHI121(T)*muB/T*muQ/T*muS/T + 1.0/2.0*CHI112(T)*muB/T*muS/T*muS/T 
            + 1.0/2.0*CHI220(T)*muB/T*muB/T*muQ/T + 1.0/2.0*CHI022(T)*muQ/T*muS/T*muS/T + 1.0/6.0*CHI310(T)*muB/T*muB/T*muB/T + 1.0/2.0*CHI130(T)*muB/T*muQ/T*muQ/T 
            + 1.0/2.0*CHI031(T)*muQ/T*muQ/T*muS/T + 1.0/6.0*CHI013(T)*muS/T*muS/T*muS/T + 1.0/6.0*CHI040(T)*muQ/T*muQ/T*muQ/T;
}


// Higher order derivatives, Taylor expanded
//Analytical expressions
// Second order
double Chi2BTaylor(double T, double muB, double muQ, double muS){
    return CHI200(T) + CHI211(T)*muQ/T*muS/T + 1.0/2.0*CHI220(T)*muQ/T*muQ/T + 1.0/2.0*CHI202(T)*muS/T*muS/T + CHI310(T)*muB/T*muQ/T + CHI301(T)*muB/T*muS/T + 1.0/2.0*CHI400(T)*muB/T*muB/T;
}
double Chi2QTaylor(double T, double muB, double muQ, double muS){
    return CHI020(T) + CHI121(T)*muB/T*muS/T + 1.0/2.0*CHI220(T)*muB/T*muB/T + 1.0/2.0*CHI022(T)*muS/T*muS/T + CHI130(T)*muB/T*muQ/T + CHI031(T)*muQ/T*muS/T + 1.0/2.0*CHI040(T)*muQ/T*muQ/T;
}
double Chi2STaylor(double T, double muB, double muQ, double muS){
    return CHI002(T) + CHI112(T)*muB/T*muQ/T + 1.0/2.0*CHI202(T)*muB/T*muB/T + 1.0/2.0*CHI022(T)*muQ/T*muQ/T + CHI103(T)*muB/T*muS/T + CHI013(T)*muQ/T*muS/T + 1.0/2.0*CHI004(T)*muS/T*muS/T;
}
double Chi11BQTaylor(double T, double muB, double muQ, double muS){
    return CHI110(T) + CHI211(T)*muB/T*muS/T + CHI121(T)*muQ/T*muS/T + 1.0/2.0*CHI112(T)*muS/T*muS/T + CHI220(T)*muB/T*muQ/T + 1.0/2.0*CHI310(T)*muB/T*muB/T + 1.0/2.0*CHI130(T)*muQ/T*muQ/T;
}
double Chi11BSTaylor(double T, double muB, double muQ, double muS){
    return CHI101(T) + CHI211(T)*muB/T*muQ/T + 1.0/2.0*CHI121(T)*muQ/T*muQ/T + CHI112(T)*muQ/T*muS/T + CHI202(T)*muB/T*muS/T + 1.0/2.0*CHI301(T)*muB/T*muB/T + 1.0/2.0*CHI103(T)*muS/T*muS/T;
}
double Chi11QSTaylor(double T, double muB, double muQ, double muS){
    return CHI011(T) + 1.0/2.0*CHI211(T)*muB/T*muB/T + CHI121(T)*muB/T*muQ/T + CHI112(T)*muB/T*muS/T + CHI022(T)*muQ/T*muS/T + 1.0/2.0*CHI031(T)*muQ/T*muQ/T + 1.0/2.0*CHI013(T)*muS/T*muS/T;
}

// Second order (one is T)
double DBarDensDTTaylor(double T, double muB, double muQ, double muS){
    return 1.0/(6.0*T*T)*(12.0*muS*T*CHI101(T) + 12.0*muQ*T*CHI110(T) + 12.0*muB*T*CHI200(T) + 6.0*muS*T*T*DCHI101DT(T) + muS*muS*muS*DCHI103DT(T) + 6.0*muQ*T*T*DCHI110DT(T) 
            + 3.0*muQ*muS*muS*DCHI112DT(T) + 3.0*muQ*muQ*muS*DCHI121DT(T) + muQ*muQ*muQ*DCHI130DT(T) + 6.0*muB*T*T*DCHI200DT(T) + 3.0*muB*muS*muS*DCHI202DT(T) 
            + 6.0*muB*muQ*muS*DCHI211DT(T) + 3.0*muB*muQ*muQ*DCHI220DT(T) + 3.0*muB*muB*muS*DCHI301DT(T) + 3.0*muB*muB*muQ*DCHI310DT(T) + muB*muB*muB*DCHI400DT(T));
}
double DStrDensDTTaylor(double T, double muB, double muQ, double muS){
    return 1.0/(6.0*T*T)*(12.0*muS*T*CHI002(T) + 12.0*muQ*T*CHI011(T) + 12.0*muB*T*CHI101(T) + 6.0*muS*T*T*DCHI002DT(T) + muS*muS*muS*DCHI004DT(T) + 6.0*muQ*T*T*DCHI011DT(T) 
            + 3.0*muQ*muS*muS*DCHI013DT(T) + 3.0*muQ*muQ*muS*DCHI022DT(T) + muQ*muQ*muQ*DCHI031DT(T) + 6.0*muB*T*T*DCHI101DT(T) + 3.0*muB*muS*muS*DCHI103DT(T) 
            + 6.0*muB*muQ*muS*DCHI112DT(T) + 3.0*muB*muQ*muQ*DCHI121DT(T) + 3.0*muB*muB*muS*DCHI202DT(T) + 3.0*muB*muB*muQ*DCHI211DT(T) + muB*muB*muB*DCHI301DT(T));
}
double DChDensDTTaylor(double T, double muB, double muQ, double muS){
    return 1.0/(6.0*T*T)*(12.0*muS*T*CHI011(T) + 12.0*muQ*T*CHI020(T) + 12.0*muB*T*CHI110(T) + 6.0*muS*T*T*DCHI011DT(T) + muS*muS*muS*DCHI013DT(T) + 6.0*muQ*T*T*DCHI020DT(T) 
            + 3.0*muQ*muS*muS*DCHI022DT(T) + 3.0*muQ*muQ*muS*DCHI031DT(T) + muQ*muQ*muQ*DCHI040DT(T) + 6.0*muB*T*T*DCHI110DT(T) + 3.0*muB*muS*muS*DCHI112DT(T) 
            + 6.0*muB*muQ*muS*DCHI121DT(T) + 3.0*muB*muQ*muQ*DCHI130DT(T) + 3.0*muB*muB*muS*DCHI211DT(T) + 3.0*muB*muB*muQ*DCHI220DT(T) + muB*muB*muB*DCHI310DT(T));
}

// Second order (both are T)
double DEntrDTTaylor(double T, double muB, double muQ, double muS){
    return (288*pow(T,2)*CHI000(T) + 24*(pow(muS,2)*CHI002(T) + 2*muQ*muS*CHI011(T) + pow(muQ,2)*CHI020(T) + 2*muB*muS*CHI101(T) + 2*muB*muQ*CHI110(T) + pow(muB,2)*CHI200(T)) 
               + 192*pow(T,3)*DCHI000DT(T) + 48*pow(muS,2)*T*DCHI002DT(T) + 96*muQ*muS*T*DCHI011DT(T) + 48*pow(muQ,2)*T*DCHI020DT(T) + 96*muB*muS*T*DCHI101DT(T) + 96*muB*muQ*T*DCHI110DT(T) 
               + 48*pow(muB,2)*T*DCHI200DT(T) + 24*pow(T,4)*D2CHI000DT2(T) + 12*pow(muS,2)*pow(T,2)*D2CHI002DT2(T) + pow(muS,4)*D2CHI004DT2(T) + 24*muQ*muS*pow(T,2)*D2CHI011DT2(T) 
               + 4*muQ*pow(muS,3)*D2CHI013DT2(T) + 12*pow(muQ,2)*pow(T,2)*D2CHI020DT2(T) + 6*pow(muQ,2)*pow(muS,2)*D2CHI022DT2(T) + 4*pow(muQ,3)*muS*D2CHI031DT2(T) + pow(muQ,4)*D2CHI040DT2(T) 
               + 24*muB*muS*pow(T,2)*D2CHI101DT2(T) + 4*muB*pow(muS,3)*D2CHI103DT2(T) + 24*muB*muQ*pow(T,2)*D2CHI110DT2(T) + 12*muB*muQ*pow(muS,2)*D2CHI112DT2(T) 
               + 12*muB*pow(muQ,2)*muS*D2CHI121DT2(T) + 4*muB*pow(muQ,3)*D2CHI130DT2(T) + 12*pow(muB,2)*pow(T,2)*D2CHI200DT2(T) + 6*pow(muB,2)*pow(muS,2)*D2CHI202DT2(T) 
               + 12*pow(muB,2)*muQ*muS*D2CHI211DT2(T) + 6*pow(muB,2)*pow(muQ,2)*D2CHI220DT2(T) + 4*pow(muB,3)*muS*D2CHI301DT2(T) + 4*pow(muB,3)*muQ*D2CHI310DT2(T) 
               + pow(muB,4)*D2CHI400DT2(T))/(24.*pow(T,2));
}


// Speed of Sound expression.
double SpSound(double T, double muB, double muQ, double muS){
   C1B = BarDensTaylor(T,muB,muS,muQ);
   C1Q = ChDensTaylor(T,muB,muS,muQ);
   C1S = StrDensTaylor(T,muB,muS,muQ);
   C1T = EntrTaylor(T,muB,muS,muQ);
   
   C2B2 = Chi2BTaylor(T,muB,muS,muQ);   
   C2Q2 = Chi2QTaylor(T,muB,muS,muQ);
   C2S2 = Chi2STaylor(T,muB,muS,muQ);
   C2BQ = Chi11BQTaylor(T,muB,muS,muQ);
   C2BS = Chi11BSTaylor(T,muB,muS,muQ);
   C2QS = Chi11QSTaylor(T,muB,muS,muQ);
   C2TB = DBarDensDTTaylor(T,muB,muS,muQ);
   C2TQ = DChDensDTTaylor(T,muB,muS,muQ);
   C2TS = DStrDensDTTaylor(T,muB,muS,muQ);
   C2T2 = DEntrDTTaylor(T,muB,muS,muQ);
      
   return T*(-(C2BQ*C2S2*C2TQ*C1B) - C2BQ*C2S2*C2TB*C1Q - pow(C2BS,2)*C2TQ*C1Q + C2B2*C2S2*C2TQ*C1Q + C2BQ*C2BS*C2TS*C1Q + C2BQ*C2BS*C2TQ*C1S - pow(C2BQ,2)*C2TS*C1S + pow(C2BQ,2)*C2S2*C1T 
            + pow(C2QS,2)*(-(C2TB*C1B) + C2B2*C1T) + C2Q2*(C2S2*C2TB*C1B - C2BS*C2TS*C1B - C2BS*C2TB*C1S + C2B2*C2TS*C1S + pow(C2BS,2)*C1T - C2B2*C2S2*C1T) + C2QS*(C2BQ*C2TS*C1B - C2B2*C2TS*C1Q 
            + C2BQ*C2TB*C1S - C2B2*C2TQ*C1S + C2BS*(C2TQ*C1B + C2TB*C1Q - 2.0*C2BQ*C1T)))
               *1.0/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) + C2Q2*C2S2*pow(C2TB,2) - 2.0*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) 
            + 2.0*C2BQ*C2QS*C2TB*C2TS - pow(C2BQ,2)*pow(C2TS,2) + 2.0*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 
            + C2S2*pow(C2TQ,2) - 2.0*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T) 
          + T*(C1B/(muB*C1B + muQ*C1Q + muS*C1S + T*C1T))*(-(C2QS*muB*(C2BQ*(C2TS*C1B + C2TB*C1S) - C2B2*(C2TS*C1Q + C2TQ*C1S) + C2BS*(C2TQ*C1B + C2TB*C1Q - 2*C2BQ*C1T))) + muB*((pow(C2BS,2) 
            - C2B2*C2S2)*C2TQ*C1Q + C2BQ*(C2S2*C2TQ*C1B + C2S2*C2TB*C1Q - C2BS*C2TS*C1Q - C2BS*C2TQ*C1S) + pow(C2BQ,2)*(C2TS*C1S - C2S2*C1T) + C2Q2*(-(C2S2*C2TB*C1B) + C2BS*C2TS*C1B 
            + C2BS*C2TB*C1S - C2B2*C2TS*C1S - pow(C2BS,2)*C1T + C2B2*C2S2*C1T)) + C2QS*(-2*C2TQ*C2TS*C1B + C2TB*C2TS*C1Q + C2TB*C2TQ*C1S - C2T2*(C2BS*C1Q + C2BQ*C1S) + C2BS*C2TQ*C1T 
            + C2BQ*C2TS*C1T)*T + (-((C2BS*C2TQ - C2BQ*C2TS)*(-(C2TS*C1Q) + C2TQ*C1S)) + C2S2*(C2BQ*C2T2*C1Q + C2TQ*(C2TQ*C1B - C2TB*C1Q - C2BQ*C1T)) + C2Q2*(C2BS*C2T2*C1S + C2TS*(C2TS*C1B 
            - C2TB*C1S - C2BS*C1T) + C2S2*(-(C2T2*C1B) + C2TB*C1T)))*T + pow(C2QS,2)*(C2TB*muB*C1B - C2B2*muB*C1T + C2T2*C1B*T - C2TB*C1T*T))/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) 
            + C2Q2*C2S2*pow(C2TB,2) - 2*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) + 2*C2BQ*C2QS*C2TB*C2TS - pow(C2BQ,2)*pow(C2TS,2) + 2*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ 
            - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 + C2S2*pow(C2TQ,2) - 2*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T)
          + T*(C1Q/(muB*C1B + muQ*C1Q + muS*C1S + T*C1T))*(pow(C2QS,2)*muQ*(C2TB*C1B - C2B2*C1T) - C2QS*muQ*(C2BQ*(C2TS*C1B + C2TB*C1S) - C2B2*(C2TS*C1Q + C2TQ*C1S) + C2BS*(C2TQ*C1B + C2TB*C1Q 
            - 2*C2BQ*C1T)) + muQ*((pow(C2BS,2) - C2B2*C2S2)*C2TQ*C1Q + C2BQ*(C2S2*C2TQ*C1B + C2S2*C2TB*C1Q - C2BS*C2TS*C1Q - C2BS*C2TQ*C1S) + pow(C2BQ,2)*(C2TS*C1S - C2S2*C1T) + C2Q2*(-(C2S2*C2TB*C1B) 
            + C2BS*C2TS*C1B + C2BS*C2TB*C1S - C2B2*C2TS*C1S - pow(C2BS,2)*C1T + C2B2*C2S2*C1T)) + C2QS*(-(C2BS*C2T2*C1B) + C2TB*C2TS*C1B + C2B2*C2T2*C1S - pow(C2TB,2)*C1S + C2BS*C2TB*C1T 
            - C2B2*C2TS*C1T)*T + (C2BS*C2TQ*C2TS*C1B + pow(C2BS,2)*C2T2*C1Q - 2*C2BS*C2TB*C2TS*C1Q + C2B2*pow(C2TS,2)*C1Q + C2BS*C2TB*C2TQ*C1S - C2B2*C2TQ*C2TS*C1S - pow(C2BS,2)*C2TQ*C1T 
            + C2S2*(-(C2TB*C2TQ*C1B) - C2B2*C2T2*C1Q + pow(C2TB,2)*C1Q + C2B2*C2TQ*C1T) + C2BQ*(-(C2BS*C2T2*C1S) + C2TS*(-(C2TS*C1B) + C2TB*C1S + C2BS*C1T) + C2S2*(C2T2*C1B - C2TB*C1T)))*T)
               *1.0/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) + C2Q2*C2S2*pow(C2TB,2) - 2*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) + 2*C2BQ*C2QS*C2TB*C2TS 
            - pow(C2BQ,2)*pow(C2TS,2) + 2*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 + C2S2*pow(C2TQ,2) 
            - 2*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T)
          + T*(C1S/(muB*C1B + muQ*C1Q + muS*C1S + T*C1T))*(pow(C2QS,2)*muS*(C2TB*C1B - C2B2*C1T) - C2QS*muS*(C2BQ*(C2TS*C1B + C2TB*C1S) - C2B2*(C2TS*C1Q + C2TQ*C1S) + C2BS*(C2TQ*C1B + C2TB*C1Q 
            - 2*C2BQ*C1T)) + muS*((pow(C2BS,2) - C2B2*C2S2)*C2TQ*C1Q + C2BQ*(C2S2*C2TQ*C1B + C2S2*C2TB*C1Q - C2BS*C2TS*C1Q - C2BS*C2TQ*C1S) + pow(C2BQ,2)*(C2TS*C1S - C2S2*C1T) 
            + C2Q2*(-(C2S2*C2TB*C1B) + C2BS*C2TS*C1B + C2BS*C2TB*C1S - C2B2*C2TS*C1S - pow(C2BS,2)*C1T + C2B2*C2S2*C1T)) + C2QS*(-(C2BQ*C2T2*C1B) + C2TB*C2TQ*C1B + C2B2*C2T2*C1Q - pow(C2TB,2)*C1Q 
            + C2BQ*C2TB*C1T - C2B2*C2TQ*C1T)*T + (C2BQ*C2TQ*C2TS*C1B + C2BQ*C2TB*C2TS*C1Q - C2B2*C2TQ*C2TS*C1Q + pow(C2BQ,2)*C2T2*C1S - 2*C2BQ*C2TB*C2TQ*C1S + C2B2*pow(C2TQ,2)*C1S 
            - pow(C2BQ,2)*C2TS*C1T + C2Q2*(-(C2TB*C2TS*C1B) - C2B2*C2T2*C1S + pow(C2TB,2)*C1S + C2B2*C2TS*C1T) + C2BS*(-(C2BQ*C2T2*C1Q) + C2TQ*(-(C2TQ*C1B) + C2TB*C1Q + C2BQ*C1T) 
            + C2Q2*(C2T2*C1B - C2TB*C1T)))*T)/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) + C2Q2*C2S2*pow(C2TB,2) - 2*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) 
            + 2*C2BQ*C2QS*C2TB*C2TS - pow(C2BQ,2)*pow(C2TS,2) + 2*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 
            + C2S2*pow(C2TQ,2) - 2*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T);
}



#undef NRANSI
