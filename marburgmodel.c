/* pomp model file: seir no seas lbv */

#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <C:/Users/David Hayman/Documents/R/win-library/3.0/pomp/include/pomp.h>

// define parameters

#define BETA        (p[parindex[0]]) // transmission rate
#define MU          (p[parindex[1]]) // adult death rate
#define DELTA       (p[parindex[2]]) // juvenile mortality rate
#define SIGMA		(p[parindex[3]]) // incubation period
#define K			(p[parindex[4]]) // carrying capacity
#define EPSILON		(p[parindex[5]]) // rate of juvenile aging
#define TAU			(p[parindex[6]]) // seroconversion
#define	KAPPA		(p[parindex[7]]) // birth rate (peak size)
#define S			(p[parindex[8]]) // synchrony
#define OMEGA		(p[parindex[9]]) // pulse per yr
#define PHI			(p[parindex[10]]) // timing during year

// define states

#define SUSJ       (x[stateindex[0]]) // number of susceptible juveniles
#define EXPJ		(x[stateindex[1]]) // number of exposed to infectious juveniles
#define INFJ       (x[stateindex[2]]) // number of infected juveniles
#define RECJ       (x[stateindex[3]]) // number of recovered juveniles
#define SUSA       (x[stateindex[4]]) // number of susceptibles adults
#define EXPA		(x[stateindex[5]]) // number of exposed to infectious adults
#define INFA       (x[stateindex[6]]) // number of infected adults
#define RECA       (x[stateindex[7]]) // number of recovered adults
#define PREVA		(x[stateindex[8]]) // adult seroprevalence
#define PREVJ		(x[stateindex[9]]) // juvenile seroprevalence
#define SPA			(x[stateindex[10]]) // adult seroprevalence
#define SPJ			(x[stateindex[11]]) // juvenile seroprevalence

// the process model:
// an SIR model with Euler-multinomial step,

void sir_euler_simulator (double *x, const double *p, 
			  const int *stateindex, const int *parindex, const int *covindex,
			  int covdim, const double *covar, 
			  double t, double dt)
{
  int nrate = 29; 			// number of rates
  double rate[nrate];		// transition rates
  double trans[nrate];		// transition numbers
  double N = x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7];		// population size
  void (*reulmult)(int,double,double*,double,double*);

  // to evaluate the basis functions and compute the transmission rate, use some of 
  // pomp's C-level eulermultinomial simulator
  reulmult = (void (*)(int,double,double*,double,double*)) R_GetCCallable("pomp","reulermultinom");

  // define pi

const double pi = 3.14159265359;

// in C --- pow(a,b) to do a^b 

  // compute the transition rates
  rate[0] = (KAPPA*(1/sqrt((1/S)*pi)*exp(-pow((cos(pi*OMEGA*t-PHI)),2)/(1/S))))*(N);		// approx delta function birth into susceptible class
  rate[1] = EPSILON;		// aging from juv to ad
  rate[2] = BETA*(INFJ+INFA)/N;			// sus to exposed- infectious route
  rate[3] = DELTA*(N/K);			// density dept mortality - Juvenile   
  rate[4] = SIGMA;		// incubation rate
  rate[5] = DELTA*(N/K);			// density dept mortality - Juvenile 
  rate[6] = EPSILON;		// aging from juv to ad
  rate[7] = TAU; // seroconversion rate 
  rate[8] = EPSILON;		// aging from juv to ad
  rate[9] = DELTA*(N/K);			// density dept mortality - Juvenile
  rate[10] = EPSILON;		// aging from juv to ad
  rate[11] = MU*(N/K);			// density dept mortality - Adult
  rate[12] = BETA*(INFJ+INFA)/N;			// sus to exposed- infectious route
  rate[13] = MU*(N/K);			// density dept mortality - Adult
  rate[14] = SIGMA;		// incubation rate
  rate[15] = MU*(N/K);			// density dept mortality - adult
  rate[16] = TAU; // seroconversion rate 
  rate[17] = MU*(N/K);			// density dept mortality - adult
  rate[18] = MU*(N/K);			// density dept mortality - adult
    
  // compute the transition numbers // in reulmult, first # is transitions, state name, rate # trans start, trans name at start
  trans[0] = rpois(rate[0]*dt);	// births are Poisson // is this correct still with ADF birth pulse?
  (*reulmult)(3,SUSJ, &rate[1],dt,&trans[1]); // euler-multinomial exits from SUSJ class - 4, death, aging, infection (2 routes, rho, 1-rho)
  (*reulmult)(3,EXPJ,&rate[4],dt,&trans[4]); // euler-multinomial exits from EJ class
  (*reulmult)(3,INFJ,&rate[7],dt,&trans[7]); // euler-multinomial exits from INFJ class
  (*reulmult)(2,RECJ,&rate[10],dt,&trans[10]); // euler-multinomial exits from RECJ class
  (*reulmult)(2,SUSA,&rate[11],dt,&trans[11]); // euler-multinomial exits from SUSA class
  (*reulmult)(2,EXPA,&rate[13],dt,&trans[13]); // euler-multinomial exits from EA class
  (*reulmult)(2,INFA,&rate[15],dt,&trans[15]); // euler-multinomial exits from INFA class
  (*reulmult)(1,RECA,&rate[17],dt,&trans[17]); // euler-multinomial exits from RECA class

  // balance the equations
  SUSJ += trans[0]-trans[1]-trans[2]-trans[3]; // IN births; OUT juv mort, aging, inf
  EXPJ += trans[4]-trans[1]-trans[13]; // IN inf from both SUSJ classes exposed; OUT incubation and juv mort
  INFJ += trans[13]-trans[17]; // IN incubation; OUT dis induced mortality
  RECJ += trans[15]-trans[18]-trans[19]; // IN seroconversion; OUT aging and juv mortality
  SUSA += trans[3]+trans[9]-trans[20]-trans[21]-trans[22]; // IN aging from both SUSJ classes; OUT adult mort, inf * 2 classes
  EIA += trans[21]-trans[23]-trans[24]; // IN from adult SUS exposed; OUT incubation and adult mortality
  INFA += trans[24]-trans[27]; // IN from incubation; OUT dis induced mortality
  RECA += trans[26]+trans[18]-trans[28]; // IN from serconverstion & aging; OUT from death
  PREVA = INFA/(SUSA+INFA+EXPA+RECA); // adult prevalence
  PREVJ = INFJ/(SUSJ+EXPJ+INFJ+RECJ); // juvenile prevalence
  SPA = RECA/(SUSA+EXPA+INFA+RECA); // adult seroprevalence
  SPJ = RECJ/(SUSJ+EXPJ+INFJ+RECJ); // juvenile seroprevalence
}
