/*
  Implements Equation 5 ( \hat\phi_{AB} ) 

  Strobeck, C and K. Morgan (1978) The effect of intragenic recombination on the 
  number of alleles in a finite population.  Genetics 88: 829-844

  Except for this comment block, this code was written by KRT in 1998
  at the University of Chicago during a rotation with Richard Hudson.
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double calc_phiab (double theta, double r) {
	double phiab;
	double numerator,denominator;
	
	numerator = 2.0*pow(theta,3.0) + 0.5*pow(theta,2.0)*r + 11.0*pow(theta,2.0) + 3.0*theta*r
				+ 0.5*pow(r,2.0) + 18.0*theta + 6.5*r +9.0;
				
	denominator = (1.0 + theta) * ( 4.0*pow(theta,3.0) + 3.0*pow(theta,2.0)*r + 0.5*pow(r,2.0)*theta
				+ 20.0*pow(theta,2.0) + 9.5*theta*r  + 0.5*pow(r,2.0) + 27.0*theta + 6.5*r + 9.0 ) ;

	phiab = numerator/denominator;
	
	return phiab;
}

int main (int argc, char *argv[]) {

	double theta,nsites;
	double *r, *r_tot, *phiab;
	int narg=0,i;
	
	r = (double *)malloc((int)(argc-3)*sizeof(double));
	r_tot = (double *)malloc((int)(argc-3)*sizeof(double));
	phiab = (double *)malloc((int)(argc-3)*sizeof(double));
	if (argc == 1) {
		printf ("usage: ./rtot  theta nsites r1 ,r2, ... , rn \nexiting...\n");
		exit (0);
	}
	theta = atof (argv[++narg]);
	nsites = atof (argv[++narg]);
	
	for ( i = 0 ; i < (argc -3); ++i)
		r[i] = atof(argv[++narg]);
		
	for ( i = 0 ; i < (argc - 3); ++i) {
		phiab[i] = calc_phiab(theta/nsites,r[i]);
	}
	
	printf ("theta/site = %f, ", theta/nsites);
	printf ("\n    r        E(phiab)    \n");
	for (i = 0 ; i < (argc - 3) ; ++i)
		printf("%5f,%11f\n",r[i],phiab[i]);
		
	return (1);
}
