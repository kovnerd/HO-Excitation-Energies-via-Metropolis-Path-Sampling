#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include "matrixmem.c"
#define TIME 10
#define N 20
#define LATTICE_SPACING 0.5 // T/N
#define MASS 1.
#define N_E 8

void metropolis(double *, double, gsl_rng*);
double G(double*, int);
double S(double*, int);
double potential(double);
void accumulateStatistics(double*, double*, double*, double*, gsl_rng*);

void metropolis(double *x, double epsilon, gsl_rng *r)
{	//lines of code in comments with *** mark alternative code
	for(int i = 0; i < N; i++)
	{
		double oldX = x[i];
		double oldS = S(x, i);
		double step_length = epsilon * (2. * gsl_rng_uniform(r) - 1.);// generates random number between -epsilon and epsilon
		
		x[i] += step_length;
		double newX = x[i];
		
		double deltaS = S(x, i) - oldS;
		double prob = exp(-deltaS);
		x[i] = (prob > gsl_rng_uniform(r)) ? newX: oldX;
	}
}

double G(double *xs, int n)
{
	double x_prod = 0.;
	for(int i = 0; i < N; i++)
	{
		//x_prod += xs[(i + n)%N] * xs[i]; 
		x_prod += xs[(i + n)%N]*xs[(i + n)%N]*xs[(i + n)%N] * xs[i]*xs[i]*xs[i];
	}
	
	return x_prod/(double)N;
}

double S(double* xs, int i)
{
	int iNext = (i == N - 1) ? 0 : i + 1;
	int iPrev = (i == 0) ? N - 1 : i - 1;
	
	return LATTICE_SPACING * potential(xs[i]) + xs[i]*(xs[i]-xs[iNext]-xs[iPrev])/LATTICE_SPACING;
}

double potential(double x)
{
	return 0.5*MASS*x*x;
}

void accumulateStatistics(double *xs, double *energy, double *prop, double *energy_error, gsl_rng *r)
{	
	int ensembleSize = 1000;//N_cf
	int N_cor = 20;
	int therm_steps = 7*N_cor;
	double epsilon = 1.4;
	double **Gs = matrix_allocate_double(ensembleSize, N);
	double *prop2 = malloc(N * sizeof(double));
	double *prop_error = malloc(N * sizeof(double));
	
	for(int i = 0; i < N; i++)
	{
		prop2[i] = 0.;
		prop_error[i] = 0.;
	}
	
	for(int t = 0; t < therm_steps; t++)//thermalize the path
		metropolis(xs, epsilon, r);
	for(int alpha = 0; alpha < ensembleSize; alpha++)//update paths via metropolis
	{
		for(int k = 0; k < N_cor; k++)
			metropolis(xs, epsilon, r);
		for(int n = 0; n < N; n++)	
			Gs[alpha][n] = G(xs, n);
	}
		
	for(int n = 0; n < N; n++)//average
	{
		for(int alpha = 0; alpha < ensembleSize; alpha++)
		{
			prop[n] += Gs[alpha][n];
			prop2[n] += Gs[alpha][n] * Gs[alpha][n];
		}
		prop[n] /= ensembleSize;
		prop2[n] /= ensembleSize;
		prop_error[n] = sqrt((prop2[n] - prop[n]*prop[n])/ensembleSize);
	}
	
	for(int n_e = 0; n_e < N_E; n_e++)
	{
		energy[n_e] = log(prop[n_e]/prop[n_e + 1]);
		//the energy error bars are calculated from the canonical error propagation formula (im hoping its right)
		energy_error[n_e] = sqrt((prop_error[n_e]/prop[n_e])*(prop_error[n_e]/prop[n_e]) + (prop_error[n_e + 1]/prop[n_e + 1])*(prop_error[n_e + 1]/prop[n_e + 1]));
		printf("%f %f %f %f\n", n_e * LATTICE_SPACING, prop[n_e], energy[n_e]/LATTICE_SPACING, energy_error[n_e]/LATTICE_SPACING);
	}
	matrix_free_double(Gs);
	free(prop2);
	free(prop_error);
				
}
