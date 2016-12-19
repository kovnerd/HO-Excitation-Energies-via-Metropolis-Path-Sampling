#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <stdlib.h>
#include <math.h>
#include "matrixmem.c"
#define N 20
#define LATTICE_SPACING 0.5 // T/N
#define MASS 1.

	
void metropolis(double *, double, gsl_rng*);
double G(double*, int);
double S(double*, int);
double potential(double);
void bootstrap_G(double**, double**, int, gsl_rng*);
void bin_G(double**, double**, int, int);
void accumulateStatistics(double*, double*, double*, double*, gsl_rng*);
void Gs_average(double**, double*, int);

void metropolis(double *x, double epsilon, gsl_rng *r)
{
	for(int i = 0; i < N; i++)
	{
		double oldX = x[i];
		double oldS = S(x, i);
		x[i] += (2 * gsl_rng_uniform(r) - 1.)*epsilon;// generates random number between -epsilon and epsilon
		double deltaS = S(x, i) - oldS;
		if(deltaS > 0 && exp(-deltaS) < gsl_rng_uniform(r))
			x[i] = oldX;
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
	
	return x_prod/N;
}

double potential(double x)
{
	return MASS*x*x/2;
}

double S(double* xs, int i)
{
	int iNext = (i+1) % N;
	int iPrev = (i-1+N) % N;
	return LATTICE_SPACING*potential(xs[i]) + xs[i]*(xs[i]-xs[iNext]-xs[iPrev])/LATTICE_SPACING;
}

void bootstrap_G(double** Gs, double** bootstrap, int n_configs, gsl_rng *r)
{
	int rand = 0;
	for(int alpha = 0; alpha < n_configs; alpha++)
	{
		rand = (int) (gsl_rng_uniform(r)*(n_configs));//choose random config from Gs
		for(int n = 0; n < N; n++)//append random config to bootstrap
			bootstrap[alpha][n] = Gs[rand][n];
	}
}

void Gs_average(double** Gs, double* avG, int n_configs)
{
	for(int n = 0; n < N; n++)//average
	{
		avG[n] = 0;
		for(int alpha = 0; alpha < n_configs; alpha++)
			avG[n] += Gs[alpha][n];
		avG[n] /= n_configs;
	}
}

void bin_G(double** Gs, double** bin, int n_configs, int bin_size)
{
	double avg;
	int i = 0;
	for(int alpha = 0; alpha < n_configs; alpha += bin_size)
	{
		avg = 0.;
		for(int n = 0; n < N; n++)
		{
			for(int b = 0; b < bin_size; b++)
				avg += Gs[alpha + b][n];
			bin[i][n] = avg/bin_size;
		}
		i++;
	}
}


void accumulateStatistics(double *xs, double *energy, double *avG, double *energy_error, gsl_rng *r)
{	
	int N_configs = 1000; //N_cf
	int N_cor = 20; 
	//int bin_size = 20;
	int N_error = 100;
	
	int therm_steps = 5*N_cor;
	double epsilon = 1.4;
	
	double **Gs = matrix_allocate_double(N_configs, N);
	double **Gs_bootstrap = matrix_allocate_double(N_configs, N);
	double *G_temp = malloc((unsigned long)N*sizeof(double));
	//double **Gs_bins = matrix_allocate_double(N_configs/bin_size, N);
	double **energy_bootstrap = matrix_allocate_double(N, N_error);
	double *temp = malloc((unsigned long)N_error*sizeof(double));
		//initialize pointers on heap memory
		for(int i = 0; i < N; i++)
		{
			xs[i] = 0;
			G_temp[i] = 0.;
		}
	for(int alpha = 0; alpha < N_configs; alpha++)
		for(int n = 0; n < N; n++)
		{
			Gs[alpha][n] = 0.;
			Gs_bootstrap[alpha][n] = 0.;
		}
	/*for(int b = 0; b < N_configs/bin_size; b++)
		for(int n = 0; n < N; n++)
			Gs_bins[b][n] = 0.;*/
	for(int n = 0; n < N; n++)
		for(int n_err = 0; n_err < N_error; n_err++)
			energy_bootstrap[n][n_err] = 0.;
			
			
	for(int t = 0; t < therm_steps; t++)//thermalize the path
		metropolis(xs, epsilon, r);
	for(int alpha = 0; alpha < N_configs; alpha++)
	{
		for(int k = 0; k < N_cor; k++)//update paths
			metropolis(xs, epsilon, r);
		for(int n = 0; n < N; n++)	
			Gs[alpha][n] = G(xs, n);
	}
	Gs_average(Gs, avG, N_configs);//compute mc average of propagator G_n
	
	
	for(int n_err = 0; n_err < N_error; n_err++)//energy error analysis
	{
		bootstrap_G(Gs, Gs_bootstrap, N_configs, r);//generate bootstrap copies of Gs
		//bin_G(Gs_bootstrap, Gs_bins, N_configs, bin_size);//bin Gs bootstrap copies
		Gs_average(Gs_bootstrap, G_temp, N_configs);//get propagator values from binned bootstrap copy of Gs
		for(int n = 0; n < N; n++)
			energy_bootstrap[n][n_err] = log(fabs(G_temp[n]/G_temp[(n + 1)%N]));//generate sample of energies
	}
	
	for(int n = 0; n < N; n++)
	{
		
		for(int n_err  = 0; n_err < N_error; n_err++)
			temp[n_err] = energy_bootstrap[n][n_err];
		energy[n] = log(fabs(avG[n]/avG[(n + 1)%N]))/LATTICE_SPACING;
		energy_error[n] = gsl_stats_sd(temp, 1, (size_t) N_error);
		printf("%f %f %f %f\n", n * LATTICE_SPACING, avG[n], energy[n], energy_error[n]);
	}
	
	
	matrix_free_double(Gs);
	matrix_free_double(Gs_bootstrap);
	//matrix_free_double(Gs_bins);	
	free(temp);
	free(G_temp);
	matrix_free_double(energy_bootstrap);
}
