#include <stdio.h>
#include "metropolis_path_sampling.c"
#include <time.h>

int main(void)
{
	double xs[N];
	double propagator[N], energy[N_E], energy_error[N_E];
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus2);
	unsigned long seed = ((unsigned long) time(NULL));
	gsl_rng_set(r, seed);
	int n = 0;
	FILE *fp = fopen("harmonic_oscillator_metropolis_x3x3.txt", "w");
	
	for(int i = 0; i < N; i++)
	{
		xs[i] = 0;
		propagator[i] = 0;
		if(i < N_E)
			energy[i] = 0;
	}
	accumulateStatistics(xs, energy, propagator, energy_error, r);
	
	for(int i = 0; i < N_E; i++)
		fprintf(fp, "%f %f %f %f\n", (n++)*LATTICE_SPACING, propagator[i], energy[i]/LATTICE_SPACING, energy_error[i]/LATTICE_SPACING);
		
	gsl_rng_free(r);
	fclose(fp);
	return 0;
}


