/* **************************************************************************
 * This code generates synthetic linear regression dataset.
 * The observation/feature matrix is generated from a gaussian distribution.
 *
 * Code by Avik Ray (avik@utexas.edu)
 *
 **************************************************************************** */

#include <cstdio>
#include <cstdlib>
#include "sgd_lib.h"

void generate_normal(config_t& config, char *filename, bool noise_on, int normalize){
	printf("\nGenerating samples from normal distribution.");
	gsamples_t gs;
	gs.generate(config,normalize,noise_on,filename);
	return;
}

int main(int argc, char *argv[]){

	printf("Generating synthetic data ...");

	// Initilization
	long num_variables = atoi(argv[1]); // Arg 1: Number of variables
	long num_obs = atoi(argv[2]); // Arg 2: Number of observations
	bool noise_on;
	if(atoi(argv[4])==1){
		noise_on = true;  // Arg 3: 1 for noise, 0 for no noise
	}
	else{
		noise_on = false;
	}

	int normalize = atoi(argv[4]); // Arg 4: 1 to normalize columns

	char filename[100];
	strcpy(filename,argv[5]); // Arg 4: filename prefix

	// Echo parameters
	printf("\nNumber of variables = %ld",num_variables);
	printf("\nNumber of observations = %ld",num_obs);
	if(noise_on){
		printf("\nNoise = ON");
	}
	else{
		printf("\nNoise = OFF");
	}
	printf("\nNormalize = %d",normalize);
	printf("\nFilename = %s",filename);
	
	// Generate Samples
	config_t config;
	config.init(num_variables,num_obs);
	
	generate_normal(config, filename, noise_on, normalize);

	printf("\nData generation complete!");
	printf("\n");

}

