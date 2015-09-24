/* **************************************************************************
 * This code generates a config file for SGD experiment.
 *
 * Code by Avik Ray (avik@utexas.edu)
 *
 **************************************************************************** */

#include <cstdio>
#include <cstdlib>
#include "sgd_lib.h"

void config_test(config_t &config, char *filename){
	config.set_default();
	config.num_steps = config.num_obs;
	config.save(filename);
	config.load(filename);
	config.print();
	return;
}

int main(int argc, char *argv[]){

	printf("Making config file ...");

	long num_variables = atoi(argv[1]);
	long num_obs = atoi(argv[2]);

	char config_filename[100];
	strcpy(config_filename,argv[3]);

	// Echo
	printf("\nNumber of variables = %ld",num_variables);
	printf("\nNumber of observations = %ld",num_obs);
	printf("\nConfig filename = %s",config_filename);

	config_t config;
	config.init(num_variables,num_obs);
	
	config_test(config, config_filename);

	printf("\nConfig file written!");

	printf("\n");
}
