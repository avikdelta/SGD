/* **************************************************************************
 * This code runs a priority SGD/Priority Kaczmarz algorithm on a linear
 * regression dataset.
 *
 * Code by Avik Ray (avik@utexas.edu)
 *
 **************************************************************************** */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "sgd_lib.h"
#include "PriorityKaczmarz.h"

int main(int argc, char *argv[]){

	char config_filename[100], filenameX[100], filenameY[100], filenameA[100];

	strcpy(config_filename,argv[1]);
	strcpy(filenameX,argv[2]);
	strcpy(filenameY,argv[3]);
	strcpy(filenameA,argv[4]);

	//long num_steps = atoi(argv[5]);

	// Initialization
	config_t config;
	config.load(config_filename);

	// load X
	config.load_x(filenameX);

	// load Y
	long num_obs = config.num_obs;
	double *Y = new double[num_obs];
	loadVec(Y,num_obs,filenameY);

	// load A
	csrmat_t A;
	long num_variables = config.num_variables;
	A.load(num_obs,num_variables,filenameA);

	// Run PriorityKaczmarz
	printf("\n----------------------------");
	printf("\nStarting PriorityKaczmarz ...");
	printf("\n----------------------------");

	//config.num_steps = num_steps;
	PriorityKaczmarz_t PK;
	PK.init(config);

	PK.Run(Y,A,config);

	printf("\n----------------------------");
	printf("\nRuntime = %lf s",PK.runtime);
	printf("\n----------------------------");

	// save results
	PK.save(config);
	
	delete [] Y;
	printf("\n");

}

