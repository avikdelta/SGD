/* **************************************************************************
 * Header file for the RandKaczmarz class
 *
 * Code by Avik Ray (avik@utexas.edu)
 *
 **************************************************************************** */

#ifndef RANDKACZMARZ_H
#define RANDKACZMARZ_H

#include <cstdio>
#include <cstdlib>
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime> 	// std::time
#include <cmath>
#include <time.h>
#include "sgd_lib.h"


class RandKaczmarz_t{
public:
	double *x;
	double *error;
	long num_variables;
	double runtime;
	long steps_till_convergence;

	RandKaczmarz_t(): x(NULL), error(NULL), num_variables(0), steps_till_convergence(0) {}

	void init(config_t &config){
		num_variables = config.num_variables;
		long num_steps = config.num_steps;

		if(x!=NULL) delete [] x;

		x = new double[num_variables];

		if(error!=NULL) delete [] error;

		error = new double[num_steps];

		return;
	}

	void clear(){
		num_variables = 0;
		if(x) delete [] x;
		x = NULL;
		if(error) delete [] error;
		error = NULL;
		return;
	}

	void save(config_t &config){
		// save x
		FILE *fp = fopen("rk_x.o","w");
		for(long i=0; i<num_variables; i++) fprintf(fp,"%lf\n",x[i]);
		fclose(fp);

		// save error
		fp = fopen("rk_e.o","w");
		for(long i=0; i<steps_till_convergence; i++) fprintf(fp,"%lf\n",error[i]);
		fclose(fp);

		return;
	}

	void Run(double y[], csrmat_t &A, config_t &config){

		clock_t begin, end;
		begin = clock();
		
		// initialize
		long num_obs = config.num_obs;
		double *x_curr = new double[num_variables];
		double *x_err = new double[num_variables];
		double *y_err = new double[num_obs];
		double *temp = new double[num_obs];
		memset(x,0,sizeof(double)*num_variables);
		memset(x_curr,0,sizeof(double)*num_variables);
		memset(error,0,sizeof(double)*config.num_steps);
		float h = config.h;

		// random permutation
		std::vector<long> row_seq_vec;

		for(long i=0; i<A.rows; i++) row_seq_vec.push_back(i);

		if(config.num_steps>A.rows){
			for(long i=A.rows; i<config.num_steps; i++) row_seq_vec.push_back(i % A.rows);
		}

		// using built-in random generator
		std::random_shuffle(row_seq_vec.begin(), row_seq_vec.end());

		// iterate
		long step_count = 0;
		long row_idx;
		double err;
		double Li;
		double step;

		for(std::vector<long>::iterator it=row_seq_vec.begin(); it!=row_seq_vec.end(); ++it)
		{	
			step_count += 1;
			if(step_count % config.display_step == 1){
				if(step_count>1){
					printf("\nSteps %ld to %ld, x_err = %lf",step_count,step_count+config.display_step-1,error[step_count-2]);
				}
				else{
					printf("\nSteps %ld to %ld, x_err = ?",step_count,step_count+config.display_step-1);
				}
			}

			row_idx = *it;

			// choose step size
			switch(config.step_type){
				case 0:
					Li = pow(A.row_norm(row_idx),2);
					step = 1/Li;
					break;
				case 1:
					step = h/step_count;
				        break;
				default:
					step = h/step_count;

			}

			err = y[row_idx] - A.row_multiply(row_idx,x_curr);
			vecAdd(x, num_variables, 1, x_curr, step*err, A.row_full(row_idx));

			// compute L2 norm error in x
			vecAdd(x_err, num_variables, 1, x, -1, config.x);
			error[step_count-1] = vecNormL2(x_err, num_variables);

			// compute L2 norm error in y
			A.multiply(temp,x);
			vecAdd(y_err, num_obs, 1, y, -1, temp);

			steps_till_convergence = step_count;
			if(vecNormL2(y_err, num_obs) < config.tolerance){
				printf("\nAlgorithm converged!");
				break;
			}
			else{
				for(long i=0; i<num_variables; i++) x_curr[i] = x[i];
			}


			if(step_count>=config.num_steps){
				printf("\nMax step over ... quiting.");
				break;
			}

		}

		end = clock();
		runtime = (double)(end - begin) / CLOCKS_PER_SEC;
		
		delete [] x_curr;
		delete [] x_err;
		delete [] y_err;
		delete [] temp;

		return;
	}

	~RandKaczmarz_t(){
		if(x) delete [] x;
		if(error) delete [] error;
	}


};



#endif /* end of RandKaczmarz.h */


