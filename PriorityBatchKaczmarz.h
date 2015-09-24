/* **************************************************************************
 * Header file for the PriorityBatchKaczmarz class
 *
 * Code by Avik Ray (avik@utexas.edu)
 *
 **************************************************************************** */

#ifndef PRIORITYBATCHKACZMARZ_H
#define PRIORITYBATCHKACZMARZ_H

#include <cstdio>
#include <cstdlib>
#include <vector>       // std::vector
#include <cmath>
#include <queue>          // std::priority_queue
#include <functional>     // std::greater, std::less
#include <algorithm>    // std::random_shuffle
#include "sgd_lib.h"

class PriorityBatchKaczmarz_t{
public:
	double *x;
	double *error;
	long num_variables;
	double runtime;
	long steps_till_convergence;
	double pscaling;

	std::priority_queue <qnode_t, std::vector<qnode_t>, std::less<qnode_t> > pq;

	PriorityBatchKaczmarz_t(): x(NULL), error(NULL), num_variables(0), steps_till_convergence(0) {
		pscaling = std::pow(10,10);
	}
	~PriorityBatchKaczmarz_t(){
		if(x) delete [] x;
		if(error) delete [] error;
	}

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

	void save(config_t& config){
		// save x
		FILE *fp = fopen("pbk_x.o","w");
		for(long i=0; i<num_variables; i++) fprintf(fp,"%lf\n",x[i]);
		fclose(fp);

		// save error
		fp = fopen("pbk_e.o","w");
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

		// initialize queue (Phase 1)
		std::vector<long> row_seq_vec;
		for(long i=0; i<A.rows; i++) row_seq_vec.push_back(i);

		long total_samples = config.pbk_init_sample_count + config.num_steps*config.pbk_sample_count;

		if(total_samples>A.rows){
			for(long i=A.rows; i<total_samples; i++) row_seq_vec.push_back(i % A.rows);
		}
		std::random_shuffle(row_seq_vec.begin(), row_seq_vec.end());


		long row_idx;
		//long *data = row_seq_vec.data();
		qnode_t node;
		for(long i=0; i<config.pbk_init_sample_count; i++){
			//row_idx = data[i];
			row_idx = row_seq_vec.back();
			row_seq_vec.pop_back();
			double priority = 0;
			switch(config.priority_val){
				case 0:
					priority = pscaling*.5*std::pow(y[row_idx]-A.row_multiply(row_idx,x_curr),2);
					break;
				case 1:
					priority = pscaling*std::pow((y[i]-A.row_multiply(i,x))*A.row_norm(i),2);
					break;
			}

			node.index = row_idx;
			node.priority = priority;

			pq.push(node); // insert in queue

		}

		// iterate (Phase 2)
		double err;
		double Li;
		double step;
		for(long it=0; it<config.num_steps; it++){

			if((it+1) % config.display_step == 1){
				if(it>0){
					printf("\nSteps %ld to %ld, x_err = %lf",it+1,it+config.display_step,error[it-1]);
				}
				else{
					printf("\nSteps %ld to %ld, x_err = ?",it+1,it+config.display_step);
				}
			}

			node = pq.top();
			pq.pop();
			row_idx = node.index;

			// choose step size
			switch(config.step_type){
				case 0:
					Li = std::pow(A.row_norm(row_idx),2);
					step = 1/Li;
					break;
				case 1:
					step = h/(it+1);
				        break;
				default:
					step = h/(it+1);

			}

			err = y[row_idx] - A.row_multiply(row_idx,x_curr);
			vecAdd(x, num_variables, 1, x_curr, step*err, A.row_full(row_idx));

			// compute L2 norm error in x
			vecAdd(x_err, num_variables, 1, x, -1, config.x);
			error[it] = vecNormL2(x_err, num_variables);

			// compute L2 norm error in y
			A.multiply(temp,x);
			vecAdd(y_err, num_obs, 1, y, -1, temp);

			steps_till_convergence = it+1;
			if(vecNormL2(y_err,num_obs) < config.tolerance){
				printf("\nAlgorithm converged!");
				break;
			}
			else{
				for(long i=0; i<num_variables; i++) x_curr[i] = x[i];
			}

			if(it==config.num_steps-1){
				printf("\nMax step over ... quiting.");
				break;
			}


			// remove old data from queue
			while(!pq.empty()){
				pq.pop();
			}

			// add new data to queue
			for(long i=0; i<config.pbk_sample_count; i++){
				//row_idx = data[i];
				row_idx = row_seq_vec.back();
				row_seq_vec.pop_back();
				double priority = 0;
				switch(config.priority_val){
					case 0:
						priority = pscaling*.5*std::pow(y[row_idx]-A.row_multiply(row_idx,x_curr),2);
						break;
					case 1:
						priority = pscaling*std::pow((y[i]-A.row_multiply(i,x))*A.row_norm(i),2);
						break;
				}

				node.index = row_idx;
				node.priority = priority;

				pq.push(node); // insert in queue

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

};

#endif /* PriorityBatchKaczmarz.h */

