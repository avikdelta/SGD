/* **************************************************************************
 * Library routines for SGD experiment.
 *
 * Code by Avik Ray (avik@utexas.edu)
 *
 **************************************************************************** */

#ifndef SGD_LIB_H
#define SGD_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <random>

/* Loads a vector from file */
void loadVec(double res[], long size, const char *filename){
	FILE *fp = fopen(filename,"r");
	std::memset(res,0,sizeof(double)*size);
	for(long i=0; i<size; i++){
		fscanf(fp,"%lf",&res[i]);
	}
	fclose(fp);
	return;
}

/* Count lines in file */
long countLine(const char *filename){
	long number_of_lines = 0;
    	FILE *fp = fopen(filename, "r");
    	int ch;

    	while(EOF != (ch=getc(fp))){
		if('\n' == ch) ++number_of_lines;
	}
	printf("\n%ld", number_of_lines);
	fclose(fp);
    	return number_of_lines;
}

/* Performs vector addition y = a1*x1 + a2*x2 */
void vecAdd(double y[], long dim, double a1, const double x1[], double a2, const double x2[]){
	
	for(long i=0; i<dim; i++){
		y[i] = a1*x1[i] + a2*x2[i];
	}
	return;
}

/* Computes L2 norm of a vector */
double vecNormL2(const double x[], long dim){

	double total = 0;
	for(long i=0; i<dim; i++){
		total += x[i]*x[i];
	}
	return(sqrt(total));
}
	
/* Random gaussian distribution */
double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }while (r == 0.0 || r > 1.0);

        double d = sqrt(-2.0*log(r)/r);
        double n1 = x*d;
        n2 = y*d;
        double result = n1*stddev + mean;
        n2_cached = 1;
        return result;        
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

/* config class */
class config_t{
public:
	long num_variables;   // number of variables
	long num_obs;         // number of observations
	long num_steps;       // number of iterative steps of SGD / Kaczmarz
	long display_step;    // number of steps between successive display output
	int step_type;        // 0 or 1, determines how value of step size is computed
	int priority_val;     // 0 or 1, determines how priority is computed
	float tolerance;      // tolerance used for testing convergence
	float noise_var;      // variance in noise
	float h;              // scaling constant for step size
	long psk_init_sample_count; // initial number of samples in queue for PriorityScheduleKaczmarz
	long psk_sample_count; // number of samples addd to priority queue in each iteration
        long pbk_init_sample_count; // initial number of samples in queue for PriorityBatchKaczmarz
	long pbk_sample_count; // number of samples addd to priority queue in each iteration	
	double *x;            // pointer to actual/oracle value of x

	config_t(): num_variables(0), num_obs(0), num_steps(0), display_step(10), step_type(0), priority_val(0), tolerance(.01), x(NULL), noise_var(.01), h(.1), psk_init_sample_count(1), psk_sample_count(1), pbk_init_sample_count(1), pbk_sample_count(1) {}

	void init(long _num_variables, long _num_obs){
		num_variables = _num_variables;
		num_obs = _num_obs;

		if(x==NULL) x = new double[num_variables];

		return;
	}

	void set_default(){
		step_type = 0;
		priority_val = 0;
		noise_var = .01;
		h = 1;
		psk_init_sample_count = 1;
		psk_sample_count = 1;
		pbk_init_sample_count = 1;
		pbk_sample_count = 1;
		return;
	}

	~config_t(){
		if(x) delete [] x;
	}

	void save(const char *filename){
		FILE *fp = fopen(filename,"w");
		fprintf(fp,"NUM_VARIABLES %ld\n",num_variables);
		fprintf(fp,"NUM_SAMPLES %ld\n",num_obs);
		fprintf(fp,"NUM_STEPS %ld\n",num_steps);
		fprintf(fp,"DISPLAY_STEP %ld\n",display_step);
		fprintf(fp,"STEP_TYPE %d\n",step_type);
		fprintf(fp,"PRIORITY %d\n",priority_val);
		fprintf(fp,"TOLERANCE %f\n",tolerance);
		fprintf(fp,"NOISE_VAR %f\n",noise_var);
		fprintf(fp,"H %f\n",h);
		fprintf(fp,"PSK_INIT_SAMPLE %ld\n",psk_init_sample_count);
		fprintf(fp,"PSK_ITER_SAMPLE %ld\n",psk_sample_count);
		fprintf(fp,"PBK_INIT_SAMPLE %ld\n",pbk_init_sample_count);
		fprintf(fp,"PBK_ITER_SAMPLE %ld\n",pbk_sample_count);
	       	fclose(fp);
		return;
	}

	void print(){
		printf("\nNUM_VARIABLES %ld",num_variables);
		printf("\nNUM_SAMPLES %ld",num_obs);
		printf("\nNUM_STEPS %ld",num_steps);
		printf("\nDISPLAY_STEP %ld",display_step);
		printf("\nSTEP_TYPE %d",step_type);
		printf("\nPRIORITY %d",priority_val);
		printf("\nTOLERANCE %f",tolerance);
		printf("\nNOISE_VAR %f",noise_var);
		printf("\nH %f",h);
		printf("\nPSK_INIT_SAMPLE %ld",psk_init_sample_count);
		printf("\nPSK_ITER_SAMPLE %ld",psk_sample_count);
		printf("\nPBK_INIT_SAMPLE %ld",pbk_init_sample_count);
		printf("\nPBK_ITER_SAMPLE %ld",pbk_sample_count);
		return;
	}

	void load(const char *filename){
		char buff[100];
		FILE *fp = fopen(filename,"r");
		fscanf(fp,"%s %ld",buff,&num_variables);
		fscanf(fp,"%s %ld",buff,&num_obs);
		fscanf(fp,"%s %ld",buff,&num_steps);
		fscanf(fp,"%s %ld",buff,&display_step);
		fscanf(fp,"%s %d",buff,&step_type);
		fscanf(fp,"%s %d",buff,&priority_val);
		fscanf(fp,"%s %f",buff,&tolerance);
		fscanf(fp,"%s %f",buff,&noise_var);
		fscanf(fp,"%s %f",buff,&h);
		fscanf(fp,"%s %ld",buff,&psk_init_sample_count);
		fscanf(fp,"%s %ld",buff,&psk_sample_count);
		fscanf(fp,"%s %ld",buff,&pbk_init_sample_count);
		fscanf(fp,"%s %ld",buff,&pbk_sample_count);
		fclose(fp);

		if(x==NULL) x = new double[num_variables];

		return;
	}

	void load_x(const char *filename){
		FILE *fp = fopen(filename,"r");
		for(long i=0; i<num_variables; i++){
			fscanf(fp,"%lf",&x[i]);
		}
		fclose(fp);
		printf("\nX loaded !");
		return;
	}

};

/* Gaussian sample class */
class gsamples_t{
public:
	double *y;
	double *x;
	double *A;

	gsamples_t(): x(NULL), y(NULL), A(NULL) {}

	void generate(config_t& config, int normalization_flag, bool noise_on, char *filename){

		FILE *fp;
		char filenameX[100], filenameY[100], filenameA[100];
		long num_variables = config.num_variables;
		long num_obs = config.num_obs;
		double noise_std = sqrt(config.noise_var);
		double x_std = 1.0;
		double x_mean = 0.0;
		double a_std = 1.0;
		double a_mean = 0.0;

		//std::default_random_engine generator;
  		//std::normal_distribution<double> xdist(x_mean,x_std);
		//std::normal_distribution<double> adist(a_mean,a_std);
		//std::normal_distribution<double> ndist(0,noise_std);
		
		

		// generate x 
		printf("\nGenerating x ...");
		strcpy(filenameX,filename);
		strcat(filenameX,"_X");
		x = new double[num_variables];
		for(long i=0; i<num_variables; i++){
			//x[i] = xdist(generator);
			x[i] = rand_normal(x_mean,x_std);
		}
		// normalize x and save
		fp = fopen(filenameX,"w");
		double xnorm = vecNormL2(x, num_variables);
		for(long i=0; i<num_variables; i++){
		       x[i] = x[i]/xnorm;
		       fprintf(fp,"%lf\n",x[i]);
		}
		fclose(fp);

		// generate samples
		printf("\nGenerating samples ...");
		y = new double[num_obs];
		A = new double[num_variables];
		strcpy(filenameY,filename);
		strcat(filenameY,"_Y");
		strcpy(filenameA,filename);
		strcat(filenameA,"_A");
		FILE *fpy = fopen(filenameY,"w");
		FILE *fpA = fopen(filenameA,"w");

		for(long sid=0; sid<num_obs; sid++){
			y[sid] = 0;
			// generate features
			for(long i=0; i<num_variables; i++){
				A[i] = 0;
				//A[i] = adist(generator);
				A[i] = rand_normal(a_mean,a_std);
				fprintf(fpA,"%ld %ld %lf\n",sid+1,i+1,A[i]);
			}
			// normalize
			switch(normalization_flag){
				case 1:
					double Anorm = vecNormL2(A, num_variables);
					for(long i=0; i<num_variables; i++){
						A[i] = A[i]/Anorm;
					}
					break;

			}
			// compute y
			for(long i=0; i<num_variables; i++) y[sid] += A[i]*x[i];

			// add noise
			//if(noise_on) y[sid] += ndist(generator);
			if(noise_on) y[sid] = rand_normal(0,noise_std);
			
			fprintf(fpy,"%lf\n",y[sid]);
		}

		fclose(fpy);
		fclose(fpA);
		printf("\nAll files saved !");

		delete [] x;
		x = NULL;
		delete [] y;
		y = NULL;
		delete [] A;
		A = NULL;

		
		return;
	}

	~gsamples_t() {
		if(x) delete [] x;
		if(y) delete [] y;
		if(A) delete [] A;
	}


};

/* CSR sparse matrix class */
class csrmat_t{
public:
	long rows, cols, nnz;
	long *csr_row_ptr;
	long *csr_col_idx;
	double *csr_val;
	double *buff;
	bool mem_alloc_by_me;

	csrmat_t(): csr_val(NULL), csr_row_ptr(NULL), csr_col_idx(NULL), buff(NULL), rows(0), cols(0), nnz(0) {}

	csrmat_t(csrmat_t& m){ *this = m; rows=0; cols=0; nnz=0; }

	void load(long _rows, long _cols, char *filename){


		// initialize
		rows = _rows;
		cols = _cols;
		nnz = countLine(filename);

		// allocate memory
		mem_alloc_by_me = true;
		csr_row_ptr = new long[rows+1];
		csr_col_idx = new long[nnz];
		csr_val = new double[nnz];

		std::memset(csr_row_ptr,0,sizeof(long)*(rows+1));

		// open data file for A matrix
		FILE *fp = fopen(filename,"r");


		// read data
		printf("\nNumber of variables = %ld",cols);
		printf("\nNumber of observations = %ld",rows);
		printf("\nReading data ...");

		for(long i=0,r,c; i<nnz; ++i){
			fscanf(fp,"%ld %ld %lf", &r, &c, &csr_val[i]);
			csr_row_ptr[r]++;
			csr_col_idx[i] = c-1;
		}
			
		fclose(fp);

		for(long r=1; r<=rows; ++r) csr_row_ptr[r] += csr_row_ptr[r-1];

		printf("\nFile Loaded !");
		return;
	}

	double* row_full(long row_index){

		if(buff==NULL) buff = new double[cols];

		std::memset(buff,0,sizeof(double)*cols);

		long start = csr_row_ptr[row_index];
		long end = csr_row_ptr[row_index+1]-1;
		
		for(long i=start; i<=end; i++){
			long col = csr_col_idx[i];
			buff[col] = csr_val[i];
		}

		return buff;
	}

	double row_multiply(long row_index, const double x[]){
		long start = csr_row_ptr[row_index];
		long end = csr_row_ptr[row_index+1]-1;
		double total = 0;
		for(long i=start; i<=end; i++)
		{
			long col = csr_col_idx[i];
			total += csr_val[i]*x[col];
		}

		return total;
	}

	void multiply(double res[], const double x[]){
		for(long i=0; i<rows; i++){
			long start = csr_row_ptr[i];
			long end = csr_row_ptr[i+1]-1;
			double total = 0;
			for(long j=start; j<=end; j++)
			{
				total += csr_val[j]*x[csr_col_idx[j]]; 
			}
			res[i] = total;
		}
		return;
	}

	double row_norm(long row_index){
		long start = csr_row_ptr[row_index];
		long end = csr_row_ptr[row_index+1]-1;
		double total = 0;
		for(long i=start; i<=end; i++)
		{
			total += pow(csr_val[i],2);
		}

		return(sqrt(total));
	}

	void print(){
		for(long i=0; i<rows; i++){
			printf("\n");
			long start = csr_row_ptr[i];
			long end = csr_row_ptr[i+1]-1;
			for(long j=start; j<=end; j++){
				printf("%lf ",csr_val[j]);
			}
		}
	}

	void clear_space() {
		if(csr_val) delete [] csr_val;
		csr_val = NULL;
		if(csr_row_ptr) delete [] csr_row_ptr;
		csr_row_ptr = NULL;
		if(csr_col_idx) delete [] csr_col_idx;
		csr_col_idx = NULL;
		if(buff) delete [] buff;
		buff = NULL;
		rows = 0; cols = 0; nnz=0;
		mem_alloc_by_me = false;
	}

	~csrmat_t(){
		if(mem_alloc_by_me) {
			if(csr_val) delete [] csr_val;
			if(csr_row_ptr) delete [] csr_row_ptr;
			if(csr_col_idx) delete [] csr_col_idx;
			if(buff) delete [] buff;
		}
	}

};

/* Queue element class */
class qnode_t{
public:
	long index;
	double priority;

	qnode_t(): index(0), priority(0) {}

	bool operator<(const qnode_t &q) const {
		return priority < q.priority;
	}

	bool operator>(const qnode_t &q) const {
		return priority > q.priority;
	}

	~qnode_t(){}
};



#endif /* sgd_lib.h */

