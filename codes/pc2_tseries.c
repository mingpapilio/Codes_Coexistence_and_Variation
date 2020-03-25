/********************************************************************************************************************
 * Species coexistence model with performance curve response to environmental variability                           *
 * This file generates time series of population dynamics                                                           *
 * File instruction:                                                                                                *
 * 1. Put "gen_beta.h" and "gen_beta.c" in the folder containing this file                                          *
 * 2. Put the dsfmt folder and the folder containing this file into the same folder                                 *
 ********************************************************************************************************************
 * Execution:
gcc pc2_tseries.c gen_beta.h gen_beta.c -lm -stdlib=libstdc++
./a.out
gnuplot
plot 'out.txt' using 3:1 title 'high temp' with lines lc rgb 'orange',\
'out.txt' using 3:2 title 'low temp' with lines lc rgb 'skyblue'
 ********************************************************************************************************************
 * Key parameters                                                                                                   *
 * mean_tmp: Average temperature                                                                                    *
 * sd_tmp_long: Amplitude of long-term variation (SD of a normal distribution)                                      *
 * sd_tmp_short: Amplitude of short-term variation (SD of a normal distribution)                                    *
 * Other crucial parameters include scaling and shaping parameters of the performance curve, see line 72-93         *
 ********************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>	
#include "gen_beta.h"

// Functions
    void message_error(char error_text[]);              // printing errors
    double *d_vector(long size);                        // creating a vector
    double **d_matrix(long size_row, long size_column); // creating a matrix
    void free_d_vector(double *x);                      // erasing a vector
    void free_d_matrix(double **x);                     // erasing a matrix
    void rk4(double p[], double k1[],int n, double t, double h, double pout[],void(*diff)(double,double[],double[]));
    void differential(double time,double in[],double out[]);
    double dx_dt(double time, double vr[]);             // Population growth equation of species 1, vr is a vector of variable
    double dy_dt(double time, double vr[]);             // Population growth equation of species 2
    double normal_dist_BM (double mean, double sd, double u1, double u2);
    double factorial (double x); 

// Global variables (parameters for the Lotka-Volterra equations)
    double a_12= 0.4;               // Intensity of interspecific interactions (species 2 to 1)
    double a_21= 0.4;               // Intensity of interspecific interactions (species 1 to 2)
    double b_1, b_2, K_1, K_2, d_1, d_2;
// Index of variabels
    int i_1= 1;                     // Index for species 1
    int i_2= 2;                     // Index for species 2

// Main function
int main (void)
{
    // Switches
        int s1= 1;                  // 1: Extinction at very low density (controlled by ext_thr)

    // Temporal variables
        int T=          3000; 	    // Duration of simulation
        double deltat=  0.005; 	    // Length of time step
        int T_short=    1;          // The length of short-term variation 
        double m=       70.0;       // The ratio of spans in sampling new long-term variation to short-term variation, notation is 'delta' in table 1
        int T_long=     T_short*m;  // The length of long-term variation

    // Basic variables
        int i,j,k;                  // For loop counters
		double t=              0.0; // Time logs
        double T_remain_long=  0.0; // Counter of remaining time to resample long-term variations
        double T_remain_short= 0.0; // Counter of remaining time to resample short-term variations
		double pp=             0.0;	// Temp space for random number
        double ext_thr= 0.5; 		// Threshold of extinction (effective when s1==1)

    // temperaure distribution (environmental factors)
        double tmp_temp;                    // mean temperature of short-term temperature distribution
        double curr_tmp;                    // current tempreature (longer-term+ shorter-term)
        double past_temp;                   // temperature of the last time-step, check for updating population growth parameters
        // Normal distribution
            double mean_tmp= 	        20.0;   // Average temperature
            double sd_tmp_long=	 	    0.5;    // Amplitude of long-term variation (SD of a normal distribution)
            double sd_tmp_short=        12.5;   // Amplitude of short-term variation (SD of a normal distribution)
            double u1, u2;

    // Thermal performance curve parameters
        double scale_K= 1E4;                    // Scaling factor of carrying capacities
        double scale_b= 0.5;                    // Scaling factor of intrinsic growth rate
        double scale_d= 0.01;                   // Scaling factor of density-independent mortality rate
        // asymmetric performance curve
            double env_opt_1=   30;             // Optimal temperature of species 1
            double env_opt_2=   17;             // Optimal temperature of species 2
            double env_max_1=   env_opt_1+ 5;   // Upper boundary of species 1's performance curve
            double env_max_2=   env_opt_2+ 6;   // Upper boundary of species 2's performance curve
            double env_sigma1=  5;              // Width parameter of species 1
            double env_sigma2=  2;              // Width parameter of species 2
            double per_shift=   0.001;          // Minimum value of fitness function, notation is 'w_base' in table 1

	// Creating temporal space
		double *p= d_vector(2);
		double *dfdt= d_vector(2);
	// Output
		FILE *out;
		out= fopen("out.txt","w");
	// Initialization
		p[i_1]= 250;// initial population of species 1
		p[i_2]= 250;// initial population of species 2
        fprintf(out,"N1\tN2\ttime\n");
        fprintf(out,"%lf\t%lf\t%lf\n",p[i_1],p[i_2],t);

    // Random number generator initialization
        double shape_tmp= 1.0;
        gen_beta_param beta_tmp;	           	    // create the random beta type variable [0:1]
        gen_beta_initialize(&beta_tmp, shape_tmp, shape_tmp); // a=b
    // Main loop for time series
    for (i=1; i<= T/deltat; i++){
        // Environment stochasticity
        if (T_remain_long<=0){
            T_remain_long= T_long;
            // Get the environment
                    u1= gen_unif(&beta_tmp);
                    u2= gen_unif(&beta_tmp);
                tmp_temp= normal_dist_BM (mean_tmp, sd_tmp_long, u1, u2);
                    u1= gen_unif(&beta_tmp);
                    u2= gen_unif(&beta_tmp);
                curr_tmp= normal_dist_BM (tmp_temp, sd_tmp_short, u1, u2);
        }
        if (T_remain_short<=0){
            T_remain_short= T_short;
            // Get the environment
                u1= gen_unif(&beta_tmp);
                u2= gen_unif(&beta_tmp);
            curr_tmp= normal_dist_BM (tmp_temp, sd_tmp_short, u1, u2);
        }
        // Calculating Population growth parameters
        if(past_tmp!= curr_tmp){
            if(curr_tmp<=env_opt_1){
                K_1= scale_K*exp(-1*((curr_tmp-env_opt_1)/2/env_sigma1)*((curr_tmp-env_opt_1)/2/env_sigma1))+scale_K*per_shift;
                b_1= scale_b*exp(-1*((curr_tmp-env_opt_1)/2/env_sigma1)*((curr_tmp-env_opt_1)/2/env_sigma1))+scale_b*per_shift;
                d_1= scale_d*(1-exp(-1*((curr_tmp-env_opt_1)/2/env_sigma1)*((curr_tmp-env_opt_1)/2/env_sigma1)));
            }
            else{
                if(curr_tmp<=env_max_1){
                    K_1= scale_K*(1-((curr_tmp-env_opt_1)/(env_opt_1-env_max_1))*((curr_tmp-env_opt_1)/(env_opt_1-env_max_1)))+scale_K*per_shift;
                    b_1= scale_b*(1-((curr_tmp-env_opt_1)/(env_opt_1-env_max_1))*((curr_tmp-env_opt_1)/(env_opt_1-env_max_1)))+scale_b*per_shift;
                    d_1= scale_d*((curr_tmp-env_opt_1)/(env_opt_1-env_max_1))*((curr_tmp-env_opt_1)/(env_opt_1-env_max_1));
                }
                else {
                    K_1= scale_K*per_shift;
                    b_1= scale_b*per_shift;
                    d_1= scale_d;
                }
            }
            if(curr_tmp<=env_opt_2){
                K_2= scale_K*exp(-1*((curr_tmp-env_opt_2)/2/env_sigma2)*((curr_tmp-env_opt_2)/2/env_sigma2))+scale_K*per_shift;
                b_2= scale_b*exp(-1*((curr_tmp-env_opt_2)/2/env_sigma2)*((curr_tmp-env_opt_2)/2/env_sigma2))+scale_b*per_shift;
                d_2= scale_d*(1-exp(-1*((curr_tmp-env_opt_2)/2/env_sigma2)*((curr_tmp-env_opt_2)/2/env_sigma2)));
            }
            else{
                if(curr_tmp<=env_max_2){
                    K_2= scale_K*(1-((curr_tmp-env_opt_2)/(env_opt_2-env_max_2))*((curr_tmp-env_opt_2)/(env_opt_2-env_max_2)))+scale_K*per_shift;
                    b_2= scale_b*(1-((curr_tmp-env_opt_2)/(env_opt_2-env_max_2))*((curr_tmp-env_opt_2)/(env_opt_2-env_max_2)))+scale_b*per_shift;
                    d_2= scale_d*((curr_tmp-env_opt_2)/(env_opt_2-env_max_2))*((curr_tmp-env_opt_2)/(env_opt_2-env_max_2));
                }
                else {
                    K_2= scale_K*per_shift;
                    b_2= scale_b*per_shift;
                    d_2= scale_d;
                }
            }
        }
        //  Calculating the next population sizes with time+deltat
            differential(t,p,dfdt);
            rk4(p,dfdt,2,t,deltat,p,differential);
        // Progression of time
            past_tmp= curr_tmp;
            T_remain_long-= deltat;
            T_remain_short-= deltat;
            t+= deltat;
        // Print
            if(i%20==0) fprintf(out,"%lf\t%lf\t%lf\t%lf\n",p[i_1],p[i_2],t, curr_tmp);
        // extinction of population
            if(s1==1){
                if(p[i_1]<ext_thr) p[i_1]=0;
                if(p[i_2]<ext_thr) p[i_2]=0;
            }
        // warning message
            if(p[i_1]<0||p[i_2]<0) {
                printf("Dynamics crashed.\n");             
            }
    }
	free_d_vector(p);
	free_d_vector(dfdt);
	fclose(out);
    return 0;    
}

// Functions
void message_error(char error_text[]) //standard error handler
{
	printf("There are some errors...\n");
	printf("%s\n",error_text);
	printf("...now existing to system...\n");
	exit(1);
}
double dx_dt(double time, double vr[])	// Species 1
{
    if(K_1>0 && b_1>0) return (b_1- b_1/K_1*vr[i_1]- b_1/K_1*a_21*vr[i_2]- d_1)*vr[i_1];
    else return -1*d_1*vr[i_1];
}
double dy_dt(double time, double vr[])	// Species 2
{
    if(K_2>0 && b_2>0) return (b_2- b_2/K_2*vr[i_2]- b_2/K_2*a_12*vr[i_1]- d_2)*vr[i_2];
    else return -1*d_2*vr[i_2];
}
void differential(double time, double in[], double out[])
{
	out[i_1]= dx_dt(time,in);
	out[i_2]= dy_dt(time,in);
}
void rk4(double p[], double k1[],int n, double t, double h, double pout[],void(*diff)(double,double[],double[]))
{
	int i;
	double tt,*k2,*k3,*k4,*pp;

	k2= d_vector(n);
	k3= d_vector(n);
	k4= d_vector(n);
	pp= d_vector(n);

	for (i=1;i<=n;i++) pp[i]= p[i]+ k1[i]*h/2;
	(*diff)(t+h/2,pp,k2);
	for (i=1;i<=n;i++) pp[i]= p[i]+ k2[i]*h/2;
	(*diff)(t+h/2,pp,k3);
	for (i=1;i<=n;i++) pp[i]= p[i]+ k3[i]*h;
	(*diff)(t+h,pp,k4);
	
	for(i=1;i<=n;i++) pout[i]= p[i]+ (k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*h/6;

	free_d_vector(k2);
	free_d_vector(k3);
	free_d_vector(k4);
	free_d_vector(pp);
}
double *d_vector(long size) 
{
	double *x;
	x= (double *) malloc((size_t)((size+1)*sizeof(double)));
	if(x==NULL) message_error((char *)"Allocation failure in d_vector()");
	return x;
}
double **d_matrix(long size_row, long size_column)
{
	double **x;
	long i;
	long size_row_P= size_row+1;
	long size_column_P= size_column+1;

	x= (double **) malloc((size_t)(size_row_P*sizeof(double *)));               //first dimension
	if (x==NULL) message_error((char *)"Allocation failure in d_vector()");
	x[0]= (double *) malloc((size_t)(size_row_P*size_column_P*sizeof(double))); //second dimension
	if (x[0]==NULL) message_error((char *)"Allocation failure in d_vector()");
	for(i=1;i<size_row_P;i++) x[i]= x[0]+ i*size_column_P;
	return x;
}
void free_d_matrix(double **x)
{
	free(x[0]);
	free(x);
}
void free_d_vector(double *x) {	free(x);}
double normal_dist_BM (double mean, double sd, double u1, double u2)
{
    /* 
     * Using Box-Muller method to generate pseudo-normal distributed numbers in [0,1]
     * Constructed in Feb, 2018
     */
    double z1;
	z1= sqrt(-2* log(u1))* cos(2* M_PI* u2);
	return z1*sd+ mean;
}
double factorial (double x)
{
	int i;
	double tmp, out;
	if(x< 0) return 0;	// not sure
	if(x<= 1) return 1;
	else{
		tmp= x;
		out= 1;
		for(i=1; i<= x; i++){
			out= out*tmp;
			tmp= tmp-1;
		}
		return out;
	}
}