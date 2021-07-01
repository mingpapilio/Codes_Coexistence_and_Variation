/********************************************************************************************************************
 * Species coexistence model with performance curve response to environmental variability                           *
 * This file gets the outcomes of coexistence and other scenarios with various mean and sd of temperature           *
 * File instruction:                                                                                                *
 * 1. Put "gen_beta.h" and "gen_beta.c" in the folder containing this file                                          *
 * 2. Put the dsfmt folder and the folder containing this file into the same folder                                 *
 ********************************************************************************************************************
 * Execution example code:
gcc rpc2_vart.c gen_beta.h gen_beta.c -lm -stdlib=libstdc++
./a.out
 ********************************************************************************************************************
 * Key parameters                                                                                                   *
 * mean_tmp: Average temperature                                                                                    *
 * sd_tmp_long: Amplitude of long-term variation (SD of a normal distribution)                                      *
 * sd_tmp_short: Amplitude of short-term variation (SD of a normal distribution)                                    *
 * The above parameters are set by three arrays: mean_temperature, short_term, and long-term (Line 93-99)           *
 *      Please adjust the setting according to requirement:                                                         *
 *      e.g. Fix long-term variation: env_long_base= desired value, env_long_span= 0.0                              *
 *      e.g. Fix short-term variation: env_short_base= desired value, env_short_span= 0.0                           *
 *      Current setting generates Fig2g (short-term= 8.5, long-term is variable)                                    *
 *      Please adjust the printed variable (Line 301) if long-term variation is fixed                               *
 *          (i.e. change sd_tmp_long to sd_tmp_short)                                                               *
 ********************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>	
#include "gen_beta.h"

// Functions
    void message_error(char error_text[]);
    double *d_vector(long size); //vector creator
    double **d_matrix(long size_row, long size_column); //matrix creater
    void free_d_vector(double *x);
    void free_d_matrix(double **x);
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
// Index of variables
    int i_1= 1;                     // Index for species 1
    int i_2= 2;                     // Index for species 2

// Main function
int main (void)
{
    // Switches
        int s1= 1;				    // 1: Extinction at very low density (controlled by ext_thr)

    // Temporal variables
        int T=          20000; 	    // Duration of simulation
        double deltat=  0.005; 	    // Length of time step
        int T_short=    1;          // The length of short-term variation 
        double m=       70.0;       // The ratio of spans in sampling new long-term variation to short-term variation, notation is 'delta' in table 1
        int T_long=     T_short*m;  // The length of long-term variation

    // Local variables
		int i,j;				    // For loop counters
        int ii,jj;                  // For loop counters (creating various combinations of averages and variations of temperature)
		double t= 0.0;			    // Time logs
        double T_remain_long= 0.0;  // Counter of remaining season length (entire episode)
        double T_remain_short=0.0;  // Counter of remaining time within one short-term condition
		double pp= 0.0;			    // Temp for random number
        double ext_thr= 0.5; 		// Threshold of extinction
        
    // Temperaure distribution (environmental factors)
        double tmp_temp;                        // mean temperature of short-term temperature distribution
        double curr_tmp;                        // current tempreature (longer-term+ shorter-term)
        double past_tmp;                        // temperature of the last time-step, check for updating population growth parameters
        // Parameter space setting (mean and variation of temperature)
            int env_mean_num= 13;               // Number of mean temperature simualted
            int env_shape_num= 5;               // Number of environmental variations simulated for each mean temperature
            double env_mean_span= 3.0;          // Increment of mean temperature
            double env_short_span= 0.0;         // Increment of short-term variation (size of SD)
            double env_long_span= 5.0;          // Increment of long-term variation (size of SD)
            double env_mean_base= 2.0;          // Base value of mean temperature
            double env_short_base= 8.5;         // Base value of short-term variation
            double env_long_base= 0.5;          // Base value of long-term variation
        // Normal dist variables
            double mean_tmp;                    // Average temperature
            double sd_tmp_long;                 // Amplitude of long-term variation (SD of a normal distribution)
            double sd_tmp_short;                // Amplitude of short-term variation (SD of a normal distribution)
            double u1, u2;
        // Arrays for mean, short- and long-term variations
            double mean_temperature[env_mean_num];
            double short_term[env_shape_num];
            double long_term[env_shape_num];
            for (i=0; i< env_mean_num; i++) mean_temperature[i]= env_mean_base+ env_mean_span*i;
            for (i=0; i< env_shape_num; i++) short_term[i]= env_short_base+ env_short_span*i;
            for (i=0; i< env_shape_num; i++) long_term[i]= env_long_base+ env_long_span*i;
            long_term[0]= 0.5; long_term[1]= 10.5; long_term[2]= 10*sqrt(2)+0.5; long_term[3]= 10*sqrt(3)+0.5; long_term[4]= 20.5;
            short_term[0]= 20.5, short_term[1]= 10*sqrt(3)+0.5; short_term[2]= 10*sqrt(2)+0.5; short_term[3]= 10.5; short_term[4]= 0.5;

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

    // episode (arithmetic mean variables)
        double K_temp1, K_temp2, b_temp1, b_temp2, d_temp1, d_temp2;
        K_temp1= K_temp2= b_temp1= b_temp2= d_temp1= d_temp2= 0.0;
        double var_temp;
    // Checking the simulation
        int restart= 1;
        int redo_counter= 0;
        int redo_limit= 100;
    // Repetition log and paramters
        int num_sim= 1E2;
        int count_coe, count_spe, count_gen, count_ext, count_fal, count_sim;
        count_coe= count_spe= count_gen= count_ext= count_fal= count_sim= 0;
        // 0= null, 1= coexist, 2= specialist(1) fix, 3= generalist(2) fix, 4= both extinct, 5= fail after 100 redo
        int result= 0;
	// Creating temporal space
		double *p= d_vector(2);
		double *dfdt= d_vector(2);
    // Random number generator initialization
        double shape_tmp= 1.0;
        gen_beta_param beta_tmp;	           	    // create the random beta type variable [0:1]
        gen_beta_initialize(&beta_tmp, shape_tmp, shape_tmp); // a=b
	// Output
		FILE *out, *summary;
		out= fopen("out.txt","w");
        fprintf(out,"mean\tsd_long\tsd_short\tresult_id\ttime\n");
        summary= fopen("summary.txt","w");
        fprintf(summary,"mean\tsd_long\tsd_short\tN_sim\tN_count\tdom_1\tdom_2\tcoexist\textinct\tff\n");
        time_t time_start, time_end;
        time_start= time(NULL);
	// Initialization
		p[i_1]= 250;// initial population of species 1
		p[i_2]= 250;// initial population of species 2

    // Simulation starts
    for (ii=0;ii<env_mean_num;ii++){
    for (jj=0;jj<env_shape_num;jj++){
        //
        mean_tmp=       mean_temperature[ii];
        sd_tmp_long=    long_term[jj];
        sd_tmp_short=   short_term[jj];
        //
        count_coe= count_spe= count_gen= count_ext= count_fal= count_sim= 0;
    // Starting of repeated simulations with the same average and variability of temperature
    for (j=1;j<= num_sim; j++){
        // Reinitialization
            p[i_1]= 250;// initial population of x
            p[i_2]= 250;// initial population of y
            restart= 1;
            t= result= redo_counter= 0;
            T_remain_long= T_remain_short= 0.0;
        // Main loop of time series
        while(restart==1&& redo_counter< redo_limit){
            for (i=1; i<=T/deltat; i++){
                if(T_remain_long<=0){			// another season
                    T_remain_long= T_long;
                    // temperature
                        u1= gen_unif(&beta_tmp);
                        u2= gen_unif(&beta_tmp);
                    tmp_temp= normal_dist_BM (mean_tmp, sd_tmp_long, u1, u2);
                        u1= gen_unif(&beta_tmp);
                        u2= gen_unif(&beta_tmp);
                    curr_tmp= normal_dist_BM (tmp_temp, sd_tmp_short, u1, u2);
                }
                // executing short term variation
                if(T_remain_short<= 0){
                    // sampling temperature
                    u1= gen_unif(&beta_tmp);
                    u2= gen_unif(&beta_tmp);
                    curr_tmp= normal_dist_BM (tmp_temp, sd_tmp_short, u1, u2);
                    T_remain_short= T_short;
                }
                // Calculating Population growth parameters
                if (past_tmp!= curr_tmp){
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
                // Calculating the next population sizes with time+deltat
                    differential(t,p,dfdt);
                    rk4(p, dfdt, 2, t, deltat, p, differential);
                // Progression of time
                    past_tmp= curr_tmp;
                    T_remain_long-= deltat;
                    T_remain_short-= deltat;
                    t+= deltat;
                // extinction of population
                    if(s1==1){
                        if(p[i_1]<ext_thr) p[i_1]=0;
                        if(p[i_2]<ext_thr) p[i_2]=0;
                    }
                // warning message
                    if(p[i_1]<0||p[i_2]<0) {
                        printf("Dynamics crashed.\n");
                        restart=1;
                        redo_counter+= 1;                
                    }
                // Condition determination (fixation and double-extinction)
                    // Both extinct
                    if(p[i_1]==0 && p[i_2]==0){
                        count_ext+=1;
                        result= 4;
                        fprintf(out,"%lf\t%lf\t%lf\t%d\t%lf\n",mean_tmp,sd_tmp_long,sd_tmp_short,result,t);
                        restart= 0;
                        t= T;
                        i= T/deltat;
                        count_sim+= 1;// Count number of simulations
                    }
                    // Fixation of generalist
                    if(p[i_1]==0 && p[i_2]>0){
                        count_gen+=1;
                        result= 3;
                        fprintf(out,"%lf\t%lf\t%lf\t%d\t%lf\n",mean_tmp,sd_tmp_long,sd_tmp_short,result,t);
                        restart= 0;
                        t= T;
                        i= T/deltat;
                        count_sim+= 1;// Count number of simulations
                    }
                    // Fixation of specialist
                    if(p[i_1]>0 && p[i_2]==0){
                        count_spe+=1;
                        result= 2;
                        fprintf(out,"%lf\t%lf\t%lf\t%d\t%lf\n",mean_tmp,sd_tmp_long,sd_tmp_short,result,t);
                        restart= 0;
                        t= T;
                        i= T/deltat;
                        count_sim+= 1;// Count number of simulations
                    }
            }
            // Coexisting
            if (p[i_1]>0 && p[i_2]>0 && t/T>0.99){
                count_coe+=1;
                result= 1;
                fprintf(out,"%lf\t%lf\t%lf\t%d\t%lf\n",mean_tmp,sd_tmp_long,sd_tmp_short,result,t);
                count_sim+= 1;// Count number of simulations
            }
            if(t/T>0.99) restart= 0;
            else {
                printf("Simulation is not completed, t is %lf.\n",t);
                restart=1;
                redo_counter+= 1;
            }
        }
        // Failing
        if (redo_counter>= redo_limit) {
            count_fal+=1;
            result= 5;
            fprintf(out,"%lf\t%lf\t%lf\t%d\t%lf\n",mean_tmp,sd_tmp_long,sd_tmp_short,result,t);
            count_sim+= 1;// Count number of simulations
        }
    }
    // IMPORTANT- co-existing processing
        if(num_sim >= count_sim){
            count_coe+= (num_sim-count_sim);
            count_sim= num_sim;
        }
    // print result 
        fprintf(summary,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",mean_tmp,sd_tmp_long,sd_tmp_short,num_sim,count_sim,count_spe,count_gen,count_coe,count_ext,count_fal);
    // printf("One simulation is completed\n");
    }}
    time_end= time(NULL);
    printf("This simulation lasted %ld seconds.\n",time_end-time_start);
	free_d_vector(p);
	free_d_vector(dfdt);
	fclose(out);
    fclose(summary);
	return 0;
}

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
	if(x==NULL) message_error("Allocation failure in d_vector()");
	return x;
}
double **d_matrix(long size_row, long size_column)
{
	double **x;
	long i;
	long size_row_P= size_row+1;
	long size_column_P= size_column+1;

	x= (double **) malloc((size_t)(size_row_P*sizeof(double *)));               //first dimension
	if (x==NULL) message_error("Allocation failure in d_vector()");
	x[0]= (double *) malloc((size_t)(size_row_P*size_column_P*sizeof(double))); //second dimension
	if (x[0]==NULL) message_error("Allocation failure in d_vector()");
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
