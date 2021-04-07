#include <math.h>
#include <float.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>

#include "fast_sin_cos.h"
#include "inverted_pendulum.hpp"
#include "read_from_file.hpp"
#include "aux_functions.hpp"

//#define test_flag
//#define COST_SQUARED

InvertedPendulum::InvertedPendulum(std::string config_file) {
	name = "Inverted Pendulum";

    load_configurations_from_file(config_file);

	// Setup Output Plot
	plot_number_of_states = 2;
    plot_number_of_controls = 1;
    plot_states_config = (int *) calloc(plot_number_of_states*Nx, sizeof(int));
    plot_control_config = (int *) calloc(plot_number_of_controls*n_U, sizeof(int));

    // State Matrix
    // r 		 1 0 0 0
    // theta 	 0 0 1 0
    plot_states_config[0] = 1;
    plot_states_config[6] = 1;

    // Control Matrix
    // u 		 1 
    plot_control_config[0] = 1;

    plot_labels = new std::string[3];
    plot_labels[0] = "r (m)";
    plot_labels[1] = "theta (rad)";
    plot_labels[2] = "u (N)";

    state_type = (int *) calloc(Nx, sizeof(int));
    state_type[0] = 0; // Position
    state_type[1] = 0; // Linear Velocity
    state_type[2] = 1; // Angle 1
    state_type[3] = 0; // Angular Velocity 1
}

void InvertedPendulum::control_from_parameters(_real * parameters, _real * previous_control, _real * control, int horizon) {

    _real t1 = exp(-Lambda*Ts)-1;
    _real t2 = exp(-q_param*Lambda*Ts)-1;
    _real t1_2 = 1/(t1-t2);

    _real d1_p1 = du_max[0]*parameters[0];
    
    _real t1_uss = t1*uss;
    _real t2_uss = t2*uss;

    _real alpha_1 =  (d1_p1 - t2*previous_control[0] + t2_uss)*t1_2;
    _real alpha_2 = -(d1_p1 - t1*previous_control[0] + t1_uss)*t1_2;
    
    _real pre_exponent;
    _real exp1;
    _real exp2;
    
    if (horizon == 1) {
        pre_exponent = 0;
        exp1 = 1;           // Transformar em TABELA
        exp2 = 1;   // Transformar em TABELA
        control[0] = std::min(std::max(uss + alpha_1 + alpha_2, u_min[0]), u_max[0]);
    }
    else{
        for (int i = 0; i < horizon; ++i) {
            pre_exponent = -Lambda*Ts*i;
            exp1 = exp(pre_exponent);           // Transformar em TABELA
            exp2 = exp(pre_exponent*q_param);   // Transformar em TABELA
            control[0*horizon+i] = std::min(std::max(uss + alpha_1 * exp1 + alpha_2 * exp2, u_min[0]), u_max[0]);
        }
    }
}

void InvertedPendulum::model(_real * state_dot, _real * state, _real * control){
	// Inverted Pendulum from:
	// Mercieca, J., & Fabri, S. G. (2011). 
	// Particle swarm optimisation for non-linear model predictive control. 
	// International Conference on Advanced Engineering Computing and Applications in Science, 5(1), 88?93.

	//_real m = 7.3;    // [Kg] Uniformily distributed mass of the pendulum
	//_real M = 14.6;   // [Kg] Cart mass
	_real g = 9.81;   // [m/s^2] Gravity
	//_real l = 1.2;    // [m] half the lenght of the pendulum
	_real b = 14.6;   // [Kg/s] surface friction
	//_real h = 0.0136;  // [Kg.m^2/s] rotation friction damping coeficient

#ifdef FAST_SINCOS
    _real c = fastcos((float) state[2]);
    _real s = fastsin((float) state[2]);
#else
    _real th = state[2];
    _real c = cos(th);
    _real s = sin(th);
#endif
    //_real c_1 = 1.0/c;

	
	_real ml = 8.76; //m*l;
	_real Mm = 21.9; //M+m;
    _real Mn_1 = 0.04566210046;
	//_real Mml43 = 35.04; //Mm*l*4/3;
	_real h_ml = 0.001552511416; // h/ml
    _real l43 = 1.6; // l = 1.2 * 4/3

    _real Mm_c = Mm/c;

	_real mls = ml*s;
	_real mlc = ml*c;
	_real u_bx1 = control[0]-b*state[1];
    _real state4_2 = state[3]*state[3]; // pow(state[3],2);

    // xdot(1) = x(2);                                                 % x dot
    // xdot(3) = x(4);                                                 % theta dot
    // xdot(4) = (u-b*x(2)+mls*x(4)^2-Mm/c*(g*s-h*x(4)/ml)) / (ml*c-(Mm)*4*l/(3*c)); % theta dot dot
    // xdot(2) = (1/(Mm))*(u-b*x(2)-ml*xdot(4)*c+mls*x(4)^2); % x dot dot


	state_dot[0] = state[1];                                                 // x dot
	state_dot[2] = state[3];                                                 // theta dot
	//state_dot[3] = (u_bx1+mls*state4_2-Mm*c_1*(g*s-(h_ml)*state[3])) / (mlc-Mml43*c_1); // theta dot dot
    state_dot[3] = (u_bx1+mls*state4_2-Mm_c*(g*s-(h_ml)*state[3])) / (mlc-Mm_c*l43); // theta dot dot
	state_dot[1] = (Mn_1)*(u_bx1-mlc*state_dot[3]+mls*state4_2); // x dot dot

#ifdef test_flag
    using namespace std;
    cout << "xdot[0] = " << state_dot[0] << endl;
	cout << "xdot[2] = " << state_dot[2] << endl;
	printf("xdot[3] = (%lf - %lf * %lf + %lf * %lf - %lf / %lf *(%lf * %lf - %lf * %lf / %lf)/(%lf * %lf - %lf *4.00* %lf /(3.00* %lf))\n", control,b,state[1],mls,state4_2,Mm,c,g,s,h,state[3],ml,ml,c,Mm,l,c); 
	printf("xdot[1] = (1/(%lf))*(%lf - %lf * %lf - %lf * %lf * %lf + %lf * %lf)\n", Mm, control[0],b,state[1],ml,state_dot[3],c,mls,state4_2);
#endif

}

_real InvertedPendulum::nmpc_cost_function(_real * control_guess, _real * xref, _real * uref, _real * xss, _real * uss){


    _real J = 0.0;
    _real * Ji ;
    _real * Jf; // = 0.0;
    _real Ju = 0.0;

    Ji = (_real *) calloc(Nx,sizeof(_real));
    Jf = (_real *) calloc(Nx,sizeof(_real));

    // Constructing the vector of guessed control actions with respect to N and Nu
    _real * uu;
    uu = (_real *) malloc(N*sizeof(_real));

    for (int i = 0; i < N; ++i) {
        if(i < Nu) {
            uu[i] = control_guess[i];
        }
        else {
            uu[i] = control_guess[Nu-1];    
        }
    }

    // Initialize x_hat vector
    _real * x_hat; //[Nx*N];
    x_hat = (_real *) malloc(Nx*N*sizeof(_real));

    for (int i = 0; i < Nx; ++i) {
        x_hat[i] = current_state[i];
    }
    for (int i = Nx; i < (Nx*N); ++i) {
        x_hat[i] = 0.0;
    }

    // Calculate cost function
    int k = 0;
    int l = 0;
    J = 0;
    //Jf = 0;
    Ju = 0;
    for (int i = 0; i < N-1; ++i) {
        k = Nx*(i+1);
        one_step_prediction(&x_hat[k], &x_hat[Nx*i], &uu[i]);
        for (int j = 0; j < Nx; ++j) {
            l = k + j;
            if(j == 2) { //} && (x_hat[l] > M_PI)) {
                //Ji[j] = Ji[j] + error_with_constrains(j,x_hat[l])*pow((M_PI*2 - x_hat[l]) - xref[l],2);
                Ji[j] = Ji[j] + error_with_constrains(j,x_hat[l])*pow(normalize_angle(x_hat[l] - xref[l]),2);
            }
            else {
                Ji[j] = Ji[j] + error_with_constrains(j,x_hat[l])*pow(x_hat[l] - xref[l],2);
            }
        }
    }
    for (int j = 0; j < Nx; j++) {
#ifdef COST_SQUARED        
        J = J + Q[j]*sqrt(Ji[j]);
#else
        J = J + Q[j]*(Ji[j]);
#endif
    }


    k = Nx*(N-1);
    for (int j = 0; j < Nx; ++j) {
        l = k + j;
        if(j == 2) {//&& (x_hat[l] > M_PI)) {
            //Jf[j] = Jf[j] + error_with_constrains(j,x_hat[l])*pow(( (M_PI*2 - x_hat[l]) - xss[j]),2);
            Jf[j] = Jf[j] + error_with_constrains(j,x_hat[l])*pow(normalize_angle(x_hat[l] - xss[j]),2);
        }
        else {
            Jf[j] = Jf[j] + error_with_constrains(j,x_hat[l])*pow(( x_hat[l] - xss[j]),2);
        }
    }
    
    for (int j = 0; j < Nx; j++) {
#ifdef COST_SQUARED    
    J = J + Qf[j]*sqrt(Jf[j]);
#else    
    J = J + Qf[j]*(Jf[j]);
#endif
    }   

    if(R[0] > 0) {
        for (int i = 1; i < N; ++i) {  
            Ju = Ju + pow(uu[i]-uu[i-1],2);
        }
    }
#ifdef COST_SQUARED                    
    J = J + R[0]*sqrt(Ju);
#else
    J = J + R[0]*Ju;
#endif

    if(x_hat != NULL) free(x_hat);
    if(uu != NULL) free(uu);
    free(Ji);
    free(Jf);

    return J;

}
