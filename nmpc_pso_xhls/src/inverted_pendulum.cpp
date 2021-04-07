#include <math.h>
#include <float.h>
#include <cmath>

#include "inverted_pendulum.hpp"

//#define test_flag
//#define COST_SQUARED

InvertedPendulum::InvertedPendulum() {
    load_configurations_from_file();

    state_type[0] = 0; // Position
    state_type[1] = 0; // Linear Velocity
    state_type[2] = 1; // Angle 1
    state_type[3] = 0; // Angular Velocity 1
}

void InvertedPendulum::model(
    _real state_dot[],
    _real state[],
    _real control[]
){
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

    _real th = state[2];
    _real c = cos(th);
    _real s = sin(th);
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


}

_real InvertedPendulum::nmpc_cost_function(
    _real control_guess[],
    _real xref[],
    _real uref[],
    _real xss[],
    _real uss[]
){
    _real J = 0.0;
    _real Ji[_Nx] = {0,0,0,0};
    _real Jf[_Nx] = {0,0,0,0}; // = 0.0;
    _real Ju = 0.0;

    // Constructing the vector of guessed control actions with respect to N and Nu
    _real uu[_N];

    for (int i = 0; i < N; ++i) {
        if(i < Nu) {
            uu[i] = control_guess[i];
        }else{
            uu[i] = control_guess[Nu-1];    
        }
    }

    // Initialize x_hat vector
    _real x_hat[Nx*N]; //[Nx*N];
    
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
            if(j == 2) { 
                //_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
            	_real tmp_err = (x_hat[l] - xref[l]);
                Ji[j] = Ji[j] + error_with_constrains(j,x_hat[l])*tmp_err*tmp_err;
            }
            else {
                _real tmp_err = x_hat[l] - xref[l];
                Ji[j] = Ji[j] + error_with_constrains(j,x_hat[l])*tmp_err*tmp_err;
            }
        }
    }
    for (int j = 0; j < Nx; j++) {
        J = J + Q[j]*(Ji[j]);
    }


    k = Nx*(N-1);
    for (int j = 0; j < Nx; ++j) {
        l = k + j;
        if(j == 2) {
            //_real tmp_err = normalize_angle(x_hat[l] - xss[j]);
            _real tmp_err = (x_hat[l] - xss[j]);
            Jf[j] = Jf[j] + error_with_constrains(j,x_hat[l])*tmp_err*tmp_err;
        }
        else {
            _real tmp_err = ( x_hat[l] - xss[j]);
            Jf[j] = Jf[j] + error_with_constrains(j,x_hat[l])*tmp_err*tmp_err;
        }
    }
    
    for (int j = 0; j < Nx; j++) {
    J = J + Qf[j]*(Jf[j]);
    }   

    if(R[0] > 0) {
        for (int i = 1; i < N; ++i) {  
            _real tmp_err = uu[i]-uu[i-1];
            Ju = Ju + tmp_err*tmp_err;
        }
    }
    J = J + R[0]*Ju;

    return J;
}
