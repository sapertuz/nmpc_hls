#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "system.hpp"
#include "config.hpp"

//#define DEBUG

void System::initializeVectors() {
}

void System::initializeParametrizationVectors(){
}

System::~System(){
}

void System::load_configurations_from_file(){
	
    uss = _uss; 
	filter_reference = _Rising_Time;
    parametrized = _Parametrization;

	for (_uchar i = 0; i < n_U; i++){
	    // Defining Input Constrains
		u_max[i] = _u_max[i];
		u_min[i] = _u_min[i];
		du_max[i]= _du_max[i];
	}

	for (_uchar i = 0; i < Nx; i++)	{
		// Definig State Constrains
		state_upper_limits[i] = _state_upper_limits[i];
		state_lower_limits[i] = _state_lower_limits[i];
		
		// Definig Cost Function Weight Matrices
		Q[i] = _Q[i];
		Qf[i] = _Qf[i];

		last_state[i] = 0.0;
		
		current_state[i] = _current_state[i];
    }
	
	for (_uchar i = 0; i < n_U; i++){
		R[i] = _R[i];
		
		// Initialize Acceleration Constrains
		ddu_max[i] = du_max[i];
	    acc_max[i] = ddu_max[i];
	    acc_min[i] = -ddu_max[i];
	}
	
#if parametrized > 0
	Lambda = _Lambda;
	q_param = _q_param;
	for (_uchar i = 0; i < parametrized; i++){
		pmax[i] = _pmax[i];
		pmin[i] = _pmin[i];
	}
#endif
#if filter_reference > 0
	for (_uchar i = 0; i < Nx; i++){
		rising_time[i] = _tr[i];
	}
#endif
}

void System::updateAcceleration(
	_real v_curr[],
	_real u_curr[],
	_real u_ant[]
){
	// Update ddu_max
    for (int i = 0; i < n_U; ++i){
        v_curr[i] = (u_curr[i] - u_ant[i]);
        acc_max[i] = v_curr[i] + ddu_max[i];
        if(acc_max[i] > du_max[i])
        	acc_max[i] = du_max[i];
        acc_min[i] = v_curr[i] - ddu_max[i];
		if(acc_min[i] < -du_max[i])
        	acc_min[i] = -du_max[i];
        //std::cout << "u_curr: " << u_curr[i] << " v_curr: " << v_curr[i] << " acc_max: " << acc_max[i] << " acc_min: " << acc_min[i] << std::endl;
    }
}


// 4th Order Runge-Kutta Method
void System::one_step_prediction(
	_real state_plus[], 
	_real state[], 
	_real control[]
){
	_real k1[Nx], k2[Nx], k3[Nx], k4[Nx];
	_real state_temp[Nx];

    //print_array("state_in", state, Nx, 0);
	model(k1,state,control);

    // print_array("k1", k1, Nx, 0);

	for (int i = 0; i < Nx; ++i) {
		state_temp[i] = state[i] + Ts*0.5*k1[i];
	}

	model(k2,state_temp,control);
    // print_array("k2", k2, Nx, 0);

	for (int i = 0; i < Nx; ++i) {
		state_temp[i] = state[i] + Ts*0.5*k2[i];
	}

	model(k3,state_temp,control);
    // print_array("k3", k3, Nx, 0);

	for (int i = 0; i < Nx; ++i) {
		state_temp[i] = state[i] + Ts*k3[i];
	}

	model(k4,state_temp,control);
    // print_array("k4", k4, Nx, 0);

	_real Ts_6 = Ts*0.1666666667; //Ts/6;

	for (int i = 0; i < Nx; ++i) {
		state_plus[i] = state[i] + Ts_6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
	}
	
}


_real System::error_with_constrains(
	int j,
	_real x
){
	const _real penality = 1e10;

	if(x > state_upper_limits[j]) {
		return (x-this->state_upper_limits[j])*penality;
	}
	else if (x < state_lower_limits[j]) {
		return (state_lower_limits[j]-x)*penality;
	}
	else {
		return 1.0;
	}
}

/*
_real System::nmpc_cost_function(
	_real control_guess[],
	_real xref[],
	_real uref[],
	_real xss[],
	_real uss[]
){
//	std::cout << "cost function system" << std::endl;

	_real J = 0.0;
	_real Jf = 0.0;

	// Constructing the vector of guessed control actions with respect to N and Nu
	for (_uchar i = 0; i < N; ++i) {
		if(i < Nu) {
			uu[i] = control_guess[i];
		}else {
			uu[i] = control_guess[Nu-1];	
		}
	}

	// Initialize x_hat vector
	_real x_hat[Nx*N]; //[Nx*N];
	
	for (_uchar i = 0; i < Nx; ++i) {
		x_hat[i] = current_state[i];
	}
	for (_uchar i = Nx; i < (Nx*N); ++i) {
		x_hat[i] = 0.0;
	}

	// Calculate cost function
	int k = 0;
	int l = 0;
    J = 0;
    Jf = 0;
	for (_uchar i = 0; i < N-1; ++i) {
		k = Nx*(i+1);
        one_step_prediction(&x_hat[k], &x_hat[Nx*i], &uu[i]);
        for (int j = 0; j < Nx; ++j) {
			l = k + j;
			if((j == 2) && (x_hat[l] > M_PI)) {
				J = J + pow(Q[j]*error_with_constrains(j,x_hat[l])*((M_PI*2 - x_hat[l]) - xref[l]),2);
			}
			else {
				J = J + pow(Q[j]*error_with_constrains(j,x_hat[l])*( x_hat[l] - xref[l]),2);
			}
		}
		J = J + R[0]*pow(uu[i],2);
	}
	k = Nx*(N-1);
	for (int j = 0; j < Nx; ++j) {
		l = k + j;
		if((j == 2) && (x_hat[l] > M_PI)) {
			Jf = Jf + Qf[j]*error_with_constrains(j,x_hat[l])*pow(( (M_PI*2 - x_hat[l]) - xss[j]),2);
		}
		else {
			Jf = Jf + Qf[j]*error_with_constrains(j,x_hat[l])*pow(( x_hat[l] - xss[j]),2);
		}
	}
	J = J + Jf + R[0]*pow(uu[N-1],2);

	return J;

}
*/

void System::getState(
	_real state[]
){
	for (_uchar i = 0; i < Nx; i++){
		state[i] = current_state[i];
	}
}
void System::setState(
	const _real state[]
){
	for (int i = 0; i < Nx; ++i){
		this->current_state[i] = state[i];
	}
}

int System::is_angle_state(
	int pos
){
    return state_type[pos];
}
