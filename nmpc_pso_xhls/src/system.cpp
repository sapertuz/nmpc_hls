#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "system.hpp"
#include "config.hpp"

#include "aux_functions.hpp"
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
	// _real k1[Nx], k2[Nx], k3[Nx], k4[Nx];
	// _real state_temp[Nx];

    // //print_array("state_in", state, Nx, 0);
	// model(k1,state,control);
	// for (int i = 0; i < Nx; ++i) {
	// 	state_temp[i] = state[i] + Ts*0.5*k1[i];
	// }

	// model(k2,state_temp,control);
	// for (int i = 0; i < Nx; ++i) {
	// 	state_temp[i] = state[i] + Ts*0.5*k2[i];
	// }

	// model(k3,state_temp,control);
	// for (int i = 0; i < Nx; ++i) {
	// 	state_temp[i] = state[i] + Ts*k3[i];
	// }

	// model(k4,state_temp,control);
	// _real Ts_6 = Ts*0.1666666667; //Ts/6;

	// for (int i = 0; i < Nx; ++i) {
	// 	state_plus[i] = state[i] + Ts_6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
	// }

	_real local_control[_n_U];
	_real local_state[_Nx];
	_real k1[_Nx], k2[_Nx], k3[_Nx], k4[_Nx];
	_real state_temp[_Nx];

	memcpy_loop_rolled(local_state, state, _Nx);
	memcpy_loop_rolled(local_control, control, _n_U);

stage_1:
//#pragma HLS DATAFLOW
        model(k1, local_state, local_control);
        update_state(state_temp, local_state, k1, Ts*.5);
//#pragma HLS DATAFLOW off

stage_2:
//#pragma HLS DATAFLOW
        model(k2,state_temp,local_control);
        update_state(state_temp, local_state, k2, Ts*.5);
//#pragma HLS DATAFLOW off

stage_3:
//#pragma HLS DATAFLOW
        model(k3,state_temp,local_control);
        update_state(state_temp, local_state, k3, Ts);
//#pragma HLS DATAFLOW off

stage_4:
//#pragma HLS DATAFLOW
        model(k4,state_temp,local_control);
        update_state_final(state_plus, local_state, k1, k2, k3, k4);
	
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

_real System::nmpc_cost_function(
	_real * control_guess, 
	_real * xref, 
	_real * uref, 
	_real * xss, 
	_real * uss){

	// New
	_real uu_1[_N*_n_U], uu_2[_N*_n_U];       
	_real x_horizon[_Nx*_N]; //[Nx*N];	
	_real Ji[_Nx]; // = 0.0;
	_real Jui[_n_U]; // = 0.0;
	_real final_x[_Nx], final_xref[_Nx];
	//memset(Jui, (const _hw_real)0.0, _n_U*sizeof(_hw_real));
	_real J_tmp1, J_tmp2;
	memset_loop(Ji, (_real)0.0, _Nx);
	memset_loop(Jui, (_real)0.0, _n_U);

	// Old
	_real J = 0.0;
	_real Jf = 0.0;
    _real Ju = 0.0;
    //_real temp;

	// // Constructing the vector of guessed control actions with respect to N and Nu
	// for (int i = 0; i < N; ++i) {
	// 	for (int j = 0; j < _n_U; ++j) {
	// 		if(i < Nu) {
	// 			uu_1[i*_n_U+j] = control_guess[j*Nu+i];
	// 		}
	// 		else {
	// 			uu_1[i*_n_U+j] = control_guess[j*Nu-1];	
	// 		}
	// 	}
	// }

	// Constructing the vector of guessed control actions with respect to N and Nu
	uu_loop(control_guess, uu_1, uu_2);
	
	// Control Error
	horizon_u_error(uu_1, Jui);
	// Predict Horizon
	horizon_step_prediction(current_state, uu_2, x_horizon);
	// State Error
	horizon_step_error(x_horizon, xref, Ji, final_x, final_xref);

	// Weigthed Errors
	Ji_error(Ji, Jui, &J_tmp1);
	Jf_error(final_x, final_xref, &J_tmp2);

	final_sum(&J, J_tmp1, J_tmp2);


//---------------------------------------------------------------------------
/*

	// Create control guess for full horizon
	// Constructing the vector of guessed control actions with respect to N and Nu
	uu_loop(control_guess, uu_1, uu_2);
	
		// Initialize x_hat vector
	for (int i = 0; i < _Nx; ++i) {
		x_hat[i] = current_state[i];
	}
	for (int i = _Nx; i < (_Nx*N); ++i) {
		x_hat[i] = 0.0;
	}

	// Calculate cost function

    // Cost of States
	for (int i = 0; i < N; ++i) {
		k = _Nx*(i+1);
        one_step_prediction(&x_hat[k], &x_hat[_Nx*i], &uu_1[_n_U*i]);
		for (int j = 0; j < _Nx; j=j+2) {
			l = k + j;
            _real tmp_err = error_with_constrains(j,x_hat[_Nx*i+j])*( x_hat[_Nx*i+j] - xref[_Nx*i+j]);
            Ji[j] = Ji[j] + tmp_err*tmp_err;
		}
	}
    for (int j = 0; j < _Nx; j=j+2) {
        J = J + Q[j]*sqrt(Ji[j]);
    }

    // Cost of States - End of Horizon
    k = _Nx*N;
    for (int j = 0; j < _Nx; j++) {
        l = k + j;
        _real tmp_err = error_with_constrains(j,x_hat[l])*( x_hat[l] - xss[j]);
        Jf = Jf + tmp_err*tmp_err;
    }   
    J = J + Qf[0]*sqrt(Jf); 

    // Cost of control
    if(R[0] > 0) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < _n_U; ++j) {
                _real tmp_err = uu_1[_n_U*i+j]-uss[j];
                Ju = Ju + tmp_err*tmp_err;
            }
        }
    }

	
*/
	return J;
}

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

void System::update_state(
	_real *state_plus, 
	_real *state, 
	_real *k, 
	_real Ts_local
	){
#pragma HLS inline
	for (unsigned i = 0; i < _Nx; ++i) {
		state_plus[i] = state[i] + Ts_local*k[i];
	}
}

void System::update_state_final(
	_real *state_plus, 
	_real *state, 
	_real *k1, 
	_real *k2, 
	_real *k3, 
	_real *k4
	){        
#pragma HLS inline
	for (unsigned i = 0; i < _Nx; ++i) {
		state_plus[i] = state[i] + Ts*0.1666666667*(k1[i] + (_real)2.0*k2[i] + (_real)2.0*k3[i] + k4[i]);
	}
}

void System::uu_loop(
	volatile _real *control_guess, 
	_real *uu_1,
	_real *uu_2
	){
// #pragma HLS inline
	_real control_val;
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < _n_U; ++j) {
			if(i < Nu) {
				control_val = control_guess[j*Nu+i];
			}
			else {
				control_val = control_guess[j*Nu-1];	
			}
			uu_1[i*_n_U+j] =  control_val;
			uu_2[i*_n_U+j] =  control_val;
		}
	}

/*
	unsigned k_cg = 0;
	unsigned k_u = 0;
	const unsigned k_cg_last = _n_U*_Nu - _n_U;
	_hw_real last_control_guess[_n_U];
	for (unsigned i = 0; i < _n_U*_N; ++i) {
		_hw_real control_val;
		if(i < _n_U*_Nu){
			control_val = control_guess[i];
			if (i >= k_cg_last){
				last_control_guess[k_cg] = control_val;
				k_cg++;
			}
		}else{
			control_val = last_control_guess[k_u];
			k_u = (k_u = _n_U-1) ? 0: k_u+1;
		}
		uu_1[i] = control_val;
		uu_2[i] = control_val;
	}
*/
}

void System::horizon_u_error(
	_real *uu,
	_real *Jui
	){
// #pragma HLS INLINE
	_hw_real Jui_buff[_n_U], Jui_buff_ant[_n_U], current_uu[_n_U];

	memset_loop(Jui_buff_ant, (const _hw_real)0.0, _n_U);
	for (unsigned i = 0; i < N; ++i) {
		memcpy_loop_rolled(current_uu, &uu[_n_U*i], _n_U);
		one_step_u_error(Jui_buff, Jui_buff_ant, current_uu);
		// one_step_u_error(Jui_buff, Jui, current_uu);
		// one_step_u_error(Jui_buff, Jui_buff, &uu[_n_U*i]);
		// memcpy_loop_enclosed<_hw_real, _hw_real, _n_U>(&uu[_n_U*i], current_uu);
		memcpy_loop_rolled(Jui_buff_ant, Jui_buff, _n_U);
	}
	memcpy_loop_rolled(Jui, Jui_buff_ant, _n_U);
}

void System::one_step_error(
	_real Ji_out[_Nx],
	volatile _real Ji_in[_Nx],
	volatile _real x_hat[_Nx],
	volatile _real xref[_Nx]
){
	_real penality;
	_real current_x_hat;
	_real current_x_ref;
//        _hw_real Ji_local[_Nx], Ji_out_local[_Nx];
//        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Ji_local, Ji_in);
	for (unsigned j = 0; j < _Nx; ++j) {
		//_hw_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
		current_x_hat = x_hat[j];
		current_x_ref = xref[j];
		_real tmp_err = (controlled_state[j] == 1) ? (current_x_hat - current_x_ref) : (_hw_real)0.0;
		penality = ((state_lower_limits[j] <= current_x_hat) && (current_x_hat <= state_upper_limits[j])) ? (_hw_real)Q[j] : (_hw_real)1e4;
//            Ji_out_local[j] = Ji_local[j] + penality*tmp_err*tmp_err;
		Ji_out[j] = Ji_in[j] + penality*tmp_err*tmp_err;
	}
//        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Ji_out, Ji_out_local);
}

void System::one_step_u_error(
	_real *Ji_out,
	volatile _real *Ji_in,
	volatile _real *uu
	// _hw_real uref[_n_U]
){
	_real penality;

	for (unsigned j = 0; j < _n_U; ++j) {
		//_hw_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
		_real current_uu = uu[j];
		_real tmp_err = (current_uu - uss);
		penality = ((u_min[j] <= current_uu) && (current_uu <= u_max[j])) ? R[j] : (_real)1e4;
		Ji_out[j] = Ji_in[j] + penality*tmp_err*tmp_err;
	}
}

void System::horizon_step_error(
	_real *x_horizon, 
	_real *xref, 
	_real *Ji,
	_real *final_x,
	_real *final_xref
	){
	_real current_x[_Nx], current_xref[_Nx]; 

	_real Ji_buff[_Nx], Ji_buff_ant[_Nx]; // = 0.0;

	memset_loop(Ji_buff_ant, (const _real)0.0, _Nx);
	for (unsigned i = 0; i < N; ++i) {
		if (i==N-1){
		split(&xref[_Nx*i], current_xref, final_x, _Nx);
		split(&x_horizon[_Nx*i], current_x, final_xref, _Nx);
		}else{
		memcpy_loop_rolled(current_xref, &xref[_Nx*i], _Nx);
		memcpy_loop_rolled(current_x, &x_horizon[_Nx*i], _Nx);
		}
		one_step_error(Ji_buff, Ji_buff_ant, current_xref, current_x);
		memcpy_loop_rolled(Ji_buff_ant, Ji_buff, _Nx);
	}
	memcpy_loop_rolled(Ji, Ji_buff_ant, _Nx);
		
}

void System::horizon_step_prediction(
	_real *x_initial,
	_real *uu, 
	_real *x_horizon
	){
	_real x_hat_out[_Nx];
	_real x_hat_in[_Nx];

	memcpy_loop_rolled(x_hat_in, x_initial, _Nx);
	for (unsigned i = 0; i < N; ++i) {
		one_step_prediction(x_hat_out, x_hat_in, &uu[_n_U*i]);
		split(x_hat_out, &x_horizon[_Nx*(i)], x_hat_in, _Nx);
	}
}

void System::final_sum(
	_real *J, 
	_real J_tmp1, 
	_real J_tmp2
	){
// #pragma HLS INLINE 
	J[0] = J_tmp1 + J_tmp2;
	
}

void System::Ji_error(
	volatile _real *Ji, 
	volatile _real *Jui,
	_real *J
	){

	_real J_local = 0.0;
	_real Ji_mul[_Nx], Jui_mul[_n_U];

	Ji_sum_loop: for (unsigned j = 0; j < _Nx; j++) {
		// J_local += Ji_mul[j];
		J_local += Ji[j];
	}

	Jui_sum_loop: for (unsigned j = 0; j < _n_U; j++) {
		// J_local += Jui_mul[j];
		J_local += Jui[j];
	}
	J[0] = J_local;
}

void System::Jf_error(
	volatile _real *final_x,
	volatile _real *final_xref,
	_real *J
	){

	// Calculate Final State Error
	_real J_local = 0.0;
	_real J_mul[_Nx];
	_real Jf; // = 0.0;
	//memset(Jf, (const _hw_real)0.0, _Nx*sizeof(_hw_real));

	Jf_mul_loop: for (unsigned j = 0; j < _Nx; ++j) {

		_real single_x_hat = final_x[j];
		_real tmp_err = (single_x_hat - final_xref[j]);
		_real penality = ((state_lower_limits[j] <= single_x_hat) && (single_x_hat <= state_upper_limits[j])) ? (_hw_real)1.0 : (_hw_real)1e4;
		Jf = penality*tmp_err*tmp_err;
		J_mul[j] = Qf[j]*Jf;
	}
	Jf_sum_loop: for (unsigned j = 0; j < _Nx; ++j) {
#pragma HLS pipeline off
		J_local += J_mul[j];
	}
	J[0] = J_local;
}