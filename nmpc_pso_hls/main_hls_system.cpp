#include <fstream>
#include <string.h>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#ifdef __SYNTHESIS__
#include "hls_math.h"
typedef half _real;
#else
typedef float _real;
#endif

#include "aux_functions.hpp"
#include "hls_system.hpp"
#include "hls_inverted_pendulum.hpp"

#define _N      3       // Prediction Horizon
#define _Nu     2       // Control Horizon
#define _Nx     4       // Number of States
#define _n_U    1       // Number of Inputs

const _real  _Ts = 0.1;
const _real  _u_max[]  = {50.0};
const _real  _u_min[]  = {-50.0};
const _real  _du_max[] = {50.0};
const _real  _uss[] = {0.0};

#define _Parametrization 0
const _real _Lambda = 10;
const _real _q_param = 10;
const _real _pmax[] = {1.0};
const _real _pmin[] = {-1.0};

const unsigned short _controlled_state[] = {1, 1, 1, 1};
const _real _state_upper_limits[] = {0.5, 1e3, 1e3, 1e3} ;
const _real _state_lower_limits[] = {-0.5, -1e3, -1e3, -1e3} ;
const _real _Q[] = {1e3, 0.0, 1e-1, 0.0};
const _real _Qf[] = {1e4, 0.0, 1e-0, 0.0};
const _real _R[] = {1e-4};

#define _Rising_Time 0
const _real _tr[] =  {0, 0, 0, 0};

float initial_state[] = {0.0, 0.0, 3.1415926536, 0.0};
// float x_ss[] = {0.4, 0.3, 0.2, 0.1};
float u_guess[] = {39, 50};
float x_ref[] = {  
    0.0, 0.0, 2.0, 0.0,
    0.0, 0.0, 3.0, 0.0,
    0.0, 0.0, 3.1415926536, 0.0
};

float cost_function_wrapper(
    const float control_guess[_n_U*_Nu],
	const float xref[_Nx*_N],
	const float current_state[_Nx]
)
{
#pragma HLS interface ap_fifo port=control_guess
#pragma HLS interface ap_fifo port=xref
#pragma HLS interface ap_fifo port=current_state

    typedef model_inverted_pendulum<_real, _Nx, _n_U> model_t;
    model_t current_model;

    typedef System<_real, model_t, _N, _Nx, _n_U, _Nu> T_system;
    T_system current_system(
        _u_max, 
        _u_min, 
        _du_max,
        _controlled_state,
        _state_upper_limits, 
        _state_lower_limits, 
        _Q, 
        _Qf, 
        _R, 
        _uss, 
        _Ts,
        current_model);

    _real local_control_guess[_n_U*_Nu];
    _real local_xref[_Nx*_N];
    _real local_current_state[_Nx];

// #pragma HLS interface ap_fifo port=local_control_guess
// #pragma HLS interface ap_fifo port=local_xref
// #pragma HLS interface ap_fifo port=local_current_state

    reg_curr_st: memcpy_loop_rolled<_real, float, _Nx>(local_current_state, (const float *)current_state);
    reg_cont_gss: memcpy_loop_rolled<_real, float, _n_U*_Nu>(local_control_guess, (const float *)control_guess);
    reg_xref: memcpy_loop_rolled<_real, float, _Nx*_N>(local_xref, (const float *)xref);
    // for (unsigned short i = 0; i < _n_U*_Nu; i++)
    //     local_control_guess[i] = control_guess[i];
    // for (unsigned short i = 0; i < _Nx*_N; i++)
    //     local_xref[i] = xref[i];

    _real cf = current_system.nmpc_cost_function(local_current_state, local_control_guess, local_xref);
    return (float) cf;
}

int main(){
    
    float j = cost_function_wrapper((const float *)u_guess, (const float *)x_ref, (const float *)initial_state);
    std::cout << "Error = " << j << std::endl;
    return 0;
}
