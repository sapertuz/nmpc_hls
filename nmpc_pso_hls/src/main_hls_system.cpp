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

#include "config.hpp"
#include "aux_functions.hpp"
//#include "hls_inverted_pendulum.hpp"
//#include "hls_sniffbot.hpp"
#include "hls_system.hpp"

#ifdef INVERTED_PENDULUM_CONFIG
    
    #define _Nh      3      // Prediction Horizon
    #define _Nu     3      // Control Horizon
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

#elif defined SNIFFBOT_CONFIG
    #define _Nh 10
    #define _Nu 10
    #define _Nx 12
    #define _n_U 4
    const _real  _Ts = 0.05;
    const _real  _u_max[] =  {100, 100, 100, 100};
    const _real  _u_min[] =  {-100, -100, -100, -100};
    const _real  _du_max[] =  {20, 20, 20, 20};
    const _real  _uss[] = {0.0, 0.0, 0.0, 0.0};

    #define  _Parametrization 0
    const _real  _Lambda = 5;
    const _real  _q_param = 2;
    const _real  _pmax[] =  {1, 1, 1, 1};
    const _real  _pmin[] =  {-1, -1, -1, -1};

    const unsigned short _controlled_state[] = {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
    const _real  _state_upper_limits[_Nx] =  {1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3};
    const _real  _state_lower_limits[_Nx] =  {-1e3, -1e3 -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3};
    const _real  _Q[] =  {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
    const _real  _Qf[] =  {10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0};
    const _real  _R[] =  {0.02, 0.0, 0.0, 0.0};

    #define  _Rising_Time 0
    const _real _tr[] =  {10, 10, 10, 0, 0, 8, 0, 0, 0, 0, 0, 0};
    const _real  initial_state[] =  {7.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    float u_guess[] = {
    0.000000e+00,   0.000000e+00,   0.000000e+00, 0.000000e+00, 
	-1.746359e-01,  -7.003172e-02,  4.187292e-02, 3.835742e-02, 
	7.178290e-01,   5.726184e-02,   3.145507e-01, 1.726252e-01, 
	-1.621489e+00,  1.954133e-01,   4.456669e+00, 2.455601e+00, 
	-9.107535e-01,  -6.533037e-01,  1.435169e+00, 5.000853e-01, 
	-4.665083e-01,  -1.498446e-02,  2.152482e+00, 4.755608e+00, 
	2.297003e+00,   7.750756e-01,   1.076399e+00, 4.698578e+00, 
	2.129837e-01,   4.823059e-01,   7.673273e+00, 3.927561e+00, 
	-9.461644e-01,  -1.957283e+00,  3.300347e+00, 2.681746e+00, 
	-1.871749e+00,  -1.185695e+00,  9.825810e-01, 7.133058e+00
    };
    float x_ref[] = {  
    7.078500e+00, 1.000100e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.141600e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
	7.078500e+00, 1.000100e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.141600e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
	7.157000e+00, 1.000500e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 6.283200e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
	7.157000e+00, 1.000500e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 6.283200e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
	7.235300e+00, 1.001100e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 9.424800e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
	7.235300e+00, 1.001100e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 9.424800e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
	7.313300e+00, 1.002000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.256600e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
	7.313300e+00, 1.002000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.256600e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
	7.391100e+00, 1.003100e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.570800e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
	7.391100e+00, 1.003100e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.570800e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
    };

#endif

float cost_function_wrapper(
    volatile float *control_guess,
	volatile float *xref,
	volatile float *current_state
)
{
#pragma HLS INTERFACE s_axilite port=return         bundle=control

#pragma HLS INTERFACE s_axilite port=control_guess  bundle=control
#pragma HLS INTERFACE s_axilite port=xref           bundle=control
#pragma HLS INTERFACE s_axilite port=current_state  bundle=control

#pragma HLS INTERFACE m_axi depth=40  port=control_guess offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=120 port=xref          offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=12  port=current_state offset=slave bundle=input


    typedef System<_real, _Nh, _Nx, _n_U, _Nu> T_system;
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
        _Ts);

    _real local_control_guess[_n_U*_Nu];
    _real local_xref[_Nx*_Nh];
    _real local_current_state[_Nx];

#pragma HLS bind_storage variable=local_control_guess type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=local_xref type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=local_current_state type=FIFO impl=LUTRAM

    reg_curr_st: memcpy_loop_rolled<_real, float, _Nx>(local_current_state, current_state);
    reg_cont_gss: memcpy_loop_rolled<_real, float, _n_U*_Nu>(local_control_guess, control_guess);
    reg_xref: memcpy_loop_rolled<_real, float, _Nx*_Nh>(local_xref, xref);

    _real cf;
    current_system.nmpc_cost_function(local_current_state, local_control_guess, local_xref, &cf);
    // current_system.nmpc_cost_function_topflow(local_current_state, local_control_guess, local_xref, &cf);
    return (float) cf;
}

int main(){
    
    float j = cost_function_wrapper((volatile float *)u_guess, (volatile float *)x_ref, (volatile float *)initial_state);
    std::cout << "J = " << j << std::endl;
    return 0;
}
