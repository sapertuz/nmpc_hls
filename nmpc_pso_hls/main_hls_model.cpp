#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#ifdef __SYNTHESIS__
#include "hls_math.h"
#include "ap_fixed.h"
typedef half _real;
//typedef ap_ufixed<16,5, AP_RND_ZERO, AP_WRAP_SM> _real;
#else
typedef float _real;
#endif

#include "aux_functions.hpp"
#include "hls_inverted_pendulum.hpp"
#include "hls_sniffbot.hpp"

#ifdef INVERTED_PENDULUM_CONFIG
    
    #define _Nh      3      // Prediction Horizon
    #define _Nu     3      // Control Horizon
    #define _Nx     4       // Number of States
    #define _n_U    1       // Number of Inputs

    const float  _Ts = 0.1;
    const float  _u_max[]  = {50.0};
    const float  _u_min[]  = {-50.0};
    const float  _du_max[] = {50.0};
    const float  _uss[] = {0.0};

    #define _Parametrization 0
    const float _Lambda = 10;
    const float _q_param = 10;
    const float _pmax[] = {1.0};
    const float _pmin[] = {-1.0};

    const unsigned short _controlled_state[] = {1, 1, 1, 1};
    const float _state_upper_limits[] = {0.5, 1e3, 1e3, 1e3} ;
    const float _state_lower_limits[] = {-0.5, -1e3, -1e3, -1e3} ;
    const float _Q[] = {1e3, 0.0, 1e-1, 0.0};
    const float _Qf[] = {1e4, 0.0, 1e-0, 0.0};
    const float _R[] = {1e-4};

    #define _Rising_Time 0
    const float _tr[] =  {0, 0, 0, 0};
    float initial_state[] = {0.0, 0.0, 3.1415926536, 0.0};
    // float x_ss[] = {0.4, 0.3, 0.2, 0.1};
    float u_guess[] = {39, 50};
    float x_ref[] = {  
        0.0, 0.0, 2.0, 0.0,
        0.0, 0.0, 3.0, 0.0,
        0.0, 0.0, 3.1415926536, 0.0
    };

#elif defined(SNIFFBOT_CONFIG)
    #define _Nh 10
    #define _Nu 10
    #define _Nx 12
    #define _n_U 4
    const float  _Ts = 0.05;
    const float  _u_max[] =  {100, 100, 100, 100};
    const float  _u_min[] =  {-100, -100, -100, -100};
    const float  _du_max[] =  {20, 20, 20, 20};
    const float  _uss[] = {0.0, 0.0, 0.0, 0.0};

    #define  _Parametrization 0
    const float  _Lambda = 5;
    const float  _q_param = 2;
    const float  _pmax[] =  {1, 1, 1, 1};
    const float  _pmin[] =  {-1, -1, -1, -1};

    const unsigned short _controlled_state[] = {1, 1, 1, 1};
    const float  _state_upper_limits[] =  {1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3};
    const float  _state_lower_limits[] =  {-1e3, -1e3 -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3};
    const float  _Q[] =  {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
    const float  _Qf[] =  {10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0};
    const float  _R[] =  {0.02, 0.0, 0.0, 0.0};

    #define  _Rising_Time 0
    const float _tr[] =  {10, 10, 10, 0, 0, 8, 0, 0, 0, 0, 0, 0};
    const float  initial_state[] =  {7.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
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

#ifdef INVERTED_PENDULUM_CONFIG
    typedef model_inverted_pendulum<_real> model_t;
#elif defined(SNIFFBOT_CONFIG)
    typedef model_sniffbot<_real> model_t;
#endif


void model_wrapper(
    float state_dot[_Nx],
	float state[_Nx],
	float u_guess[_n_U]
)
{

#pragma HLS INTERFACE s_axilite port=return     bundle=control

#pragma HLS INTERFACE s_axilite port=u_guess    bundle=control
#pragma HLS INTERFACE s_axilite port=state      bundle=control
#pragma HLS INTERFACE s_axilite port=state_dot  bundle=control

#pragma HLS INTERFACE m_axi depth=4 port=u_guess offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=12 port=state offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=12 port=state_dot offset=slave bundle=output

    model_t current_model;

    _real local_state_dot[_Nx];
    _real local_state[_Nx];
    _real local_u_guess[_n_U];

    reg_local_st: memcpy_loop_rolled<_real, float, _Nx>(local_state, (const float *)state);
    reg_local_u_guess: memcpy_loop_rolled<_real, float, _n_U>(local_u_guess, (const float *)u_guess);

    current_model.model(local_state_dot, local_state, local_u_guess);

    reg_st_dot: memcpy_loop_rolled<float, _real, _Nx>(state_dot, (const _real *)local_state_dot);
    
}

int main(){
    float state_dot[_Nx];
    float state_dot_ant[_Nx];
    memcpy_loop_unrolled<float, float, _Nx>(state_dot_ant, (const float *)initial_state);
    for (unsigned i = 0; i < 10; i++)
    {
        model_wrapper((float *)state_dot, (float *)state_dot_ant, (float *)&u_guess[i*_n_U]);
        //std::cout << "state = " << initial_state[0] << " , " << initial_state[1] << " , "  << initial_state[2] << " , "  << initial_state[3] << std::endl;
        std::cout << "State = \t";
        print_formatted_float_array(state_dot, _Nx, 4, 7);
        std::cout << std::endl;
        memcpy_loop_unrolled<float, float, _Nx>(state_dot_ant, (const float *)state_dot);
    }
    
    return 0;
}
