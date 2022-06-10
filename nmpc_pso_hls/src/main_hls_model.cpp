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
#include "config.hpp"
#include "hls_system.hpp"

#ifdef INVERTED_PENDULUM_CONFIG
    
    #define _Nh      3      // Prediction Horizon
    #define _Nu     3      // Control Horizon
    #define _Nx     4       // Number of States
    #define _n_U    1       // Number of Inputs

    const float _Ts = 0.1;
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
    // float u_guess[] = {
    // 0.000000e+00,   0.000000e+00,   0.000000e+00, 0.000000e+00, 
	// -1.746359e-01,  -7.003172e-02,  4.187292e-02, 3.835742e-02, 
	// 7.178290e-01,   5.726184e-02,   3.145507e-01, 1.726252e-01, 
	// -1.621489e+00,  1.954133e-01,   4.456669e+00, 2.455601e+00, 
	// -9.107535e-01,  -6.533037e-01,  1.435169e+00, 5.000853e-01, 
	// -4.665083e-01,  -1.498446e-02,  2.152482e+00, 4.755608e+00, 
	// 2.297003e+00,   7.750756e-01,   1.076399e+00, 4.698578e+00, 
	// 2.129837e-01,   4.823059e-01,   7.673273e+00, 3.927561e+00, 
	// -9.461644e-01,  -1.957283e+00,  3.300347e+00, 2.681746e+00, 
	// -1.871749e+00,  -1.185695e+00,  9.825810e-01, 7.133058e+00
    // };

    float u_guess [] = {
    -1.683541e+01,   -1.204010e+01,   -8.516966e+00,   -1.808013e+01 ,
	1.906993e+00,    -5.580946e+00,   1.148303e+01,    1.919867e+00 ,
	4.614624e+00,    8.722816e-01,    9.207100e+00,    2.855330e+00 ,
	-5.692396e+00,   -1.055326e+01,   -1.079290e+01,   7.825756e+00 ,
	1.402380e+01,    9.446742e+00,    -1.744006e+01,   -1.217363e+01 ,
	-4.229003e-01,   4.586183e+00,    1.230374e+00,    9.478782e-01 ,
	5.715226e-02,    -1.143759e+00,   1.346339e+01,    4.769753e+00 ,
	5.770670e+00,    9.241539e+00,    6.024680e+00,    -3.689350e+00 ,
	-1.397318e+00,   1.040769e+01,    1.519665e+01,    -5.472915e-01 ,
	-5.020665e+00,   1.624450e+00,    -2.605072e+00,   4.887021e+00 ,
	-1.434748e+00,   -7.306583e+00,   1.435703e+00,    2.127525e+01 ,
	-5.170852e+00,   -1.273662e+01,   6.583066e+00,    1.275253e+00 ,
	-3.892238e+00,   -1.035981e+01,   1.873818e+01,    -2.092604e+00 ,
	-1.018937e+01,   9.640190e+00,    1.691805e+01,    -2.209081e+01 ,
	9.810628e+00,    2.665881e+00,    1.887086e+01,    -1.538159e+01 ,
	-1.018937e+01,   1.961893e+01,    -1.129141e+00,   -8.202435e+00 ,
	9.810628e+00,    -3.810730e-01,   8.161014e+00,    1.179756e+01 ,
	-1.018937e+01,   1.913718e+01,    8.255074e+00,    -8.193270e+00 ,
	9.810628e+00,    -8.628197e-01,   2.825507e+01,    -2.819327e+01 ,
	-1.018937e+01,   8.523878e-01,    8.255074e+00,    -8.193270e+00 ,
	7.047241e+00,    4.492685e+00,    2.499290e+01,    1.180673e+01 ,
	-1.294143e+01,   -1.410734e+01,   4.992903e+00,    -8.193270e+00 ,
	7.058572e+00,    -2.583131e+01,   -1.500710e+01,   -1.476354e+01 ,
	-1.449816e-01,   -5.831314e+00,   -2.568273e+01,   -1.042127e+01 ,
	1.985502e+01,    1.416869e+01,    -2.539221e+01,   -3.042127e+01 ,
	-1.449814e-01,   1.813846e+01,    -4.539122e+01,   -1.062083e+01 ,
	1.985502e+01,    6.231643e+00,    -3.964610e+01,   3.104842e+00 ,
	1.248140e+01,    4.336839e+00,    -5.144139e+01,   -1.689516e+01 ,
	6.354776e+00,    -1.565768e+01,   -3.144139e+01,   3.104843e+00 
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

// #ifdef INVERTED_PENDULUM_CONFIG
//     typedef model_inverted_pendulum<_real> model_t;
// #elif defined(SNIFFBOT_CONFIG)
//     typedef model_sniffbot<_real> model_t;
// #endif


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

    reg_local_st: memcpy_loop_unrolled<_real, float, _Nx>(local_state, (const float *)state);
    reg_local_u_guess: memcpy_loop_unrolled<_real, float, _n_U>(local_u_guess, (const float *)u_guess);

    current_model.model(local_state_dot, local_state, local_u_guess);

    reg_st_dot: memcpy_loop_unrolled<float, _real, _Nx>(state_dot, (const _real*)local_state_dot);
    
}

void one_step_pred_wrapper(
    float state_dot[_Nx],
	float state[_Nx],
	float u_guess[_n_U]
)
{
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

    _real local_state_dot[_Nx];
    _real local_state[_Nx];
    _real local_u_guess[_n_U];

    reg_local_st: memcpy_loop_unrolled<_real, float, _Nx>(local_state, (const float *)state);
    reg_local_u_guess: memcpy_loop_unrolled<_real, float, _n_U>(local_u_guess, (const float *)u_guess);

    current_system.one_step_prediction(local_state_dot, local_state, local_u_guess);

    reg_st_dot: memcpy_loop_unrolled<float, _real, _Nx>(state_dot, (const _real*)local_state_dot);
}

int main(){
    float state_dt_dot[_Nx];
    float state_ant[_Nx];
    float state_dot[_Nx];

    float state_rk4_ant[_Nx];
    float state_rk4_dot[_Nx];

    memcpy_loop_unrolled<float, float, _Nx>(state_ant, (const float *)initial_state);
    memcpy_loop_unrolled<float, float, _Nx>(state_rk4_ant, (const float *)initial_state);
    int n_steps = sizeof(u_guess)/(sizeof(u_guess[0])*_n_U);
    printf("\n %d \n", n_steps);
    for (unsigned i = 0; i < n_steps; i++)
    {
        // Normal Prediction
        model_wrapper((float *)state_dt_dot, (float *)state_ant, (float *)&u_guess[i*_n_U]);
        //std::cout << "state = " << initial_state[0] << " , " << initial_state[1] << " , "  << initial_state[2] << " , "  << initial_state[3] << std::endl;
        for (int j = 0; j < _Nx; j++){
            state_dot[j] = state_dt_dot[j]*_Ts + state_ant[j];
        }
        // Ruhge Kutta
        one_step_pred_wrapper((float *)state_rk4_dot, (float *)state_rk4_ant, (float *)&u_guess[i*_n_U]);
        
        // std::cout << "State dt = \t";
        // print_formatted_float_array(state_dt_dot, _Nx, 4, 7);
        // std::cout << std::endl;
        std::cout << "State Normal[ " << i <<" ] = \t";
        print_formatted_float_array<float>(state_dot, _Nx, 4, 7);
        std::cout << std::endl;
        std::cout << "State RK4[ " << i <<" ] = \t";
        print_formatted_float_array<float>(state_rk4_dot, _Nx, 4, 7);
        std::cout << std::endl;
        // std::cout << "--------------------------------------------------------------------------------------------"<<std::endl;
        memcpy_loop_unrolled<float, float, _Nx>(state_ant, (const float *)state_dot);
        memcpy_loop_unrolled<float, float, _Nx>(state_rk4_ant, (const float *)state_rk4_dot);
    }
    
    return 0;
}
