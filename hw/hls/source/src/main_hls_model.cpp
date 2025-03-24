#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "aux_functions.hpp"
#include "config.hpp"
#include "hls_model.hpp"

// #ifdef __SYNTHESIS__
// #include "hls_math.h"
// #include "ap_fixed.h"
// typedef half _model_real;
// //typedef ap_ufixed<16,5, AP_RND_ZERO, AP_WRAP_SM> _real;
// #else
// typedef float _model_real;
// #endif

typedef _hw_top_real _model_real;


#ifdef INVERTED_PENDULUM_CONFIG
float initial_state[] = {0.0, 0.0, 3.1415926536, 0.0};
// float x_ss[] = {0.4, 0.3, 0.2, 0.1};
float u_guess[] = {39, 50};
float x_ref[] = {  
    0.0, 0.0, 2.0, 0.0,
    0.0, 0.0, 3.0, 0.0,
    0.0, 0.0, 3.1415926536, 0.0
};

#elif defined(SNIFFBOT_CONFIG)
const float  initial_state[] =  {7.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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

void model_wrapper(
    volatile float state_dot[_Nx],
	volatile float state[_Nx],
	volatile float u_guess[_n_U]
)
{

#pragma HLS INTERFACE s_axilite port=return     bundle=control

#pragma HLS INTERFACE s_axilite port=u_guess    bundle=control
#pragma HLS INTERFACE s_axilite port=state      bundle=control
#pragma HLS INTERFACE s_axilite port=state_dot  bundle=control

#pragma HLS INTERFACE m_axi depth=4  port=u_guess   offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=12 port=state     offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=12 port=state_dot offset=slave bundle=output

    // top_model_t current_model;

    _model_real local_state_dot[_Nx];
    _model_real local_state[_Nx];
    _model_real local_u_guess[_n_U];

    reg_local_st: memcpy_loop_rolled<_model_real, float, _Nx>(local_state, (float *)state);
    reg_local_u_guess: memcpy_loop_rolled<_model_real, float, _n_U>(local_u_guess, (float *)u_guess);

    __model<_model_real>(local_state_dot, local_state, local_u_guess);

    reg_st_dot: memcpy_loop_rolled<float, _model_real, _Nx>((float*)state_dot, local_state_dot);
    
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

        std::cout << "State Normal[ " << i <<" ] = \t";
        print_formatted_float_array<float>(state_dot, _Nx, 4, 7);
        std::cout << std::endl;

        memcpy_loop_unrolled<float, float, _Nx>(state_ant, (const float *)state_dot);

    }
    
    return 0;
}
