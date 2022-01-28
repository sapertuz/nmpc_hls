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


void model_wrapper(
    float state_dot[_Nx],
	const float state[_Nx],
	const float u_guess[_n_U]
)
{
    typedef model_inverted_pendulum<_real, _Nx, _n_U> model_t;
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
    float state_dot[4];
    model_wrapper((float *)state_dot, (const float *)initial_state, (const float *)&u_guess[0]);
    std::cout << "state = " << initial_state[0] << " , " << initial_state[1] << " , "  << initial_state[2] << " , "  << initial_state[3] << std::endl;
    std::cout << "state_dot = " << state_dot[0] << " , "  << state_dot[1] << " , "  << state_dot[2] << " , "  << state_dot[3] << std::endl;
    return 0;
}
