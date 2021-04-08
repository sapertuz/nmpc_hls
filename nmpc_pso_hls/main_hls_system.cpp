#include <fstream>
#include <string.h>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "hls_system.hpp"

#define _N      3       // Prediction Horizon
#define _Nu     2       // Control Horizon
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

const unsigned char _controlled_state[] = {1, 1, 1, 1};
const float _state_upper_limits[] = {0.5, 1e10, 1e10, 1e10} ;
const float _state_lower_limits[] = {-0.5, -1e10, -1e10, -1e10} ;
const float _Q[] = {1e4, 1.0, 1e4, 1.0};
const float _Qf[] = {1e3, 1e3, 1e3, 1e3};
const float _R[] = {1.0};

#define _Rising_Time 0
const float _tr[] =  {0, 0, 0, 0};

float initial_state[] = {0.0, 0.0, 3.1415926536, 0.0};
// float x_ss[] = {0.4, 0.3, 0.2, 0.1};
float u_guess[] = {50, 50};
float x_ref[] = {  
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 3.1415926536, 0.0,
    0.0, 0.0, 3.1415926536, 0.0
};
float u_ref[] = {0, 0};

float cost_function_wrapper(
    float control_guess[_n_U*_Nu],
    float xref[_Nx*_N],
    
    float current_state[_Nx]
)
{
    typedef System<float, _N, _Nx, _n_U, _Nu> T_system;
    T_system current_system(
        _u_max, 
        _u_min, 
        _du_max, 
        _controlled_state,
        _state_upper_limits, 
        _state_lower_limits, 
        _Q, 
        _Qf, 
        // current_state, 
        _R, 
        _uss, 
        _Ts);

    float cf = current_system.nmpc_cost_function(current_state, control_guess, xref);
    return cf;
}

int main(){
    
    float j = cost_function_wrapper(u_guess, x_ref, initial_state);
    std::cout << "Error = " << j << std::endl;
    return 0;
}