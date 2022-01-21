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

float pso_wrapper(
)
{


}

int main(){
    
    float j = cost_function_wrapper((const float *)u_guess, (const float *)x_ref, (const float *)initial_state);
    std::cout << "Error = " << j << std::endl;
    return 0;
}
