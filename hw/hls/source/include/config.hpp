#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <time.h>
#include <sys/time.h>

// #ifdef __SYNTHESIS__
// #include "hls_math.h"
// typedef half _real_model;
// #else
// typedef float _real_model;
// #endif
// Definitions

#define N_DOFS 7
#define BILLION 1E9

#ifdef __SYNTHESIS__
#include "hls_math.h"
#include "hls_stream.h"
#include "ap_fixed.h"

// typedef half _hw_top_real;
typedef ap_fixed<32,17,AP_RND_INF,AP_SAT_SYM> _hw_top_real; // round to +/-inf and simmetrical saturation
typedef half _hw_model_real;
// typedef ap_fixed<32,16,AP_RND_INF,AP_SAT_SYM> _hw_model_real;
typedef hls::stream<_hw_top_real> _rand_real_stream;
#else
typedef float _hw_top_real;
typedef float _hw_model_real;
#endif

// #ifdef PSO_CONFIG
#define _KPSO 1
#define _stable_zero 1
#define _n_S 10
#define _maxiter 100
#define _max_v 30
#define _w0 0.9
#define _wf 0.1 // Drone
// #define _wf 0.4 // Inverted Pendulum
#define _c1 2.1 
#define _c2 1.0
#define _threshold 1e-2
#define _stop_criteria 0
#define _slope (_wf-_w0)/_maxiter
const _hw_top_real rand_min = 0.0;
const _hw_top_real rand_max = 1.0;
// const int rand_seed[] = {
    // 84680,
    // 577726,
    // 273600,
    // 804402,
    // 747952
// };
// #endif

#ifdef INVERTED_PENDULUM_CONFIG

    #define _Nh     25       // Prediction Horizon
    #define _Nu     25       // Control Horizon
    #define _Nx     4         // Number of States
    #define _n_U    1         // Number of Inputs

    const _hw_top_real  _Ts = 0.1;
    const _hw_top_real  _u_max[]  = {50.0};
    const _hw_top_real  _u_min[]  = {-50.0};
    const _hw_top_real  _du_max[] = {50.0};
    const _hw_top_real  _uss[] = {0.0};

    #define _Parametrization 0
    const _hw_top_real _Lambda = 10;
    const _hw_top_real _q_param = 10;
    const _hw_top_real _pmax[] = {1.0};
    const _hw_top_real _pmin[] = {-1.0};

    const unsigned short _controlled_state[] = {1, 1, 2, 0};
    const _hw_top_real _state_upper_limits[] = {1, 1e3, 1e3, 1e3} ;
    const _hw_top_real _state_lower_limits[] = {-1, -1e3, -1e3, -1e3} ;
    const _hw_top_real _Q[] = { 2  , 0.015, 1  , 0.0};
    const _hw_top_real _Qf[] = {.1, 0.015, 0.1, 0.0};
    const _hw_top_real _R[] = {.001};

    #define _Rising_Time 0
    const _hw_top_real _tr[] =  {0, 0, 0, 0};
    // float x_ss[] = {0.4, 0.3, 0.2, 0.1};
#elif defined(SNIFFBOT_CONFIG)

    #define _Nh 25
    #define _Nu 25
    #define _Nx 12
    #define _n_U 4
    const _hw_top_real  _Ts = 0.05;
    const _hw_top_real  _u_max[] =  {100, 100, 100, 100};
    const _hw_top_real  _u_min[] =  {-100, -100, -100, -100};
    const _hw_top_real  _du_max[] =  {
        20, 
        20, 
        20, 
        20
    };
    const _hw_top_real  _uss[] = {0.0, 0.0, 0.0, 0.0};

    #define  _Parametrization 0
    const _hw_top_real  _Lambda = 5;
    const _hw_top_real  _q_param = 2;
    const _hw_top_real  _pmax[] =  {1, 1, 1, 1};
    const _hw_top_real  _pmin[] =  {-1, -1, -1, -1};

    const _hw_top_real drone_pi = 3.14159265359;
    const unsigned short _controlled_state[] = {1, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0};
    const _hw_top_real  _state_upper_limits[_Nx] =  { 
         1e3,  
         1e3,  
         1e3,  
         drone_pi*(_hw_top_real).25,  
         drone_pi*(_hw_top_real).25,  
         1e3,  
         1e3, 1e3, 1e3, 1e3, 1e3, 1e3
    };
    const _hw_top_real  _state_lower_limits[_Nx] =  {
        -1e3, 
        -1e3, 
        -1e3, 
        -drone_pi*(_hw_top_real).25, 
        -drone_pi*(_hw_top_real).25, 
        -1e3, 
        -1e3, -1e3, -1e3, -1e3, -1e3, -1e3};
    const _hw_top_real  _Q[] =  {
        1, 
        1, 
        1, 
        3, 
        3, 
        3, 
        0, 0, 0, 0, 0, 0};
    const _hw_top_real  _Qf[] =  {
        10, 
        10, 
        10, 
        30, 
        30, 
        30, 
        0, 0, 0, 0, 0, 0};
    const _hw_top_real  _R[] =  {
        0.0005, 
        0.0005, 
        0.0005, 
        0.0005
    };

    #define  _Rising_Time 0
    const _hw_top_real _tr[] =  {10, 10, 10, 0, 0, 8, 0, 0, 0, 0, 0, 0};
#endif

#if defined(_POSIX_MONOTONIC_CLOCK)
/*  The identifier for the system-wide monotonic clock, which is defined
 *  as a clock whose value cannot be set via clock_settime() and which
 *  cannot have backward clock jumps. */

    #define CLOCK_ID CLOCK_MONOTONIC
#else
    #define CLOCK_ID CLOCK_REALTIME
#endif

const int _pragma_Nh = _Nh;
const int _pragma_Nu = _Nu;
const int _pragma_Nx = _Nx;
const int _pragma_n_U = _n_U;

#ifndef __SYNTHESIS__
const _hw_top_real H_MAX = 3.4028234664e30;
#else
const _hw_top_real H_MAX = 60000.0;
#endif

#endif
