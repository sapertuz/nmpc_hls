#ifndef __ARITH_AUX_HPP__
#define __ARITH_AUX_HPP__

#ifdef __SYNTHESIS__
#include "fast_sin_cos.hpp"
using namespace hls;
//typedef ap_ufixed<32,8, AP_RND_ZERO, AP_WRAP_SM> _model_real;
#else
#include "math.h"
#endif

template <
    typename _data_type
>_data_type value_map(
    _data_type input, 
    _data_type in_min, 
    _data_type in_max, 
    _data_type out_min, 
    _data_type out_max
){
#pragma HLS inline
    _data_type output = (input - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
    return output;
}

template <
    typename _data_type
>_data_type local_sin(_data_type angle){
#if defined(__SYNTHESIS__)
#ifdef USE_FAST_SIN_COS
#pragma HLS inline
    _data_type new_angle = angle * (_data_type)HALF_MAX_CIRCLE_ANGLE_PI;
    return fastsin<_data_type>(new_angle);
    // return half_sin(angle);
#else
    return half_sin(angle);
#endif
#else        
    return sin(angle);
#endif
}

template <
    typename _data_type
>_data_type local_cos(_data_type angle){
#pragma HLS inline
#if defined(__SYNTHESIS__)
#ifdef USE_FAST_SIN_COS
	_data_type new_angle = angle * (_data_type)HALF_MAX_CIRCLE_ANGLE_PI;
    return fastcos<_data_type>(new_angle);
    // return half_cos(angle);
#else
    return half_cos(angle);
#endif
#else        
    return cos(angle);
#endif
}

#endif