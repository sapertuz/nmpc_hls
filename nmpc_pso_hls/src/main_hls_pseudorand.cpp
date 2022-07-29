#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "aux_functions.hpp"
#include "hls_pseudorand.hpp"

//#define __VITIS__

#if (defined(__SYNTHESIS__) || (defined(__VITIS__)))
#include "hls_math.h"
#include "hls_stream.h"
#include "ap_fixed.h"
typedef float _rand_real;
// typedef ap_fixed<48,5, AP_RND, AP_SAT> _rand_real;
typedef hls::stream< _rand_real> _rand_real_stream;
#else
typedef float _rand_real;
typedef _rand_real _rand_real_stream;
#endif



const _rand_real rand_min = 0.0;
const _rand_real rand_max = 1.0;
typedef pseudoRand_gen<_rand_real> _randCore_t;
_randCore_t _hw_rand_core;
// Nonlinear PSO Solver
// _system_t *_hw_system_ptr = &_hw_system;
// _randCore_t *_hw_rand_core_ptr = &_hw_rand_core;


void rand_wrapper(
	// bool enable,
    _rand_real_stream &rand_out
)
{
// #pragma HLS INTERFACE mode=axis 	port=rand_out
// #pragma HLS interface mode=ap_none  port=enable
// #pragma HLS interface mode=axis     port=rand_out
// #pragma HLS interface mode=ap_ctrl_none port=return

#pragma HLS stream variable=rand_out type=fifo depth=16

#if  (defined(__SYNTHESIS__) || (defined(__VITIS__)))
    _rand_real rand_num_out;
	while(1){
#pragma HLS pipeline II=12
	if (!rand_out.full()){
		_hw_rand_core.rand_num(rand_num_out);
		rand_out.write(rand_num_out);
	}
    }
#else
    _rand_real rand_num_out;
    _hw_rand_core.rand_num(rand_num_out);
    rand_out = rand_num_out;
#endif
}

// _rand_real rand_wrapper()
// {
// #pragma HLS interface s_axilite port=return
//     _rand_real rand_num_out; 
//     _hw_rand_core.rand_num(rand_num_out);
//     return rand_num_out;
// }

int main(){
    std::cout << "Test Pseudorand Core" << std::endl;
    std::cout << "--------------------" << std::endl;
    _rand_real_stream num;
    float num_f;
    
    for (int k = 0; k < 10; k++){
    for (int i = 0; i < 10; i++){
        rand_wrapper(num);
#if (defined(__SYNTHESIS__) || (defined(__VITIS__)))
        num.read(num_f);
#else
        num_f = num;
#endif
        std::cout << std::fixed << std::setprecision(3) << num_f << ", \t";
    }
    std::cout << std::endl;
    }
    rand_wrapper(num);
    
    std::cout << "--------------------" << std::endl;
    std::cout << "End Test" << std::endl;
    
    return 0;
}
