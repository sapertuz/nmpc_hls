#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "aux_functions.hpp"
#include "hls_pseudorand.hpp"

#ifdef __SYNTHESIS__
#include "hls_math.h"
#include "hls_stream.h"
#include "ap_fixed.h"
typedef float _rand_real;
// typedef ap_fixed<48,5, AP_RND, AP_SAT> _rand_real;
#else
typedef float _rand_real;
#endif

const _rand_real rand_min = 0.0;
const _rand_real rand_max = 1.0;
typedef pseudoRand_gen<_rand_real> _randCore_t;
_randCore_t _hw_rand_core(
    (const _rand_real)rand_min,
    (const _rand_real)rand_max
);
// Nonlinear PSO Solver
// _system_t *_hw_system_ptr = &_hw_system;
_randCore_t *_hw_rand_core_ptr = &_hw_rand_core;


void rand_wrapper(_rand_real *rand_num_arr)
{
// #pragma HLS interface s_axilite port=return
    _rand_real rand_num_out;
    for (int i = 0; i < 20; i++){   
// #ifndef __SYNTHESIS__
    _hw_rand_core.rand_num(rand_num_out);
// #else
//     hls::stream< _rand_real > rand_num_stream;
// #pragma HLS STREAM variable=rand_num_stream depth=10 type=fifo
//     _hw_rand_core.rand_num(rand_num_stream);
//     rand_num_stream.read(rand_num_out);
// #endif
    rand_num_arr[i] = rand_num_out;
    }
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
    float num[20];
    
    for (int k = 0; k < 10; k++){
    rand_wrapper(num);
    for (int i = 0; i < 10; i++){
        std::cout << std::fixed << std::setprecision(3) << num[i] << ", \t";
    }
    std::cout << std::endl;
    }

    std::cout << "--------------------" << std::endl;
    std::cout << "End Test" << std::endl;
    return 0;
}
