#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "aux_functions.hpp"
#include "hls_pseudorand.hpp"

#ifdef __SYNTHESIS__
#include "hls_math.h"
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


float rand_wrapper()
{
#pragma HLS interface s_axilite port=return
    float rand_num = _hw_rand_core.rand_num();
    return rand_num;
}

int main(){
    std::cout << "Test Pseudorand Core" << std::endl;
    std::cout << "--------------------" << std::endl;

    for (int k = 0; k < 10; k++){
    for (int i = 0; i < 5; i++){
        float num = rand_wrapper();
        std::cout << num << ", \t";
    }
    std::cout << std::endl;
    }

    std::cout << "--------------------" << std::endl;
    std::cout << "End Test" << std::endl;
    return 0;
}
