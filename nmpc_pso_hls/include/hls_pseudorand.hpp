#ifndef HLS_PSEUDORAND
#define HLS_PSEUDORAND

#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>


#ifdef __SYNTHESIS__
#include "hls_math.h"
#include "hls_stream.h"
#include "ap_fixed.h"
#endif

#define __pseudoRand_mem_size 16

template<
    class _rand_hw_real
>class pseudoRand_gen{
    protected:
#ifdef __SYNTHESIS__
    typedef ap_fixed<64,32, AP_RND, AP_SAT> _rand_local_real;
#else
    typedef float _rand_local_real;
#endif    
    
    _rand_hw_real __rnd_num[__pseudoRand_mem_size];
    // {
        // 0.0911807,
        // 0.6093829,
        // 0.6900185,
        // 0.9341996,
        // 0.7065054,
        // 0.1883367,
        // 0.3572521,
        // 0.5649552
        // ,
        // 0.3093766,
        // 0.5622700
    // };
    
    uint32_t __rnd_seed_co[__pseudoRand_mem_size]={
        17279329,
        722157382,
        8550184,
        14541370,
        13613342,
        18664584,
        387505572,
        12121103,
        129057043,
        18431762,
        317028252,
        73549602,
        27304155,
        448127162,
        528520216,
        166120444
    };
    
    const _rand_local_real RAND_INT_MIN = 0.0;           // Minimum value for a variable of type int.
    const _rand_local_real RAND_INT_MAX = 2147483647.0;  // Maximum value for a variable of type int.

    int32_t core_index = 0;
    char flag = 0;

    const _rand_local_real out_min = 0.0;
    const _rand_local_real out_max = 1.0;
    const _rand_local_real divider = 1.0/2147483647.0;
public:

pseudoRand_gen()
{}

//---------------------------------------------------------------------
void rand_num (
    _rand_hw_real &rand_out
){
#pragma HLS array_partition variable=__rnd_seed_co type=complete
#pragma HLS pipeline II=10

// #pragma HLS dependence variable=__rnd_seed_co class=array type=inter \
//     direction=WAR distance=8 true 

    _rand_hw_real random_number;
    this->aux_rand_num();
    rand_out = this->__rnd_num[this->core_index];
    this->core_index = (this->core_index < __pseudoRand_mem_size-1)?\
        this->core_index + 1 : 0;

};
//---------------------------------------------------------------------

private:

void aux_rand_num()
{
#pragma HLS inline
    uint32_t hi,lo;
    uint32_t local_seed = this->__rnd_seed_co[this->core_index];
    hi = 16807 * (local_seed >> 16);
    lo = 16807 * (local_seed & 0xFFFF);
    lo += (hi & 0x7FFF) << 16;
    lo += hi >> 15;
    lo = (lo > 2147483647)? lo - 2147483647 : lo;
    this->__rnd_seed_co[this->core_index] = lo;
    
    _rand_local_real output_fp = (_rand_local_real)lo * divider; //+ out_min;
    this->__rnd_num[this->core_index] = output_fp;
    // std::cout << " -> " << (_rand_local_real)this->__rnd_seed_co[core] << " * " << divider << " = " << output_f << std::endl;
    // return output_f;    
}

};
// }
#endif