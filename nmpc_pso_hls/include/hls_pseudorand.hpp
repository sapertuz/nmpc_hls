#ifndef HLS_PSEUDORAND
#define HLS_PSEUDORAND

#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>


#if (defined(__SYNTHESIS__) || (defined(__VITIS__)))
#include "hls_math.h"
#include "hls_stream.h"
#include "ap_fixed.h"
#endif

#define __pseudoRand_mem_size 8

template<
    class _rand_hw_real
>class pseudoRand_gen{
    protected:
#if (defined(__SYNTHESIS__) || (defined(__VITIS__)))
    typedef ap_ufixed<64,32> _rand_local_real;
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
        12121103
        // 129057043,
        // 18431762,
        // 317028252,
        // 73549602,
        // 27304155,
        // 448127162,
        // 528520216,
        // 166120444
    };
    uint32_t __rnd_seed_co_2[__pseudoRand_mem_size];
    
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
#pragma HLS array_partition variable=__rnd_seed_co_2 type=complete
// #pragma HLS pipeline II=10
#pragma HLS inline
// #pragma HLS dependence variable=__rnd_seed_co type=inter \
//     direction=WAR distance=16*2 true 
// #pragma HLS dependence variable=__rnd_seed_co_2 type=inter \
//     direction=WAR distance=16*2 true 
    
    volatile uint32_t output_num_int = this->aux_rand_num();
    if (this->core_index < __pseudoRand_mem_size-1){
        this->core_index++; 
    }else{
        this->core_index = 0;
        flag = !flag;
    }

    _rand_local_real lo_fp = (_rand_local_real)output_num_int;
    _rand_local_real output_fp = lo_fp * divider; //+ out_min;
    // volatile _rand_local_real output_fp = (_rand_local_real)lo_fp / RAND_INT_MAX; //+ out_min;
    _rand_hw_real output_f = (_rand_local_real)output_fp;
    rand_out = (_rand_hw_real)output_f;

};
//---------------------------------------------------------------------

private:

uint32_t aux_rand_num()
{
#pragma HLS inline

    uint32_t hi,hi1,hi2,lo1,lo2,lo3,lo4;
    // volatile _rand_local_real lo_fp, output_fp;
    uint32_t local_seed;

    if (!flag)
        local_seed = this->__rnd_seed_co[this->core_index];
    else
        local_seed = this->__rnd_seed_co_2[this->core_index];

    hi = 16807 * (local_seed >> 16);
    hi1 = ((hi & 0x7FFF) << 16);//((hi & 0x7FFF) << 16);
    hi2 = (hi >> 15);

    lo1 = 16807 * (local_seed & 0xFFFF);
    lo2 = lo1 + hi1;
    
    lo3 = lo2 + hi2;

    lo4 = (lo3 > 2147483647)? lo3 - 2147483647 : lo3;
    
    if (!flag)
        this->__rnd_seed_co_2[this->core_index] = lo4;
    else
        this->__rnd_seed_co[this->core_index] = lo4;

    // lo_fp = (_rand_local_real)lo4;
    // output_fp = (_rand_local_real)lo_fp * divider; //+ out_min;
    // output_f = (_rand_local_real)output_fp;
    // this->__rnd_num[this->core_index] = output_fp;
    // std::cout << " -> " << (_rand_local_real)this->__rnd_seed_co[core] << " * " << divider << " = " << output_f << std::endl;
    // return output_f;    
    return lo4;
}

};
// }
#endif