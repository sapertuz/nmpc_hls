#ifndef HLS_PSEUDORAND
#define HLS_PSEUDORAND

#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>


#ifdef __SYNTHESIS__
#include "hls_math.h"
#include "ap_fixed.h"
#endif

template<
    class _rand_hw_real
>class pseudoRand_gen{
protected:
#ifdef __SYNTHESIS__
    typedef ap_fixed<64,32, AP_RND, AP_SAT> _rand_local_real;
#else
    typedef float _rand_local_real;
#endif    

    _rand_hw_real __rnd_num[10]={
        0.0911807,
        0.6093829,
        0.6900185,
        0.9341996,
        0.7065054,
        0.1883367,
        0.3572521,
        0.5649552,
        0.3093766,
        0.5622700
    };
    
    int __rnd_seed_co[10]={
        17279329,
        334905,
        8550184,
        14541370,
        13613342,
        18664584,
        3189472,
        12121103,
        23714,
        18431762
    };

    const _rand_local_real RAND_INT_MIN = 0.0;           // Minimum value for a variable of type int.
    const _rand_local_real RAND_INT_MAX = 2147483647.0;  // Maximum value for a variable of type int.

    char core_index = 0;

    const _rand_local_real out_min = 0.0;
    const _rand_local_real out_max = 1.0;
    const _rand_local_real divider = 1.0/2147483647.0;
public:

constexpr pseudoRand_gen(
    const _rand_hw_real _out_min,
    const _rand_hw_real _out_max
)
{
}

//---------------------------------------------------------------------
_rand_hw_real rand_num (void)
{
#pragma PIPELINE II=7
    _rand_hw_real random_number = this->__rnd_num[this->core_index];

    this->__rnd_num[this->core_index] = this->aux_rand_num(this->core_index);
    this->core_index = (this->core_index < 9)? this->core_index + 1 : 0;

    return random_number;
};
//---------------------------------------------------------------------

private:

_rand_hw_real aux_rand_num(int core)
{
#pragma HLS inline
    unsigned int hi,lo;
    hi = 16807 * (this->__rnd_seed_co[core] >> 16);
    lo = 16807 * (this->__rnd_seed_co[core] & 0xFFFF);
    lo += (hi & 0x7FFF) << 16;
    lo += hi >> 15;
    if (lo > 2147483647)
        lo -= 2147483647;
    this->__rnd_seed_co[core] = lo;
    
    _rand_hw_real output_f = (_rand_local_real)this->__rnd_seed_co[core] * divider; //+ out_min;
    // std::cout << " -> " << (_rand_local_real)this->__rnd_seed_co[core] << " * " << divider << " = " << output_f << std::endl;
    return output_f;    
}

};
// }
#endif