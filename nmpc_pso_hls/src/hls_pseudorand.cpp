#include "hls_pseudorand.hpp"

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

//---------------------------------------------------------------------
float pseudoRand_gen::rand_num (void)
{
//#pragma HLS inline
#pragma HLS ALLOCATION type=operation instances=mul limit=1
#pragma HLS ALLOCATION type=operation instances=sub limit=1
#pragma HLS ALLOCATION type=operation instances=add limit=1
#pragma HLS ALLOCATION type=operation instances=fmul limit=1
#pragma HLS ALLOCATION type=operation instances=fsub limit=1
#pragma HLS ALLOCATION type=operation instances=fadd limit=1

    return this->aux_rand_num(&__rnd_seed_co[0]);
};
//---------------------------------------------------------------------
/*
_hw_real pseudoRand_gen::rand_num (unsigned core)
{
//#pragma HLS inline
#pragma HLS ALLOCATION type=operation instances=mul limit=1
#pragma HLS ALLOCATION type=operation instances=sub limit=1
#pragma HLS ALLOCATION type=operation instances=add limit=1
#pragma HLS ALLOCATION type=operation instances=fmul limit=1
#pragma HLS ALLOCATION type=operation instances=fsub limit=1
#pragma HLS ALLOCATION type=operation instances=fadd limit=1

    return aux_rand_num(&rnd_seed_co[core]);
};
*/
//---------------------------------------------------------------------
float pseudoRand_gen::aux_rand_num(int *_rnd_seed)
{
#pragma HLS inline
    unsigned int hi,lo;
    hi = 16807 * (_rnd_seed[0] >> 16);
    lo = 16807 * (_rnd_seed[0] & 0xFFFF);
    lo += (hi & 0x7FFF) << 16;
    lo += hi >> 15;
    if (lo > 2147483647)
        lo -= 2147483647;
    _rnd_seed[0] = lo;
    //printf("%.1f ", ((float)rnd_seed[0]));

    // int k1;
    // int ix = rnd_seed;
	
    // k1 = ix / 127773;
    // ix = 16807 * (ix - k1 * 127773) - k1 * 2836;
    // if (ix < 0)
    //     ix += 2147483647;
    // rnd_seed = ix;
    
    float output_f = ((float)_rnd_seed[0] - RAND_INT_MIN) * (out_max - out_min) * divider + out_min;
    return output_f;    
}