#ifndef HLS_PSEUDORAND
#define HLS_PSEUDORAND

#define RAND_INT_MIN 0.0f           // Minimum value for a variable of type int.
#define RAND_INT_MAX 2147483647.0f  // Maximum value for a variable of type int.

template<
    class _hw_real
>class pseudoRand_gen{
protected:
    int rnd_seed;
    const float out_min;
    const float out_max;
    const float divider;

public:
constexpr pseudoRand_gen(
    int new_seed,
    const float _out_min,
    const float _out_max
) : rnd_seed(new_seed), out_min(_out_min), out_max(_out_max),
    divider(1.0f/(RAND_INT_MAX - RAND_INT_MIN))
{
}

_hw_real rand_num (void)
{
#pragma HLS inline
#pragma HLS ALLOCATION instances=mul limit=1 operation
#pragma HLS ALLOCATION instances=fmul limit=1 operation
#pragma HLS ALLOCATION instances=fsub limit=1 operation
#pragma HLS ALLOCATION instances=fadd limit=1 operation

    unsigned int hi,lo;
    hi = 16807 * (rnd_seed >> 16);
    lo = 16807 * (rnd_seed & 0xFFFF);
    lo += (hi & 0x7FFF) << 16;
    lo += hi >> 15;
    if (lo > 2147483647)
        lo -= 2147483647;
    rnd_seed = lo;

    // int k1;
    // int ix = rnd_seed;
	
    // k1 = ix / 127773;
    // ix = 16807 * (ix - k1 * 127773) - k1 * 2836;
    // if (ix < 0)
    //     ix += 2147483647;
    // rnd_seed = ix;
    
    float output = ((float)rnd_seed - RAND_INT_MIN) * (out_max - out_min) * divider + out_min;;

    return (_hw_real)output;
}
};

#endif