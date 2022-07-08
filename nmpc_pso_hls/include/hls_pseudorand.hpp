#ifndef HLS_PSEUDORAND
#define HLS_PSEUDORAND

#define RAND_INT_MIN 0.0f           // Minimum value for a variable of type int.
#define RAND_INT_MAX 2147483647.0f  // Maximum value for a variable of type int.

// namespace pseudo_rand{
class pseudoRand_gen{
protected:
    // int __rnd_seed_co[10]={
    //     17279329,
    //     334905,
    //     8550184,
    //     14541370,
    //     13613342,
    //     18664584,
    //     3189472,
    //     12121103,
    //     23714,
    //     18431762
    // };
    const float out_min;
    const float out_max;
    const float divider;
public:
constexpr pseudoRand_gen(
    const float _out_min,
    const float _out_max
) : out_min(_out_min), 
    out_max(_out_max),
    divider(1.0f/(RAND_INT_MAX - RAND_INT_MIN))
{
}
    float rand_num (void);
    float rand_num (unsigned core);

private:
    float aux_rand_num(int *_rnd_seed);
};
// }
#endif