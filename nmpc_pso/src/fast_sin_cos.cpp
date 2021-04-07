#include "fast_sin_cos.h"

float fast_cossin_table[MAX_CIRCLE_ANGLE];           // Declare table of fast cosinus and sinus

void initSinCosTable(){
    // Build cossin table
    for (int i = 0 ; i < MAX_CIRCLE_ANGLE ; i++)
    {
        fast_cossin_table[i] = (float)sin((double)i * PI / HALF_MAX_CIRCLE_ANGLE);
    }
}

float fastcos(float n)
{
    float f = n * HALF_MAX_CIRCLE_ANGLE / PI;
    int i;

    i = (int) f;
    if (i < 0) {
        return fast_cossin_table[((-i) + QUARTER_MAX_CIRCLE_ANGLE)&MASK_MAX_CIRCLE_ANGLE];
    }
    else {
        return fast_cossin_table[(i + QUARTER_MAX_CIRCLE_ANGLE)&MASK_MAX_CIRCLE_ANGLE];
    }
    assert(0);
}

float fastsin(float n)
{
    float f = n * HALF_MAX_CIRCLE_ANGLE / PI;
    int i;
    
    i = (int) f;
    if (i < 0) {
        return fast_cossin_table[(-((-i)&MASK_MAX_CIRCLE_ANGLE)) + MAX_CIRCLE_ANGLE];
    }
    else {
        return fast_cossin_table[i&MASK_MAX_CIRCLE_ANGLE];
    }
    assert(0);
}