#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "fast_sin_cos.hpp"
#include "aux_functions.hpp"

typedef fast_sin_cos<float> fast_sin_cos_t;
fast_sin_cos_t fast_sin_cos_p;

int main(){
    float data[3];
    // Build cossin table
    for (int i = 0 ; i < HALF_MAX_CIRCLE_ANGLE ; i++)
    {   
        data[0] = (float)i * LOCAL_PI / (float)HALF_MAX_CIRCLE_ANGLE;
        data[1] = (float)sin(data[0]);
        data[2] = fast_sin_cos_p.fastsin(data[0]);
        
        print_formatted_float_array(data, 3, 7, 10);
        std::cout << std::endl;
    }
    return 0;
}