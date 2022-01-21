#include <fstream>
#include <string.h>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "hls_pseudorand.hpp"

int main(int argc, char *argv[]) {
    // printf() displays the string inside quotation
    if (argc == 5)
    {
        int seed = atoi(argv[1]);
        float min = atof(argv[2]);
        float max = atof(argv[3]);
        int count = atoi(argv[4]);

        typedef pseudoRand_gen<float> rand_core_t;
        rand_core_t rand_core(seed, min, max);
        printf("Random Numbers : ");
        for (int i = 0; i < count; i++)
        {
            printf("%.3f, ", rand_core.rand_num());
        }
        printf("\n");
        return 0;
    }else{
        printf ("Invalid input arguments -> [1] int seed [2] float out_min [3] float out_max [4] int quantity to print \n");
        return -1;
    }
}