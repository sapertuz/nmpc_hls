#include <fstream>
#include <string.h>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "hls_pseudorand.hpp"

#define _ncores 2

int main(int argc, char *argv[]) {
    // printf() displays the string inside quotation
    if (argc == 5)
    {
        int seed[_ncores];
        seed[0] = atoi(argv[1]); 
        float min = atof(argv[2]);
        float max = atof(argv[3]);
        int count = atoi(argv[4]);

        typedef pseudoRand_gen<float, _ncores> rand_core_t;
        rand_core_t rand_core((int *)seed, (const float)min, (const float)max);
        printf("Random Numbers : ");
        for (int i = 0; i < count; i++)
        {
            printf("| %.3f |", rand_core.rand_num());
        }
        printf("\n");
        return 0;
    }else{
        printf ("Invalid input arguments -> [1] int seed [2] float out_min [3] float out_max [4] int quantity to print \n");
        return -1;
    }
}