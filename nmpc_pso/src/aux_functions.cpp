#include "aux_functions.hpp"
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#ifdef __MACH__
#include <sys/time.h>
//#define CLOCK_REALTIME 0
//#define CLOCK_MONOTONIC 0
//clock_gettime is not implemented on OSX
int clock_gettime(int clock, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}
#endif

_real min_array(_real * array, int * pos, int size) {
	//[bestfitness(k), p] = min(f_ind);
	_real min = array[0];
	int pos_temp = 0;
	for (int i = 1; i < size; ++i) {
        if(array[i] < min) {
			min = array[i];
			pos_temp = i;
        }
	}
	*pos = pos_temp;
	return min;
}

_real rand_real(){
	return (_real)rand()/(_real)RAND_MAX;
}

_real ** alloc_matrix(int n, int m){
    _real ** matrix;
    if ((matrix = (_real **) malloc(n*sizeof(_real *))) == NULL) {
        std::cout << "Erro de alocação" << std::endl;
        return NULL;
    }
    else {
        for (int i = 0; i < n; i++) {
            if((matrix[i] = (_real *) calloc(m, sizeof(_real))) == NULL) {
                std::cout << "Erro de alocação 2" << std::endl;
                return NULL;
            }
        }
    }
    return matrix;
}

void free_matrix(_real ** ptr, int n) {
    for (int i = 0; i < n; ++i) {
        free(ptr[i]);
    }
    free(ptr);
}

void print_array(std::string name, _real * array, int size, int format){
    std::cout << name << ":";
    for (int i=0; i < size; i++) {
        //std::cout << " " << array[i];
        if(format == 1)
            printf(" %.4e", array[i]);
        else
            printf(" %.4f", array[i]);
    }
    std::cout << std::endl;
}

void print_array_int(std::string name, int * array, int size){
    std::cout << name << ":";
    for (int i=0; i < size; i++) {
        printf(" %d", array[i]);
    }
    std::cout << std::endl;
}

_real normalize_angle(_real angle) {
    _real out_angle = fmod(angle,2*M_PI);
    if(out_angle < 0)
        out_angle += 2*M_PI;
    if(out_angle>M_PI)
        out_angle -= 2*M_PI;
    return out_angle;
}
void print_state_of_cost_function(_real * state, _real * xref, _real * xss, _real * control, int N, int Nx){
    printf("initial_state\n");
    for (int i = 0; i < Nx; ++i) {
        printf("%f ", state[i]);
    }
    printf("\nreference at steady state xss\n");
    for (int i = 0; i < Nx; ++i) {
        printf("%f ", xss[i]);
    }
    printf("\nxref0 xref1 xref2 xref3 u1\n");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < Nx; ++j) {
            printf("%f ", xref[i*Nx+j]);
        }
        printf("%f \n", control[i]);
        
    }
    printf("\n");
}

_real maximum(_real a, _real b){
    if(a > b)
        return a;
    else
        return b;
}
