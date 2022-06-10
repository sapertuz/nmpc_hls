#ifndef AUX_FUNCTIONS_HPP
#define AUX_FUNCTIONS_HPP

#include <string>
#include "config.hpp"

_real ** alloc_matrix(int n, int m);
void free_matrix(_real ** ptr, int n);
_real min_array(_real * array, int * pos, int size);
_real rand_real();
void print_array(std::string name, _real * array, int size, int format);
void print_array_int(std::string name, int * array, int size);
_real normalize_angle(_real angle);
void print_state_of_cost_function(_real * state, _real * xref, _real * xss, _real * control, int N, int Nx);
_real maximum(_real a, _real b);

void memcpy_loop_rolled(_real *dest, _real *src, unsigned n);
void memset_loop(_real *array, const _real data, unsigned n);
void split (volatile _real *in, _real *out1, _real *out2, unsigned n);

#endif
