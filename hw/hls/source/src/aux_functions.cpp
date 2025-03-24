#include "aux_functions.hpp"
#include <iostream>
#include <iomanip>

void print_formatted_float_array(const float * array, const unsigned short size, const unsigned short precision, const unsigned short width){
    for (unsigned short i = 0; i < size; i++){
        std::cout << std::fixed << std::setprecision(precision) << std::right << std::setw(width) << array[i];
        if (i<size-1){
            std::cout << "\t";
        }
    }
}
float ** alloc_matrix(int n, int m){
    float ** matrix;
    if ((matrix = (float **) malloc(n*sizeof(float *))) == NULL) {
        std::cout << "Erro de alocação" << std::endl;
        return NULL;
    }
    else {
        for (int i = 0; i < n; i++) {
            if((matrix[i] = (float *) calloc(m, sizeof(float))) == NULL) {
                std::cout << "Erro de alocação 2" << std::endl;
                return NULL;
            }
        }
    }
    return matrix;
}
void free_matrix(float ** ptr, int n) {
    for (int i = 0; i < n; ++i) {
        free(ptr[i]);
    }
    free(ptr);
}
