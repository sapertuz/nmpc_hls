#ifndef AUX_FUNC_HPP
#define AUX_FUNC_HPP

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

template<typename float_type> void print_formatted_float_array(const float_type * array, const unsigned short size, const unsigned short precision, const unsigned short width){
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

template<typename array_type_dest, typename array_type_src, unsigned n> void memcpy_loop_unrolled
    (array_type_dest dest[n], const array_type_src src[n]){
#pragma HLS inline
    // Copy contents of src[] to dest[]
    for (unsigned short i=0; i<n; i++){
#pragma HLS unroll
        dest[i] = src[i];
    }
}

template<typename array_type_dest, typename array_type_src, unsigned n> void memcpy_loop_rolled
    (volatile array_type_dest dest[n], volatile array_type_src src[n]){
#pragma HLS inline
    // Copy contents of src[] to dest[]
    memcpy_loop: for (unsigned short i=0; i<n; i++){
#pragma HLS pipeline off
        dest[i] = src[i];
    }
}

template<typename array_type_dest, typename array_type_src, unsigned n> void memcpy_loop_enclosed
    (array_type_dest dest[n], const array_type_src src[n]){
    // Copy contents of src[] to dest[]
    for (unsigned short i=0; i<n; i++){
#pragma HLS unroll
        dest[i] = src[i];
    }
}

template<typename data_type> void memset_loop
    (data_type *array, const data_type data, unsigned n){
    // Copy contents of src[] to dest[]
    for (unsigned short i=0; i<n; i++){
#pragma HLS unroll
        array[i] = data;
    }
}

template<typename array_type_dest, typename array_type_src, unsigned n> void split 
    (volatile array_type_src *in, array_type_dest *out1, array_type_dest *out2) {
#pragma HLS inline
        split: for(int i=0; i<n; i++) {
#pragma HLS pipeline off
            array_type_dest a = in[i];
            out1[i] = a;
            out2[i] = a;
        }
    }
    
template<typename array_type_dest, typename array_type_src, unsigned n> void split_4 
    (array_type_src in[n], 
    array_type_dest out1[n], array_type_dest out2[n], array_type_dest out3[n], array_type_dest out4[n]) {
#pragma HLS inline
        split_4: for(int i=0; i<n; i++) {
            array_type_dest a = in[i];
            out1[i] = a;
            out2[i] = a;
            out3[i] = a;
            out4[i] = a;
        }
    }

template<typename array_type_dest, typename array_type_src, unsigned n> void split_5 
    (array_type_src in[n], 
    array_type_dest out1[n], array_type_dest out2[n], array_type_dest out3[n], array_type_dest out4[n], array_type_dest out5[n]) {
#pragma HLS inline
        split_5 :for(int i=0; i<n; i++) {
            array_type_dest a = in[i];
            out1[i] = a;
            out2[i] = a;
            out3[i] = a;
            out4[i] = a;
            out5[i] = a;
        }
    }

#endif
