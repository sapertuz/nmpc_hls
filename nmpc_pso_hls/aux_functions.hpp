#ifndef AUX_FUNC_HPP
#define AUX_FUNC_HPP

void print_formatted_float_array(const float * array, const unsigned short size, const unsigned short precision, const unsigned short width){
    for (unsigned short i = 0; i < size; i++){
        std::cout << std::fixed << std::setprecision(precision) << std::right << std::setw(width) << array[i] << "\t";
    }
}

template<typename array_type_dest, typename array_type_src, unsigned n> void memcpy_loop
    (array_type_dest dest[n], const array_type_src src[n]){
#pragma HLS inline
    // Copy contents of src[] to dest[]
    for (unsigned short i=0; i<n; i++){
#pragma HLS unroll
        dest[i] = src[i];
    }
}

template<typename array_type_dest, typename array_type_src, unsigned n> void memcpy_loop_rolled
    (array_type_dest dest[n], const array_type_src src[n]){
#pragma HLS inline
    // Copy contents of src[] to dest[]
    for (unsigned short i=0; i<n; i++){
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

#endif
