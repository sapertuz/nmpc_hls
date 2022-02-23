#ifndef AUX_FUNC_HPP
#define AUX_FUNC_HPP

void print_formatted_float_array(const float * array, const unsigned short size, const unsigned short precision, const unsigned short width){
    for (unsigned short i = 0; i < size; i++){
        std::cout << std::fixed << std::setprecision(precision) << std::right << std::setw(width) << array[i] << "\t";
    }
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

template<typename array_type_dest, typename array_type_src, unsigned n> void memcpy_loop_feedback
    (array_type_dest dest[n], array_type_src src[n]){
    // Copy contents of src[] to dest[]
#pragma HLS inline
    for (unsigned short i=0; i<n; i++){
        array_type_dest a = src[i]; 
        dest[i] = a;
        array_type_src b = a;
        src[i] = b;
    }
}

template<typename array_type_dest, typename array_type_src, unsigned n> void memcpy_loop_rolled_2dest
    (array_type_dest dest_1[n], array_type_dest dest_2[n], const array_type_src src[n]){
    // Copy contents of src[] to dest[]
#pragma HLS inline
    for (unsigned short i=0; i<n; i++){
        array_type_dest a = src[i]; 
        dest_1[i] = a;
        dest_2[i] = a;
    }
}

#endif
