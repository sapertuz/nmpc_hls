#include <stdio.h>
#include <iostream>
#include <cstring>

template <
    typename _hw_real,
    unsigned _Nx
>class test_class {
protected:
    _hw_real current_state[_Nx];
public:
    explicit test_class(const void *_current_state){
        std::memcpy(current_state, (const _hw_real*)_current_state, _Nx*sizeof(_hw_real));
    }
    _hw_real get_current_state(unsigned pos){
        return (current_state[pos]);
    }
};

int main(int argc, char const *argv[])
{
    const float test_array[4] = {1, 2.3, 3.4, 5};
    test_class<float,4> test_template_class(test_array);
    unsigned pos = 1;
    std::cout << "content[" << pos << "]= " << test_template_class.get_current_state(pos) << std::endl; 
    return 0;
}
