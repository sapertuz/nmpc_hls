#include "random_number.hpp"

RandomNumber::RandomNumber(int test) {
	count = 0;
	this->test = test;
	if(test) {
		random_vector = fopen("rand_vector.txt", "r");
	    if (random_vector == NULL) {
		    throw(1);
		}
	}
}

RandomNumber::~RandomNumber() {
	if(test) {
		fclose(random_vector);
	}
}

_real RandomNumber::read(){
	if(test) {
		double random_number;
		fscanf(random_vector, "%lf\n", &random_number);
		return (_real) random_number;
	}
	else
		return (_real)rand()/(_real)RAND_MAX;
}