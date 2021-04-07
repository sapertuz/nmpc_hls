#ifndef RANDOM_NUMBER_HPP
#define RANDOM_NUMBER_HPP

#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

class RandomNumber {

private:
    int count;
    int test;
    FILE * random_vector;

public:
	RandomNumber(int test);
    ~RandomNumber();
    _real read();
};

#endif