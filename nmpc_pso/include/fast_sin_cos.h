// Fast Sin and Cosine Library
#ifndef FAST_SIN_COS_HPP
#define FAST_SIN_COS_HPP

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#define MAX_CIRCLE_ANGLE 2048     //512
#define HALF_MAX_CIRCLE_ANGLE (MAX_CIRCLE_ANGLE/2)
#define QUARTER_MAX_CIRCLE_ANGLE (MAX_CIRCLE_ANGLE/4)
#define MASK_MAX_CIRCLE_ANGLE (MAX_CIRCLE_ANGLE - 1)
#define PI 3.14159265358979323846f

void initSinCosTable();
float fastcos(float n);
float fastsin(float n);

#endif