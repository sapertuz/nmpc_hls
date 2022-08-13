#ifndef HLS_INVERTED_PENDULUM_HPP
#define HLS_INVERTED_PENDULUM_HPP

#include "hls_arith_aux_functions.hpp"
#ifdef __SYNTHESIS__
#include "hls_math.h"
#include "ap_fixed.h"
using namespace hls;
//typedef ap_ufixed<32,8, AP_RND_ZERO, AP_WRAP_SM> _model_real;
#else
#include "math.h"
#endif

template <
    typename _model_real
>void __model(
    volatile _model_real *state_dot,
    volatile _model_real *state,
    volatile _model_real *control
){
#pragma HLS interface mode=ap_fifo port=state_dot   depth=12
#pragma HLS interface mode=ap_fifo port=state       depth=12
#pragma HLS interface mode=ap_fifo port=control     depth=4

#pragma HLS allocation operation instances=hmul limit=2
#pragma HLS allocation operation instances=hdiv limit=1
#pragma HLS allocation operation instances=hadd limit=1
#pragma HLS allocation operation instances=hsub limit=1

#pragma HLS BIND_OP variable=state_dot op=hadd impl=fulldsp latency=4
#pragma HLS BIND_OP variable=state_dot op=hsub impl=fulldsp latency=4
#pragma HLS BIND_OP variable=state_dot op=hmul impl=fulldsp latency=4

    // Inverted Pendulum from:
    // Mercieca, J., & Fabri, S. G. (2011). 
    // Particle swarm optimisation for non-linear model predictive control. 
    // International Conference on Advanced Engineering Computing and Applications in Science, 5(1), 88?93.
    
#ifdef DEBUG_SYSTEM
// std::cout << "\t states: "  << std::fixed << std::setprecision(2) << std::right << std::setw(6) << x0 << "\t";
// std::cout                   << std::fixed << std::setprecision(2) << std::right << std::setw(6) << x1 << "\t";
// std::cout                   << std::fixed << std::setprecision(2) << std::right << std::setw(6) << th << "\t";
// std::cout                   << std::fixed << std::setprecision(2) << std::right << std::setw(6) << th1 << std::endl;
// std::cout << "states: " << x0 << " " << x1 << " " << th << " " << th1 << std::endl;
#endif
    // Inverted Pendulum from:
	// Mercieca, J., & Fabri, S. G. (2011). 
	// Particle swarm optimisation for non-linear model predictive control. 
	// International Conference on Advanced Engineering Computing and Applications in Science, 5(1), 88?93.

	const _model_real m = 7.3;    // [Kg] Uniformily distributed mass of the pendulum
	const _model_real M = 14.6;   // [Kg] Cart mass
	const _model_real g = 9.81;   // [m/s^2] Gravity
	const _model_real l = 1.2;    // [m] half the lenght of the pendulum
	const _model_real b = 14.6;   // [Kg/s] surface friction
	const _model_real h = 0.0136;  // [Kg.m^2/s] rotation friction damping coeficient

	
	const _model_real ml = m*l;
	const _model_real Mm = M+m;
    const _model_real Mn_1 = 1/(M+m);//0.04566210046;
	//_model_real Mml43 = 35.04; //Mm*l*4/3;
	const _model_real h_ml = h/(m*l);
    const _model_real l43 = l * (_model_real)4.0/(_model_real)3.0;

    _model_real th = state[2];
    _model_real c = cos(th);
    _model_real s = sin(th);
    //_model_real c_1 = 1.0/c;

    _model_real Mm_c = Mm/c;

	_model_real mls = ml*s;
	_model_real mlc = ml*c;
	_model_real u_bx1 = control[0]-b*state[1];
    _model_real state4_2 = state[3]*state[3]; // pow(state[3],2);

    // xdot(1) = x(2);                                                 % x dot
    // xdot(3) = x(4);                                                 % theta dot
    // xdot(4) = (u-b*x(2)+mls*x(4)^2-Mm/c*(g*s-h*x(4)/ml)) / (ml*c-(Mm)*4*l/(3*c)); % theta dot dot
    // xdot(2) = (1/(Mm))*(u-b*x(2)-ml*xdot(4)*c+mls*x(4)^2); % x dot dot


	state_dot[0] = state[1];                                                 // x dot
	state_dot[2] = state[3];                                                 // theta dot
	//state_dot[3] = (u_bx1+mls*state4_2-Mm*c_1*(g*s-(h_ml)*state[3])) / (mlc-Mml43*c_1); // theta dot dot
    state_dot[3] = (u_bx1+mls*state4_2-Mm_c*(g*s-(h_ml)*state[3])) / (mlc-Mm_c*l43); // theta dot dot
	state_dot[1] = (Mn_1)*(u_bx1-mlc*state_dot[3]+mls*state4_2); // x dot dot

}
#endif
