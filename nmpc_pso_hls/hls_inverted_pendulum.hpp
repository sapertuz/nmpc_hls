#ifndef HLS_INVERTED_PENDULUM_HPP
#define HLS_INVERTED_PENDULUM_HPP

#ifdef __SYNTHESIS__
#include "hls_math.h"
using namespace hls;
#endif

template <
    typename _hw_real
>class model_inverted_pendulum{
public:
    model_inverted_pendulum(){}

    void model(
        _hw_real *state_dot,
        _hw_real *state,
        _hw_real *control
    ){
// #pragma HLS interface ap_fifo  port=state_dot
// #pragma HLS interface ap_fifo  port=state
// #pragma HLS interface ap_fifo  port=control

#pragma HLS allocation operation instances=hmul limit=2
#pragma HLS allocation operation instances=hdiv limit=1
#pragma HLS allocation operation instances=hadd limit=1
#pragma HLS allocation operation instances=hsub limit=1

        // Inverted Pendulum from:
        // Mercieca, J., & Fabri, S. G. (2011). 
        // Particle swarm optimisation for non-linear model predictive control. 
        // International Conference on Advanced Engineering Computing and Applications in Science, 5(1), 88?93.
        //_hw_real m = 7.3;    // [Kg] Uniformily distributed mass of the pendulum
        //_hw_real M = 14.6;   // [Kg] Cart mass
        _hw_real g = 9.81;   // [m/s^2] Gravity
        //_hw_real l = 1.2;    // [m] half the lenght of the pendulum
        _hw_real b = 14.6;   // [Kg/s] surface friction
        //_hw_real h = 0.0136;  // [Kg.m^2/s] rotation friction damping coeficient

        _hw_real u = control[0];
        _hw_real x0 = state[0];
        _hw_real x1 = state[1];
        _hw_real th = state[2];
        _hw_real th1 = state[3];
        
#ifdef DEBUG_SYSTEM
    // std::cout << "\t states: "  << std::fixed << std::setprecision(2) << std::right << std::setw(6) << x0 << "\t";
    // std::cout                   << std::fixed << std::setprecision(2) << std::right << std::setw(6) << x1 << "\t";
    // std::cout                   << std::fixed << std::setprecision(2) << std::right << std::setw(6) << th << "\t";
    // std::cout                   << std::fixed << std::setprecision(2) << std::right << std::setw(6) << th1 << std::endl;
    // std::cout << "states: " << x0 << " " << x1 << " " << th << " " << th1 << std::endl;
#endif
        _hw_real c = cos(th);
        _hw_real s = sin(th);
        //_hw_real c_1 = 1.0/c;

        
        _hw_real ml = 8.76; //m*l;
        _hw_real Mm = 21.9; //M+m;
        _hw_real Mn_1 = 0.04566210046;
        //_hw_real Mml43 = 35.04; //Mm*l*4/3;
        _hw_real h_ml = 0.001552511416; // h/ml
        _hw_real l43 = 1.6; // l = 1.2 * 4/3

        _hw_real Mm_c = Mm/c;

        _hw_real mls = ml*s;
        _hw_real mlc = ml*c;
        _hw_real u_bx1 = u-b*x1;
        _hw_real state4_2 = th1*th1; // pow(state[3],2);

        // xdot(1) = x(2);                                                 % x dot
        // xdot(3) = x(4);                                                 % theta dot
        // xdot(4) = (u-b*x(2)+mls*x(4)^2-Mm/c*(g*s-h*x(4)/ml)) / (ml*c-(Mm)*4*l/(3*c)); % theta dot dot
        // xdot(2) = (1/(Mm))*(u-b*x(2)-ml*xdot(4)*c+mls*x(4)^2); % x dot dot

        _hw_real state_dot_3 = (u_bx1+mls*state4_2-Mm_c*(g*s-(h_ml)*th1)) / (mlc-Mm_c*l43);

        state_dot[0] = x1;                                                 // x dot
        state_dot[1] = (Mn_1)*(u_bx1-mlc*state_dot_3+mls*state4_2); // x dot dot
        state_dot[2] = th1;                                                 // theta dot
        //state_dot[3] = (u_bx1+mls*state4_2-Mm*c_1*(g*s-(h_ml)*state[3])) / (mlc-Mml43*c_1); // theta dot dot
        state_dot[3] = state_dot_3;// theta dot dot
        
    }
};
#endif
