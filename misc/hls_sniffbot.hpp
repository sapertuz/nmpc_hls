#ifndef HLS_SNIFFBOT_HPP
#define HLS_SNIFFBOT_HPP

#include "fast_sin_cos.hpp"

#ifdef __SYNTHESIS__
#include "hls_math.h"
#include "ap_fixed.h"
using namespace hls;
//typedef ap_ufixed<32,8, AP_RND_ZERO, AP_WRAP_SM> _hw_real;
#endif


/*
#ifndef hls
namespace hls{
    float sinf(float angle){
        return sinf(angle);
    }
    float cosf(float angle){
        return cosf(angle);
    }
}
#endif
*/
template <
    typename _hw_real
>class model_sniffbot{

protected:
    const _hw_real u_max[4] =  {100, 100, 100, 100};
    const _hw_real u_min[4] =  {-100, -100, -100, -100};
    const _hw_real Ixx = 1.2;    // Moment of Inertia (Ixx Iyy Izz)[kg.m^2]
    const _hw_real Iyy = 1.2;
    const _hw_real Izz = 2.3;
    const _hw_real l  = .25;     // Arm length [m]
    const _hw_real kf = 1;       // Thrust (lift) coefficient [Ns^2]
    const _hw_real kM = 0.2;     // Moment (drag) coefficient [Nms^2]
    const _hw_real m  = 2;       // Mass of quadcopter [kg]
    const _hw_real g  = 9.81;    // Gravity [m/s^2]
    const _hw_real b  = 4.9050;  // Speed offset [rad^2/sec^2]4.9050

    const _hw_real M_mma[16] = { 
        1,  1,  1,  1,  // thrust (T)
        0, -1,  0,  1,  // roll  (tx)
        -1,  0,  1,  0,  // pitch (ty)
        -1,  1, -1,  1}; // yaw   (tz)
    const _hw_real M_mma_inv[16] = { 
        0.2500,         0,   -0.5000,   -0.2500,
        0.2500,   -0.5000,         0,    0.2500,
        0.2500,         0,    0.5000,   -0.2500,
        0.2500,    0.5000,         0,    0.2500};
    // Input limits
    const _hw_real phys_u_min[4] = {-15,  -3,  -3,  -3};
    const _hw_real phys_u_max[4] = {15,   3,   3,   3};

    const _hw_real local_pi = 3.14159265358979323846;
    const _hw_real local_pi_half = local_pi*.5;
#ifdef USE_FAST_SIN_COS
    typedef fast_sin_cos<_hw_real> fast_sin_cos_t;
    fast_sin_cos_t fast_sin_cos_p;
#endif

public:
    model_sniffbot(){}

    void model(
        _hw_real state_dot[12],
        _hw_real state[12],
        _hw_real control[4]
    ){
// #pragma HLS interface ap_fifo  port=state_dot
// #pragma HLS interface ap_fifo  port=state
// #pragma HLS interface ap_fifo  port=control

//#pragma HLS INLINE
//#pragma HLS pipeline II=4


#pragma HLS allocation operation instances=hmul limit=4
#pragma HLS allocation operation instances=hdiv limit=1
#pragma HLS allocation operation instances=hadd limit=2
#pragma HLS allocation operation instances=hsub limit=2

#ifdef __VITIS_HLS__
#pragma HLS allocation function instances=half_sincos limit=1
#else
#pragma HLS allocation function instances=sin_or_cos limit=1
#endif

#pragma HLS allocation function instances=local_sin limit=1
#pragma HLS allocation function instances=local_cos limit=1
// #pragma HLS allocation function instances=hls::sinf limit=1
// #pragma HLS allocation function instances=hls::cosf limit=1

        // _real _state_dot[_Nx];
        // _real _state[_Nx];
        // _real _control[_n_U];

        // for (short i=0; i<_Nx; i++){
        //     _state[i] = state[i];
        // }
        // for (short i=0; i<_n_U; i++){
        //     _control[i] = control[i];
        // }

        // int i;
        // int j;
        _hw_real w[4] = {0.0, 0.0, 0.0, 0.0};
        _hw_real param[4], u_param[4];
        //---------------------------------------------------
        _hw_real t2, t3, t4, t5, t6, t7, t8, t9, t12, t21, t22, t23, t24,t13;
        _hw_real t14, t16, t17, t18, t19, t20, t25, t15, t27, t10, t11, t26;
        _hw_real x, y, z, xdot, ydot, zdot, phi, phidot, psidot, psi, thetadot, theta;
        _hw_real u1, u2, u3, u4;
        //---------------------------------------------------
        _hw_real phys_control[4];
        for (short i = 0; i < 4; i++){
            phys_control[i] = value_map(control[i], u_min[i], u_max[i], phys_u_min[i], phys_u_max[i]);
        }

        param[0] = (_hw_real)1.0 / kf;
        t3 = (_hw_real)1.0 / (l * kf);
        param[1] = t3;
        param[2] = t3;
        param[3] = (_hw_real)1.0 / kM;
        for (unsigned i = 0; i < 4; i++) {
//            w[i] = (_hw_real)0.0;
#pragma HLS unroll
            u_param[i] = phys_control[i] * param[i];
        }

        for (unsigned i = 0; i < 4; i++) {
            for (unsigned j = 0; j < 4; j++) {
                w[i] += M_mma_inv[j + (i << 2)] * u_param[j];
            }
            w[i] += b;
        }

        x = state[0];
        y = state[1];
        z = state[2];
        phi = state[3];
        theta = state[4];
        psi = state[5];
        xdot = state[6];
        ydot = state[7];
        zdot = state[8];
        phidot = state[9];
        psidot = state[11];
        thetadot = state[10];

        u1 = w[0];
        u2 = w[1];
        u3 = w[2];
        u4 = w[3];
/*
        t2 = hls::cosf(phi);
        t4 = hls::cosf(theta);
        t3 = hls::cosf(psi);
        t5 = hls::sinf(phi);
        t7 = hls::sinf(theta);
        t6 = hls::sinf(psi);
*/
/*
        t2 = local_cos(phi);
        t4 = local_cos(theta);
        t3 = local_cos(psi);
*/
        t2 = local_sin(phi + (_hw_real)local_pi_half);
        t4 = local_sin(theta + (_hw_real)local_pi_half);
        t3 = local_sin(psi + (_hw_real)local_pi_half);
        t5 = local_sin(phi);        
        t7 = local_sin(theta);
        t6 = local_sin(psi);
        t8 = Iyy*Iyy;
        t9 = Izz*Izz;
        t10 = phi*(_hw_real)2.0;
        t11 = theta*(_hw_real)2.0;
/*
        t17 = hls::sinf(t10);
        t19 = hls::sinf(t11);
*/
        t17 = local_sin(t10);
        t19 = local_sin(t11);

        t12 = thetadot*thetadot;
        t21 = (_hw_real)1.0/Iyy;
        t22 = (_hw_real)1.0/Izz;
        t23 = (_hw_real)1.0/m;
        t24 = u1+u2+u3+u4;
        t13 = t2*t2;
        t14 = t2*t2*t2;
        t16 = t4*t4;
        t18 = t5*t5;
        t20 = t7*t7;
        t25 = (_hw_real)1.0/t4;
        t15 = t13*t13;
        t26 = t20-(_hw_real)1.0;
        t27 = (_hw_real)1.0/t26;

        state_dot[0] = xdot;
        state_dot[1] = ydot;
        state_dot[2] = zdot;
        state_dot[3] = phidot;
        state_dot[4] = thetadot;
        state_dot[5] = psidot;

        state_dot[6] = kf*t23*t24*(t5*t6+t2*t3*t7);

        state_dot[7] = -kf*t23*t24*(t3*t5-t2*t6*t7);

        state_dot[8] = -g+kf*t2*t4*t23*t24;

        state_dot[9] = -t7*t21*t22*t27*(Iyy-Iyy*t18+Izz*t18)*(-kM*u1+kM*u2-kM*u3+kM*u4+
            (Ixx*phidot*t4*thetadot)*(_hw_real)0.5+(Iyy*phidot*t4*thetadot)*(_hw_real)0.5-
            (Izz*phidot*t4*thetadot)*(_hw_real)0.5-(Ixx*psidot*t19*thetadot)*(_hw_real)0.5+
            (Iyy*psidot*t19*thetadot)*(_hw_real)0.5-Iyy*phidot*t4*t13*thetadot+
            Izz*phidot*t4*t13*thetadot+(Iyy*t2*t5*t7*t12)*(_hw_real)0.5-
            (Izz*t2*t5*t7*t12)*(_hw_real)0.5-Iyy*phidot*psidot*t2*t5*t16+
            Izz*phidot*psidot*t2*t5*t16-Iyy*psidot*t4*t7*t13*thetadot+
            Izz*psidot*t4*t7*t13*thetadot)-(t21*t22*t27*(-kf*l*u2+kf*l*u4+
            (Ixx*psidot*t4*thetadot)*(_hw_real)0.5)*(Iyy*Izz+Ixx*Iyy*t20-Iyy*Izz*t20-
            Ixx*Iyy*t18*t20+Ixx*Izz*t18*t20))/Ixx+t7*t17*t21*t22*t25*
            (Iyy*(_hw_real)0.5-Izz*(_hw_real)0.5)*(kf*l*u1-kf*l*u3-(Iyy*phidot*t17*thetadot)*(_hw_real)0.5+
            (Izz*phidot*t17*thetadot)*(_hw_real)0.5-(Iyy*phidot*psidot*t4)*(_hw_real)0.5+
            (Izz*phidot*psidot*t4)*(_hw_real)0.5+Iyy*phidot*psidot*t4*t13-
            Izz*phidot*psidot*t4*t13-(Iyy*psidot*t2*t5*t7*thetadot)*(_hw_real)0.5+
            (Izz*psidot*t2*t5*t7*thetadot)*(_hw_real)0.5);

        state_dot[10] = (t21*t22*t25*(phidot*psidot*t8*t16-t7*t8*t12*t13-t7*t9*t12*t13+
            t7*t8*t12*t15+t7*t9*t12*t15+Iyy*kM*t17*u1-Iyy*kM*t17*u2+Iyy*kM*t17*u3-
            Iyy*kM*t17*u4-Izz*kM*t17*u1+Izz*kM*t17*u2-Izz*kM*t17*u3+Izz*kM*t17*u4-
            Iyy*Izz*phidot*psidot*t16+Iyy*Izz*t7*t12*t13*(_hw_real)2.0-
            Iyy*Izz*t7*t12*t15*(_hw_real)2.0-Iyy*kf*l*t4*u1*(_hw_real)2.0+Iyy*kf*l*t4*u3*(_hw_real)2.0-
            phidot*psidot*t8*t13*t16+phidot*psidot*t9*t13*t16+
            phidot*t2*t4*t5*t8*thetadot-phidot*t2*t4*t5*t9*thetadot+
            Iyy*kf*l*t4*t13*u1*(_hw_real)2.0-Iyy*kf*l*t4*t13*u3*(_hw_real)2.0-Izz*kf*l*t4*t13*u1*(_hw_real)2.0+
            Izz*kf*l*t4*t13*u3*(_hw_real)2.0-Ixx*Iyy*phidot*t2*t4*t5*thetadot+
            Ixx*Izz*phidot*t2*t4*t5*thetadot+Iyy*kf*l*t2*t5*t7*u2*(_hw_real)2.0-
            Iyy*kf*l*t2*t5*t7*u4*(_hw_real)2.0-Izz*kf*l*t2*t5*t7*u2*(_hw_real)2.0+
            Izz*kf*l*t2*t5*t7*u4*(_hw_real)2.0-psidot*t2*t4*t5*t7*t8*thetadot+
            psidot*t4*t5*t7*t8*t14*thetadot+psidot*t4*t5*t7*t9*t14*thetadot+
            Ixx*Iyy*psidot*t2*t4*t5*t7*thetadot-Ixx*Izz*psidot*t2*t4*t5*t7*thetadot
            +Iyy*Izz*psidot*t2*t4*t5*t7*thetadot-
            Iyy*Izz*psidot*t4*t5*t7*t14*thetadot*(_hw_real)2.0))*(_hw_real)0.5;

        state_dot[11] = (t21*t22*(Izz*kM*u1*-(_hw_real)2.0+Izz*kM*u2*(_hw_real)2.0-Izz*kM*u3*(_hw_real)2.0+Izz*kM*u4*(_hw_real)2.0-
            phidot*t4*t9*thetadot-Iyy*kM*t13*u1*(_hw_real)2.0+Iyy*kM*t13*u2*(_hw_real)2.0-
            Iyy*kM*t13*u3*(_hw_real)2.0+Iyy*kM*t13*u4*(_hw_real)2.0+Izz*kM*t13*u1*(_hw_real)2.0-Izz*kM*t13*u2*(_hw_real)2.0+
            Izz*kM*t13*u3*(_hw_real)2.0-Izz*kM*t13*u4*(_hw_real)2.0+Ixx*Izz*phidot*t4*thetadot+
            Iyy*Izz*phidot*t4*thetadot-(Ixx*Izz*psidot*t19*thetadot)*(_hw_real)0.5+
            Iyy*Izz*psidot*t19*thetadot-Izz*kf*l*t7*u2*(_hw_real)2.0+Izz*kf*l*t7*u4*(_hw_real)2.0-
            phidot*t4*t8*t13*thetadot+phidot*t4*t9*t13*thetadot-t2*t5*t7*t9*t12+
            t5*t7*t8*t12*t14+t5*t7*t9*t12*t14-phidot*psidot*t2*t5*t8*t16+
            phidot*psidot*t2*t5*t9*t16+psidot*t4*t7*t8*t13*thetadot+
            psidot*t4*t7*t9*t13*thetadot-psidot*t4*t7*t8*t15*thetadot-
            psidot*t4*t7*t9*t15*thetadot+Ixx*Iyy*phidot*t4*t13*thetadot-
            Ixx*Izz*phidot*t4*t13*thetadot+Iyy*Izz*t2*t5*t7*t12-
            Iyy*Izz*t5*t7*t12*t14*(_hw_real)2.0-Iyy*kf*l*t7*t13*u2*(_hw_real)2.0+Iyy*kf*l*t7*t13*u4*(_hw_real)2.0+
            Izz*kf*l*t7*t13*u2*(_hw_real)2.0-Izz*kf*l*t7*t13*u4*(_hw_real)2.0-
            Ixx*Iyy*psidot*t4*t7*t13*thetadot+Ixx*Izz*psidot*t4*t7*t13*thetadot-
            Iyy*Izz*psidot*t4*t7*t13*thetadot*(_hw_real)2.0+
            Iyy*Izz*psidot*t4*t7*t15*thetadot*(_hw_real)2.0+Iyy*kf*l*t2*t4*t5*u1*(_hw_real)2.0-
            Iyy*kf*l*t2*t4*t5*u3*(_hw_real)2.0-Izz*kf*l*t2*t4*t5*u1*(_hw_real)2.0+
            Izz*kf*l*t2*t4*t5*u3*(_hw_real)2.0))/(t16*(_hw_real)2.0);

        // for (short i=0; i<_Nx; i++){
        //     state_dot[i] = _state_dot[i];
        // }
        
    }
    
private:
    _hw_real value_map(_hw_real input, _hw_real in_min, _hw_real in_max, _hw_real out_min, _hw_real out_max){
#pragma HLS inline
    _hw_real output = (input - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
        return output;
    }

    _hw_real local_sin(_hw_real angle){
//#pragma HLS inline
#if defined(__VITIS_HLS__) && defined(__SYNTHESIS__)
#ifndef USE_FAST_SIN_COS
        return half_sin(angle);
#else
        return fast_sin_cos_p.fastsin(angle);
#endif
#else
        return sin(angle);
#endif
    }
    _hw_real local_cos(_hw_real angle){
//#pragma HLS inline
#if defined(__VITIS_HLS__) && defined(__SYNTHESIS__)
#ifndef USE_FAST_SIN_COS
        return half_cos(angle);
#else
        return fast_sin_cos_p.fastcos(angle);
#endif    
#else
        return cos(angle);
#endif
    }

};
#endif
