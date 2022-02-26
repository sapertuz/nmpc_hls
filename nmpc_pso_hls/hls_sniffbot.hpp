#ifndef HLS_SNIFFBOT_HPP
#define HLS_SNIFFBOT_HPP

#ifdef __SYNTHESIS__
#include "fast_sin_cos.hpp"
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

public:
    model_sniffbot(){}

    void model(
        _hw_real *state_dot, // [12]
        volatile _hw_real *state,     // [12]
        volatile _hw_real *control    // [4]
    ){
// #pragma HLS interface ap_fifo  port=state_dot
// #pragma HLS interface ap_fifo  port=state
// #pragma HLS interface ap_fifo  port=control
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
    const _hw_real factor[4] = {0.15, 0.03, 0.03, 0.03}; // factor = (out_max[i] - out_min[i]) / (in_max[i] - in_min[i])
    
    const _hw_real local_pi = 3.14159265358979323846;
    const _hw_real local_pi_half = 1.57079632679;

#pragma HLS INLINE
//#pragma HLS pipeline II=4
//#pragma HLS expression_balance on

#pragma HLS allocation operation instances=hmul limit=4
#pragma HLS allocation operation instances=hdiv limit=1
#pragma HLS allocation operation instances=hadd limit=2
#pragma HLS allocation operation instances=hsub limit=2

#pragma HLS allocation operation instances=mul limit=4
#pragma HLS allocation operation instances=div limit=1
#pragma HLS allocation operation instances=add limit=2
#pragma HLS allocation operation instances=sub limit=2

#ifdef __VITIS_HLS__
#pragma HLS allocation function instances=half_sincos limit=1
#else
#pragma HLS allocation function instances=sin_or_cos limit=1
#endif

#ifndef USE_FAST_SIN_COS
#pragma HLS allocation function instances=local_sin limit=1
#pragma HLS allocation function instances=local_cos limit=1
#else
#pragma HLS allocation function instances=fastsin<_hw_real> limit=1
#endif
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
        _hw_real w[4];
        _hw_real param[4], u_param[4];
        //---------------------------------------------------
        _hw_real t[28];
#ifdef __VITIS_HLS__
//#pragma HLS bind_storage variable=t type=RAM_2P impl=AUTO
#endif
        // _hw_real t2, t3, t[4], t5, t6, t7, t8, t9, t[12], t21, t22, t23, t24,t[13];
        // _hw_real t[14], t[16], t17, t18, t19, t20, t25, t[15], t27, t10, t11, t26;
        _hw_real x, y, z, xdot, ydot, zdot, phi, phidot, psidot, psi, thetadot, theta;
        _hw_real u1, u2, u3, u4;
        //---------------------------------------------------
        _hw_real phys_control[4];
        map_loop: for (short i = 0; i < 4; i++){
#pragma HLS pipeline off
            phys_control[i] = (control[i] - u_min[i]) * factor[i] + phys_u_min[i];
//            phys_control[i] = value_map(control[i], u_min[i], u_max[i], phys_u_min[i], phys_u_max[i]);
        }

        param[0] = (_hw_real)1.0 / kf;
        t[3] = (_hw_real)1.0 / (l * kf);
        param[1] = t[3];
        param[2] = t[3];
        param[3] = (_hw_real)1.0 / kM;
        u_param_loop: for (unsigned i = 0; i < 4; i++) {
//            w[i] = (_hw_real)0.0;
#pragma HLS unroll
            u_param[i] = phys_control[i] * param[i];
        }

        w_loop: for (unsigned i = 0; i < 4; i++) {
#pragma HLS pipeline off
            w[i] = b;
            for (unsigned j = 0; j < 4; j++) {
#pragma HLS pipeline off
                w[i] += M_mma_inv[j + (i << 2)] * u_param[j];
            }
        }

        for (unsigned i = 0; i < 12; i++){
#pragma HLS unroll factor=1
#pragma HLS pipeline off 
              if (i==0){ x = state[0];
        }else if (i==1){ y = state[1];
        }else if (i==2){ z = state[2];
        }else if (i==3){ phi = state[3];
        }else if (i==4){ theta = state[4];
        }else if (i==5){ psi = state[5];
        }else if (i==6){ xdot = state[6];
        }else if (i==7){ ydot = state[7];
        }else if (i==8){ zdot = state[8];
        }else if (i==9){ phidot = state[9];
        }else if (i==10){ psidot = state[10];
        }else if (i==11){ thetadot = state[11];}
        }       

        u1 = w[0];
        u2 = w[1];
        u3 = w[2];
        u4 = w[3];

        t[2] = local_sin(phi + (_hw_real)local_pi_half);
        t[4] = local_sin(theta + (_hw_real)local_pi_half);
        t[3] = local_sin(psi + (_hw_real)local_pi_half);
        t[5] = local_sin(phi);        
        t[7] = local_sin(theta);
        t[6] = local_sin(psi);
        t[8] = Iyy*Iyy;
        t[9] = Izz*Izz;
        t[10] = phi*(_hw_real)2.0;
        t[11] = theta*(_hw_real)2.0;

        t[17] = local_sin(t[10]);
        t[19] = local_sin(t[11]);

        t[12] = thetadot*thetadot;
        t[21] = (_hw_real)1.0/Iyy;
        t[22] = (_hw_real)1.0/Izz;
        t[23] = (_hw_real)1.0/m;
        t[24] = u1+u2+u3+u4;
        t[13] = t[2]*t[2];
        t[14] = t[2]*t[2]*t[2];
        t[16] = t[4]*t[4];
        t[18] = t[5]*t[5];
        t[20] = t[7]*t[7];
        t[25] = (_hw_real)1.0/t[4];
        t[15] = t[13]*t[13];
        t[26] = t[20]-(_hw_real)1.0;
        t[27] = (_hw_real)1.0/t[26];

        for (unsigned i = 0; i < 12; i++){
#pragma HLS unroll factor=1
#pragma HLS pipeline off 
        if (i==0){
        state_dot[0] = xdot;
        }else if (i==1){
        state_dot[1] = ydot;

        }else if (i==2){
        state_dot[2] = zdot;
        
        }else if (i==3){
        state_dot[3] = phidot;
        
        }else if (i==4){
        state_dot[4] = thetadot;
        
        }else if (i==5){
        state_dot[5] = psidot;
        
        }else if (i==6){
        state_dot[6] = kf*t[23]*t[24]*(t[5]*t[6]+t[2]*t[3]*t[7]);
        
        }else if (i==7){
        state_dot[7] = -kf*t[23]*t[24]*(t[3]*t[5]-t[2]*t[6]*t[7]);

        }else if (i==8){
        state_dot[8] = -g+kf*t[2]*t[4]*t[23]*t[24];
        
        }else if (i==9){
        state_dot[9] = -t[7]*t[21]*t[22]*t[27]*(Iyy-Iyy*t[18]+Izz*t[18])*(-kM*u1+kM*u2-kM*u3+kM*u4+
            (Ixx*phidot*t[4]*thetadot)*(_hw_real)0.5+(Iyy*phidot*t[4]*thetadot)*(_hw_real)0.5-
            (Izz*phidot*t[4]*thetadot)*(_hw_real)0.5-(Ixx*psidot*t[19]*thetadot)*(_hw_real)0.5+
            (Iyy*psidot*t[19]*thetadot)*(_hw_real)0.5-Iyy*phidot*t[4]*t[13]*thetadot+
            Izz*phidot*t[4]*t[13]*thetadot+(Iyy*t[2]*t[5]*t[7]*t[12])*(_hw_real)0.5-
            (Izz*t[2]*t[5]*t[7]*t[12])*(_hw_real)0.5-Iyy*phidot*psidot*t[2]*t[5]*t[16]+
            Izz*phidot*psidot*t[2]*t[5]*t[16]-Iyy*psidot*t[4]*t[7]*t[13]*thetadot+
            Izz*psidot*t[4]*t[7]*t[13]*thetadot)-(t[21]*t[22]*t[27]*(-kf*l*u2+kf*l*u4+
            (Ixx*psidot*t[4]*thetadot)*(_hw_real)0.5)*(Iyy*Izz+Ixx*Iyy*t[20]-Iyy*Izz*t[20]-
            Ixx*Iyy*t[18]*t[20]+Ixx*Izz*t[18]*t[20]))/Ixx+t[7]*t[17]*t[21]*t[22]*t[25]*
            (Iyy*(_hw_real)0.5-Izz*(_hw_real)0.5)*(kf*l*u1-kf*l*u3-(Iyy*phidot*t[17]*thetadot)*(_hw_real)0.5+
            (Izz*phidot*t[17]*thetadot)*(_hw_real)0.5-(Iyy*phidot*psidot*t[4])*(_hw_real)0.5+
            (Izz*phidot*psidot*t[4])*(_hw_real)0.5+Iyy*phidot*psidot*t[4]*t[13]-
            Izz*phidot*psidot*t[4]*t[13]-(Iyy*psidot*t[2]*t[5]*t[7]*thetadot)*(_hw_real)0.5+
            (Izz*psidot*t[2]*t[5]*t[7]*thetadot)*(_hw_real)0.5);

        
        }else if (i==10){
        state_dot[10] = (t[21]*t[22]*t[25]*(phidot*psidot*t[8]*t[16]-t[7]*t[8]*t[12]*t[13]-t[7]*t[9]*t[12]*t[13]+
            t[7]*t[8]*t[12]*t[15]+t[7]*t[9]*t[12]*t[15]+Iyy*kM*t[17]*u1-Iyy*kM*t[17]*u2+Iyy*kM*t[17]*u3-
            Iyy*kM*t[17]*u4-Izz*kM*t[17]*u1+Izz*kM*t[17]*u2-Izz*kM*t[17]*u3+Izz*kM*t[17]*u4-
            Iyy*Izz*phidot*psidot*t[16]+Iyy*Izz*t[7]*t[12]*t[13]*(_hw_real)2.0-
            Iyy*Izz*t[7]*t[12]*t[15]*(_hw_real)2.0-Iyy*kf*l*t[4]*u1*(_hw_real)2.0+Iyy*kf*l*t[4]*u3*(_hw_real)2.0-
            phidot*psidot*t[8]*t[13]*t[16]+phidot*psidot*t[9]*t[13]*t[16]+
            phidot*t[2]*t[4]*t[5]*t[8]*thetadot-phidot*t[2]*t[4]*t[5]*t[9]*thetadot+
            Iyy*kf*l*t[4]*t[13]*u1*(_hw_real)2.0-Iyy*kf*l*t[4]*t[13]*u3*(_hw_real)2.0-Izz*kf*l*t[4]*t[13]*u1*(_hw_real)2.0+
            Izz*kf*l*t[4]*t[13]*u3*(_hw_real)2.0-Ixx*Iyy*phidot*t[2]*t[4]*t[5]*thetadot+
            Ixx*Izz*phidot*t[2]*t[4]*t[5]*thetadot+Iyy*kf*l*t[2]*t[5]*t[7]*u2*(_hw_real)2.0-
            Iyy*kf*l*t[2]*t[5]*t[7]*u4*(_hw_real)2.0-Izz*kf*l*t[2]*t[5]*t[7]*u2*(_hw_real)2.0+
            Izz*kf*l*t[2]*t[5]*t[7]*u4*(_hw_real)2.0-psidot*t[2]*t[4]*t[5]*t[7]*t[8]*thetadot+
            psidot*t[4]*t[5]*t[7]*t[8]*t[14]*thetadot+psidot*t[4]*t[5]*t[7]*t[9]*t[14]*thetadot+
            Ixx*Iyy*psidot*t[2]*t[4]*t[5]*t[7]*thetadot-Ixx*Izz*psidot*t[2]*t[4]*t[5]*t[7]*thetadot
            +Iyy*Izz*psidot*t[2]*t[4]*t[5]*t[7]*thetadot-
            Iyy*Izz*psidot*t[4]*t[5]*t[7]*t[14]*thetadot*(_hw_real)2.0))*(_hw_real)0.5;
        
        }else if (i==11){
        state_dot[11] = (t[21]*t[22]*(Izz*kM*u1*-(_hw_real)2.0+Izz*kM*u2*(_hw_real)2.0-Izz*kM*u3*(_hw_real)2.0+Izz*kM*u4*(_hw_real)2.0-
            phidot*t[4]*t[9]*thetadot-Iyy*kM*t[13]*u1*(_hw_real)2.0+Iyy*kM*t[13]*u2*(_hw_real)2.0-
            Iyy*kM*t[13]*u3*(_hw_real)2.0+Iyy*kM*t[13]*u4*(_hw_real)2.0+Izz*kM*t[13]*u1*(_hw_real)2.0-Izz*kM*t[13]*u2*(_hw_real)2.0+
            Izz*kM*t[13]*u3*(_hw_real)2.0-Izz*kM*t[13]*u4*(_hw_real)2.0+Ixx*Izz*phidot*t[4]*thetadot+
            Iyy*Izz*phidot*t[4]*thetadot-(Ixx*Izz*psidot*t[19]*thetadot)*(_hw_real)0.5+
            Iyy*Izz*psidot*t[19]*thetadot-Izz*kf*l*t[7]*u2*(_hw_real)2.0+Izz*kf*l*t[7]*u4*(_hw_real)2.0-
            phidot*t[4]*t[8]*t[13]*thetadot+phidot*t[4]*t[9]*t[13]*thetadot-t[2]*t[5]*t[7]*t[9]*t[12]+
            t[5]*t[7]*t[8]*t[12]*t[14]+t[5]*t[7]*t[9]*t[12]*t[14]-phidot*psidot*t[2]*t[5]*t[8]*t[16]+
            phidot*psidot*t[2]*t[5]*t[9]*t[16]+psidot*t[4]*t[7]*t[8]*t[13]*thetadot+
            psidot*t[4]*t[7]*t[9]*t[13]*thetadot-psidot*t[4]*t[7]*t[8]*t[15]*thetadot-
            psidot*t[4]*t[7]*t[9]*t[15]*thetadot+Ixx*Iyy*phidot*t[4]*t[13]*thetadot-
            Ixx*Izz*phidot*t[4]*t[13]*thetadot+Iyy*Izz*t[2]*t[5]*t[7]*t[12]-
            Iyy*Izz*t[5]*t[7]*t[12]*t[14]*(_hw_real)2.0-Iyy*kf*l*t[7]*t[13]*u2*(_hw_real)2.0+Iyy*kf*l*t[7]*t[13]*u4*(_hw_real)2.0+
            Izz*kf*l*t[7]*t[13]*u2*(_hw_real)2.0-Izz*kf*l*t[7]*t[13]*u4*(_hw_real)2.0-
            Ixx*Iyy*psidot*t[4]*t[7]*t[13]*thetadot+Ixx*Izz*psidot*t[4]*t[7]*t[13]*thetadot-
            Iyy*Izz*psidot*t[4]*t[7]*t[13]*thetadot*(_hw_real)2.0+
            Iyy*Izz*psidot*t[4]*t[7]*t[15]*thetadot*(_hw_real)2.0+Iyy*kf*l*t[2]*t[4]*t[5]*u1*(_hw_real)2.0-
            Iyy*kf*l*t[2]*t[4]*t[5]*u3*(_hw_real)2.0-Izz*kf*l*t[2]*t[4]*t[5]*u1*(_hw_real)2.0+
            Izz*kf*l*t[2]*t[4]*t[5]*u3*(_hw_real)2.0))/(t[16]*(_hw_real)2.0);        
        }
        }

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
#if defined(__VITIS_HLS__) && defined(__SYNTHESIS__)
#ifdef USE_FAST_SIN_COS
#pragma HLS inline
    _hw_real new_angle = angle * (_real)HALF_MAX_CIRCLE_ANGLE_PI;
    return fastsin<_hw_real>(new_angle);
    // return half_sin(angle);
#else
    return half_sin(angle);
#endif
#else        
    return sin(angle);
#endif
        //_hw_real new_angle = angle * (_hw_real)HALF_MAX_CIRCLE_ANGLE_PI;
    }
    _hw_real local_cos(_hw_real angle){
//#pragma HLS inline
#if defined(__VITIS_HLS__) && defined(__SYNTHESIS__)
#ifdef USE_FAST_SIN_COS
    return fastcos<_hw_real>(angle);
    // return half_cos(angle);
#else
    return half_cos(angle);
#endif
#else        
    return cos(angle);
#endif
    }

};
#endif
