#include <math.h>
#include <float.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <algorithm>    // std::max std::min

#include "sniffbot.hpp"
#include "read_from_file.hpp"

#ifdef DEBUG
#include "aux_functions.hpp"
#endif

Sniffbot::Sniffbot(std::string config_file){
  name = "Sniffbot";

  load_configurations_from_file(config_file); 
    /*
    NonLinear model predivtive control for Quadrotor
                (y)
              
                m_1
                |
                |
        m_4 ---- o ---- m_2 (x)
                |
                |
                m_3
    */ 

    // Control Matrix
    // [T    tx   ty   tz]
    uu = (_real *) malloc(N*n_U*sizeof(_real));
    // State Matrix
    // % 12 states [x,y,z,phi,theta,psi,... diff ant]
    x_hat = (_real *) malloc(Nx*(N+1)*sizeof(_real));

    state_type[0]  = 0; // 
    state_type[1]  = 0; // 
    state_type[2]  = 0; // 
    state_type[3]  = 1; // 
    state_type[4]  = 1; // 
    state_type[5]  = 1; // 
    state_type[6]  = 0; // 
    state_type[7]  = 0; // 
    state_type[8]  = 0; // 
    state_type[9]  = 1; // 
    state_type[10] = 1; // 
    state_type[11] = 1; // 
}

void Sniffbot::model(_real * state_dot, _real * state, _real * control){

  _real Ixx = 1.2;    // Moment of Inertia (Ixx Iyy Izz)[kg.m^2] 
  _real Iyy = 1.2; 
  _real Izz = 2.3;	
  _real l  = .25;     // Arm length [m]
  _real kf = 1;       // Thrust (lift) coefficient [Ns^2]
  _real kM = 0.2;     // Moment (drag) coefficient [Nms^2]
  _real m  = 2;       // Mass of quadcopter [kg]
  _real g  = 9.81;    // Gravity [m/s^2]
  _real b  = 4.9050;  // Speed offset [rad^2/sec^2]4.9050
  
  _real M_mma[16] ={ 1,  1,  1,  1,  // thrust (T)
                     0, -1,  0,  1,  // roll  (tx)
                    -1,  0,  1,  0,  // pitch (ty)
                    -1,  1, -1,  1}; // yaw   (tz)
  _real M_mma_inv[16] = { 0.2500,         0,   -0.5000,   -0.2500,
                          0.2500,   -0.5000,         0,    0.2500,
                          0.2500,         0,    0.5000,   -0.2500,
                          0.2500,    0.5000,         0,    0.2500};
  // Input limits
  _real phys_u_min[4] = {-15,  -3,  -3,  -3};
  _real phys_u_max[4] = {15,   3,   3,   3};

  int i;
  _real w[4];
  _real param[4], u_param[4];
  //---------------------------------------------------
  _real t3;
  int j;
  _real t2, t4, t5, t6, t7, t8, t9, t12, t21, t22, t23, t24,t13;
  _real t14, t16, t17, t18, t19, t20, t25, t15, t27, t10, t11, t26;
  _real xdot, ydot, zdot, phi, phidot, psidot, psi, thetadot, theta;
  _real u1, u2, u3, u4;
  //---------------------------------------------------
  _real phys_control[4];
  //phys_control = (_real *) calloc(n_U,sizeof(_real));
  for (uint8_t i = 0; i < n_U; i++)
  {
    phys_control[i] = value_map(control[i], u_min[i], u_max[i], phys_u_min[i], phys_u_max[i]);
  }
  
  param[0] = 1.0 / kf;
  t3 = 1.0 / (l * kf);
  param[1] = t3;
  param[2] = t3;
  param[3] = 1.0 / kM;
  for (i = 0; i < n_U; i++) {
    w[i] = 0.0;
    u_param[i] = phys_control[i] * param[i];
  }

  for (i = 0; i < n_U; i++) {
    for (j = 0; j < n_U; j++) {
      w[i] += M_mma_inv[j + (i << 2)] * u_param[j];
    }
    w[i] += b;
  }

  phi = state[3];
  phidot = state[9];
  psidot = state[11];
  psi = state[5];
  thetadot = state[10];
  theta = state[4];
  u1 = w[0];
  u2 = w[1];
  u3 = w[2];
  u4 = w[3];
  xdot = state[6];
  ydot = state[7];
  zdot = state[8];
  t2 = cos(phi);
  t3 = cos(psi);
  t4 = cos(theta);
  t5 = sin(phi);
  t6 = sin(psi);
  t7 = sin(theta);
  t8 = pow(Iyy,2);
  t9 = pow(Izz,2);
  t10 = phi*2.0;
  t11 = theta*2.0;
  t12 = pow(thetadot,2);
  t21 = 1.0/Iyy;
  t22 = 1.0/Izz;
  t23 = 1.0/m;
  t24 = u1+u2+u3+u4;
  t13 = pow(t2,2);
  t14 = pow(t2,3);
  t16 = pow(t4,2);
  t17 = sin(t10);
  t18 = pow(t5,2);
  t19 = sin(t11);
  t20 = pow(t7,2);
  t25 = 1.0/t4;
  t15 = pow(t13,2);
  t26 = t20-1.0;
  t27 = 1.0/t26;

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
      (Ixx*phidot*t4*thetadot)/2.0+(Iyy*phidot*t4*thetadot)/2.0-
      (Izz*phidot*t4*thetadot)/2.0-(Ixx*psidot*t19*thetadot)/2.0+
      (Iyy*psidot*t19*thetadot)/2.0-Iyy*phidot*t4*t13*thetadot+
      Izz*phidot*t4*t13*thetadot+(Iyy*t2*t5*t7*t12)/2.0-
      (Izz*t2*t5*t7*t12)/2.0-Iyy*phidot*psidot*t2*t5*t16+
      Izz*phidot*psidot*t2*t5*t16-Iyy*psidot*t4*t7*t13*thetadot+
      Izz*psidot*t4*t7*t13*thetadot)-(t21*t22*t27*(-kf*l*u2+kf*l*u4+
      (Ixx*psidot*t4*thetadot)/2.0)*(Iyy*Izz+Ixx*Iyy*t20-Iyy*Izz*t20-
      Ixx*Iyy*t18*t20+Ixx*Izz*t18*t20))/Ixx+t7*t17*t21*t22*t25*
      (Iyy/2.0-Izz/2.0)*(kf*l*u1-kf*l*u3-(Iyy*phidot*t17*thetadot)/2.0+
      (Izz*phidot*t17*thetadot)/2.0-(Iyy*phidot*psidot*t4)/2.0+
      (Izz*phidot*psidot*t4)/2.0+Iyy*phidot*psidot*t4*t13-
      Izz*phidot*psidot*t4*t13-(Iyy*psidot*t2*t5*t7*thetadot)/2.0+
      (Izz*psidot*t2*t5*t7*thetadot)/2.0);

  state_dot[10] = (t21*t22*t25*(phidot*psidot*t8*t16-t7*t8*t12*t13-t7*t9*t12*t13+
      t7*t8*t12*t15+t7*t9*t12*t15+Iyy*kM*t17*u1-Iyy*kM*t17*u2+Iyy*kM*t17*u3-
      Iyy*kM*t17*u4-Izz*kM*t17*u1+Izz*kM*t17*u2-Izz*kM*t17*u3+Izz*kM*t17*u4-
      Iyy*Izz*phidot*psidot*t16+Iyy*Izz*t7*t12*t13*2.0-
      Iyy*Izz*t7*t12*t15*2.0-Iyy*kf*l*t4*u1*2.0+Iyy*kf*l*t4*u3*2.0-
      phidot*psidot*t8*t13*t16+phidot*psidot*t9*t13*t16+
      phidot*t2*t4*t5*t8*thetadot-phidot*t2*t4*t5*t9*thetadot+
      Iyy*kf*l*t4*t13*u1*2.0-Iyy*kf*l*t4*t13*u3*2.0-Izz*kf*l*t4*t13*u1*2.0+
      Izz*kf*l*t4*t13*u3*2.0-Ixx*Iyy*phidot*t2*t4*t5*thetadot+
      Ixx*Izz*phidot*t2*t4*t5*thetadot+Iyy*kf*l*t2*t5*t7*u2*2.0-
      Iyy*kf*l*t2*t5*t7*u4*2.0-Izz*kf*l*t2*t5*t7*u2*2.0+
      Izz*kf*l*t2*t5*t7*u4*2.0-psidot*t2*t4*t5*t7*t8*thetadot+
      psidot*t4*t5*t7*t8*t14*thetadot+psidot*t4*t5*t7*t9*t14*thetadot+
      Ixx*Iyy*psidot*t2*t4*t5*t7*thetadot-Ixx*Izz*psidot*t2*t4*t5*t7*thetadot
      +Iyy*Izz*psidot*t2*t4*t5*t7*thetadot-
      Iyy*Izz*psidot*t4*t5*t7*t14*thetadot*2.0))/2.0;

  state_dot[11] = (t21*t22*(Izz*kM*u1*-2.0+Izz*kM*u2*2.0-Izz*kM*u3*2.0+Izz*kM*u4*2.0-
      phidot*t4*t9*thetadot-Iyy*kM*t13*u1*2.0+Iyy*kM*t13*u2*2.0-
      Iyy*kM*t13*u3*2.0+Iyy*kM*t13*u4*2.0+Izz*kM*t13*u1*2.0-Izz*kM*t13*u2*2.0+
      Izz*kM*t13*u3*2.0-Izz*kM*t13*u4*2.0+Ixx*Izz*phidot*t4*thetadot+
      Iyy*Izz*phidot*t4*thetadot-(Ixx*Izz*psidot*t19*thetadot)/2.0+
      Iyy*Izz*psidot*t19*thetadot-Izz*kf*l*t7*u2*2.0+Izz*kf*l*t7*u4*2.0-
      phidot*t4*t8*t13*thetadot+phidot*t4*t9*t13*thetadot-t2*t5*t7*t9*t12+
      t5*t7*t8*t12*t14+t5*t7*t9*t12*t14-phidot*psidot*t2*t5*t8*t16+
      phidot*psidot*t2*t5*t9*t16+psidot*t4*t7*t8*t13*thetadot+
      psidot*t4*t7*t9*t13*thetadot-psidot*t4*t7*t8*t15*thetadot-
      psidot*t4*t7*t9*t15*thetadot+Ixx*Iyy*phidot*t4*t13*thetadot-
      Ixx*Izz*phidot*t4*t13*thetadot+Iyy*Izz*t2*t5*t7*t12-
      Iyy*Izz*t5*t7*t12*t14*2.0-Iyy*kf*l*t7*t13*u2*2.0+Iyy*kf*l*t7*t13*u4*2.0+
      Izz*kf*l*t7*t13*u2*2.0-Izz*kf*l*t7*t13*u4*2.0-
      Ixx*Iyy*psidot*t4*t7*t13*thetadot+Ixx*Izz*psidot*t4*t7*t13*thetadot-
      Iyy*Izz*psidot*t4*t7*t13*thetadot*2.0+
      Iyy*Izz*psidot*t4*t7*t15*thetadot*2.0+Iyy*kf*l*t2*t4*t5*u1*2.0-
      Iyy*kf*l*t2*t4*t5*u3*2.0-Izz*kf*l*t2*t4*t5*u1*2.0+
      Izz*kf*l*t2*t4*t5*u3*2.0))/(t16*2.0);

}

void Sniffbot::control_from_parameters(_real * parameters, _real * previous_control, _real * control, int horizon) {

    _real t1 = exp(-Lambda*Ts)-1;
    _real t2 = exp(-q_param*Lambda*Ts)-1;
    _real t1_2 = 1/(t1-t2);

    _real d1_p1 = du_max[0]*parameters[0];
    _real d2_p2 = du_max[1]*parameters[1];
    _real d3_p3 = du_max[2]*parameters[2];
    _real d4_p4 = du_max[3]*parameters[3];
    
    _real t1_uss = t1*uss;
    _real t2_uss = t2*uss;
    
    _real alpha_1 =  (d1_p1 - t2*previous_control[0] + t2_uss)*t1_2;
    _real alpha_2 = -(d1_p1 - t1*previous_control[0] + t1_uss)*t1_2;
    _real alpha_3 =  (d2_p2 - t2*previous_control[1] + t2_uss)*t1_2;
    _real alpha_4 = -(d2_p2 - t1*previous_control[1] + t1_uss)*t1_2;
    _real alpha_5 =  (d3_p3 - t2*previous_control[2] + t2_uss)*t1_2;
    _real alpha_6 = -(d3_p3 - t1*previous_control[2] + t1_uss)*t1_2;
    _real alpha_7 =  (d4_p4 - t2*previous_control[3] + t2_uss)*t1_2;
    _real alpha_8 = -(d4_p4 - t1*previous_control[3] + t1_uss)*t1_2;
    
    _real pre_exponent;
    _real exp1;
    _real exp2;
    
    if (horizon == 1) {
        pre_exponent = 0;
        exp1 = 1;           // Transformar em TABELA
        exp2 = 1;   // Transformar em TABELA
        control[0] = std::min(std::max(uss + alpha_1 + alpha_2, u_min[0]), u_max[0]);
        control[1] = std::min(std::max(uss + alpha_3 + alpha_4, u_min[1]), u_max[1]);
        control[2] = std::min(std::max(uss + alpha_5 + alpha_6, u_min[2]), u_max[2]);
        control[3] = std::min(std::max(uss + alpha_7 + alpha_8, u_min[3]), u_max[3]);
    }
    else{
        for (int i = 0; i < horizon; ++i) {
            pre_exponent = -Lambda*Ts*i;
            exp1 = exp(pre_exponent);           // Transformar em TABELA
            exp2 = exp(pre_exponent*q_param);   // Transformar em TABELA
            control[0*horizon+i] = std::min(std::max(uss + alpha_1 * exp1 + alpha_2 * exp2, u_min[0]), u_max[0]);
            control[1*horizon+i] = std::min(std::max(uss + alpha_3 * exp1 + alpha_4 * exp2, u_min[1]), u_max[1]);
            control[2*horizon+i] = std::min(std::max(uss + alpha_5 * exp1 + alpha_6 * exp2, u_min[2]), u_max[2]);
            control[3*horizon+i] = std::min(std::max(uss + alpha_7 * exp1 + alpha_8 * exp2, u_min[3]), u_max[3]);

        }
    }
}


_real Sniffbot::nmpc_cost_function(_real * control_guess, _real * xref, _real * uref, _real * xss, _real * uss){
	
	_real J = 0.0;
    _real * Ji ;
	_real Jf = 0.0;
    _real Ju = 0.0;
    //_real temp;

    Ji = (_real *) calloc(Nx,sizeof(_real));

	// Constructing the vector of guessed control actions with respect to N and Nu
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < n_U; ++j) {
			if(i < Nu) {
				uu[i*n_U+j] = control_guess[j*Nu+i];
			}
			else {
				uu[i*n_U+j] = control_guess[j*Nu-1];	
			}
		}
	}
    
//#ifdef DEBUG
//    print_array("xss", xss, Nx, 0);
//#endif
    
	// Initialize x_hat vector
	for (int i = 0; i < Nx; ++i) {
		x_hat[i] = current_state[i];
	}
	for (int i = Nx; i < (Nx*N); ++i) {
		x_hat[i] = 0.0;
	}

	// Calculate cost function
	int k = 0;
	int l = 0;
    J = 0;
    Jf = 0;
    Ju = 0;
    // Cost of States
	for (int i = 0; i < N; ++i) {
		k = Nx*(i+1);
//#ifdef DEBUG
//      print_array("x_hat", &x_hat[Nx*i], Nx, 0);
//#endif
        one_step_prediction(&x_hat[k], &x_hat[Nx*i], &uu[n_U*i]);
		for (int j = 0; j < Nx; j=j+2) {
			l = k + j;
            //printf("Ji[%d] = error * (x_hat[%d] - xref[%d]) = %f * ( %f - %f ) \t ^2 = %f\n", j, Nx*i+j, Nx*i+j, error_with_constrains(j,x_hat[Nx*i+j]), x_hat[Nx*i+j], xref[Nx*i+j], pow(error_with_constrains(j,x_hat[Nx*i+j])*( x_hat[Nx*i+j] - xref[Nx*i+j]),2));
            Ji[j] = Ji[j] + pow(error_with_constrains(j,x_hat[Nx*i+j])*( x_hat[Nx*i+j] - xref[Nx*i+j]),2);
		}
	}
    for (int j = 0; j < Nx; j=j+2) {
//#ifdef DEBUG
//        printf("J[%d]: %f\n", j, Q[j]*sqrt(Ji[j]));
//#endif
        J = J + Q[j]*sqrt(Ji[j]);
    }

    // Cost of States - End of Horizon
    k = Nx*N;
    for (int j = 0; j < Nx; j++) {
        l = k + j;
        Jf = Jf + pow(error_with_constrains(j,x_hat[l])*( x_hat[l] - xss[j]),2);
    }   
    J = J + Qf[0]*sqrt(Jf); 
//#ifdef DEBUG
//      printf("Jf: %f\n", Qf[0]*sqrt(Jf));
//#endif
    // Cost of control
    if(R[0] > 0) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < n_U; ++j) {
                //printf("Ju = %f |\t uu[%d]-uss[%d] = %f | ^2 = %f\n", Ju, n_U*i+j, j, uu[n_U*i+j]-uss[j], pow(uu[n_U*i+j]-uss[j],2));
                Ju = Ju + pow(uu[n_U*i+j]-uss[j],2);
            }
        }
    }
//#ifdef DEBUG
//    printf("Ju: %f sqrt(Ju)*R:  %f\n", Ju, R[0]*sqrt(Ju));
//#endif
    J = J + R[0]*sqrt(Ju);
    
#ifdef DEBUG
    if(isnan(J) || isinf(J)){
        for (int i=0; i < N; i++) {
            print_array("uu", &uu[n_U*i], n_U, 0);
        }
        for (int i = 0; i < N; i++) {
            print_array("x_hat", &x_hat[Nx*i], Nx, 1);
        }
    }
#endif
    
	return J;
}



_real Sniffbot::value_map(_real input, _real in_min, _real in_max, _real out_min, _real out_max){
  _real output = (input - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
  return output;
}