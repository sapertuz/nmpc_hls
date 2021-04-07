#include <math.h>
#include <float.h>
#include <cmath>

#include "sniffbot.hpp"

Sniffbot::Sniffbot(){
  load_configurations_from_file(); 
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
    // State Matrix
    // % 12 states [x,y,z,phi,theta,psi,... diff ant]
    state_type[0]  = 0; // 
    state_type[1]  = 0; // 
    state_type[2]  = 0; // 
    state_type[3]  = 0; // 
    state_type[4]  = 0; // 
    state_type[5]  = 0; // 
    state_type[6]  = 0; // 
    state_type[7]  = 0; // 
    state_type[8]  = 0; // 
    state_type[9]  = 0; // 
    state_type[10] = 0; // 
    state_type[11] = 0; // 
}

void Sniffbot::model(
    _real state_dot[],
    _real state[],
    _real control[]
){
	_real _state_dot[_Nx];
	_real _state[_Nx];
	_real _control[_n_U];

	for (short i=0; i<_Nx; i++){
		_state[i] = state[i];
	}
	for (short i=0; i<_n_U; i++){
		_control[i] = control[i];
	}

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
    _real phys_u_min[] = {-15,  -3,  -3,  -3};
    _real phys_u_max[] = {15,   3,   3,   3};

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
    for (short i = 0; i < _n_U; i++)
    {
    phys_control[i] = value_map(_control[i], u_min[i], u_max[i], phys_u_min[i], phys_u_max[i]);
    }

    param[0] = (_real)1.0 / kf;
    t3 = (_real)1.0 / (l * kf);
    param[1] = t3;
    param[2] = t3;
    param[3] = (_real)1.0 / kM;
    for (i = 0; i < _n_U; i++) {
        w[i] = 0.0;
        u_param[i] = phys_control[i] * param[i];
    }

    for (i = 0; i < _n_U; i++) {
        for (j = 0; j < _n_U; j++) {
            w[i] += M_mma_inv[j + (i << 2)] * u_param[j];
        }
        w[i] += b;
    }

    phi = _state[3];
    theta = _state[4];
    psi = _state[5];
    xdot = _state[6];
    ydot = _state[7];
    zdot = _state[8];
    phidot = _state[9];
    psidot = _state[11];
    thetadot = _state[10];

    u1 = w[0];
    u2 = w[1];
    u3 = w[2];
    u4 = w[3];

    t2 = cos(phi);
    t4 = cos(theta);
    t3 = cos(psi);
    t5 = sin(phi);
    t7 = sin(theta);
    t6 = sin(psi);
    t8 = Iyy*Iyy;
    t9 = Izz*Izz;
    t10 = phi*(_real)2.0;
    t11 = theta*(_real)2.0;
    t12 = thetadot*thetadot;
    t21 = (_real)1.0/Iyy;
    t22 = (_real)1.0/Izz;
    t23 = (_real)1.0/m;
    t24 = u1+u2+u3+u4;
    t13 = t2*t2;
    t14 = t2*t2*t2;
    t16 = t4*t4;
    t17 = sin(t10);
    t18 = t5*t5;
    t19 = sin(t11);
    t20 = t7*t7;
    t25 = (_real)1.0/t4;
    t15 = t13*t13;
    t26 = t20-(_real)1.0;
    t27 = (_real)1.0/t26;

    _state_dot[0] = xdot;
    _state_dot[1] = ydot;
    _state_dot[2] = zdot;
    _state_dot[3] = phidot;
    _state_dot[4] = thetadot;
    _state_dot[5] = psidot;

    _state_dot[6] = kf*t23*t24*(t5*t6+t2*t3*t7);

    _state_dot[7] = -kf*t23*t24*(t3*t5-t2*t6*t7);

    _state_dot[8] = -g+kf*t2*t4*t23*t24;

    _state_dot[9] = -t7*t21*t22*t27*(Iyy-Iyy*t18+Izz*t18)*(-kM*u1+kM*u2-kM*u3+kM*u4+
        (Ixx*phidot*t4*thetadot)*(_real)0.5+(Iyy*phidot*t4*thetadot)*(_real)0.5-
        (Izz*phidot*t4*thetadot)*(_real)0.5-(Ixx*psidot*t19*thetadot)*(_real)0.5+
        (Iyy*psidot*t19*thetadot)*(_real)0.5-Iyy*phidot*t4*t13*thetadot+
        Izz*phidot*t4*t13*thetadot+(Iyy*t2*t5*t7*t12)*(_real)0.5-
        (Izz*t2*t5*t7*t12)*(_real)0.5-Iyy*phidot*psidot*t2*t5*t16+
        Izz*phidot*psidot*t2*t5*t16-Iyy*psidot*t4*t7*t13*thetadot+
        Izz*psidot*t4*t7*t13*thetadot)-(t21*t22*t27*(-kf*l*u2+kf*l*u4+
        (Ixx*psidot*t4*thetadot)*(_real)0.5)*(Iyy*Izz+Ixx*Iyy*t20-Iyy*Izz*t20-
        Ixx*Iyy*t18*t20+Ixx*Izz*t18*t20))/Ixx+t7*t17*t21*t22*t25*
        (Iyy*(_real)0.5-Izz*(_real)0.5)*(kf*l*u1-kf*l*u3-(Iyy*phidot*t17*thetadot)*(_real)0.5+
        (Izz*phidot*t17*thetadot)*(_real)0.5-(Iyy*phidot*psidot*t4)*(_real)0.5+
        (Izz*phidot*psidot*t4)*(_real)0.5+Iyy*phidot*psidot*t4*t13-
        Izz*phidot*psidot*t4*t13-(Iyy*psidot*t2*t5*t7*thetadot)*(_real)0.5+
        (Izz*psidot*t2*t5*t7*thetadot)*(_real)0.5);

    _state_dot[10] = (t21*t22*t25*(phidot*psidot*t8*t16-t7*t8*t12*t13-t7*t9*t12*t13+
        t7*t8*t12*t15+t7*t9*t12*t15+Iyy*kM*t17*u1-Iyy*kM*t17*u2+Iyy*kM*t17*u3-
        Iyy*kM*t17*u4-Izz*kM*t17*u1+Izz*kM*t17*u2-Izz*kM*t17*u3+Izz*kM*t17*u4-
        Iyy*Izz*phidot*psidot*t16+Iyy*Izz*t7*t12*t13*(_real)2.0-
        Iyy*Izz*t7*t12*t15*(_real)2.0-Iyy*kf*l*t4*u1*(_real)2.0+Iyy*kf*l*t4*u3*(_real)2.0-
        phidot*psidot*t8*t13*t16+phidot*psidot*t9*t13*t16+
        phidot*t2*t4*t5*t8*thetadot-phidot*t2*t4*t5*t9*thetadot+
        Iyy*kf*l*t4*t13*u1*(_real)2.0-Iyy*kf*l*t4*t13*u3*(_real)2.0-Izz*kf*l*t4*t13*u1*(_real)2.0+
        Izz*kf*l*t4*t13*u3*(_real)2.0-Ixx*Iyy*phidot*t2*t4*t5*thetadot+
        Ixx*Izz*phidot*t2*t4*t5*thetadot+Iyy*kf*l*t2*t5*t7*u2*(_real)2.0-
        Iyy*kf*l*t2*t5*t7*u4*(_real)2.0-Izz*kf*l*t2*t5*t7*u2*(_real)2.0+
        Izz*kf*l*t2*t5*t7*u4*(_real)2.0-psidot*t2*t4*t5*t7*t8*thetadot+
        psidot*t4*t5*t7*t8*t14*thetadot+psidot*t4*t5*t7*t9*t14*thetadot+
        Ixx*Iyy*psidot*t2*t4*t5*t7*thetadot-Ixx*Izz*psidot*t2*t4*t5*t7*thetadot
        +Iyy*Izz*psidot*t2*t4*t5*t7*thetadot-
        Iyy*Izz*psidot*t4*t5*t7*t14*thetadot*(_real)2.0))*(_real)0.5;

    _state_dot[11] = (t21*t22*(Izz*kM*u1*-(_real)2.0+Izz*kM*u2*(_real)2.0-Izz*kM*u3*(_real)2.0+Izz*kM*u4*(_real)2.0-
        phidot*t4*t9*thetadot-Iyy*kM*t13*u1*(_real)2.0+Iyy*kM*t13*u2*(_real)2.0-
        Iyy*kM*t13*u3*(_real)2.0+Iyy*kM*t13*u4*(_real)2.0+Izz*kM*t13*u1*(_real)2.0-Izz*kM*t13*u2*(_real)2.0+
        Izz*kM*t13*u3*(_real)2.0-Izz*kM*t13*u4*(_real)2.0+Ixx*Izz*phidot*t4*thetadot+
        Iyy*Izz*phidot*t4*thetadot-(Ixx*Izz*psidot*t19*thetadot)*(_real)0.5+
        Iyy*Izz*psidot*t19*thetadot-Izz*kf*l*t7*u2*(_real)2.0+Izz*kf*l*t7*u4*(_real)2.0-
        phidot*t4*t8*t13*thetadot+phidot*t4*t9*t13*thetadot-t2*t5*t7*t9*t12+
        t5*t7*t8*t12*t14+t5*t7*t9*t12*t14-phidot*psidot*t2*t5*t8*t16+
        phidot*psidot*t2*t5*t9*t16+psidot*t4*t7*t8*t13*thetadot+
        psidot*t4*t7*t9*t13*thetadot-psidot*t4*t7*t8*t15*thetadot-
        psidot*t4*t7*t9*t15*thetadot+Ixx*Iyy*phidot*t4*t13*thetadot-
        Ixx*Izz*phidot*t4*t13*thetadot+Iyy*Izz*t2*t5*t7*t12-
        Iyy*Izz*t5*t7*t12*t14*(_real)2.0-Iyy*kf*l*t7*t13*u2*(_real)2.0+Iyy*kf*l*t7*t13*u4*(_real)2.0+
        Izz*kf*l*t7*t13*u2*(_real)2.0-Izz*kf*l*t7*t13*u4*(_real)2.0-
        Ixx*Iyy*psidot*t4*t7*t13*thetadot+Ixx*Izz*psidot*t4*t7*t13*thetadot-
        Iyy*Izz*psidot*t4*t7*t13*thetadot*(_real)2.0+
        Iyy*Izz*psidot*t4*t7*t15*thetadot*(_real)2.0+Iyy*kf*l*t2*t4*t5*u1*(_real)2.0-
        Iyy*kf*l*t2*t4*t5*u3*(_real)2.0-Izz*kf*l*t2*t4*t5*u1*(_real)2.0+
        Izz*kf*l*t2*t4*t5*u3*(_real)2.0))/(t16*(_real)2.0);

	for (short i=0; i<_Nx; i++){
		state_dot[i] = _state_dot[i];
	}

}

_real Sniffbot::nmpc_cost_function(_real * control_guess, _real * xref, _real * uref, _real * xss, _real * uss){
	
	_real J = 0.0;
    _real Ji[] = {0,0,0,0,0,0,0,0,0,0,0,0};
	_real Jf = 0.0;
    _real Ju = 0.0;
    //_real temp;

	// Constructing the vector of guessed control actions with respect to N and Nu
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < _n_U; ++j) {
			if(i < Nu) {
				uu[i*_n_U+j] = control_guess[j*Nu+i];
			}
			else {
				uu[i*_n_U+j] = control_guess[j*Nu-1];	
			}
		}
	}
    
	// Initialize x_hat vector
	for (int i = 0; i < _Nx; ++i) {
		x_hat[i] = current_state[i];
	}
	for (int i = _Nx; i < (_Nx*N); ++i) {
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
		k = _Nx*(i+1);
        one_step_prediction(&x_hat[k], &x_hat[_Nx*i], &uu[_n_U*i]);
		for (int j = 0; j < _Nx; j=j+2) {
			l = k + j;
            _real tmp_err = error_with_constrains(j,x_hat[_Nx*i+j])*( x_hat[_Nx*i+j] - xref[_Nx*i+j]);
            Ji[j] = Ji[j] + tmp_err*tmp_err;
		}
	}
    for (int j = 0; j < _Nx; j=j+2) {
        J = J + Q[j]*sqrt(Ji[j]);
    }

    // Cost of States - End of Horizon
    k = _Nx*N;
    for (int j = 0; j < _Nx; j++) {
        l = k + j;
        _real tmp_err = error_with_constrains(j,x_hat[l])*( x_hat[l] - xss[j]);
        Jf = Jf + tmp_err*tmp_err;
    }   
    J = J + Qf[0]*sqrt(Jf); 

    // Cost of control
    if(R[0] > 0) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < _n_U; ++j) {
                _real tmp_err = uu[_n_U*i+j]-uss[j];
                Ju = Ju + tmp_err*tmp_err;
            }
        }
    }

	return J;
}



_real Sniffbot::value_map(_real input, _real in_min, _real in_max, _real out_min, _real out_max){
  _real output = (input - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
  return output;
}