#include <math.h>
#include "mex.h"

typedef double _real;

/* Input Arguments */
#define INPUT_N 2
#define T_IN prhs[0]
#define Y_IN prhs[1]

/* Output Arguments */
#define OUTPUT_N 1
#define YP_OUT plhs[0]

const uint16_t n_U = 4;

static _real value_map(_real input, _real in_min, _real in_max, _real out_min, _real out_max){
  _real output = (input - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
  return output;
}

static void model(_real state_dot[], _real state[], _real control[]){
//Test variables
_real u_min[4] = {-100,-100,-100,-100};
_real u_max[4] = {100,100,100,100};
_real Ts = 0.1;
uint16_t Nx = 12;
uint16_t n_U = 4;

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

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    _real* state_dot;
    _real *state, *control;
    size_t m, n;

    /* Check for proper number of arguments */

    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:yprime:invalidNumInputs", "Two input arguments required.");
    } else if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:yprime:maxlhs", "Too many output arguments.");
    }
    
    /* Create a matrix for the return argument */
    YP_OUT = mxCreateDoubleMatrix((mwSize)1, (mwSize)12, mxREAL);

    /* Assign pointers to the various parameters */
    state_dot = mxGetPr(YP_OUT);
    state = mxGetPr(T_IN);
    control = mxGetPr(Y_IN);
    /* Do the actual computations in a subroutine */
    model(state_dot,state,control);
    return;
}