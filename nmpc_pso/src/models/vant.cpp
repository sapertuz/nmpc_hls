#include <math.h>
#include <float.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <algorithm>    // std::max std::min

#include "vant.hpp"
#include "read_from_file.hpp"

#ifdef DEBUG
#include "aux_functions.hpp"
#endif

Vant::Vant(std::string config_file) {
	name = "Vant";

    load_configurations_from_file(config_file);    

	// Setup Output Plot
	plot_number_of_states = 3;
    plot_number_of_controls = 1;
    plot_states_config = (int *) calloc(plot_number_of_states*Nx, sizeof(int));
    plot_control_config = (int *) calloc(plot_number_of_controls*n_U, sizeof(int));

    // State Matrix
    // xyz   	      1 0 1 0 1 0 0 0 0 0 0 0
    // x.y.z.         0 1 0 1 0 1 0 0 0 0 0 0
    // phi theta psi  0 0 0 0 0 0 1 0 1 0 1 0
    plot_states_config[0] = 1;
    plot_states_config[2] = 1;
    plot_states_config[4] = 1;
    plot_states_config[1] = 1;
    plot_states_config[3] = 1;
    plot_states_config[5] = 1;
    plot_states_config[7] = 1;
    plot_states_config[9] = 1;
    plot_states_config[11] = 1;

    // Control Matrix
    // u 		 1 
    plot_control_config[0] = 1;
    plot_control_config[1] = 1;
    plot_control_config[2] = 1;
    plot_control_config[3] = 1;

    plot_labels = new std::string[4];
    plot_labels[0] = "position (m)";
    plot_labels[0] = "velocity (m/s)";
    plot_labels[1] = "angles (rad/s)";
    plot_labels[2] = "u (N)";

    state_type = (int *) calloc(Nx, sizeof(int));
    state_type[0] = 0; // Pos 1
    state_type[1] = 0; // Linear Velocity 1
    state_type[2] = 0; // Pos 2
    state_type[3] = 0; // Linear Velocity 2
    state_type[4] = 0; // Pos 3
    state_type[5] = 0; // Linear Velocity 3
    state_type[6] = 1; // Angle 1
    state_type[7] = 0; // Angular Velocity 1
    state_type[8] = 1; // Angle 2
    state_type[9] = 0; // Angular Velocity 2
    state_type[10] = 1; // Angle 3
    state_type[11] = 0; // Angular Velocity 3


    uu = (_real *) malloc(N*n_U*sizeof(_real));
    x_hat = (_real *) malloc(Nx*(N+1)*sizeof(_real));

}

void Vant::model(_real * state_dot, _real * state, _real * control){
    // VANT qudrotor
    // States = x 
    //          x_dot 
    //          y 
    //          y_dot 
    //          z 
    //          z_dot 
    //          Phi 
    //          Phi_dot 
    //          Theta 
    //          Theta_dot 
    //          Psi 
    //          Psi_dot

    _real b = 2.64e-4;    // Coeficiente de empuxo dos motores [N.(s^2)]
    _real L = 0.5;        // Distância dos motores ao centro do quadrirotor (meia envergadura)[m]
    _real d = 7.5e-7;     // Coeficiente de arrasto [N.m.(s^2)]            
    _real m = 4;          // Massa total [Kg]
    _real g = 9.81;       // Gravidade local

    _real Ixx = 0.033;    // Momentos de inércia em torno dos eixos de rolagem, arfagem e guinada, respectivamente. [Kg.(m^2)]
    _real Iyy = 0.033;
    _real Izz = 0.066;
    
    // _real u1 = control[0];
    // _real u2 = control[1];
    // _real u3 = control[2];
    // _real u4 = control[3];

//    printf("u = %f\t%f\t%f\t%f\n", u1, u2, u3, u4);
    _real u1_2 = control[0]*control[0];
    _real u2_2 = control[1]*control[1];
    _real u3_2 = control[2]*control[2];
    _real u4_2 = control[3]*control[3];

    _real U=b*(u1_2+u2_2+u3_2+u4_2);
    _real taux=b*L*(u4_2-u2_2);
    _real tauy=b*L*(u3_2-u1_2);
    _real tauz=d*(u1_2-u2_2+u3_2-u4_2);

//    printf("U = %f\ntaux%f\ntauy%f\nyauz%f\n", U, taux, tauy, tauz);
    // state_dot = zeros(12,1);
    // state_dot(1)  = state[2);
    // state_dot(2)  = (cos(state[11))*sin(state[9))*cos(state[7))+sin(state[11))*sin(state[7)))*(U/m);
    // state_dot(3)  = state[4);
    // state_dot(4)  = (sin(state[11))*sin(state[9))*cos(state[7))-sin(state[7))*cos(state[11)))*(U/m);
    // state_dot(5)  = state[6);
    // state_dot(6)  = -g+(cos(state[9))*cos(state[7)))*(U/m);
    // state_dot(7)  = state[8);
    // state_dot(8)  = -(((state[12)*state[10)*cos(state[9))*(Ixx+(Iyy-Izz)*(2*(cos(state[7)))^2 -1)))/(Ixx)-0.5*(state[10)^2*sin(2*state[7))*(Iyy-Izz))/(Ixx)+(0.5*state[12)^2*sin(2*state[7))*(cos(state[9)))^2*(Iyy-Izz))/(Ixx)+(taux/Ixx))*sin(state[9))*Ixx*(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2)-0.5*(-0.5*(state[12))^2*sin(2*state[9))*(-Ixx+Iyy*(sin(state[7)))^2+Izz*(cos(state[7)))^2)+state[10)*state[8)*sin(2*state[7))*(Izz-Iyy)+state[12)*state[8)*cos(state[9))*(cos(2*state[7))*(Iyy-Izz)+Ixx)+tauy)*sin(2*state[7))*cos(state[9))*(Iyy-Izz)+(-state[10)*state[12)*sin(2*state[9))*(Ixx-Izz*(cos(state[7))^2)+Iyy*(sin(state[7)))^2)+state[12)*state[8)*sin(2*state[7))*(cos(state[9)))^2*(Iyy-Izz)-state[10)*state[8)*cos(state[9))*(Ixx+(Iyy-Izz)*(2*(cos(state[7)))^2-1))+0.5*(state[10))^2*sin(2*state[7))*sin(state[9))*(Iyy-Izz)+tauz)*(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2))*sin(state[9))/((sin(state[9)))^2*Ixx*(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2)+0.25*(sin(2*state[7)))^2*(cos(state[9)))^2*(Iyy-Izz)^2-(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2)*((cos(state[9)))^2*Izz*(cos(state[7)))^2+Iyy*(sin(state[7)))^2+(sin(state[9)))^2*Ixx))+(state[12)*state[10)*cos(state[9))*(Ixx+(Iyy-Izz)*(2*(cos(state[7)))^2-1)))/(Ixx)-0.5*(state[10)^2*sin(2*state[7))*(Iyy-Izz))/(Ixx)+(0.5*state[12)^2*sin(2*state[7))*(cos(state[9)))^2*(Iyy-Izz))/(Ixx)+(taux)/(Ixx);
    // state_dot(9)  = state[10);
    // state_dot(10) = (0.5*(((state[12)*state[10)*cos(state[9))*(Ixx+(Iyy-Izz)*(2*(cos(state[7)))^2-1)))/(Ixx)-0.5*(state[10)^2*sin(2*state[7))*(Iyy-Izz))/(Ixx)+(0.5*state[12)^2*sin(2*state[7))*(cos(state[9)))^2*(Iyy-Izz))/(Ixx)+(taux)/(Ixx))*sin(state[9))*Ixx*(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2)-0.5*(-0.5*(state[12)^2)*sin(2*state[9))*(-Ixx+Iyy*sin(state[7))^2+Izz*(cos(state[7)))^2)+state[10)*state[8)*sin(2*state[7))*(Izz-Iyy)+state[12)*state[8)*cos(state[9))*(cos(2*state[7))*(Iyy-Izz)+Ixx)+tauy)*sin(2*state[7))*cos(state[9))*(Iyy-Izz)+(-state[10)*state[12)*sin(2*state[9))*(Ixx-Izz*(cos(state[7)))^2+Iyy*(sin(state[7)))^2)+state[12)*state[8)*sin(2*state[7))*(cos(state[9)))^2*(Iyy-Izz)-state[10)*state[8)*cos(state[9))*(Ixx+(Iyy-Izz)*(2*(cos(state[7)))^2-1))+0.5*state[10)^2*sin(2*state[7))*sin(state[9))*(Iyy-Izz)+tauz)*(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2))*(sin(2*state[7))*cos(state[9))*(Iyy-Izz))/((sin(state[9)))^2*Ixx*(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2)+0.25*(sin(2*state[7)))^2*(cos(state[9)))^2*(Iyy-Izz)^2-(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2)*((cos(state[9)))^2*Izz*(cos(state[7)))^2+Iyy*(sin(state[7)))^2+(sin(state[9)))^2*Ixx))-0.5*(state[12)^2)*sin(2*state[9))*(-Ixx+Iyy*(sin(state[7)))^2+Izz*(cos(state[7)))^2)+state[10)*state[8)*sin(2*state[7))*(Izz-Iyy)+state[12)*state[8)*cos(state[9))*(cos(2*state[7))*(Iyy-Izz)+Ixx)+tauy)/(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2);
    // state_dot(11) = state[12);
    // state_dot(12) = -(((state[12)*state[10)*cos(state[9))*(Ixx+(Iyy-Izz)*(2*(cos(state[7)))^2-1)))/(Ixx)-0.5*((state[10)^2)*sin(2*state[7))*(Iyy-Izz))/(Ixx)+(0.5*state[12)*sin(2*state[7))*(cos(state[9)))^2*(Iyy-Izz))/(Ixx)+(taux)/(Ixx))*sin(state[9))*Ixx*(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2)-0.5*(-0.5*state[12)^2*sin(2*state[9))*(-Ixx+Iyy*(sin(state[7)))^2+Izz*(cos(state[7)))^2)+state[10)*state[8)*sin(2*state[7))*(Izz-Iyy)+state[12)*state[8)*cos(state[9))*(cos(2*state[7))*(Iyy-Izz)+Ixx)+tauy)*sin(2*state[7))*cos(state[9))*(Iyy-Izz)+(-state[10)*state[12)*sin(2*state[9))*(Ixx-Izz*(cos(state[7)))^2+Iyy*(sin(state[7)))^2)+state[12)*state[8)*sin(2*state[7))*(cos(state[9)))^2*(Iyy-Izz)-state[10)*state[8)*cos(state[9))*(Ixx+(Iyy-Izz)*(2*cos((state[7)))^2-1))+0.5*(state[10)^2)*sin(2*state[7))*sin(state[9))*(Iyy-Izz)+tauz)*(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2))/((sin(state[9)))^2*Ixx*(Iyy*(cos(state[7)))^2+Izz*(sin(state[7)))^2)+0.25*(sin(2*state[7)))^2*(cos(state[9)))^2*(Iyy-Izz)^2-(Iyy*cos(state[7))^2+Izz*(sin(state[7)))^2)*((cos(state[9)))^2*Izz*(cos(state[7)))^2+Iyy*(sin(state[7)))^2+(sin(state[9)))^2*Ixx));

    _real s7 = sin(state[6]);
    _real s9 = sin(state[8]);
    _real s11 = sin(state[10]);

    _real c7 = cos(state[6]);
    _real c9 = cos(state[8]);
    _real c11 = cos(state[10]);

//    printf("s7 = %f\t s9 = %f\t s11 = %f\t c7 = %f\t c9 = %f\t c11 = %f\n", s7,s9,s11,c7,c9,c11);
    _real s2_7 = sin(2*state[6]);
    _real s2_9 = sin(2*state[8]);

    _real c2_7 = cos(2*state[6]);

    _real s7_2 = s7*s7; //(s7)^2;
    _real s9_2 = s9*s9; //(s9)^2;

    _real c9_2 = c9*c9; //(c9)^2;
    _real c7_2 = c7*c7; //(c7)^2;

    _real s2_7_c9_2 = s2_7*c9_2;

    _real x10_x8 = state[9]*state[7];
    _real x12_x8 = state[11]*state[7];
    _real x10_2 = state[9]*state[9];
    _real x12_2 = state[11]*state[11];

    _real Iyy_Izz = (Iyy-Izz);
    //_real Ixx_Iyy_Izz = Ixx+Iyy_Izz;

    _real U_m = U/m;
    _real taux_Ixx = taux/Ixx;

    _real x12_x11_init = (state[11]*state[9]*c9*(Ixx+Iyy_Izz*(2*c7_2-1)));
    _real x12_x11_init_Ixx = x12_x11_init/(Ixx);

    _real oper0 = (Iyy*c7_2+Izz*s7_2);
    _real oper01 = (c9_2*Izz*c7_2+Iyy*s7_2+s9_2*Ixx);
    _real oper1 = x12_x11_init_Ixx-0.5*(x10_2*s2_7*Iyy_Izz)/(Ixx);
    _real oper2 = (0.5*x12_2*s2_7_c9_2*Iyy_Izz);
    _real oper3 = s9*Ixx*oper0;
    _real oper4 = 0.5*(-0.5*x12_2*s2_9*(-Ixx+Iyy*s7_2+Izz*c7_2)+x10_x8*s2_7*(Izz-Iyy)+x12_x8*c9*(c2_7*Iyy_Izz+Ixx)+tauy)*s2_7*c9*Iyy_Izz;
    _real oper5 = (-state[9]*state[11]*s2_9*(Ixx-Izz*c7_2+Iyy*s7_2)+x12_x8*s2_7_c9_2*Iyy_Izz-x10_x8*c9*(Ixx+Iyy_Izz*(2*c7_2-1))+0.5*x10_2*s2_7*s9*Iyy_Izz+tauz)*oper0;
    _real oper6 = 0.25*s2_7*s2_7*c9_2*Iyy_Izz*Iyy_Izz;
   
    state_dot[0]  = state[1];
    state_dot[1]  = (c11*s9*c7+s11*s7)*(U_m);
    state_dot[2]  = state[3];
    state_dot[3]  = (s11*s9*c7-s7*c11)*(U_m);
    state_dot[4]  = state[5];
    state_dot[5]  = -g+(c9*c7)*(U_m);
    state_dot[6]  = state[7];
    state_dot[7]  =     -((oper1+                        oper2/(Ixx)+taux_Ixx)*oper3-oper4+oper5)*s9/(s9_2*Ixx*oper0+oper6-oper0*oper01)+oper1+oper2/(Ixx)+taux_Ixx;
    state_dot[8]  = state[9];
    state_dot[9] = (0.5*((oper1+                        oper2/(Ixx)+taux_Ixx)*oper3-oper4+oper5)*(s2_7*c9*Iyy_Izz)/(s9_2*Ixx*oper0+oper6-oper0*oper01)-0.5*(x12_2)*s2_9*(-Ixx+Iyy*s7_2+Izz*c7_2)+x10_x8*s2_7*(Izz-Iyy)+x12_x8*c9*(c2_7*Iyy_Izz+Ixx)+tauy)/oper0;
    state_dot[10] = state[11];
    state_dot[11] =     -((oper1+(0.5*state[11]*s2_7_c9_2*Iyy_Izz)/(Ixx)+taux_Ixx)*oper3-oper4+oper5)/(s9_2*Ixx*oper0+oper6-(Iyy*c7_2+Izz*s7_2)*oper01);
    
}

void Vant::control_from_parameters(_real * parameters, _real * previous_control, _real * control, int horizon) {

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


_real Vant::nmpc_cost_function(_real * control_guess, _real * xref, _real * uref, _real * xss, _real * uss){
	
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
