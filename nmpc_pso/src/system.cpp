#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include "aux_functions.hpp"
#include "read_from_file.hpp"
#include "system.hpp"

//#define DEBUG

void System::initializeVectors() {
	Q = (_real *) malloc(Nx*sizeof(_real));
	Qf = (_real *) malloc(Nx*sizeof(_real));
	R = (_real *) malloc(n_U*sizeof(_real));
	u_max = (_real *) malloc(n_U*sizeof(_real));
	u_min = (_real *) malloc(n_U*sizeof(_real));
	du_max = (_real *) malloc(n_U*sizeof(_real));
	state_upper_limits = (_real *) malloc(Nx*sizeof(_real));
	state_lower_limits = (_real *) malloc(Nx*sizeof(_real));
	current_state = (_real *) malloc(Nx*sizeof(_real));
	last_state = (_real *) malloc(Nx*sizeof(_real));
    rising_time = (_real *) malloc(Nx*sizeof(_real));
	
	if((Q == NULL) || (Qf == NULL) || (state_upper_limits == NULL) || (state_lower_limits == NULL) || (current_state == NULL) || (last_state == NULL)) {
		throw(1);
	}
}

void System::initializeParametrizationVectors(){
    pmax = (_real *) malloc(parametrized*sizeof(_real));
    pmin = (_real *) malloc(parametrized*sizeof(_real));
}

System::~System(){
/*	if(Q != NULL) free(Q);
	if(Qf != NULL) free(Qf);
    if(R != NULL) free(R);
    if(u_max != NULL) free(u_max);
    if(u_min != NULL) free(u_min);
    if(du_max != NULL) free(du_max);
    if(pmax != NULL) free(pmax);
    if(pmin != NULL) free(pmin);
	if(state_upper_limits != NULL) free(state_upper_limits);
	if(state_lower_limits != NULL) free(state_lower_limits);
	if(current_state != NULL) free(current_state);
	if(last_state != NULL) free(last_state);
    if(rising_time != NULL) free(rising_time);
 */
}

void System::load_configurations_from_file(std::string config_file){
	// Accessing config file
	std::ifstream sim_config;
    sim_config.open(config_file, std::ios::in);

    std::string error_message = name + " config file not found.";
    if(!sim_config.is_open())
    	throw std::runtime_error(error_message.c_str());
    
    N = read_int(&sim_config, "N");		// Prediction Horizon
    Nu = read_int(&sim_config, "Nu");	// Control Horizon
    Nx = read_int(&sim_config, "Nx");	// Number of States
    n_U = read_int(&sim_config, "n_U");	// Number of Inputs
    Ts = read_real(&sim_config, "Ts");	// Sampling Rate
    
    initializeVectors();

    // Defining Input Constrains
    read_real_vector(&sim_config, "u_max", u_max, n_U);
    read_real_vector(&sim_config, "u_min", u_min, n_U);
    read_real_vector(&sim_config, "du_max", du_max, n_U);      

    uss = read_real(&sim_config, "uss"); 
    
	// Definig State Constrains
    read_real_vector(&sim_config, "state_upper_limits", state_upper_limits, Nx);
    read_real_vector(&sim_config, "state_lower_limits", state_lower_limits, Nx);
    
	// Definig Cost Function Weight Matrices
    read_real_vector(&sim_config, "Q", Q, Nx);
    read_real_vector(&sim_config, "Qf", Qf, Nx);
    read_real_vector(&sim_config, "R", R, n_U);

	read_real_vector(&sim_config, "current_state", current_state, Nx);
	
	for (int i = 0; i < Nx; ++i) {
		last_state[i] = 0.0;	
	}
    
    parametrized = read_int(&sim_config, "Parametrization");
    if(parametrized) {
        initializeParametrizationVectors();
        Lambda = read_real(&sim_config, "Lambda");
        q_param = read_real(&sim_config, "q_param");
        read_real_vector(&sim_config, "pmax", pmax, parametrized);
        read_real_vector(&sim_config, "pmin", pmin, parametrized);
    }

    find_token(sim_config, "Rising_Time");
    filter_reference = read_int(&sim_config, "Rising_Time"); 
    
//    filter_reference = read_int(&sim_config, "Rising_Time");
    if(filter_reference){
        read_real_vector(&sim_config, "tr", rising_time, Nx);
    }

	sim_config.close();

    // Initialize Acceleration Constrains
    ddu_max = (_real *) malloc(n_U*sizeof(_real));
    acc_max = (_real *) malloc(n_U*sizeof(_real));
    acc_min = (_real *) malloc(n_U*sizeof(_real));
    
	state_type = (int *) malloc(Nx*sizeof(int));

	for (int i = 0; i < n_U; ++i) {
	    ddu_max[i] = du_max[i];
	    acc_max[i] = ddu_max[i];
	    acc_min[i] = -ddu_max[i];
	}
}

// Get/Set Methods
_real * System::getU_max(){return u_max;}
_real * System::getU_min(){return u_min;}
_real * System::getDU_max(){return du_max;}
int System::getN(){return N;}
int System::getNu(){return Nu;}
int System::getn_U(){return n_U;}

std::string System::getName(){
	return name;
}

_real * System::getStateUpperLimits() {
	return state_upper_limits;
}
_real * System::getStateLowerLimits() {
	return state_lower_limits;
}

void System::setN(int N){this->N = N;}
void System::setNu(int Nu){this->Nu = Nu;}
void System::setTs(_real Ts){this->Ts = Ts;}

int System::getNx(){return Nx;}
_real System::getTs(){return Ts;}
_real System::getState(int num){return current_state[num];}
void System::setState(_real * state){
	for (int i = 0; i < Nx; ++i){
		this->current_state[i] = state[i];
	}
}

int System::getParametrized(){
	return parametrized;
}

_real * System::get_pmin(){return pmin;}
_real * System::get_pmax(){return pmax;}

_real System::get_Q_index(int position){
	if(position >= Nx)
		return -1.0;
	else
		return Q[position];
}
_real System::get_Qf_index(int position){
	if(position >= Nx)
		return -1.0;
	else
		return Qf[position];
}
_real System::get_R_index(int position){
	if(position >= n_U)
		return -1.0;
	else
		return R[position];
}
void System::setCostMatrix(_real * Q, _real * Qf, _real * R){
	for (int i = 0; i < getNx(); ++i) {
		this->Q[i] = Q[i];
		this->Qf[i] = Qf[i];
	}
    for (int i = 0; i < getn_U(); ++i) {
        this->R[i] = R[i];
    } 
}
void System::setRMatrix(_real * R){
    for (int i = 0; i < getn_U(); ++i) {
        this->R[i] = R[i];
    }   
}

void System::setQfMatrix(_real * Qf){
	for (int i = 0; i < getNx(); ++i) {
		this->Qf[i] = Qf[i];
	}

}

int System::getFilterReference(){
    return filter_reference;
}

_real System::getRisingTime(int pos){
    return rising_time[pos];
}

int System::is_angle_state(int pos){
    return state_type[pos];
}

void System::updateAcceleration(_real * v_curr, _real * u_curr, _real * u_ant){
	// Update ddu_max
    for (int i = 0; i < n_U; ++i)
    {
        v_curr[i] = (u_curr[i] - u_ant[i]);
        acc_max[i] = v_curr[i] + ddu_max[i];
        if(acc_max[i] > du_max[i])
        	acc_max[i] = du_max[i];
        acc_min[i] = v_curr[i] - ddu_max[i];
		if(acc_min[i] < -du_max[i])
        	acc_min[i] = -du_max[i];
        //std::cout << "u_curr: " << u_curr[i] << " v_curr: " << v_curr[i] << " acc_max: " << acc_max[i] << " acc_min: " << acc_min[i] << std::endl;
    }
}


// 4th Order Runge-Kutta Method
void System::one_step_prediction(_real * state_plus, _real * state, _real * control){
	_real k1[Nx], k2[Nx], k3[Nx], k4[Nx];
	_real state_temp[Nx];

    //print_array("state_in", state, Nx, 0);
	model(k1,state,control);

    // print_array("k1", k1, Nx, 0);

	for (int i = 0; i < Nx; ++i) {
		state_temp[i] = state[i] + Ts*0.5*k1[i];
	}

	model(k2,state_temp,control);
    // print_array("k2", k2, Nx, 0);

	for (int i = 0; i < Nx; ++i) {
		state_temp[i] = state[i] + Ts*0.5*k2[i];
	}

	model(k3,state_temp,control);
    // print_array("k3", k3, Nx, 0);

	for (int i = 0; i < Nx; ++i) {
		state_temp[i] = state[i] + Ts*k3[i];
	}

	model(k4,state_temp,control);
    // print_array("k4", k4, Nx, 0);

	_real Ts_6 = Ts*0.1666666667; //Ts/6;

	for (int i = 0; i < Nx; ++i) {
		state_plus[i] = state[i] + Ts_6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
	}
	
}

_real System::error_with_constrains(int j, _real x) {
	const _real penality = 1e10;

	if(x > state_upper_limits[j]) {
		return (x-this->state_upper_limits[j])*penality;
	}
	else if (x < state_lower_limits[j]) {
		return (state_lower_limits[j]-x)*penality;
	}
	else {
		return 1.0;
	}
}

_real System::nmpc_cost_function(_real * control_guess, _real * xref, _real * uref, _real * xss, _real * uss){

	_real J = 0.0;
	_real Jf = 0.0;

#ifdef DEBUG2
	int debug = 1;
	float J_ant = 0;
#else
	int debug = 0;
	float J_ant = 0;
#endif

	// Constructing the vector of guessed control actions with respect to N and Nu
	//_real * uu;
	//uu = (_real *) malloc(N*sizeof(_real));
	_real uu[N];

	for (int i = 0; i < N; ++i) {				
		if(i < Nu) {
			uu[i] = control_guess[i];
		}else {
			uu[i] = control_guess[Nu-1];	
		}
	}
    if(debug) print_array("control_guess", uu, N, 0);

	// Initialize x_hat vector
	_real * x_hat; //[Nx*N];
	x_hat = (_real *) malloc(Nx*N*sizeof(_real));

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
	for (int i = 0; i < N-1; ++i) {
		k = Nx*(i+1);
        one_step_prediction(&x_hat[k], &x_hat[Nx*i], &uu[i]);
        if(debug) printf("=== %d ===\nJ = J + pow(Q[j]*( x_hat[l] - xref[l]),2);\n",i+1);
		for (int j = 0; j < Nx; ++j) {
			l = k + j;
			if(debug) J_ant = J;
			if((j == 2) && (x_hat[l] > M_PI)) {
				J = J + pow(Q[j]*error_with_constrains(j,x_hat[l])*((M_PI*2 - x_hat[l]) - xref[l]),2);
			}
			else {
				J = J + pow(Q[j]*error_with_constrains(j,x_hat[l])*( x_hat[l] - xref[l]),2);
			}
			if(debug) printf("l = %d | %f = %f + pow(%f * (%f - %f), 2)\n", l, J, J_ant, Q[j], x_hat[l], xref[l]);
		}
		if(debug) printf("J = J + R*pow(uu[i],2);\n");
		if(debug) J_ant = J;
		J = J + R[0]*pow(uu[i],2);
		if(debug) printf("%f = %f + %f * pow(%f)^2\n", J, J_ant, R[0], uu[i]);
	}
	k = Nx*(N-1);
	if(debug) printf("== END OF J \n");
	if(debug) printf("Jf = Jf + pow(( x_hat[l] - xref[l]),2);\n");
	for (int j = 0; j < Nx; ++j) {
		l = k + j;
		if(debug) J_ant = Jf;
		if((j == 2) && (x_hat[l] > M_PI)) {
			Jf = Jf + Qf[j]*error_with_constrains(j,x_hat[l])*pow(( (M_PI*2 - x_hat[l]) - xss[j]),2);
		}
		else {
			Jf = Jf + Qf[j]*error_with_constrains(j,x_hat[l])*pow(( x_hat[l] - xss[j]),2);
		}
		if(debug) printf("l = %d | %f = %f + %f * pow(%f - %f)^2\n", l, Jf, J_ant, Qf[j], x_hat[l], xref[l]);
	}
	J = J + Jf + R[0]*pow(uu[N-1],2);
	if(debug) printf("=== FINAL ===\n %f = %f + %f * pow(%f)^2\n", J, Jf, R[0], uu[N-1]);

	if(x_hat != NULL) free(x_hat);
//	if(uu != NULL) free(uu);
	return J;

}
