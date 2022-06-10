#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <string>
#include <fstream>
#include "config.hpp"

class System {

protected:
//	std::string name;
	// System Properties
	_real current_state[_Nx];
	_real last_state[_Nx];
	// NMPC Properties
    const int N = _N;       // Prediction Horizon
    const int Nu = _Nu;     // Control Horizon
    const int n_U = _n_U;   // Number of Inputs
    const int Nx = _Nx;     // Number of States
    
    const int controlled_state[12] = {1,1,1,1,1,1,0,0,0,0,0,0};
    
    // Simulation Properties
    _real Ts = _Ts; // Sampling rate (sec)
    // Constrains
    _real u_max[_n_U];
    _real u_min[_n_U];
    _real du_max[_n_U];

    _real uss; // Control Signal at Steady State

    _real state_upper_limits[_Nx];
    _real state_lower_limits[_Nx];
    // Weight Matrices
    _real Q[_Nx]; 
    _real Qf[_Nx];
    _real R[_n_U];

    // Temporary variables for cost functions computations
    _real uu[_N];
    _real x_hat[_Nx*_N];

    // Parametrization variables
    int parametrized;
    _real Lambda;
    _real q_param;
    _real pmax[_Parametrization];
    _real pmin[_Parametrization];

    // Array to define each system state as an angle or not for relative MSE computing compensating for angles above 360º 
    int state_type[_Nx];

public:
    // Acceleration Control
    _real ddu_max[_n_U];
    _real acc_max[_n_U];
    _real acc_min[_n_U];

    // Rising time for Reference Filtering
    int filter_reference;
    _real rising_time[_Nx];

protected:
    void initializeVectors();
    void initializeParametrizationVectors();
    void load_configurations_from_file();
    _real error_with_constrains(
        int j, 
        _real x
    );
    virtual void model(
        _real state_dot[_Nx], 
        _real state[_Nx],
        _real control[_n_U]
    ) = 0;
void update_state(
	_real *state_plus, 
	_real *state, 
	_real *k, 
	_real Ts_local
	);
void update_state_final(
	_real *state_plus, 
	_real *state, 
	_real *k1, 
	_real *k2, 
	_real *k3, 
	_real *k4
	);

void one_step_error(
	_real Ji_out[_Nx],
	volatile _real Ji_in[_Nx],
	volatile _real x_hat[_Nx],
	volatile _real xref[_Nx]
);
void one_step_u_error(
	_real *Ji_out,
	volatile _real *Ji_in,
	volatile _real *uu
	// _hw_real uref[_n_U]
);

void uu_loop(
	volatile _real *control_guess, 
	_real *uu_1,
	_real *uu_2
	);
void horizon_step_prediction(
	_hw_real *x_initial,
	_hw_real *uu, 
	_hw_real *x_horizon
	);
void horizon_step_error(
	_real *x_horizon, 
	_real *xref, 
	_real *Ji,
	_real *final_x,
	_real *final_xref
	);
void horizon_u_error(
	_real *uu,
	_real *Jui
    );

void Jf_error(
	volatile _real *final_x,
	volatile _real *final_xref,
	_real *J
	);
void Ji_error(
	volatile _real *Ji, 
	volatile _real *Jui,
	_real *J
	);
void final_sum(
	_real *J, 
	_real J_tmp1, 
	_real J_tmp2
	);

public:
    virtual ~System();
	void one_step_prediction(
        _real state_plus[_Nx], 
        _real state[_Nx], 
        _real control[_n_U]
    );

	// virtual 
    _real nmpc_cost_function(
        _real * control_guess,
        _real * xref, 
        _real * uref, 
        _real * xss, 
        _real * uss
    ) 
    // = 0
    ;

    void updateAcceleration(
        _real v_curr[_n_U],
        _real u_curr[_n_U],
        _real u_ant[_n_U]
    );

    // Get/Set Methods
    void getState(
        _real state[_Nx]
    );
    void setState(
        const _real state[_Nx]
    );
    int is_angle_state(
        int pos
    );
};

#endif
