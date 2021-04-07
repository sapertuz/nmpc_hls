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

public:
    virtual ~System();
	void one_step_prediction(
        _real state_plus[_Nx], 
        _real state[_Nx], 
        _real control[_n_U]
    );

	virtual _real nmpc_cost_function(
        _real * control_guess,
        _real * xref, 
        _real * uref, 
        _real * xss, 
        _real * uss
    ) = 0;

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
