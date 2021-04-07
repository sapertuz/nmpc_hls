#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <string>
#include <fstream>
#include "config.h"

class System {

protected:
	std::string name;
	// System Properties
	_real * current_state;
	_real * last_state;
	// NMPC Properties
    int N;   // Prediction Horizon
    int Nu;  // Control Horizon
    int n_U; // Number of Inputs
    int Nx;  // Number of States
    // Simulation Properties
    _real Ts; // Sampling rate (sec)
    // Constrains
    _real * u_max;
    _real * u_min;
    _real * du_max;

    _real uss; // Control Signal at Steady State

    _real * state_upper_limits;
    _real * state_lower_limits;
    // Weight Matrices
    _real * Q; 
    _real * Qf;
    _real * R;

    // Temporary variables for cost functions computations
    _real * uu;
    _real * x_hat;

    // Parametrization variables
    int parametrized;
    _real Lambda;
    _real q_param;
    _real * pmax;
    _real * pmin;

    // Array to define each system state as an angle or not for relative MSE computing compensating for angles above 360º 
    int * state_type;
public:
    // Acceleration Control
    _real * ddu_max;
    _real * acc_max;
    _real * acc_min;

    // Plot Configuration
    int plot_number_of_states;
    int plot_number_of_controls;
    int * plot_states_config;
    int * plot_control_config;
    std::string * plot_labels;

    // Rising time for Reference Filtering
    int filter_reference;
    _real * rising_time;

    // This variable is on the first iteration of the PSO to inform the hardware (VHDL) modules of reference variables just once
    int pso_init;
    int fpga = 0;

protected:
    void initializeVectors();
    void initializeParametrizationVectors();
    void load_configurations_from_file(std::string config_file);
    _real error_with_constrains(int j, _real x);
    virtual void model(_real * state_dot, _real * state, _real * control) = 0;

public:
    virtual ~System();
    virtual void control_from_parameters(_real * parameters, _real * previous_control, _real * control, int horizon) = 0;
	void one_step_prediction(_real * state_plus, _real * state, _real * control);
	virtual _real nmpc_cost_function(_real * control_guess, _real * xref, _real * uref, _real * xss, _real * uss);

    void updateAcceleration(_real * v_curr, _real * u_curr, _real * u_ant);

    // Get/Set Methods
    _real * getU_max();
    _real * getU_min();
    _real * getDU_max();
    _real * getStateUpperLimits();
    _real * getStateLowerLimits();
    std::string getName();

    int getN();
    int getNu();
    int getn_U();
    void setN(int N);
    void setNu(int Nu);
    void setTs(_real Ts);
    int getNx();
    _real getTs();
    _real getState(int num);

    int getParametrized();
    _real * get_pmin();
    _real * get_pmax();


    _real get_Q_index(int position);
    _real get_Qf_index(int position);
    _real get_R_index(int position);
    void setState(_real * state);
    void setCostMatrix(_real * Q, _real * Qf, _real * R);
    void setQfMatrix(_real * Qf);
    void setRMatrix(_real * R);

    int getFilterReference();
    _real getRisingTime(int pos);

    int is_angle_state(int pos);
};

#endif
