#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#define DEBUG_SYSTEM

#include "hls_inverted_pendulum.hpp"

template <
    typename _hw_real,
    unsigned _N,
    unsigned _Nx,
    unsigned _n_U,
    unsigned _Nu
>class System {
protected:
	// System Properties
	_hw_real current_state[_Nx];
	// _hw_real last_state[_Nx];
	
    // NMPC Properties
    // const unsigned N = _N;       // Prediction Horizon
    // const unsigned Nu = _Nu;     // Control Horizon
    // const unsigned n_U = _n_U;   // Number of Inputs
    // const unsigned Nx = _Nx;     // Number of States

    // Simulation Properties
    _hw_real Ts; // Sampling rate (sec)
    _hw_real Ts_2;
    _hw_real Ts_6;
    
    // Constrains
    _hw_real u_max[_n_U];
    _hw_real u_min[_n_U];
    _hw_real du_max[_n_U];

    _hw_real uss[_n_U]; // Control Signal at Steady State

    _hw_real controlled_state[_Nx];
    _hw_real state_upper_limits[_Nx];
    _hw_real state_lower_limits[_Nx];
    // Weight Matrices
    _hw_real Q[_Nx]; 
    _hw_real Qf[_Nx];
    _hw_real R[_n_U];

    // Temporary variables for cost functions computations
    _hw_real uu[_N];
    _hw_real x_hat[_Nx*_N];

    // Array to define each system state as an angle or not for relative MSE computing compensating for angles above 360º 
    int state_type[_Nx];

    // Acceleration Control
    _hw_real ddu_max[_n_U];
    _hw_real acc_max[_n_U];
    _hw_real acc_min[_n_U];

public:
    System(
        const _hw_real _u_max[_n_U],
        const _hw_real _u_min[_n_U],
        const _hw_real _du_max[_n_U],
        const unsigned char _controlled_state[_Nx],
        const _hw_real _state_upper_limits[_Nx],
        const _hw_real _state_lower_limits[_Nx],
        const _hw_real _Q[_Nx],
        const _hw_real _Qf[_Nx],
        // const _hw_real _current_state[_Nx],
        const _hw_real _R[_n_U],
        
        const _hw_real _uss[_n_U],
        const _hw_real _Ts
    ){
        Ts  = _Ts;
        Ts_2 = Ts*(_hw_real)0.5; //Ts/6;
        Ts_6 = Ts*(_hw_real)0.1666666667; //Ts/6;
        // uss = _uss;
        memcpy(uss,                 (_hw_real *)_uss, _n_U*sizeof(_hw_real));
        memcpy(u_max,               (_hw_real *)_u_max, _n_U*sizeof(_hw_real));
        memcpy(u_min,               (_hw_real *)_u_min, _n_U*sizeof(_hw_real));
        memcpy(du_max,              (_hw_real *)_du_max, _n_U*sizeof(_hw_real));
        memcpy(controlled_state,  (_hw_real *)_controlled_state, _Nx*sizeof(_hw_real));
        memcpy(state_upper_limits,  (_hw_real *)_state_upper_limits, _Nx*sizeof(_hw_real));
        memcpy(state_lower_limits,  (_hw_real *)_state_lower_limits, _Nx*sizeof(_hw_real));
        memcpy(Q,                   (_hw_real *)_Q, _Nx*sizeof(_hw_real));
        memcpy(Qf,                  (_hw_real *)_Qf, _Nx*sizeof(_hw_real));
        // memcpy(last_state,          (_hw_real *)_last_state, _Nx*sizeof(_hw_real));
        // memcpy(current_state,       (_hw_real *)_current_state, _Nx*sizeof(_hw_real));
        memcpy(R,                   (_hw_real *)_R, _n_U*sizeof(_hw_real));
        for (unsigned i = 0; i < _n_U; i++){
            ddu_max[i] = du_max[i];
            acc_max[i] = ddu_max[i];
            acc_min[i] = -ddu_max[i];
	    }
    };
// ---------------------------------------------------
	_hw_real nmpc_cost_function(
        _hw_real current_state[_Nx],

        _hw_real control_guess[_n_U*_Nu],
        _hw_real xref[_Nx*_N]
    ){
#pragma HLS interface ap_fifo port=control_guess
#pragma HLS interface ap_fifo port=xref

#pragma HLS interface ap_fifo port=uref
#pragma HLS interface ap_fifo port=xss
#pragma HLS interface ap_fifo port=uss

#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation

#pragma HLS ALLOCATION instances=one_step_prediction limit=1 function
        // Register inputs locally
        _hw_real local_xref[_Nx*_N];
        memcpy(local_xref, (_hw_real *)xref, _Nx*_N*sizeof(_hw_real));

        // Constructing the vector of guessed control actions with respect to N and Nu
        _hw_real uu[_n_U*_N];
        _hw_real last_control_guess[_n_U];
        unsigned k_cg = 0;
        const unsigned k_cg_last = _n_U*_Nu - _n_U;
        for (unsigned i = 0; i < _n_U*_Nu; ++i) {
            _hw_real control_guess_val = control_guess[i];
            uu[i] = control_guess_val;
            if (i == k_cg_last){
                last_control_guess[k_cg] = control_guess_val;
                k_cg = (k_cg<_n_U) ? k_cg+1 : 0;
            }
        }
#if _Nu < _N
        unsigned k_u = 0;
        for (unsigned i = _n_U*_Nu; i < _n_U*_N; ++i) {    
            uu[i] = last_control_guess[k_u];
            k_u = (k_u<_n_U) ? k_u+1 : 0;
        }
#endif

        // Initialize x_hat vector
        _hw_real x_hat[_Nx*_N]; //[Nx*N];
        for (unsigned i = 0; i < _Nx; ++i) {
            x_hat[i] = current_state[i];
        }
        for (unsigned i = _Nx; i < (_Nx*_N); ++i) {
            x_hat[i] = (_hw_real)0.0;
        }

        // Calculate State error
        _hw_real J = 0.0;
        _hw_real Ji[_Nx], Ji_buff[_Nx]; // = 0.0;
        _hw_real Jui[_n_U], Jui_buff[_n_U]; // = 0.0;
        _hw_real current_x_hat[_Nx];
        _hw_real current_uu[_n_U];
        _hw_real current_xref[_Nx];
        memset(Ji, (const _hw_real)0.0, _Nx*sizeof(_hw_real));
        memset(Jui, (const _hw_real)0.0, _n_U*sizeof(_hw_real));

        unsigned k = 0;
        unsigned l = 0;
        _hw_real x_buffer[_Nx];

#pragma HLS ARRAY_RESHAPE variable=x_hat block factor=2 dim=1

        one_step_error_loop: for (unsigned i = 0; i < _N-1; ++i) {
#ifdef DEBUG_SYSTEM
    std::cout << "Horizon[" << i << "]"<< std::endl;
#endif
            k = _Nx*(i+1);
#pragma HLS dataflow
            memcpy(current_x_hat, (_hw_real *)&x_hat[_Nx*i], _Nx*sizeof(_hw_real));
            memcpy(current_uu, (_hw_real *)&uu[_n_U*i], _n_U*sizeof(_hw_real));
            
            one_step_prediction(x_buffer, current_x_hat, current_uu);
            
            memcpy(&x_hat[k], (_hw_real *)x_buffer, _Nx*sizeof(_hw_real));
            memcpy(current_xref, (_hw_real *)&local_xref[_Nx*i], _Nx*sizeof(_hw_real));
            
            one_step_error(Ji_buff, Ji, x_buffer, current_xref);
            one_step_u_error(Jui_buff, Jui, current_uu, uss);

            memcpy(Ji, (_hw_real *)Ji_buff, _Nx*sizeof(_hw_real));
            memcpy(Jui, (_hw_real *)Jui_buff, _n_U*sizeof(_hw_real));
        }
        for (unsigned j = 0; j < _Nx; j++) {
            J = J + Q[j]*(Ji[j]);
        }
        for (unsigned j = 0; j < _n_U; j++) {
            J = J + R[j]*(Jui[j]);
        }

        // Calculate Final State Error
        _hw_real Jf[_Nx]; // = 0.0;
        memset(Jf, (const _hw_real)0.0, _Nx*sizeof(_hw_real));

        for (unsigned j = 0; j < _Nx; ++j) {
            _hw_real single_x_hat = x_buffer[j];
            _hw_real tmp_err = (single_x_hat - current_xref[j]);
            _hw_real penality = ((state_lower_limits[j] < single_x_hat) && (single_x_hat < state_upper_limits[j])) ? (_hw_real)1.0 : (_hw_real)1e10;
            Jf[j] = penality*tmp_err*tmp_err;
        }
        for (unsigned j = 0; j < _Nx; j++) {
            J = J + Qf[j]*(Jf[j]);
        }

        // Squared Error
        // if(R[0] > 0) {
        //     for (unsigned i = 1; i < _N; ++i) {  
        //         _hw_real tmp_err = uu[i]-uu[i-1];
        //         Ju = Ju + tmp_err*tmp_err;
        //     }
        // }
        // J = J + R[0]*Ju;

        return J;
    }

protected:
    void model(
        _hw_real state_dot[_Nx],
        _hw_real state[_Nx],
        _hw_real control[_n_U]
    ){
        model_inverted_pendulum<_hw_real, _Nx, _n_U>(state_dot, state, control);
    }

    void one_step_error(
        _hw_real Ji_out[_Nx],
        _hw_real Ji_in[_Nx],
        _hw_real x_hat[_Nx],
        _hw_real xref[_Nx]
    ){
#pragma HLS inline
        _hw_real penality;
        // _hw_real Ji_local[_Nx];
        for (unsigned j = 0; j < _Nx; ++j) {
            //_hw_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
            _hw_real current_x_hat = x_hat[j];
            _hw_real tmp_err = (controlled_state[j] == 1) ? (current_x_hat - xref[j]) : (_hw_real)0.0;
            penality = ((state_lower_limits[j] < current_x_hat) && (current_x_hat < state_upper_limits[j])) ? (_hw_real)1.0 : (_hw_real)1e10;
            Ji_out[j] = Ji_in[j] + penality*tmp_err*tmp_err;
        }
    }

    void one_step_u_error(
        _hw_real Ji_out[_n_U],
        _hw_real Ji_in[_n_U],
        _hw_real uu[_n_U],
        _hw_real uref[_n_U]
    ){
#pragma HLS inline
        _hw_real penality;
        // _hw_real Ji_local[_n_U];
        for (unsigned j = 0; j < _n_U; ++j) {
            //_hw_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
            _hw_real current_uu = uu[j];
            _hw_real tmp_err = (current_uu - uref[j]);
            penality = ((u_min[j] < current_uu) && (current_uu < u_max[j])) ? (_hw_real)1.0 : (_hw_real)1e10;
            Ji_out[j] = Ji_in[j] + penality*tmp_err*tmp_err;
        }
    }

	void one_step_prediction(
        _hw_real state_plus[_Nx],
        _hw_real state[_Nx],
        _hw_real control[_n_U]
    ){
#pragma HLS interface ap_fifo  port=state_plus
#pragma HLS interface ap_fifo  port=state
#pragma HLS interface ap_fifo  port=control

#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation

#pragma HLS ALLOCATION instances=model_inverted_pendulum limit=1 function

        _hw_real local_control[_n_U];
        _hw_real local_state[_Nx];
        memcpy(local_control, (_hw_real *)control, _n_U*sizeof(_hw_real));
        memcpy(local_state, (_hw_real *)state, _Nx*sizeof(_hw_real));

        _hw_real k1[_Nx], k2[_Nx], k3[_Nx], k4[_Nx];
        _hw_real state_temp1[_Nx];
        _hw_real state_temp2[_Nx];
        _hw_real state_temp3[_Nx];

        model(k1,local_state,local_control);
        for (unsigned i = 0; i < _Nx; ++i) {
            state_temp1[i] = local_state[i] + Ts_2*k1[i];
        }

        model(k2,state_temp1,local_control);
        for (unsigned i = 0; i < _Nx; ++i) {
            state_temp2[i] = local_state[i] + Ts_2*k2[i];
        }

        model(k3,state_temp2,local_control);
        for (unsigned i = 0; i < _Nx; ++i) {
            state_temp3[i] = local_state[i] + Ts*k3[i];
        }

        model(k4,state_temp3,local_control);
        for (unsigned i = 0; i < _Nx; ++i) {
            state_plus[i] = local_state[i] + Ts_6*(k1[i] + (_hw_real)2.0*k2[i] + (_hw_real)2.0*k3[i] + k4[i]);
        }
    };


/*
    void updateAcceleration(
        _hw_real v_curr[_n_U],
        _hw_real u_curr[_n_U],
        _hw_real u_ant[_n_U]
    );

    // Get/Set Methods
    void getState(
        _hw_real state[_Nx]
    );
    void setState(
        const _hw_real state[_Nx]
    );
    int is_angle_state(
        int pos
    );
*/

};

#endif
