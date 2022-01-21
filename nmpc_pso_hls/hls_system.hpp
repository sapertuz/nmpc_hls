#ifndef SYSTEM_HPP
#define SYSTEM_HPP

// #define DEBUG_SYSTEM
#define def_N 3

//#include "hls_inverted_pendulum.hpp"
#include "aux_functions.hpp"

template <
    class _hw_real,
    class _model_t,
    unsigned _N,
    unsigned _Nx,
    unsigned _n_U,
    unsigned _Nu
>class System {
protected:
	// System Properties
	// _hw_real *current_state;//[_Nx];
	// _hw_real last_state;//[_Nx];
	
    // NMPC Properties
    // const unsigned N = _N;       // Prediction Horizon
    // const unsigned Nu = _Nu;     // Control Horizon
    // const unsigned n_U = _n_U;   // Number of Inputs
    // const unsigned Nx = _Nx;     // Number of States

    // Simulation Properties
    const _hw_real Ts; // Sampling rate (sec)
    const _hw_real Ts_2;
    const _hw_real Ts_6;
    
    // Constrains
    const _hw_real *u_max;//[_n_U];
    const _hw_real *u_min;//[_n_U];
    const _hw_real *du_max;//[_n_U];

    const _hw_real *uss;//[_n_U]; // Control Signal at Steady State

    const unsigned short *controlled_state;//[_Nx];
    const _hw_real *state_upper_limits;//[_Nx];
    const _hw_real *state_lower_limits;//[_Nx];
    // Weight Matrices
    const _hw_real *Q;//[_Nx]; 
    const _hw_real *Qf;//[_Nx];
    const _hw_real *R;//[_n_U];

    const unsigned N;
    // Temporary variables for cost functions computations
    //const _hw_real *uu;//[_N];
    //_hw_real *x_hat;//[_Nx*_N];

    // Array to define each system state as an angle or not for relative MSE computing compensating for angles above 360º 
    //int *state_type;//[_Nx];

    // Acceleration Control
    // _hw_real ddu_max[_n_U];
    // _hw_real acc_max[_n_U];
    // _hw_real acc_min[_n_U];
    _model_t model_ptr;
public:
    constexpr System(
        const _hw_real _u_max[_n_U],
        const _hw_real _u_min[_n_U],
        const _hw_real _du_max[_n_U],
        const unsigned short _controlled_state[_Nx],
        const _hw_real _state_upper_limits[_Nx],
        const _hw_real _state_lower_limits[_Nx],
        const _hw_real _Q[_Nx],
        const _hw_real _Qf[_Nx],
        const _hw_real _R[_n_U],
        const _hw_real _uss[_n_U],
        const _hw_real _Ts,

        _model_t _model_ptr
    ) : Ts(_Ts), Ts_2(_Ts*(_hw_real)0.5), Ts_6(_Ts*(_hw_real)0.1666666667),
        uss(_uss), u_max(_u_max), u_min(_u_min), du_max(_du_max),
        controlled_state(_controlled_state), 
        state_upper_limits(_state_upper_limits), state_lower_limits(_state_lower_limits),
        Q(_Q), Qf(_Qf), R(_R),
		N(_N),
        model_ptr(_model_ptr)
    {
        // Ts  = _Ts;
        // Ts_2 = Ts*(_hw_real)0.5; //Ts/6;
        // Ts_6 = Ts*(_hw_real)0.1666666667; //Ts/6;
        // memcpy_loop<_hw_real, _hw_real, _n_U>(uss,                 (const _hw_real *)_uss);
        // memcpy_loop<_hw_real, _hw_real, _n_U>(u_max,               (const _hw_real *)_u_max);
        // memcpy_loop<_hw_real, _hw_real, _n_U>(u_min,               (const _hw_real *)_u_min);
        // memcpy_loop<_hw_real, _hw_real, _n_U>(du_max,              (const _hw_real *)_du_max);
        // memcpy_loop<unsigned short, _Nx>(controlled_state,  (const unsigned short *)_controlled_state);
        // memcpy_loop<_hw_real, _hw_real, _Nx>(state_upper_limits,  (const _hw_real *)_state_upper_limits);
        // memcpy_loop<_hw_real, _hw_real, _Nx>(state_lower_limits,  (const _hw_real *)_state_lower_limits);
        // memcpy_loop<_hw_real, _hw_real, _Nx>(Q,                   (const _hw_real *)_Q);
        // memcpy_loop<_hw_real, _hw_real, _Nx>(Qf,                  (const _hw_real *)_Qf);
        // memcpy_loop<_hw_real, _hw_real, _n_U>(R,                  (const _hw_real *)_R);
    };
// ---------------------------------------------------
	_hw_real nmpc_cost_function(
        _hw_real current_state[_Nx],

        _hw_real control_guess[_n_U*_Nu],
        _hw_real xref[_Nx*_N]
    ){
// #pragma HLS interface ap_fifo port=current_state
// #pragma HLS interface ap_fifo port=control_guess
// #pragma HLS interface ap_fifo port=xref

// #pragma HLS ALLOCATION instances=hmul limit=1 operation
// #pragma HLS ALLOCATION instances=hadd limit=1 operation
// #pragma HLS ALLOCATION instances=hsub limit=1 operation

#pragma HLS ALLOCATION instances=one_step_prediction limit=1 function

        // Register inputs locally
        // _hw_real local_xref[_Nx*_N];
        // reg_xref: memcpy_loop<_hw_real, _hw_real,_Nx*_N>(local_xref, (const _hw_real *)xref);

        // Constructing the vector of guessed control actions with respect to N and Nu
        _hw_real uu[_n_U*_N];
        _hw_real last_control_guess[_n_U];
        unsigned k_cg = 0;
        unsigned k_u = 0;
        const unsigned k_cg_last = _n_U*_Nu - _n_U;
        /*
        uu_loop1: for (unsigned short i = 0; i < _n_U*_Nu; i = (i+1)*_n_U )
        {
            memcpy_loop<_hw_real, _hw_real, _n_U>(last_control_guess,&control_guess[i]);
            memcpy_loop<_hw_real, _hw_real, _n_U>(&uu[i],last_control_guess);
        }
        uu_loop2: for (unsigned short i = _n_U*_Nu; i < _n_U*_N; i = (i+1)*_n_U)
        {
            memcpy_loop<_hw_real, _hw_real, _n_U>(&uu[i],last_control_guess);
        }
        */      
        
        uu_loop: for (unsigned i = 0; i < _n_U*_N; ++i) {
            _hw_real control_val;
            if(i < _n_U*_Nu){
                control_val = control_guess[i];
                if (i >= k_cg_last){
                    last_control_guess[k_cg] = control_val;
                    k_cg++;
                }
            }else{
                control_val = last_control_guess[k_u];
                k_u = (k_u<_n_U) ? k_u+1 : 0;
            }
            uu[i] = control_val;
        }
        
        /*
        uu_loop1: for (unsigned i = 0; i < _n_U*_Nu; ++i) {
            _hw_real control_val;
            control_val = control_guess[i];
            if (i >= k_cg_last){
                last_control_guess[k_cg] = control_val;
                k_cg++;
            }
            uu[i] = control_val;
        }
        uu_loop2: for (unsigned i = _n_U*_Nu; i < _n_U*_N; ++i) {
            _hw_real control_val;
            control_val = last_control_guess[k_u];
            k_u = (k_u<_n_U) ? k_u+1 : 0;

            uu[i] = control_val;
        }
        */

        // Initialize x_hat vector
        _hw_real x_hat[_Nx*_N]; //[Nx*N];

        // _hw_real x_value;
        // for (unsigned i = 0; i < _Nx*_N; ++i) {
        //     if (i < _Nx){
        //         x_value = current_state[i];
        //     }else{
        //         x_value = (_hw_real)0.0;
        //     }
        //     x_hat[i] = x_value;
        // }

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
        _hw_real prev_x_hat[_Nx]; //[Nx*N];
        reg_prev_x: memcpy_loop<_hw_real, _hw_real,_Nx>(prev_x_hat, (const _hw_real *)current_state);
        // memcpy_loop<_hw_real, _hw_real,>(&x_hat[0], prev_x_hat, _Nx*sizeof(_hw_real));

        // one_step_error_loop: for (unsigned i = 0; i < _N-1; ++i) {
        one_step_error_loop: for (unsigned i = 0; i < N; ++i) {
            // k = _Nx*(i+1);
            // k = _Nx*(i);
// #pragma HLS array_partition variable=x_hat block factor=_N
// #pragma HLS array_partition variable=uu block factor=_N
#pragma HLS array_partition variable=xref block factor=_N
#pragma HLS array_partition variable=control_guess block factor=_Nu
#pragma HLS dataflow

//#pragma HLS dependence variable=prev_x_hat inter RAW distance=1 true
//#pragma HLS dependence variable=current_uu inter RAW distance=1 true
//#pragma HLS dependence variable=x_buffer inter RAW distance=2 true
//#pragma HLS dependence variable=Ji inter RAW distance=1 true
//#pragma HLS dependence variable=Ji_buff inter RAW distance=1 true
//#pragma HLS dependence variable=Jui inter RAW distance=1 true
//#pragma HLS dependence variable=Jui_buff inter RAW distance=1 true

        	// memcpy_loop<_hw_real, _hw_real,>(current_x_hat, (const _hw_real *)&x_hat[_Nx*i], _Nx*sizeof(_hw_real));
			// memcpy_loop<_hw_real, _hw_real, _Nx>(current_x_hat, (const _hw_real *)prev_x_hat);
        	memcpy_loop_enclosed<_hw_real, _hw_real, _n_U>(current_uu, (const _hw_real *)&uu[_n_U*i]);
        	memcpy_loop_enclosed<_hw_real, _hw_real, _Nx>(current_xref, (const _hw_real *)&xref[_Nx*i]);

			// _hw_real Jui_buff[_n_U]; // = 0.0;
			one_step_u_error(Jui_buff, Jui, current_uu);
			memcpy_loop_enclosed<_hw_real, _hw_real, _n_U>(Jui, (const _hw_real *)Jui_buff);

			// _hw_real Ji_buff[_Nx]; // = 0.0;
			one_step_error(Ji_buff, Ji, prev_x_hat, current_xref);
			memcpy_loop_enclosed<_hw_real, _hw_real, _Nx>(Ji, (const _hw_real *)Ji_buff);

			// one_step_prediction(x_buffer, current_x_hat, current_uu);
			one_step_prediction(x_buffer, prev_x_hat, current_uu);

			//memcpy_loop<_hw_real, _hw_real, _Nx>(&x_hat[_Nx*(i)], (const _hw_real *)x_buffer);

			memcpy_loop_enclosed<_hw_real, _hw_real, _Nx>(prev_x_hat, (const _hw_real *)x_buffer);
            
#ifdef DEBUG_SYSTEM
    std::cout << "Horizon[" << i << "]"<< std::endl;
    std::cout << "\t State"; print_formatted_float_array(prev_x_hat, _Nx, 2, 6); 
    std::cout << "\t Control"; print_formatted_float_array(current_uu, _n_U, 2, 6); std::cout << std::endl;
    std::cout << "\t Set  "; print_formatted_float_array(current_xref, _Nx, 2, 6); std::cout << std::endl;
#endif
        }
        Ji_loop: for (unsigned j = 0; j < _Nx; j++) {
            J = J + Q[j]*(Ji[j]);
        }
        Jui_loop: for (unsigned j = 0; j < _n_U; j++) {
            J = J + R[j]*(Jui[j]);
        }

        // Calculate Final State Error
        _hw_real Jf[_Nx]; // = 0.0;
        memset(Jf, (const _hw_real)0.0, _Nx*sizeof(_hw_real));

        Jf_error_loop: for (unsigned j = 0; j < _Nx; ++j) {
            _hw_real single_x_hat = x_buffer[j];
            _hw_real tmp_err = (single_x_hat - current_xref[j]);
            _hw_real penality = ((state_lower_limits[j] <= single_x_hat) && (single_x_hat <= state_upper_limits[j])) ? (_hw_real)1.0 : (_hw_real)1e10;
            Jf[j] = penality*tmp_err*tmp_err;
        }
        Jf_loop: for (unsigned j = 0; j < _Nx; j++) {
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
    void memcpy_loop_checklastu(
        unsigned _i,
        _hw_real _current_uu[_n_U], 
        _hw_real _control_guess[_n_U]
    ){
#pragma HLS inline
        if (_n_U*_i < _n_U*_Nu){
            memcpy_loop<_hw_real, _hw_real, _n_U>(_current_uu, (const _hw_real *)_control_guess);
        }else{
            memcpy_loop<_hw_real, _hw_real, _n_U>(_current_uu, (const _hw_real *)_current_uu);
        }
    }

    void model(
        _hw_real state_dot[_Nx],
        _hw_real state[_Nx],
        _hw_real control[_n_U]
    ){
        //model_inverted_pendulum<_hw_real, _Nx, _n_U>(state_dot, state, control);
        model_ptr.model(state_dot, state, control);
    }

    void one_step_error(
        _hw_real Ji_out[_Nx],
        _hw_real Ji_in[_Nx],
        _hw_real x_hat[_Nx],
        _hw_real xref[_Nx]
    ){
// #pragma HLS inline
        _hw_real penality;
        // _hw_real Ji_local[_Nx];
        for (unsigned j = 0; j < _Nx; ++j) {
            //_hw_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
            _hw_real current_x_hat = x_hat[j];
            _hw_real tmp_err = (controlled_state[j] == 1) ? (current_x_hat - xref[j]) : (_hw_real)0.0;
            penality = ((state_lower_limits[j] <= current_x_hat) && (current_x_hat <= state_upper_limits[j])) ? (_hw_real)1.0 : (_hw_real)1e10;
            Ji_out[j] = Ji_in[j] + penality*tmp_err*tmp_err;
        }
    }

    void one_step_u_error(
        _hw_real Ji_out[_n_U],
        _hw_real Ji_in[_n_U],
        _hw_real uu[_n_U]
        // _hw_real uref[_n_U]
    ){
// #pragma HLS inline
        _hw_real penality;
        // _hw_real Ji_local[_n_U];
        for (unsigned j = 0; j < _n_U; ++j) {
            //_hw_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
            _hw_real current_uu = uu[j];
            _hw_real tmp_err = (current_uu - uss[j]);
            penality = ((u_min[j] <= current_uu) && (current_uu <= u_max[j])) ? (_hw_real)1.0 : (_hw_real)1e10;
            Ji_out[j] = Ji_in[j] + penality*tmp_err*tmp_err;
        }
    }

	void one_step_prediction(
        _hw_real state_plus[_Nx],
        _hw_real state[_Nx],
        _hw_real control[_n_U]
    ){
// #pragma HLS interface ap_fifo  port=state_plus
// #pragma HLS interface ap_fifo  port=state
// #pragma HLS interface ap_fifo  port=control

//#pragma HLS inline

#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation

#pragma HLS ALLOCATION instances=model_inverted_pendulum limit=1 function

        _hw_real local_control[_n_U];
        _hw_real local_state[_Nx];
        memcpy_loop<_hw_real, _hw_real, _n_U>(local_control, (const _hw_real *)control);
        memcpy_loop<_hw_real, _hw_real, _Nx>(local_state, (const _hw_real *)state);

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

//    void update_model(_hw_real state_out[_Nx], _hw_real state_prevp[], _hw_real )
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
