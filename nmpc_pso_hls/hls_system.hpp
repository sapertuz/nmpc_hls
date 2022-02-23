#ifndef SYSTEM_HPP
#define SYSTEM_HPP

// #ifdef __SYNTHESIS__
// #include "ap_int.h"
// #include "hls_stream.h"
// #endif

#include "aux_functions.hpp"
//#include "hls_inverted_pendulum.hpp"

// #define DEBUG_SYSTEM
#define def_N 3

// #ifndef __VITIS_HLS__
// #define __VITIS_HLS__
// #endif

template <
    class _hw_real,
    class _model_t,
    unsigned _Nh,
    unsigned _Nx,
    unsigned _n_U,
    unsigned _Nu
>class System {
protected:
	// System Properties
	// _hw_real *current_state;//[_Nx];
	// _hw_real last_state;//[_Nx];
	
    // NMPC Properties
    // const unsigned N = _Nh;       // Prediction Horizon
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
    //const _hw_real *uu;//[_Nh];
    //_hw_real *x_hat;//[_Nx*_Nh];

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
		N(_Nh),
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
	void nmpc_cost_function(
        _hw_real current_state[_Nx],

        _hw_real control_guess[_n_U*_Nu],
        _hw_real xref[_Nx*_Nh],
        _hw_real *J
    ){
// #pragma HLS interface ap_fifo port=current_state
// #pragma HLS interface ap_fifo port=control_guess
// #pragma HLS interface ap_fifo port=xref
// #pragma HLS interface ap_vld port=J

//#pragma HLS INLINE
#pragma HLS DATAFLOW
//         _hw_real uu[_n_U*_Nh];
// #pragma HLS STREAM variable=uu depth=10
// //#pragma HLS data_pack variable=uu
// #pragma HLS RESOURCE variable=uu core=FIFO_BRAM

    _hw_real uu_1[_n_U*_Nh], uu_2[_n_U*_Nh];       
#pragma HLS stream variable=uu_1 depth=_n_U type=pipo
#pragma HLS stream variable=uu_2 depth=_n_U type=pipo
#if defined(__VITIS_HLS__)
#pragma HLS bind_storage variable=uu_1 type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=uu_2 type=FIFO impl=LUTRAM
#else
#pragma HLS data_pack variable=uu_1
#pragma HLS RESOURCE variable=uu_1 core=FIFO_LUTRAM
#pragma HLS data_pack variable=uu_2
#pragma HLS RESOURCE variable=uu_2 core=FIFO_LUTRAM
#endif


        _hw_real x_horizon[_Nx*_Nh]; //[Nx*N];
#pragma HLS STREAM variable=x_horizon depth=_Nx type=pipo
#if defined(__VITIS_HLS__)
#pragma HLS bind_storage variable=x_horizon type=FIFO impl=LUTRAM
#else
#pragma HLS data_pack variable=x_horizon
#pragma HLS RESOURCE variable=x_horizon core=FIFO_LUTRAM
#endif
        _hw_real J_tmp1, J_tmp2;
        // _hw_real J = 0.0;
        _hw_real Ji[_Nx]; // = 0.0;
#pragma HLS STREAM variable=Ji depth=_Nx
#if defined(__VITIS_HLS__)
#pragma HLS bind_storage variable=Ji type=FIFO impl=LUTRAM
#else
#pragma HLS data_pack variable=Ji
#pragma HLS RESOURCE variable=Ji core=FIFO_LUTRAM
#endif
        _hw_real Jui[_n_U]; // = 0.0;
#pragma HLS STREAM variable=Jui depth=_n_U
#if defined(__VITIS_HLS__)
#pragma HLS bind_storage variable=Jui type=FIFO impl=LUTRAM
#else
#pragma HLS RESOURCE variable=Jui core=FIFO_LUTRAM
#pragma HLS data_pack variable=Jui
#endif
        _hw_real final_x[_Nx], final_xref[_Nx];
#pragma HLS STREAM variable=final_x depth=_Nx
#pragma HLS STREAM variable=final_xref depth=_Nx
#if defined(__VITIS_HLS__)
#pragma HLS bind_storage variable=final_x type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=final_xref type=FIFO impl=LUTRAM
#else
#pragma HLS data_pack variable=final_x
#pragma HLS RESOURCE variable=final_x core=FIFO_LUTRAM
#pragma HLS data_pack variable=final_xref
#pragma HLS RESOURCE variable=final_xref core=FIFO_LUTRAM
#endif
        //memset(Jui, (const _hw_real)0.0, _n_U*sizeof(_hw_real));

        // Create control guess for full horizon
        // Constructing the vector of guessed control actions with respect to N and Nu
        uu_loop(control_guess, uu_1, uu_2);
        // Control Error
        horizon_u_error(uu_1, Jui);
        // Predict Horizon
        horizon_step_prediction(current_state, uu_2, x_horizon);
        // State Error
        horizon_step_error(x_horizon, xref, Ji, final_x, final_xref);
        // memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(final_x, (const _hw_real *)&x_horizon[(N-1)*_Nx]);
        // memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(final_xref, (const _hw_real *)&xref[(N-1)*_Nx]);
        // Weigthed Errors
        Ji_error(Ji, Jui, &J_tmp1);
        Jf_error(final_x, final_xref, &J_tmp2);
#ifdef DEBUG_SYSTEM
        std::cout << "Ji = " << J_tmp1 << std::endl;
        std::cout << "Jf = " << J_tmp2 << std::endl;
#endif
        final_sum(J, J_tmp1, J_tmp2);
        // return J;
    }

protected:
    void final_sum(_hw_real *J, _hw_real J_tmp1, _hw_real J_tmp2){
        J[0] = J_tmp1 + J_tmp2;
    }
    void memcpy_loop_checklastu(
        unsigned _i,
        _hw_real _current_uu[_n_U], 
        _hw_real _control_guess[_n_U]
    ){
#pragma HLS inline
        if (_n_U*_i < _n_U*_Nu){
            memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(_current_uu, (const _hw_real *)_control_guess);
        }else{
            memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(_current_uu, (const _hw_real *)_current_uu);
        }
    }

    void model(
        _hw_real state_dot[_Nx],
        _hw_real state[_Nx],
        _hw_real control[_n_U]
    ){
#pragma HLS inline
        //model_inverted_pendulum<_hw_real, _Nx, _n_U>(state_dot, state, control);
        model_ptr.model(state_dot, state, control);
    }

    void one_step_error(
        _hw_real Ji_out[_Nx],
        _hw_real Ji_in[_Nx],
        _hw_real x_hat[_Nx],
        _hw_real xref[_Nx]
    ){
#pragma HLS INLINE2048
        _hw_real penality;
        _hw_real current_x_hat;
        _hw_real Ji_local[_Nx];
        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Ji_local, Ji_in);
        for (unsigned j = 0; j < _Nx; ++j) {
            //_hw_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
            current_x_hat = x_hat[j];
            _hw_real tmp_err = (controlled_state[j] == 1) ? (current_x_hat - xref[j]) : (_hw_real)0.0;
            penality = ((state_lower_limits[j] <= current_x_hat) && (current_x_hat <= state_upper_limits[j])) ? (_hw_real)1.0 : (_hw_real)1e10;
            Ji_out[j] = Ji_local[j] + penality*tmp_err*tmp_err;
        }
    }

    void one_step_u_error(
        _hw_real Ji_out[_n_U],
        _hw_real Ji_in[_n_U],
        _hw_real uu[_n_U]
        // _hw_real uref[_n_U]
    ){
#pragma HLS INLINE
        _hw_real penality;
        _hw_real Ji_local[_n_U];
        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Ji_local, Ji_in);
        for (unsigned j = 0; j < _n_U; ++j) {
            //_hw_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
            _hw_real current_uu = uu[j];
            _hw_real tmp_err = (current_uu - uss[j]);
            penality = ((u_min[j] <= current_uu) && (current_uu <= u_max[j])) ? (_hw_real)1.0 : (_hw_real)1e10;
            Ji_out[j] = Ji_local[j] + penality*tmp_err*tmp_err;
        }
    }

	void one_step_prediction(
        _hw_real state_plus[_Nx],
        _hw_real state[_Nx],
        _hw_real control[_n_U]
    ){
#pragma HLS INLINE

#pragma HLS ALLOCATION operation instances=hmul limit=1
#pragma HLS ALLOCATION operation instances=hadd limit=1
#pragma HLS ALLOCATION operation instances=hsub limit=1

#pragma HLS ALLOCATION function instances=model limit=1

        _hw_real local_control[_n_U];
        _hw_real local_state[_Nx];
        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(local_control, (const _hw_real *)control);
        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(local_state, (const _hw_real *)state);

#if defined(__VITIS_HLS__)
#pragma HLS bind_storage variable=local_control type=RAM_1P impl=auto
#pragma HLS bind_storage variable=local_state type=FIFO impl=memory
#endif

        _hw_real k1[_Nx], k2[_Nx], k3[_Nx], k4[_Nx];
        _hw_real state_temp[_Nx];
        
        model(k1,local_state,local_control);
        for (unsigned i = 0; i < _Nx; ++i) {
            state_temp[i] = local_state[i] + Ts_2*k1[i];
        }

        model(k2,state_temp,local_control);
        for (unsigned i = 0; i < _Nx; ++i) {
            state_temp[i] = local_state[i] + Ts_2*k2[i];
        }

        model(k3,state_temp,local_control);
        for (unsigned i = 0; i < _Nx; ++i) {
            state_temp[i] = local_state[i] + Ts*k3[i];
        }

        model(k4,state_temp,local_control);
        for (unsigned i = 0; i < _Nx; ++i) {
            state_plus[i] = local_state[i] + Ts_6*(k1[i] + (_hw_real)2.0*k2[i] + (_hw_real)2.0*k3[i] + k4[i]);
        }
    };

    void Ji_error(
        _hw_real Ji[_Nx], 
        _hw_real Jui[_n_U],
        _hw_real *J
        ){
#pragma HLS INLINE

#pragma HLS ALLOCATION operation instances=hmul limit=1
#pragma HLS ALLOCATION operation instances=hadd limit=1
#pragma HLS ALLOCATION operation instances=hsub limit=1

        _hw_real J_local = 0.0;
        _hw_real J_local_ant = 0.0;
        Ji_loop: for (unsigned j = 0; j < _Nx; j++) {
            J_local_ant = J_local;
            J_local = J_local_ant + Q[j]*(Ji[j]);
        }
        Jui_loop: for (unsigned j = 0; j < _n_U; j++) {
            J_local_ant = J_local;
            J_local = J_local_ant + R[j]*(Jui[j]);
        }
        J[0] = J_local;
    }

    void Jf_error(
        _hw_real final_x[_Nx],
        _hw_real final_xref[_Nx],
        _hw_real *J
        ){
#pragma HLS INLINE

#pragma HLS ALLOCATION operation instances=hmul limit=1
#pragma HLS ALLOCATION operation instances=hadd limit=1
#pragma HLS ALLOCATION operation instances=hsub limit=1
        // Calculate Final State Error
        _hw_real J_local = 0.0;
        _hw_real J_local_ant = 0.0;
        _hw_real Jf = 0.0; // = 0.0;
        //memset(Jf, (const _hw_real)0.0, _Nx*sizeof(_hw_real));

        Jf_error_loop: for (unsigned j = 0; j < _Nx; ++j) {
            J_local_ant = J_local;
            _hw_real single_x_hat = final_x[j];
            _hw_real tmp_err = (single_x_hat - final_xref[j]);
            _hw_real penality = ((state_lower_limits[j] <= single_x_hat) && (single_x_hat <= state_upper_limits[j])) ? (_hw_real)1.0 : (_hw_real)1e10;
            Jf = penality*tmp_err*tmp_err;
            J_local = J_local_ant + Qf[j]*Jf;
        }
        J[0] = J_local;
    }

    void horizon_u_error(
        _hw_real uu[_n_U*_Nh],
        _hw_real Jui[_n_U]
        ){
#pragma HLS INLINE
        _hw_real Jui_buff[_n_U], current_uu[_n_U];
        memset(Jui_buff, (const _hw_real)0.0, _n_U*sizeof(_hw_real));
        for (unsigned i = 0; i < N; ++i) {
        	memcpy_loop_feedback<_hw_real, _hw_real, _n_U>(current_uu, &uu[_n_U*i]);
			one_step_u_error(Jui_buff, Jui_buff, current_uu);
            // one_step_u_error(Jui_buff, Jui, current_uu);
            // one_step_u_error(Jui_buff, Jui_buff, &uu[_n_U*i]);
            // memcpy_loop_enclosed<_hw_real, _hw_real, _n_U>(&uu[_n_U*i], (const _hw_real *)current_uu);
            // memcpy_loop_enclosed<_hw_real, _hw_real, _n_U>(Jui, (const _hw_real *)Jui_buff);
        }
        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Jui, (const _hw_real *)Jui_buff);
    }

    void horizon_step_error(
        _hw_real x_horizon[_Nx*_Nh], 
        _hw_real xref[_Nx*_Nh], 
        _hw_real Ji[_Nx],
        _hw_real final_x[_Nx],
        _hw_real final_xref[_Nx]
        ){
#pragma HLS INLINE
#pragma HLS allocation function instances=one_step_error limit=1 
        _hw_real current_x[_Nx], current_xref[_Nx]; 
// #pragma HLS STREAM variable=current_xref depth=8
// //#pragma HLS data_pack variable=current_xref
// #pragma HLS RESOURCE variable=current_xref core=FIFO_LUTRAM
        _hw_real Ji_buff[_Nx]; // = 0.0;
// #pragma HLS STREAM variable=Ji_buff depth=8
// //#pragma HLS data_pack variable=Ji_buff
// #pragma HLS RESOURCE variable=Ji_buff core=FIFO_LUTRAM
        memset(Ji_buff, (const _hw_real)0.0, _Nx*sizeof(_hw_real));
        for (unsigned i = 0; i < N; ++i) {
        	memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(current_xref, (const _hw_real *)&xref[_Nx*i]);
        	memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(current_x, (const _hw_real *)&x_horizon[_Nx*i]);
			one_step_error(Ji_buff, Ji_buff, current_xref, current_x);
#ifdef DEBUG_SYSTEM
            std::cout << "Horizon[" << i << "]"<< std::endl;
            std::cout << "\t State"; print_formatted_float_array(current_x, _Nx, 2, 6);
            std::cout << std::endl;
            std::cout << "\t Set  "; print_formatted_float_array(current_xref, _Nx, 2, 6); 
            std::cout << std::endl;
            std::cout << "\t Ji   "; print_formatted_float_array(Ji_buff, _Nx, 2, 6); 
            std::cout << std::endl;
            //std::cout << "\t Control"; print_formatted_float_array(&uu[_n_U*i], _n_U, 2, 6); std::cout << std::endl;
#endif
        }
        memcpy_loop_enclosed<_hw_real, _hw_real, _Nx>(Ji, (const _hw_real *)Ji_buff);
        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(final_x, (const _hw_real *)current_x);
        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(final_xref, (const _hw_real *)current_xref);
        	
    }

    void horizon_step_prediction(
        _hw_real x_initial[_Nx],
        _hw_real uu[_n_U*_Nh], 
        _hw_real x_horizon[_Nx*_Nh]
        ){
        _hw_real x_hat_out[_Nx];
        _hw_real x_hat_in[_Nx];
// #pragma HLS STREAM variable=x_hat depth=8
// //#pragma HLS data_pack variable=x_hat
// #pragma HLS RESOURCE variable=x_hat core=FIFO_LUTRAM
        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(x_hat_in, (const _hw_real *)x_initial);
        for (unsigned i = 0; i < N; ++i) {
			// one_step_prediction(x_buffer, current_x_hat, current_uu);
            one_step_prediction(x_hat_out, x_hat_in, &uu[_n_U*i]);
            // memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(&x_horizon[_Nx*(i)], (const _hw_real *)x_hat);
            memcpy_loop_rolled_2dest<_hw_real, _hw_real, _Nx>(&x_horizon[_Nx*(i)], x_hat_in, (const _hw_real *)x_hat_out);
        }
    }

    void uu_loop(
        _hw_real control_guess[_Nx*_Nu], 
        _hw_real uu[_n_U*_Nh]
        ){
#pragma HLS inline
        unsigned k_cg = 0;
        unsigned k_u = 0;
        const unsigned k_cg_last = _n_U*_Nu - _n_U;
        _hw_real last_control_guess[_n_U];
        for (unsigned i = 0; i < _n_U*_Nh; ++i) {
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
    }

    void uu_loop(
        _hw_real control_guess[_Nx*_Nu], 
        _hw_real uu_1[_n_U*_Nh],
        _hw_real uu_2[_n_U*_Nh]
        ){
#pragma HLS inline
        unsigned k_cg = 0;
        unsigned k_u = 0;
        const unsigned k_cg_last = _n_U*_Nu - _n_U;
        _hw_real last_control_guess[_n_U];
        for (unsigned i = 0; i < _n_U*_Nh; ++i) {
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
            uu_1[i] = control_val;
            uu_2[i] = control_val;
        }
    }
};

#endif
