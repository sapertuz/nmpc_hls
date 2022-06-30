#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#ifdef __SYNTHESIS__
#include "ap_int.h"
#include "hls_stream.h"
#endif

// #include "config.hpp"
#include "aux_functions.hpp"
//#include "hls_inverted_pendulum.hpp"

// #define DEBUG_SYSTEM
#define def_N 3

// #ifndef __VITIS_HLS__
// #define __VITIS_HLS__
// #endif

template <
    class _system_hw_real,
    class _system_model_t,
    unsigned _system_Nh,
    unsigned _system_Nx,
    unsigned _system_n_U,
    unsigned _system_Nu
>class System {
protected:
	// System Properties
	// _system_hw_real *current_state;//[_system_Nx];
	// _system_hw_real last_state;//[_system_Nx];
	
    // NMPC Properties
    // const unsigned N = _system_Nh;       // Prediction Horizon
    // const unsigned Nu = _system_Nu;     // Control Horizon
    // const unsigned n_U = _system_n_U;   // Number of Inputs
    // const unsigned Nx = _system_Nx;     // Number of States

    // Simulation Properties
    const _system_hw_real Ts; // Sampling rate (sec)
    const _system_hw_real Ts_2;
    const _system_hw_real Ts_6;
    
    // Constrains
    const _system_hw_real *u_max;//[_system_n_U];
    const _system_hw_real *u_min;//[_system_n_U];
    const _system_hw_real *du_max;//[_system_n_U];

    const _system_hw_real *uss;//[_system_n_U]; // Control Signal at Steady State

    const unsigned short *controlled_state;//[_system_Nx];
    const _system_hw_real *state_upper_limits;//[_system_Nx];
    const _system_hw_real *state_lower_limits;//[_system_Nx];
    // Weight Matrices
    const _system_hw_real *Q;//[_system_Nx]; 
    const _system_hw_real *Qf;//[_system_Nx];
    const _system_hw_real *R;//[_system_n_U];

    const unsigned N = _system_Nh;
    // Temporary variables for cost functions computations
    //const _system_hw_real *uu;//[_system_Nh];
    //_system_hw_real *x_hat;//[_system_Nx*_system_Nh];

    // Array to define each system state as an angle or not for relative MSE computing compensating for angles above 360º 
    //int *state_type;//[_system_Nx];

    // Acceleration Control
    // _system_hw_real ddu_max[_system_n_U];
    // _system_hw_real acc_max[_system_n_U];
    // _system_hw_real acc_min[_system_n_U];
    _system_model_t *model_ptr;
public:
    constexpr System(
        const _system_hw_real __u_max[_system_n_U],
        const _system_hw_real __u_min[_system_n_U],
        const _system_hw_real __du_max[_system_n_U],
        const unsigned short __controlled_state[_system_Nx],
        const _system_hw_real __state_upper_limits[_system_Nx],
        const _system_hw_real __state_lower_limits[_system_Nx],
        const _system_hw_real __Q[_system_Nx],
        const _system_hw_real __Qf[_system_Nx],
        const _system_hw_real __R[_system_n_U],
        const _system_hw_real __uss[_system_n_U],
        const _system_hw_real __Ts
        ,
        _system_model_t* __model_ptr
    ) : Ts(__Ts), Ts_2(__Ts*(_system_hw_real)0.5), Ts_6(__Ts*(_system_hw_real)0.1666666667),
        uss(__uss), u_max(__u_max), u_min(__u_min), du_max(__du_max),
        controlled_state(__controlled_state), 
        state_upper_limits(__state_upper_limits), state_lower_limits(__state_lower_limits),
        Q(__Q), Qf(__Qf), R(__R),
        model_ptr(__model_ptr)
    {};


// ---------------------------------------------------
	void nmpc_cost_function_lowflow(
        volatile _system_hw_real *current_state,
        volatile _system_hw_real *control_guess,
        volatile _system_hw_real *xref,

        volatile _system_hw_real *Jui_in,
        volatile _system_hw_real *Ji_in,
        
        unsigned int i_step,

        volatile _system_hw_real *next_state,
        volatile _system_hw_real *Jui_out,
        volatile _system_hw_real *Ji_out
    ){
//#pragma HLS INLINE
#pragma HLS DATAFLOW

    _system_hw_real uu_1[_system_n_U], uu_2[_system_n_U];       
#pragma HLS stream          variable=uu_1   type=FIFO depth=_system_n_U
#pragma HLS stream          variable=uu_2   type=FIFO depth=_system_n_U
#pragma HLS bind_storage    variable=uu_1   type=FIFO impl=LUTRAM
#pragma HLS bind_storage    variable=uu_2   type=FIFO impl=LUTRAM
    _system_hw_real local_next_state[_system_Nx], xx_1[_system_Nx];
#pragma HLS stream          variable=local_next_state   type=FIFO depth=_system_Nx
#pragma HLS stream          variable=xx_1               type=FIFO depth=_system_Nx
#pragma HLS bind_storage    variable=local_next_state   type=FIFO impl=LUTRAM
#pragma HLS bind_storage    variable=xx_1               type=FIFO impl=LUTRAM

    // Create control guess for full horizon
    // Constructing the vector of guessed control actions with respect to N and Nu
    split<_system_hw_real, _system_hw_real, _system_n_U>(control_guess, uu_1, uu_2);
    // Control Error
    one_step_u_error((_system_hw_real *)Jui_out, (_system_hw_real *)Jui_in, (_system_hw_real *)uu_1);
    // Predict Horizon
    one_step_prediction(local_next_state, current_state, uu_2);
    split<_system_hw_real, _system_hw_real, _system_Nx>(local_next_state, xx_1, (_system_hw_real *)next_state);
    // State Error
    one_step_error((_system_hw_real *)Ji_out, i_step, (_system_hw_real *)Ji_in, xx_1, (_system_hw_real *)xref);
}
// ---------------------------------------------------
	void nmpc_cost_function_topflow(
        volatile _system_hw_real *current_state,
        volatile _system_hw_real *control_guess,
        volatile _system_hw_real *xref,
        _system_hw_real *J
    ){
#pragma HLS interface mode=ap_fifo port=current_state   depth=_system_Nx
#pragma HLS interface mode=ap_fifo port=control_guess   depth=_system_n_U
#pragma HLS interface mode=ap_fifo port=xref            depth=_system_Nx

    _system_hw_real Ji_in[_system_Nx], Ji_out[_system_Nx]; // = 0.0;
#pragma HLS stream       variable=Ji_in   type=FIFO depth=_system_Nx
#pragma HLS stream       variable=Ji_out  type=FIFO depth=_system_Nx
#pragma HLS bind_storage variable=Ji_in   type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=Ji_out  type=FIFO impl=LUTRAM

    _system_hw_real Jui_in[_system_n_U], Jui_out[_system_n_U]; // = 0.0;
#pragma HLS stream       variable=Jui_in  type=FIFO depth=_system_n_U
#pragma HLS stream       variable=Jui_out type=FIFO depth=_system_n_U
#pragma HLS bind_storage variable=Jui_in  type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=Jui_out type=FIFO impl=LUTRAM

    _system_hw_real x_current[_system_Nx], x_next[_system_Nx];
#pragma HLS stream       variable=x_current type=FIFO depth=_system_Nx
#pragma HLS stream       variable=x_next    type=FIFO depth=_system_Nx
#pragma HLS bind_storage variable=x_current type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=x_next    type=FIFO impl=LUTRAM

    _system_hw_real uu_last_ram[_system_n_U];
#pragma HLS bind_storage variable=uu_last_ram type=RAM_1P impl=LUTRAM
    _system_hw_real uu_fifo[_system_n_U];
#pragma HLS stream       variable=uu_fifo type=fifo depth=_system_n_U
#pragma HLS bind_storage variable=uu_fifo type=FIFO impl=LUTRAM

    // unsigned k_xref = 0;
    // unsigned k_u = 0;
    // const unsigned k_cg_last = _system_n_U*_system_Nu - _system_n_U;
    memset_loop<_system_hw_real>((_system_hw_real *)Jui_in, (const _system_hw_real)0.0, _system_n_U);
    memset_loop<_system_hw_real>((_system_hw_real *)Ji_in, (const _system_hw_real)0.0, _system_Nx);
    memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(x_current, (_system_hw_real *)current_state);
    memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>(uu_fifo, (_system_hw_real *)&control_guess[0]);
        
    // volatile _system_hw_real *x_ref_ptr, *u_current_ptr;
    step_loop: for (unsigned i = 0; i < N; ++i) {
#ifdef DEBUG_SYSTEM
        std::cout << "Horizon[" << i << "]"<< std::endl;
#endif
#pragma HLS PIPELINE
        // unsigned k_u_aux = i*_system_n_U;
        // k_u = (k_u == k_cg_last)? k_cg_last : k_u_aux;
        
        nmpc_cost_function_lowflow(
            x_current, 
            uu_fifo, 
            &xref[i*_system_Nx], 
            Jui_in, 
            Ji_in, 
            i,
            x_next,
            Jui_out, 
            Ji_out
        );

        // x_ref_ptr = &xref[k_xref];
        get_next_uu(i, control_guess, uu_fifo, uu_last_ram);
        memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>(Jui_in, Jui_out);
        memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(Ji_in, Ji_out);
        memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(x_current, x_next);
    }
    // J_error(Jui_in, Ji_in, &xref[(N-1)*_system_Nx], x_current, J);
    J_error(Jui_in, Ji_in, J);
}

    void get_next_uu(
        unsigned local_i,
        volatile _system_hw_real *control_guess,
        volatile _system_hw_real *uu_fifo,
        volatile _system_hw_real *uu_last_ram
    ){
#pragma HLS inline
        const unsigned k_cg_last = _system_n_U*_system_Nu - _system_n_U;
        unsigned k_u = (local_i+1)*_system_n_U;
        if (k_u < k_cg_last){
            // u_current_ptr = &control_guess[k_u];
            memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>((_system_hw_real *)uu_fifo, (_system_hw_real *)&control_guess[k_u]);
        }else if (k_u >= k_cg_last){
            if(k_u == k_cg_last)
                memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>((_system_hw_real *)uu_last_ram, (_system_hw_real *)&control_guess[k_u]);
            memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>((_system_hw_real *)uu_fifo, (_system_hw_real *)uu_last_ram);
            // u_current_ptr = uu_last_fifo;
        }

    }

// ---------------------------------------------------
	void nmpc_cost_function(
        volatile _system_hw_real *current_state,
        volatile _system_hw_real *control_guess,
        volatile _system_hw_real *xref,
        _system_hw_real *J
    ){

#pragma HLS DATAFLOW

#pragma HLS interface mode=ap_fifo port=current_state   depth=_system_Nx
#pragma HLS interface mode=ap_fifo port=control_guess   depth=_system_Nx
#pragma HLS interface mode=ap_fifo port=xref            depth=_system_Nx
#pragma HLS interface mode=ap_vld  port=J

    _system_hw_real uu_1[_system_n_U*_system_Nh], uu_2[_system_n_U*_system_Nh];       
#pragma HLS stream variable=uu_1 type=fifo depth=_system_n_U
#pragma HLS stream variable=uu_2 type=fifo depth=_system_n_U
#pragma HLS bind_storage variable=uu_1 type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=uu_2 type=FIFO impl=LUTRAM

    _system_hw_real x_horizon[_system_Nx*_system_Nh]; //[Nx*N];
#pragma HLS STREAM variable=x_horizon type=fifo depth=_system_Nx
#pragma HLS bind_storage variable=x_horizon type=FIFO impl=LUTRAM

    _system_hw_real Ji[_system_Nx]; // = 0.0;
#pragma HLS STREAM variable=Ji depth=_system_Nx
#pragma HLS bind_storage variable=Ji type=FIFO impl=LUTRAM

    _system_hw_real Jui[_system_n_U]; // = 0.0;
#pragma HLS STREAM variable=Jui depth=_system_n_U
#pragma HLS bind_storage variable=Jui type=FIFO impl=LUTRAM

    _system_hw_real final_x[_system_Nx], final_xref[_system_Nx];
#pragma HLS STREAM variable=final_x depth=_system_Nx
#pragma HLS STREAM variable=final_xref depth=_system_Nx
#pragma HLS bind_storage variable=final_x type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=final_xref type=FIFO impl=LUTRAM

        //memset(Jui, (const _system_hw_real)0.0, _system_n_U*sizeof(_system_hw_real));
        _system_hw_real J_tmp1, J_tmp2;
		// memset_loop<_system_hw_real>(Ji, (_system_hw_real)0.0, _system_Nx);
		// memset_loop<_system_hw_real>(Jui, (_system_hw_real)0.0, _system_n_U);

        // Create control guess for full horizon
        // Constructing the vector of guessed control actions with respect to N and Nu
        uu_loop(control_guess, uu_1, uu_2);
        // Control Error
        horizon_u_error(uu_1, Jui);
        // Predict Horizon
        horizon_step_prediction(current_state, uu_2, x_horizon);
        // State Error
        horizon_step_error(x_horizon, xref, Ji, final_x, final_xref);

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

// ------------------------------------------------------
	void one_step_prediction(
        _system_hw_real *state_plus,
        volatile _system_hw_real *state,
        volatile _system_hw_real *control
    ){
//#pragma HLS INLINE
//#pragma HLS DATAFLOW
//#pragma HLS pipeline off

#pragma HLS ALLOCATION operation instances=hmul limit=1
#pragma HLS ALLOCATION operation instances=hadd limit=1

#pragma HLS ALLOCATION function instances=model limit=1

        _system_hw_real local_control[_system_n_U];
        _system_hw_real local_state[_system_Nx];
        _system_hw_real k1[_system_Nx], k2[_system_Nx], k3[_system_Nx], k4[_system_Nx];
        _system_hw_real state_temp[_system_Nx];

#pragma HLS STREAM variable=local_state depth=_system_Nx*2 type=shared
#pragma HLS STREAM variable=local_control depth=_system_n_U*2 type=shared

#pragma HLS STREAM variable=k1 depth=_system_Nx*2 type=shared
#pragma HLS STREAM variable=k2 depth=_system_Nx type=shared
#pragma HLS STREAM variable=k3 depth=_system_Nx type=shared
#pragma HLS STREAM variable=k4 depth=_system_Nx type=shared

#pragma HLS STREAM variable=state_temp depth=_system_Nx type=shared

        memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(local_state, (_system_hw_real *)state);
        memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>(local_control, (_system_hw_real *)control);

stage_1:
//#pragma HLS DATAFLOW
        model(k1, local_state, local_control);
        update_state(state_temp, local_state, k1, Ts_2);
//#pragma HLS DATAFLOW off

stage_2:
//#pragma HLS DATAFLOW
        model(k2,state_temp,local_control);
        update_state(state_temp, local_state, k2, Ts_2);
//#pragma HLS DATAFLOW off

stage_3:
//#pragma HLS DATAFLOW
        model(k3,state_temp,local_control);
        update_state(state_temp, local_state, k3, Ts);
//#pragma HLS DATAFLOW off

stage_4:
//#pragma HLS DATAFLOW
        model(k4,state_temp,local_control);
        update_state_final(state_plus, local_state, k1, k2, k3, k4);
//#pragma HLS DATAFLOW off
    };

protected:
    void memcpy_loop_checklastu(
        unsigned _i,
        _system_hw_real *_current_uu, 
        _system_hw_real *_control_guess
    ){
#pragma HLS inline
        if (_system_n_U*_i < _system_n_U*_system_Nu){
            memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>(_current_uu, _control_guess);
        }else{
            memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>(_current_uu, _current_uu);
        }
    }

    void model(
        volatile _system_hw_real *state_dot,
        volatile _system_hw_real *state,
        volatile _system_hw_real *control
    ){
//#pragma HLS inline
// #pragma HLS interface mode=ap_fifo port=state_dot   depth=_system_Nx*2
// #pragma HLS interface mode=ap_fifo port=state       depth=_system_Nx*2
// #pragma HLS interface mode=ap_fifo port=control     depth=_system_n_U*2
// #pragma HLS pipeline off
        //model_inverted_pendulum<_system_hw_real, _system_Nx, _system_n_U>(state_dot, state, control);
        model_ptr->model(state_dot, state, control);
    }

    void one_step_error(
        _system_hw_real *Ji_out,
        unsigned int i_step,
        _system_hw_real *Ji_in,
        _system_hw_real *x_hat,
        _system_hw_real *xref
    ){
// #pragma HLS INLINE
        _system_hw_real penality;
        _system_hw_real current_x_hat;
        _system_hw_real current_x_ref;
        for (unsigned j = 0; j < _system_Nx; ++j) {
#pragma HLS PIPELINE II=11
            current_x_hat = x_hat[j];
            current_x_ref = xref[j];
            _system_hw_real tmp_err = (controlled_state[j] == 1) ? (current_x_hat - current_x_ref) : (_system_hw_real)0.0;
            _system_hw_real Q_local = (i_step == N-1)? Qf[j] : Q[j];
            penality = ((state_lower_limits[j] <= current_x_hat) && (current_x_hat <= state_upper_limits[j])) ? (_system_hw_real)Q_local : (_system_hw_real)1e4;
            Ji_out[j] = Ji_in[j] + penality*tmp_err*tmp_err;
        }
#ifdef DEBUG_SYSTEM
        std::cout << " State"; print_formatted_float_array(x_hat, _system_Nx, 2, 6);
        std::cout << std::endl;
        std::cout << " Set  "; print_formatted_float_array(xref, _system_Nx, 2, 6); 
        std::cout << std::endl;
        std::cout << " Ji   "; print_formatted_float_array(Ji_out, _system_Nx, 2, 6); 
        std::cout << std::endl;
        //std::cout << "\t Control"; print_formatted_float_array(&uu[_system_n_U*i], _system_n_U, 2, 6); std::cout << std::endl;
#endif

    }

// ---------------------------------------------------

    void one_step_u_error(
        _system_hw_real *Ji_out,
        _system_hw_real *Ji_in,
        _system_hw_real *uu
        // _system_hw_real uref[_system_n_U]
    ){
// #pragma HLS INLINE
        _system_hw_real penality;
        for (unsigned j = 0; j < _system_n_U; ++j) {
#pragma HLS pipeline II=9
            _system_hw_real current_uu = uu[j];
            _system_hw_real tmp_err = (current_uu - uss[j]);
            penality = ((u_min[j] <= current_uu) && (current_uu <= u_max[j])) ? R[j] : (_system_hw_real)1e4;
            Ji_out[j] = Ji_in[j] + penality*tmp_err*tmp_err;
        }
    }

    void final_sum(
        _system_hw_real *J, 
        _system_hw_real J_tmp1, 
        _system_hw_real J_tmp2
        ){
#pragma HLS INLINE 
        J[0] = J_tmp1 + J_tmp2;
        
    }

void J_error(
        volatile _system_hw_real *local_Jui,
        volatile _system_hw_real *local_Ji,
        
        _system_hw_real *J
    ){
    // Weigthed Errors    
#pragma HLS DATAFLOW
// #pragma HLS loop_merge 
// #pragma HLS ALLOCATION operation instances=hmul limit=1
#pragma HLS ALLOCATION operation instances=hadd limit=2

    _system_hw_real Ji_acum = 0.0;
    _system_hw_real Jui_acum = 0.0;

    Ji_sum_loop: for (unsigned j = 0; j < _system_Nx; j++) {
#pragma HLS pipeline off
        // J_local += Ji_mul[j];
        Ji_acum += local_Ji[j];
    }

    Jui_sum_loop: for (unsigned j = 0; j < _system_n_U; j++) {
#pragma HLS pipeline off
        // J_local += Jui_mul[j];
        Jui_acum += local_Jui[j];
    }
    J[0] = Ji_acum + Jui_acum;
#ifdef DEBUG_SYSTEM
    std::cout << "Ji = " << Ji_acum << std::endl;
    std::cout << "Ju = " << Jui_acum << std::endl;

#endif
}

    void Ji_error(
        volatile _system_hw_real *Ji, 
        volatile _system_hw_real *Jui,
        _system_hw_real *J
        ){
#pragma HLS INLINE
// #pragma HLS dataflow

// #pragma HLS ALLOCATION operation instances=hmul limit=1
// #pragma HLS ALLOCATION operation instances=hadd limit=1

        _system_hw_real Ji_acum = 0.0;
        _system_hw_real Jui_acum = 0.0;

        Ji_sum_loop: for (unsigned j = 0; j < _system_Nx; j++) {
#pragma HLS pipeline off
            // J_local += Ji_mul[j];
            Ji_acum += Ji[j];;
        }

        Jui_sum_loop: for (unsigned j = 0; j < _system_n_U; j++) {
#pragma HLS pipeline off
            // J_local += Jui_mul[j];
            Jui_acum += Jui[j];
        }
        J[0] = Ji_acum + Jui_acum;
    }

    void Jf_error(
        volatile _system_hw_real *final_x,
        volatile _system_hw_real *final_xref,
        _system_hw_real *J
        ){
#pragma HLS INLINE

        // Calculate Final State Error
        _system_hw_real J_local = 0.0;
        // _system_hw_real J_mul[_system_Nx];
        _system_hw_real Jf; // = 0.0;

        Jf_mul_loop: for (unsigned j = 0; j < _system_Nx; ++j) {
#pragma HLS pipeline II=6
            _system_hw_real single_x_hat = final_x[j];
            _system_hw_real tmp_err = (single_x_hat + (-final_xref[j]));
            _system_hw_real penality = ((state_lower_limits[j] <= single_x_hat) && (single_x_hat <= state_upper_limits[j])) ? (_system_hw_real)1.0 : (_system_hw_real)1e4;
            Jf = penality*tmp_err*tmp_err;
            _system_hw_real J_mul = Qf[j]*Jf;
            J_local += J_mul;
        }
        J[0] = J_local;
    }

    void horizon_u_error(
        volatile _system_hw_real *uu,
        _system_hw_real *Jui
        ){
// #pragma HLS INLINE
        _system_hw_real Jui_buff[_system_n_U], Jui_buff_ant[_system_n_U], current_uu[_system_n_U];
#ifdef __VITIS_HLS__
#pragma HLS bind_storage variable=Jui_buff type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=Jui_buff_ant type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=current_uu type=FIFO impl=LUTRAM
#endif
        memset_loop<_system_hw_real>(Jui_buff_ant, (const _system_hw_real)0.0, _system_n_U);
        for (unsigned i = 0; i < N; ++i) {
        	memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>(current_uu, (_system_hw_real *)&uu[_system_n_U*i]);
			one_step_u_error(Jui_buff, Jui_buff_ant, current_uu);
            memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>(Jui_buff_ant, Jui_buff);
        }
        memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_n_U>(Jui, Jui_buff_ant);
    }

    void horizon_step_error(
        volatile _system_hw_real *x_horizon, 
        volatile _system_hw_real *xref, 
        _system_hw_real *Ji,
        _system_hw_real *final_x,
        _system_hw_real *final_xref
        ){
#pragma HLS allocation function instances=one_step_error limit=1 
        _system_hw_real current_x[_system_Nx], current_xref[_system_Nx]; 
        _system_hw_real Ji_buff[_system_Nx], Ji_buff_ant[_system_Nx]; // = 0.0;
#ifdef __VITIS_HLS__
#pragma HLS bind_storage variable=Ji_buff type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=Ji_buff_ant type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=current_x type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=current_xref type=FIFO impl=LUTRAM
#endif

        memset_loop<_system_hw_real>(Ji_buff_ant, (const _system_hw_real)0.0, _system_Nx);
        for (unsigned i = 0; i < N-1; ++i) {
#pragma HLS PIPELINE off
//#pragma HLS loop_merge
#pragma HLS allocation operation instances=hmul limit=1
#pragma HLS allocation operation instances=hadd limit=1
            memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(current_xref, (_system_hw_real *)&xref[_system_Nx*i]);
        	memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(current_x, (_system_hw_real *)&x_horizon[_system_Nx*i]);
#ifdef DEBUG_SYSTEM
            std::cout << "Horizon[" << i << "]"<< std::endl;
#endif
			one_step_error(Ji_buff, i, Ji_buff_ant, current_x, current_xref);
            memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(Ji_buff_ant, Ji_buff);
        }

        memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(final_x, (_system_hw_real *)&xref[_system_Nx*(N-1)]);
        memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(final_xref, (_system_hw_real *)&x_horizon[_system_Nx*(N-1)]);
        memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(Ji, Ji_buff_ant);
        	
    }

    void horizon_step_prediction(
        volatile _system_hw_real *x_initial,
        volatile _system_hw_real *uu, 
        _system_hw_real *x_horizon
        ){
        _system_hw_real x_hat_out[_system_Nx];
        _system_hw_real x_hat_in[_system_Nx];

#pragma HLS bind_storage variable=x_hat_out type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=x_hat_in type=FIFO impl=LUTRAM
        memcpy_loop_rolled<_system_hw_real, _system_hw_real, _system_Nx>(x_hat_in, (_system_hw_real *)x_initial);
        for (unsigned i = 0; i < N; ++i) {
#pragma HLS pipeline off
            one_step_prediction(x_hat_out, x_hat_in, &uu[_system_n_U*i]);
            split<_system_hw_real, _system_hw_real, _system_Nx>(x_hat_out, &x_horizon[_system_Nx*(i)], x_hat_in);
        }
    }

    void uu_loop(
        _system_hw_real volatile *control_guess, 
        _system_hw_real *uu
        ){
#pragma HLS inline
        unsigned k_cg = 0;
        unsigned k_u = 0;
        const unsigned k_cg_last = _system_n_U*_system_Nu - _system_n_U;
        _system_hw_real last_control_guess[_system_n_U];
        for (unsigned i = 0; i < _system_n_U*_system_Nh; ++i) {
            _system_hw_real control_val;
            if(i < _system_n_U*_system_Nu){
                control_val = control_guess[i];
                if (i >= k_cg_last){
                    last_control_guess[k_cg] = control_val;
                    k_cg++;
                }
            }else{
                control_val = last_control_guess[k_u];
                k_u = (k_u<_system_n_U) ? k_u+1 : 0;
            }
            uu[i] = control_val;
        }
    }

    void uu_loop(
        volatile _system_hw_real *control_guess, 
        volatile _system_hw_real *uu_1,
        volatile _system_hw_real *uu_2
        ){
// #pragma HLS inline
        _system_hw_real control_val;
        _system_hw_real control_last[_system_n_U];
        for (int i = 0; i < _system_Nu; ++i) {
            for (int j = 0; j < _system_n_U; ++j) {
                control_val = control_guess[i*_system_n_U+j];
                if (i == _system_Nu-1)
                    control_last[j] = control_val;
                uu_1[i*_system_n_U+j] =  control_val;
                uu_2[i*_system_n_U+j] =  control_val;
            }
        }
        for (int i = _system_Nu; i < _system_Nh; ++i) {
            for (int j = 0; j < _system_n_U; ++j) {
                control_val = control_last[j];
                uu_1[i*_system_n_U+j] =  control_val;
                uu_2[i*_system_n_U+j] =  control_val;
            }
        }
    }

    void update_state(
        _system_hw_real *state_plus, 
        _system_hw_real *state, 
        _system_hw_real *k, 
        _system_hw_real Ts_local
        ){        
#pragma HLS inline
        for (unsigned i = 0; i < _system_Nx; ++i) {
            state_plus[i] = state[i] + Ts_local*k[i];
        }
    }

    void update_state_final(
        _system_hw_real *state_plus, 
        _system_hw_real *state, 
        _system_hw_real *k1, 
        _system_hw_real *k2, 
        _system_hw_real *k3, 
        _system_hw_real *k4
        ){        
#pragma HLS inline
        for (unsigned i = 0; i < _system_Nx; ++i) {
            state_plus[i] = state[i] + Ts_6*(k1[i] + (_system_hw_real)2.0*k2[i] + (_system_hw_real)2.0*k3[i] + k4[i]);
        }
    }

};

    
#endif
