#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#ifdef __SYNTHESIS__
#include "ap_int.h"
#include "hls_stream.h"
#endif

#include "config.hpp"
#include "aux_functions.hpp"
//#include "hls_inverted_pendulum.hpp"

// #define DEBUG_SYSTEM
#define def_N 3

// #ifndef __VITIS_HLS__
// #define __VITIS_HLS__
// #endif

template <
    class _hw_real,
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

    const unsigned N = _Nh;
    // Temporary variables for cost functions computations
    //const _hw_real *uu;//[_Nh];
    //_hw_real *x_hat;//[_Nx*_Nh];

    // Array to define each system state as an angle or not for relative MSE computing compensating for angles above 360º 
    //int *state_type;//[_Nx];

    // Acceleration Control
    // _hw_real ddu_max[_n_U];
    // _hw_real acc_max[_n_U];
    // _hw_real acc_min[_n_U];
    model_t model_ptr;
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
        const _hw_real _Ts
    ) : Ts(_Ts), Ts_2(_Ts*(_hw_real)0.5), Ts_6(_Ts*(_hw_real)0.1666666667),
        uss(_uss), u_max(_u_max), u_min(_u_min), du_max(_du_max),
        controlled_state(_controlled_state), 
        state_upper_limits(_state_upper_limits), state_lower_limits(_state_lower_limits),
        Q(_Q), Qf(_Qf), R(_R)
    {};

// ---------------------------------------------------
	void nmpc_cost_function_lowflow(
        volatile _hw_real *current_state,
        volatile _hw_real *control_guess,
        volatile _hw_real *xref,

        volatile _hw_real *Jui_in,
        volatile _hw_real *Ji_in,
        
        _hw_real *next_state,
        _hw_real *Jui_out,
        _hw_real *Ji_out
    ){
#pragma HLS DATAFLOW

// #pragma HLS interface mode=ap_fifo port=current_state   depth=_Nx
// #pragma HLS interface mode=ap_fifo port=control_guess   depth=_n_U
// #pragma HLS interface mode=ap_fifo port=xref            depth=_Nx

// #pragma HLS interface mode=ap_fifo port=Ji_in           depth=_Nx
// #pragma HLS interface mode=ap_fifo port=Jui_in          depth=_n_U

// #pragma HLS interface mode=ap_fifo port=next_state      depth=_Nx
// #pragma HLS interface mode=ap_fifo port=Ji_out          depth=_Nx
// #pragma HLS interface mode=ap_fifo port=Jui_out         depth=_n_U

    _hw_real uu_1[_n_U], uu_2[_n_U];       
#pragma HLS stream          variable=uu_1   type=FIFO depth=_n_U
#pragma HLS stream          variable=uu_2   type=FIFO depth=_n_U
#pragma HLS bind_storage    variable=uu_1   type=FIFO impl=LUTRAM
#pragma HLS bind_storage    variable=uu_2   type=FIFO impl=LUTRAM
    _hw_real local_next_state[_Nx], xx_1[_Nx];
#pragma HLS stream          variable=local_next_state   type=FIFO depth=_Nx
#pragma HLS stream          variable=xx_1               type=FIFO depth=_Nx
#pragma HLS bind_storage    variable=local_next_state   type=FIFO impl=LUTRAM
#pragma HLS bind_storage    variable=xx_1               type=FIFO impl=LUTRAM

    // Create control guess for full horizon
    // Constructing the vector of guessed control actions with respect to N and Nu
    split<_hw_real, _hw_real, _n_U>(control_guess, uu_1, uu_2);
    // Control Error
    one_step_u_error(Jui_out, Jui_in, uu_1);
    // Predict Horizon
    one_step_prediction(local_next_state, current_state, uu_2);
    split<_hw_real, _hw_real, _Nx>(local_next_state, xx_1, next_state);
    // State Error
    one_step_error(Ji_out, Ji_in, xx_1, xref);
}
// ---------------------------------------------------
	void nmpc_cost_function_topflow(
        volatile _hw_real current_state[_Nx],
        volatile _hw_real control_guess[_Nu*_n_U],
        volatile _hw_real xref[_Nx*_Nh],
        _hw_real *J
    ){
// #pragma HLS DATAFLOW

#pragma HLS interface mode=ap_fifo port=current_state   depth=_Nx
#pragma HLS interface mode=ap_fifo port=control_guess   depth=_n_U
#pragma HLS interface mode=ap_fifo port=xref            depth=_Nx

    _hw_real Ji_in[_Nx], Ji_out[_Nx]; // = 0.0;
#pragma HLS stream       variable=Ji_in   type=FIFO depth=_Nx
#pragma HLS stream       variable=Ji_out  type=FIFO depth=_Nx
#pragma HLS bind_storage variable=Ji_in   type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=Ji_out  type=FIFO impl=LUTRAM

    _hw_real Jui_in[_n_U], Jui_out[_n_U]; // = 0.0;
#pragma HLS stream       variable=Jui_in  type=FIFO depth=_n_U
#pragma HLS stream       variable=Jui_out type=FIFO depth=_n_U
#pragma HLS bind_storage variable=Jui_in  type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=Jui_out type=FIFO impl=LUTRAM

    _hw_real x_current[_Nx], x_next[_Nx];
#pragma HLS stream       variable=x_current type=FIFO depth=_Nx
#pragma HLS stream       variable=x_next    type=FIFO depth=_Nx
#pragma HLS bind_storage variable=x_current type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=x_next    type=FIFO impl=LUTRAM

    _hw_real uu_last_ram[_n_U];
#pragma HLS bind_storage variable=uu_last_ram type=RAM_1P impl=LUTRAM
    _hw_real uu_last_fifo[_n_U];
#pragma HLS stream       variable=uu_last_fifo type=fifo depth=_Nx
#pragma HLS bind_storage variable=uu_last_fifo type=FIFO impl=LUTRAM

    _hw_real J_tmp1, J_tmp2;
    // memset_loop<_hw_real>(Ji, (_hw_real)0.0, _Nx);
    // memset_loop<_hw_real>(Jui, (_hw_real)0.0, _n_U);

    unsigned k_xref = 0;
    unsigned k_u = 0;
    const unsigned k_cg_last = _n_U*_Nu - _n_U;
    memset_loop<_hw_real>(Jui_in, (const _hw_real)0.0, _n_U);
    memset_loop<_hw_real>(Ji_in, (const _hw_real)0.0, _Nx);
    memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(x_current, current_state);

    // volatile _hw_real *x_ref_ptr, *u_current_ptr;
    for (unsigned i = 0; i < N-1; ++i) {
#ifdef DEBUG_SYSTEM
        std::cout << "Horizon[" << i << "]"<< std::endl;
#endif
        unsigned k_u_aux = i*_n_U;
        k_u = (k_u == k_cg_last)? k_cg_last : k_u_aux;
        k_xref = i*_Nx;

        if (k_u_aux < k_cg_last){
            // u_current_ptr = &control_guess[k_u];
            memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(uu_last_fifo, &control_guess[k_u]);
        }else if (k_u_aux >= k_cg_last){
            if(k_u_aux == k_cg_last)
                memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(uu_last_ram, &control_guess[k_u]);
            memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(uu_last_fifo, uu_last_ram);
            // u_current_ptr = uu_last_fifo;
        }
        // x_ref_ptr = &xref[k_xref];

        nmpc_cost_function_lowflow(
            x_current, 
            uu_last_fifo, 
            &xref[k_xref], 
            Jui_in, 
            Ji_in, 
            x_next,
            Jui_out, 
            Ji_out
        );
        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Jui_in, Jui_out);
        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(Ji_in, Ji_out);
        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(x_current, x_next);
    }
    // x_ref_ptr = &xref[(N-1)*_Nx];
    // Weigthed Errors
    Ji_error(Ji_out, Jui_out, &J_tmp1);
    Jf_error(x_current, &xref[(N-1)*_Nx], &J_tmp2);
#ifdef DEBUG_SYSTEM
    std::cout << "Ji = " << J_tmp1 << std::endl;
    std::cout << "Jf = " << J_tmp2 << std::endl;
#endif
    final_sum(J, J_tmp1, J_tmp2);
    // return J;
}
// ---------------------------------------------------
	void nmpc_cost_function(
        volatile _hw_real *current_state,
        volatile _hw_real *control_guess,
        volatile _hw_real *xref,
        _hw_real *J
    ){

#pragma HLS DATAFLOW

#pragma HLS interface mode=ap_fifo port=current_state   depth=_Nx
#pragma HLS interface mode=ap_fifo port=control_guess   depth=_Nx
#pragma HLS interface mode=ap_fifo port=xref            depth=_Nx
#pragma HLS interface mode=ap_vld  port=J

    _hw_real uu_1[_n_U*_Nh], uu_2[_n_U*_Nh];       
#pragma HLS stream variable=uu_1 type=fifo depth=_n_U
#pragma HLS stream variable=uu_2 type=fifo depth=_n_U
#pragma HLS bind_storage variable=uu_1 type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=uu_2 type=FIFO impl=LUTRAM

    _hw_real x_horizon[_Nx*_Nh]; //[Nx*N];
#pragma HLS STREAM variable=x_horizon type=fifo depth=_Nx
#pragma HLS bind_storage variable=x_horizon type=FIFO impl=LUTRAM

    _hw_real Ji[_Nx]; // = 0.0;
#pragma HLS STREAM variable=Ji depth=_Nx
#pragma HLS bind_storage variable=Ji type=FIFO impl=LUTRAM

    _hw_real Jui[_n_U]; // = 0.0;
#pragma HLS STREAM variable=Jui depth=_n_U
#pragma HLS bind_storage variable=Jui type=FIFO impl=LUTRAM

    _hw_real final_x[_Nx], final_xref[_Nx];
#pragma HLS STREAM variable=final_x depth=_Nx
#pragma HLS STREAM variable=final_xref depth=_Nx
#pragma HLS bind_storage variable=final_x type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=final_xref type=FIFO impl=LUTRAM

        //memset(Jui, (const _hw_real)0.0, _n_U*sizeof(_hw_real));
        _hw_real J_tmp1, J_tmp2;
		// memset_loop<_hw_real>(Ji, (_hw_real)0.0, _Nx);
		// memset_loop<_hw_real>(Jui, (_hw_real)0.0, _n_U);

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
        _hw_real *state_plus,
        volatile _hw_real *state,
        volatile _hw_real *control
    ){
//#pragma HLS INLINE
//#pragma HLS DATAFLOW
//#pragma HLS pipeline off

#pragma HLS ALLOCATION operation instances=hmul limit=1
#pragma HLS ALLOCATION operation instances=hadd limit=1

#pragma HLS ALLOCATION function instances=model limit=1

        _hw_real local_control[_n_U];
        _hw_real local_state[_Nx];
        _hw_real k1[_Nx], k2[_Nx], k3[_Nx], k4[_Nx];
        _hw_real state_temp[_Nx];

#pragma HLS STREAM variable=local_state depth=_Nx type=shared
#pragma HLS STREAM variable=local_control depth=_Nx type=shared

#pragma HLS STREAM variable=k1 depth=_Nx type=shared
#pragma HLS STREAM variable=k2 depth=_Nx type=shared
#pragma HLS STREAM variable=k3 depth=_Nx type=shared
#pragma HLS STREAM variable=k4 depth=_Nx type=shared

#pragma HLS STREAM variable=state_temp depth=_Nx type=shared

        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(local_state, state);
        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(local_control, control);

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
        _hw_real _current_uu[_n_U], 
        _hw_real _control_guess[_n_U]
    ){
#pragma HLS inline
        if (_n_U*_i < _n_U*_Nu){
            memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(_current_uu, _control_guess);
        }else{
            memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(_current_uu, _current_uu);
        }
    }

    void model(
        _hw_real *state_dot,
        volatile _hw_real *state,
        volatile _hw_real *control
    ){
//#pragma HLS inline
#pragma HLS interface mode=ap_fifo port=state_dot   depth=_Nx*2
#pragma HLS interface mode=ap_fifo port=state       depth=_Nx*2
#pragma HLS interface mode=ap_fifo port=control     depth=_n_U*2
#pragma HLS pipeline off
        //model_inverted_pendulum<_hw_real, _Nx, _n_U>(state_dot, state, control);
        model_ptr.model(state_dot, state, control);
    }

    void one_step_error(
        _hw_real Ji_out[_Nx],
        volatile _hw_real Ji_in[_Nx],
        volatile _hw_real x_hat[_Nx],
        volatile _hw_real xref[_Nx]
    ){
#pragma HLS INLINE
        _hw_real penality;
        _hw_real current_x_hat;
        _hw_real current_x_ref;
//        _hw_real Ji_local[_Nx], Ji_out_local[_Nx];
//        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Ji_local, Ji_in);
        for (unsigned j = 0; j < _Nx; ++j) {
#pragma HLS PIPELINE II=11
            //_hw_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
            current_x_hat = x_hat[j];
            current_x_ref = xref[j];
            _hw_real tmp_err = (controlled_state[j] == 1) ? (current_x_hat - current_x_ref) : (_hw_real)0.0;
            penality = ((state_lower_limits[j] <= current_x_hat) && (current_x_hat <= state_upper_limits[j])) ? (_hw_real)Q[j] : (_hw_real)1e4;
//            Ji_out_local[j] = Ji_local[j] + penality*tmp_err*tmp_err;
            Ji_out[j] = Ji_in[j] + penality*tmp_err*tmp_err;
        }
//        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Ji_out, Ji_out_local);
#ifdef DEBUG_SYSTEM
        std::cout << " State"; print_formatted_float_array(x_hat, _Nx, 2, 6);
        std::cout << std::endl;
        std::cout << " Set  "; print_formatted_float_array(xref, _Nx, 2, 6); 
        std::cout << std::endl;
        std::cout << " Ji   "; print_formatted_float_array(Ji_out, _Nx, 2, 6); 
        std::cout << std::endl;
        //std::cout << "\t Control"; print_formatted_float_array(&uu[_n_U*i], _n_U, 2, 6); std::cout << std::endl;
#endif

    }

    void one_step_u_error(
        _hw_real *Ji_out,
        volatile _hw_real *Ji_in,
        volatile _hw_real *uu
        // _hw_real uref[_n_U]
    ){
#pragma HLS INLINE
        _hw_real penality;
//        _hw_real Ji_local[_n_U];
//        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Ji_local, Ji_in);
        for (unsigned j = 0; j < _n_U; ++j) {
#pragma HLS pipeline II=9
            //_hw_real tmp_err = normalize_angle(x_hat[l] - xref[l]);
            _hw_real current_uu = uu[j];
            _hw_real tmp_err = (current_uu - uss[j]);
            penality = ((u_min[j] <= current_uu) && (current_uu <= u_max[j])) ? R[j] : (_hw_real)1e4;
            Ji_out[j] = Ji_in[j] + penality*tmp_err*tmp_err;
        }
    }

    void final_sum(
        _hw_real *J, 
        _hw_real J_tmp1, 
        _hw_real J_tmp2
        ){
// #pragma HLS INLINE 
        J[0] = J_tmp1 + J_tmp2;
        
    }

    void Ji_error(
        volatile _hw_real *Ji, 
        volatile _hw_real *Jui,
        _hw_real *J
        ){
// #pragma HLS INLINE
// #pragma HLS dataflow

#pragma HLS ALLOCATION operation instances=hmul limit=1
#pragma HLS ALLOCATION operation instances=hadd limit=1

        _hw_real J_local = 0.0;
        _hw_real Ji_mul[_Nx], Jui_mul[_n_U];
//         Ji_mul_loop: for (unsigned j = 0; j < _Nx; j++) {
// #pragma HLS pipeline II=6
//             Ji_mul[j] = Q[j]*(Ji[j]);
//         }
        Ji_sum_loop: for (unsigned j = 0; j < _Nx; j++) {
#pragma HLS pipeline off
            // J_local += Ji_mul[j];
            J_local += Ji[j];
        }
//         Jui_mul_loop: for (unsigned j = 0; j < _n_U; j++) {
// #pragma HLS pipeline II=6
//             Jui_mul[j] = R[j]*(Jui[j]);
//         }
        Jui_sum_loop: for (unsigned j = 0; j < _n_U; j++) {
#pragma HLS pipeline off
            // J_local += Jui_mul[j];
            J_local += Jui[j];
        }
        J[0] = J_local;
    }

    void Jf_error(
        volatile _hw_real *final_x,
        volatile _hw_real *final_xref,
        _hw_real *J
        ){
// #pragma HLS INLINE

#pragma HLS ALLOCATION operation instances=hmul limit=1
#pragma HLS ALLOCATION operation instances=hadd limit=1
        // Calculate Final State Error
        _hw_real J_local = 0.0;
        _hw_real J_mul[_Nx];
        _hw_real Jf; // = 0.0;
        //memset(Jf, (const _hw_real)0.0, _Nx*sizeof(_hw_real));

        Jf_mul_loop: for (unsigned j = 0; j < _Nx; ++j) {
#pragma HLS pipeline II=6
            _hw_real single_x_hat = final_x[j];
            _hw_real tmp_err = (single_x_hat - final_xref[j]);
            _hw_real penality = ((state_lower_limits[j] <= single_x_hat) && (single_x_hat <= state_upper_limits[j])) ? (_hw_real)1.0 : (_hw_real)1e4;
            Jf = penality*tmp_err*tmp_err;
            J_mul[j] = Qf[j]*Jf;
        }
        Jf_sum_loop: for (unsigned j = 0; j < _Nx; ++j) {
#pragma HLS pipeline off
            J_local += J_mul[j];
        }
        J[0] = J_local;
    }

    void horizon_u_error(
        volatile _hw_real *uu,
        _hw_real *Jui
        ){
// #pragma HLS INLINE
        _hw_real Jui_buff[_n_U], Jui_buff_ant[_n_U], current_uu[_n_U];
#ifdef __VITIS_HLS__
#pragma HLS bind_storage variable=Jui_buff type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=Jui_buff_ant type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=current_uu type=FIFO impl=LUTRAM
#endif
        memset_loop<_hw_real>(Jui_buff_ant, (const _hw_real)0.0, _n_U);
        for (unsigned i = 0; i < N; ++i) {
        	memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(current_uu, &uu[_n_U*i]);
			one_step_u_error(Jui_buff, Jui_buff_ant, current_uu);
            // one_step_u_error(Jui_buff, Jui, current_uu);
            // one_step_u_error(Jui_buff, Jui_buff, &uu[_n_U*i]);
            // memcpy_loop_enclosed<_hw_real, _hw_real, _n_U>(&uu[_n_U*i], current_uu);
            memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Jui_buff_ant, Jui_buff);
        }
        memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(Jui, Jui_buff_ant);
    }

    void horizon_step_error(
        volatile _hw_real *x_horizon, 
        volatile _hw_real *xref, 
        _hw_real *Ji,
        _hw_real *final_x,
        _hw_real *final_xref
        ){
// #pragma HLS INLINE
#pragma HLS allocation function instances=one_step_error limit=1 
        _hw_real current_x[_Nx], current_xref[_Nx]; 
// #pragma HLS STREAM variable=current_xref depth=8
// //#pragma HLS data_pack variable=current_xref
// #pragma HLS RESOURCE variable=current_xref core=FIFO_LUTRAM
        _hw_real Ji_buff[_Nx], Ji_buff_ant[_Nx]; // = 0.0;
#ifdef __VITIS_HLS__
#pragma HLS bind_storage variable=Ji_buff type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=Ji_buff_ant type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=current_x type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=current_xref type=FIFO impl=LUTRAM
#endif

// #pragma HLS STREAM variable=Ji_buff depth=8
// //#pragma HLS data_pack variable=Ji_buff
// #pragma HLS RESOURCE variable=Ji_buff core=FIFO_LUTRAM
        memset_loop<_hw_real>(Ji_buff_ant, (const _hw_real)0.0, _Nx);
        for (unsigned i = 0; i < N; ++i) {
#pragma HLS PIPELINE off
//#pragma HLS loop_merge
#pragma HLS allocation operation instances=hmul limit=1
#pragma HLS allocation operation instances=hadd limit=1
            if (i==N-1){
        	split<_hw_real, _hw_real, _Nx>(&xref[_Nx*i], current_xref, final_x);
        	split<_hw_real, _hw_real, _Nx>(&x_horizon[_Nx*i], current_x, final_xref);
            }else{
        	memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(current_xref, &xref[_Nx*i]);
        	memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(current_x, &x_horizon[_Nx*i]);
            }
#ifdef DEBUG_SYSTEM
            std::cout << "Horizon[" << i << "]"<< std::endl;
#endif
			one_step_error(Ji_buff, Ji_buff_ant, current_x, current_xref);
            memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(Ji_buff_ant, Ji_buff);
        }
        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(Ji, Ji_buff_ant);
//        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(final_x, current_x);
//        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(final_xref, current_xref);
        	
    }

    void horizon_step_prediction(
        volatile _hw_real *x_initial,
        volatile _hw_real *uu, 
        _hw_real *x_horizon
        ){
        _hw_real x_hat_out[_Nx];
        _hw_real x_hat_in[_Nx];
// #pragma HLS STREAM variable=x_hat depth=8
// //#pragma HLS data_pack variable=x_hat
#pragma HLS bind_storage variable=x_hat_out type=FIFO impl=LUTRAM
#pragma HLS bind_storage variable=x_hat_in type=FIFO impl=LUTRAM
        memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(x_hat_in, x_initial);
        for (unsigned i = 0; i < N; ++i) {
#pragma HLS pipeline off
			// one_step_prediction(x_buffer, current_x_hat, current_uu);
            one_step_prediction(x_hat_out, x_hat_in, &uu[_n_U*i]);
            // memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(&x_horizon[_Nx*(i)], x_hat);
            split<_hw_real, _hw_real, _Nx>(x_hat_out, &x_horizon[_Nx*(i)], x_hat_in);
        }
    }

    void uu_loop(
        _hw_real volatile *control_guess, 
        _hw_real *uu
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
        volatile _hw_real *control_guess, 
        _hw_real *uu_1,
        _hw_real *uu_2
        ){
// #pragma HLS inline
        _real control_val;
        _real control_last[_n_U];
        for (int i = 0; i < _Nu; ++i) {
            for (int j = 0; j < _n_U; ++j) {
                control_val = control_guess[i*_n_U+j];
                if (i == _Nu-1)
                    control_last[j] = control_val;
                uu_1[i*_n_U+j] =  control_val;
                uu_2[i*_n_U+j] =  control_val;
            }
        }
        for (int i = _Nu; i < _Nh; ++i) {
            for (int j = 0; j < _n_U; ++j) {
                control_val = control_last[j];
                uu_1[i*_n_U+j] =  control_val;
                uu_2[i*_n_U+j] =  control_val;
            }
        }
    }

    void update_state(
        _hw_real *state_plus, 
        _hw_real *state, 
        _hw_real *k, 
        _hw_real Ts_local
        ){        
#pragma HLS inline
        for (unsigned i = 0; i < _Nx; ++i) {
            state_plus[i] = state[i] + Ts_local*k[i];
        }
    }

    void update_state_final(
        _hw_real *state_plus, 
        _hw_real *state, 
        _hw_real *k1, 
        _hw_real *k2, 
        _hw_real *k3, 
        _hw_real *k4
        ){        
#pragma HLS inline
        for (unsigned i = 0; i < _Nx; ++i) {
            state_plus[i] = state[i] + Ts_6*(k1[i] + (_hw_real)2.0*k2[i] + (_hw_real)2.0*k3[i] + k4[i]);
        }
    }

};

    
#endif
