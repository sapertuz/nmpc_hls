############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2019 Xilinx, Inc. All Rights Reserved.
############################################################
set ip_path "/home/chello/Documents/Vivado_WS/vitis_ip_repo"

open_project nmpc_hls_model
set_top model_wrapper
add_files ../aux_functions.hpp
add_files ../fast_sin_cos.hpp
add_files ../hls_inverted_pendulum.hpp
add_files ../hls_sniffbot.hpp
add_files ../main_hls_model.cpp -cflags "-DSNIFFBOT_CONFIG -DUSE_FAST_SIN_COS"
add_files -tb ../main_hls_model.cpp -cflags "-DSNIFFBOT_CONFIG -DUSE_FAST_SIN_COS -Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"

open_solution "solution_model" -flow_target vivado

set_part {xc7z020-clg400-1}
create_clock -period 10 -name default
config_compile -pipeline_loops 6
config_export -display_name model_sniffbot -format ip_catalog -output $ip_path/model_sniffbot_fast_sincos.zip -rtl verilog -vendor tu-dresden -version 1.0

# config_export -vivado_phys_opt all -vivado_optimization_level 3
# config_schedule -effort high
# config_bind -effort high
config_core DSP48 -latency 4

csim_design -clean -O -profile
csynth_design
cosim_design -O -rtl vhdl
export_design -rtl verilog -format ip_catalog -vendor "tu-dresden" -display_name "model_sniffbot" -output $ip_path/model_wrapper_fast_sincos.zip

