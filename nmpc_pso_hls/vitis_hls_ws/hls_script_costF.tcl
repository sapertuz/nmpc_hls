############################################################
## This file is generated automatically by Vitis HLS.
## Please DO NOT edit it.
## Copyright 1986-2021 Xilinx, Inc. All Rights Reserved.
############################################################
set ip_path "/home/chello/Documents/Vivado_WS/vitis_ip_repo"

open_project nmpc_hls_costF
set_top cost_function_wrapper
add_files ../main_hls_system.cpp -cflags "-DSNIFFBOT_CONFIG -std=c++11"
add_files ../hls_system.hpp
add_files ../hls_sniffbot.hpp
add_files ../hls_inverted_pendulum.hpp
add_files ../fast_sin_cos.hpp
add_files ../aux_functions.hpp
add_files -tb ../main_hls_system.cpp -cflags "-DDEBUG_SYSTEM -DSNIFFBOT_CONFIG -std=c++11 -Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"

open_solution "solution_system" -flow_target vivado

set_part {xc7z020-clg400-1}
create_clock -period 10 -name default
config_compile -pipeline_loops 6
config_interface -m_axi_addr64=0
config_export -display_name sniffbot_costF -format ip_catalog -output $ip_path/sniffbot_costF.zip -rtl verilog -vendor tu-dresden -version 1.0

config_core DSP48 -latency 4

csim_design -clean -O -profile
csynth_design
# cosim_design -O -rtl vhdl
export_design -format ip_catalog
