############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2019 Xilinx, Inc. All Rights Reserved.
############################################################
set ip_path "/home/chello/Documents/Vivado_WS/vitis_ip_repo"

set prj_name pseudorand_stream
set prj_top rand_wrapper

set workspace [pwd]
set workspace [file dirname $workspace]

set src_path ${workspace}/src
set incl_path ${workspace}/include
set main_name "main_hls_pseudorand" 

set c_flags "-I${incl_path} -D__VITIS__ -std=c++11 -Wno-unknown-pragmas"
set csim_tb_flags "-I${incl_path} -D__VITIS__ -std=c++11 -Wno-unknown-pragmas"

open_project $prj_name
set_top $prj_top

add_files $src_path/$main_name.cpp -cflags $c_flags -csimflags $csim_tb_flags
add_files -tb $src_path/$main_name.cpp -cflags $c_flags -csimflags $csim_tb_flags

open_solution "solution_pseudorand" -flow_target vivado

set_part {xczu3eg-sbva484-1-i}
create_clock -period 10 -name default
config_export -display_name $prj_name -format ip_catalog -output $ip_path/$prj_name.zip -rtl verilog -vendor tu-dresden -version 1.0
config_interface -m_axi_addr64=0 -s_axilite_auto_restart_counter 1
config_rtl -reset state

# config_core DSP48 -latency 4

# csim_design -clean -O -profile
csynth_design
# cosim_design -O -rtl vhdl
export_design -format ip_catalog
