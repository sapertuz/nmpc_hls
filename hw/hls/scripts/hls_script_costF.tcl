############################################################
## This file is generated automatically by Vitis HLS.
## Please DO NOT edit it.
## Copyright 1986-2021 Xilinx, Inc. All Rights Reserved.
############################################################
set ip_path "${workspace}/vitis_ip_repo"

set prj_name test_costF
set prj_top cost_function_wrapper

set workspace [pwd]
set workspace [file dirname $workspace]

set src_path ${workspace}/src
set incl_path ${workspace}/include
set main_name "main_hls_system"

set c_flags "-DSNIFFBOT_CONFIG -DUSE_FAST_SIN_COS -I${incl_path} -I${incl_path}/models -std=c++11 -Wno-unknown-pragmas"
set csim_tb_flags "-DSNIFFBOT_CONFIG -I${incl_path} -I${incl_path}/models -std=c++11 -DDEBUG_SYSTEM -Wno-unknown-pragmas"

open_project $prj_name
set_top $prj_top

add_files $src_path/$main_name.cpp -cflags $c_flags -csimflags $csim_tb_flags
# foreach file [glob -dir $incl_path *.hpp] {
#     add_files $file
# }

add_files -tb ${src_path}/${main_name}.cpp -cflags $c_flags -csimflags $csim_tb_flags
add_files -tb ${src_path}/aux_functions.cpp -cflags $c_flags -csimflags $csim_tb_flags

open_solution "solution_system" -flow_target vivado

set_part {xczu3eg-sbva484-1-i}
create_clock -period 10 -name default
config_compile -pipeline_loops 6
config_interface -m_axi_addr64=0
config_export -display_name sniffbot_costF -format ip_catalog -output $ip_path/sniffbot_costF.zip -rtl verilog -vendor tu-dresden -version 1.0

# config_core DSP48 -latency 4

csim_design -clean -O -profile
csynth_design
# cosim_design -O -rtl vhdl
# export_design -format ip_catalog
