############################################################
## This file is generated automatically by Vitis HLS.
## Please DO NOT edit it.
## Copyright 1986-2021 Xilinx, Inc. All Rights Reserved.
############################################################
set prj_name nmpc_hls_pso
set prj_top nonlinear_solver_wrapper

set ip_path "/home/chello/Documents/Vivado_WS/vitis_ip_repo"
set workspace [pwd]
set workspace [file dirname $workspace]

set src_path ${workspace}/src
set incl_path ${workspace}/include
set module_file "hls_nonlinear_solver" 
set main_name "main_hls_pso" 

set arg_str "${workspace}/config/sniffbot/project_config.txt ${workspace}/config/sniffbot/simulation_config_ring.txt"
set c_flags "-DSNIFFBOT_CONFIG -DPSO_CONFIG -DUSE_FAST_SIN_COS -I${incl_path} -I${incl_path}/models -std=c++11 -Wno-unknown-pragmas"
set csim_tb_flags "-DSNIFFBOT_CONFIG -DPSO_CONFIG -I${incl_path} -I${incl_path}/models -std=c++11  -DDEBUG_FILE -DPRINT_TO_TERMINAL -Wno-unknown-pragmas"

open_project $prj_name

set_top $prj_top
# foreach file [glob -dir $incl_path *.hpp] {
#     add_files $file
# }

add_files ${incl_path}/${module_file}.hpp -cflags $c_flags -csimflags $csim_tb_flags
add_files ${src_path}/${module_file}.cpp -cflags $c_flags -csimflags $csim_tb_flags
add_files -tb ${src_path}/${main_name}.cpp -cflags $csim_tb_flags -csimflags $csim_tb_flags
add_files -tb ${src_path}/aux_functions.cpp -cflags $csim_tb_flags -csimflags $csim_tb_flags

open_solution "solution_system" -flow_target vivado

set_part {xczu3eg-sbva484-1-i}
create_clock -period 10 -name default
# config_compile -pipeline_loops 6
config_interface -m_axi_addr64=0
config_export -display_name sniffbot_nmpc -format ip_catalog -output $ip_path/sniffbot_nmpc.zip -rtl verilog -vendor tu-dresden -version 1.0

# config_core DSP48 -latency 4

csim_design -argv $arg_str -clean -O -profile
csynth_design
# cosim_design -O -rtl vhdl
# export_design -format ip_catalog
