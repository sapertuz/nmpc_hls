array unset init_s_ports
array set init_s_ports {
    u_curr_local
    x_curr_local
    xref_local
    f_ind_local
    local_du_min
    local_du_max
    local_u_min
    local_u_max
}

array set costF_ports {
    y
    local_x_curr
    local_xref
    local_f_ind
}

# set accel_block initializeParticles_0
set accel_block evaluateFitnessAndDe_0
foreach port [array get port_name] {
   puts $port
 }
foreach port [array get port_name] {
    set bram_name ${port}_bram
    puts $bram_name
    create_bd_cell -type ip -vlnv xilinx.com:ip:blk_mem_gen:8.4 $bram_name
    # set_property -dict [list CONFIG.Memory_Type {Single_Port_RAM} CONFIG.Enable_B {Always_Enabled} CONFIG.Use_RSTB_Pin {false} CONFIG.Port_B_Clock {0} CONFIG.Port_B_Write_Rate {0} CONFIG.Port_B_Enable_Rate {0}] [get_bd_cells $bram_name]
    set_property -dict [list CONFIG.Memory_Type {True_Dual_Port_RAM} CONFIG.Enable_B {Use_ENB_Pin} CONFIG.Use_RSTB_Pin {true} CONFIG.Port_B_Clock {100} CONFIG.Port_B_Write_Rate {50} CONFIG.Port_B_Enable_Rate {100}] [get_bd_cells $bram_name]
    connect_bd_intf_net [get_bd_intf_pins $bram_name/BRAM_PORTA] [get_bd_intf_pins $accel_block/${port}_PORTA]
}

foreach port [array get init_s_ports] {
    set const_name ${port}_addr
    puts $const_name
    startgroup
    create_bd_cell -type ip -vlnv xilinx.com:ip:xlconstant:1.1 nmpc_solver/$const_name
    set_property -dict [list CONFIG.CONST_WIDTH {32} CONFIG.CONST_VAL {0}] [get_bd_cells nmpc_solver/$const_name]
    connect_bd_net [get_bd_pins nmpc_solver/$const_name/dout] [get_bd_pins nmpc_solver/initializeParticles_0/$port]
    endgroup
}
