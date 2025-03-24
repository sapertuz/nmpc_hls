namespace eval _tcl {
proc get_script_folder {} {
   set script_path [file normalize [info script]]
   set script_folder [file dirname $script_path]
   return $script_folder
}
}
variable script_folder
set script_folder [_tcl::get_script_folder]

proc lshift listVar {
  upvar 1 $listVar L
  set r [lindex $L 0]
  set L [lreplace $L [set L 0] 0]
  return $r
}

 #-------------------------------------------------------
 # Process command line arguments
 #------------------------------------------------------- 
set error 2
set help 0
set synth_flag 0
set force_flag 0
set gui_flag 0

set project_name ""
set design_file ""
set constraint_file ""
# set hw_dir "./hw_files"

set source_files ""
set top_module ""

set proj_dir "."
set ip_repo_dir ""

set board_name ""
set part_name ""

# if {[llength $argc] == 0} { incr help }; # Uncomment if necessary

if {$argc > 0} {
    while {[llength $argv]} {
        set flag [lshift argv]
        switch -- $flag {
        -d -
        -design {
            set design_file [lshift argv]
            incr error -1
        }
        -n -
        -name {
            set project_name [lshift argv]
            incr error -1
        }
        -pt -
        -part {
            set part_name [lshift argv]
        }
        -bd -
        -board {
            set board_name [lshift argv]
        }
        -ip -
        -iprepo {
            set ip_repo_dir [lshift argv]
        }
        -c -
        -const_file {
            set constraint_file [lshift argv]
        }
        -top{
            set top_module [lshift argv]
        }
        -s -
        -source-files {
            while {[llength $argv] > 0 && ![string match -* [lindex $argv 0]]} {
                lappend source_files [lshift argv]
            }
        }
        -g -
        -gui {
            incr gui_flag
        }
        -synth {
            incr synth_flag
        }
        -force {
            incr force_flag
        }
        -h -
        -help {
            incr help
        }
        default {
            if {[string match "-*" $flag]} {
            puts " ERROR - option '$flag' is not a valid option."
            incr error
            } else {
            puts "ERROR - option '$flag' is not a valid option."
            incr error
            }
        }
        }
        # puts $flag
    }
} else  {
    incr help
}

if {$help} {
    # set callerflag [lindex [info level [expr [info level] -1]] 0]
    # <-- HELP
    puts [format "
Usage: create_project.tcl
    \[-design|-d <path to tcl file with design>\]
    \[-name|-n <project name>\]
    \[-board|-bd <board name>\]
    \[-part|-pt <part name>\]
    \[-iprepo|-ip <path to ip reposotory folder>\]
    \[-c|-const_file <path to constraint file>\]
    \[-source-files <path to source files>\]
    \[-top-module <name of top module, in case default is not correct>\]
    \[-g|-gui <flag to open the gui or not>\]
    \[-synth <flag that synthesizes and generate xsa>\]
    \[-force <flag to delete folder with same name>\]
    \[-help|-h <this menu>\]                          

Description: . 
    With this file you can create a vivado design having the block design tcl scripts and other sources.

Example:
vivado -mode tcl -notrace -source create_project.tcl -tclargs -n <project-name> -d scripts/<script-name>
    " ]
    # HELP -->
    return -code ok {}
}

# Check validity of arguments. Increment $error to generate an error

if {$error > 0} {
return -code error {Oops, something is not correct, remember that the parameter -name|-n and -design|-d is mandatory}
}

# fake vivado version to make it backwards compatible
set scripts_vivado_version [version -short]
if {![regexp {^2019\.1$|^2019\.2$|^2020\.1$} $scripts_vivado_version]} {
    return -code error "ERROR: Unsupported Vivado version ${scripts_vivado_version}. Supported versions are 2019.1, 2019.2, and 2020.1."
}
puts "INFO: Vivado version is ${scripts_vivado_version}"

# Do something
set bd_design_name "design_${project_name}"

# Create project 
if {$gui_flag > 0} {
start_gui
}


## Create project
if {$force_flag} {
    create_project $project_name $proj_dir/$project_name -force     
} else {
    create_project $project_name $proj_dir/$project_name
}

## Set project properties
set obj [current_project]
if {[string length $part_name] > 0} {
    set_property part -value $part_name -objects $obj
}
if {[string length $board_name] > 0} {
    set_property board_part -value $board_name -objects $obj
}

# Add repo (if included)
if {[string length $ip_repo_dir] > 0} {
    set_property ip_repo_paths $ip_repo_dir $obj
}

# Create and construct system
update_ip_catalog
create_bd_design $bd_design_name

# Add VHDL sources
if {[string length $source_files] > 0} {
    add_files -norecurse ${source_files}
}

# create design
source $design_file

# Add constraints files (if included)
if {[string length $constraint_file] > 0} {
    add_files -fileset constrs_1 -norecurse $constraint_file
}

# Validate design
validate_bd_design
regenerate_bd_layout

# Wrap and save design
make_wrapper -files [get_files $bd_design_name.bd] -top
add_files -norecurse $proj_dir/$project_name/$project_name.gen/sources_1/bd/$bd_design_name/hdl/${bd_design_name}_wrapper.v
update_compile_order -fileset sources_1
update_compile_order -fileset sim_1
save_bd_design

# Set Top
if {[string length $top_module] == 0} {
    set top_module "design_${project_name}_wrapper"
}
set_property top $top_module [current_fileset]
update_compile_order -fileset sources_1

if {$synth_flag} {
    
    set synth_name synth_1
    set impl_name impl_1
    
    # Synthesize
    launch_runs $synth_name -j 4
    ## Wait run
    wait_on_run $synth_name

    # Implement
    launch_runs $impl_name  -jobs 4 -to_step write_bitstream
    ## Wait run
    wait_on_run $impl_name
    
    ## Export HW
    update_compile_order -fileset sources_1
    write_hw_platform -fixed -include_bit -force -file ./$project_name/$project_name.xsa 
    
}

return -code ok {}