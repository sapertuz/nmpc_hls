# Vivado Project Generation

## Overview

This folder contains TCL scripts for generating Vivado projects from predefined design blocks. 
Follow the instructions below to create a Vivado project based on the desired design.

## Prerequisites

Make sure you have sourced Vivado in the terminal:

```bash
source <path-to-xilinx>/Vivado/<version>/settings64.sh
```
**IMPORTANT**: The curent version of vivado used for all projects is either `2019.1`, `2019.2`, or `2020.1`

## Instructions

1. Change directory to the Vivado project folder:

    ```bash
    $ cd <path-to-git>/hw/hdl
    ```

2. Run the following command to create a Vivado project:

    ```bash
    vivado -mode batch -notrace -source create_project.tcl -tclargs -n orb_extract \
    -d scripts/designs/design_1.tcl -ip ../hls/vitis_ip_repo -pt xczu3eg-sbva484-1-i \
    -bd avnet-tria:ultra96v2:part0:1.3 -force -synth
    ```
## Some Device Names

| Name          | Part Number              | Board Name                                  |
|---------------|--------------------------|---------------------------------------------|
| Ultra96-v2    | `xczu3eg-sbva484-1-i`    | `avnet.com:ultra96v2:part0:1.2`             |
| Ultra96-v2    | `xczu3eg-sbva484-1-i`    | `avnet-tria:ultra96v2:part0:1.3`            |
| PYNQ-Z1       | `xc7z020clg400-1`        | `www.digilentinc.com:pynq-z1:part0:1.0`     |
