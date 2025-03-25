# NMPC Controller HLS Build Instructions

This instruction file explains how to build and extend the Vivado HLS implementation of the NMPC that uses Particle Swarm Optimization (PSO) as its backend. The HLS project supports multiple models, such as the inverted pendulum and the sniffbot (quadrotor drone).

---

## 🛠 Requirements

- **Vivado HLS**: Version 2019.1 to 2020.1 (with the Y2K22 patch applied) is required to build the project.
- **Operating System**: Ubuntu LTS (up to version 2022) is recommended for compatibility.
- **MATLAB**: Optional, used for running simulations and visualizations.
- **Cmake**: Optional, used to build standalone version.

## Folder Structure Overview (Relevant Parts Only)

```
/hw/hls
├── scripts/                  # TCL scripts to drive Vivado HLS
│   └── hls_script_pso.tcl   # Main build script
└── source/
    ├── config/              # Simulation and project config per model
    ├── include/             # Header files, including models and configuration
    ├── src/                 # C++ sources for HLS
    └── matlab/              # MATLAB scripts for simulation and output visualization
```

---

## How to Build the HLS IP Core

```bash
source /opt/Xilinx/Vivado/2020.1/settings64.sh
cd <path-to-repo>/hw/hls
vivado_hls -f scripts/hls_script_pso.tcl
```

### Selecting a Model to Build
The model is selected via a flag defined in the TCL script. Open [`scripts/hls_script_pso.tcl`](scripts/hls_script_pso.tcl) and set one of the following:

```tcl
set flag "-DINVERTED_PENDULUM_CONFIG"
# or
set flag "-DSNIFFBOT_CONFIG"
```

These flags determine which system model is compiled into hardware. Each configuration adjusts the controller size, input/output dimensions, bounds, and cost parameters.

---

## Adding a New Model Configuration

To extend the HLS project with a new system model:

1. **Create a new model header** in `include/models/` (e.g. `hls_myrobot.hpp`).
2. **Add a preprocessor flag** in `config.hpp`:

```cpp
#elif defined(MYROBOT_CONFIG)
#include "models/hls_myrobot.hpp"
```

### Configuration Parameters

As described in the published [original paper](https://dx.doi.org/10.1109/ASAP57973.2023.00038), the following parameters are associated with each model:

| Parameter        | Inverted Pendulum      | Quadrotor Drone       |
|------------------|------------------------|------------------------|
| S                | 10                     | 5                      |
| Dim              | 1                      | 4                      |
| u<sub>min</sub>  | [-50]                  | [-100, -100, -100, -100] |
| u<sub>max</sub>  | [50]                   | [100, 100, 100, 100]     |
| k<sub>max</sub>  | 100                    | 25                     |
| c₁, c₂           | [2.1, 1.0]             | [2.1, 1.0]             |
| ω₀, ωf     | [0.9, 0.1]             | [0.9, 0.1]             |
| Δu             | [50]                   | [20, 20, 20, 20]       |
| T<sub>s</sub>     | 0.1                    | 0.05                   |
| x<sub>min</sub>   | [-2, -, -, -]         | [-, 0, -, -π/4, -π/4, ...] |
| x<sub>max</sub>   | [2, -, -, -]          | [-, 0, -, π/4, π/4, ...]   |
| u<sub>ss</sub>    | [0]                    | [0, 0, 0, 0]           |
| Q                | [2, -, 1, -]           | [1.0, 1.0, ..., 3.0]   |
| Q<sub>f</sub>     | [0.1, 0.015, 0.1, -]   | [10, 10, ..., 30]      |
| R                | [0]                    | [5e-4, 5e-4, 5e-4, 5e-4]|

**Note**: `-` represents irrelevant entries for that specific model.

---

## 📁 MATLAB Scripts

The folder `source/matlab/` contains scripts to:
- Run controller simulations
- Visualize results
- Compare PSO vs KPSO

Examples:
- `inverted_pendulum_kpso.m`
- `drone_pso.m`

---
