# Standalone HLS Project Build and Memory Testing Guide

## Build HLS Project Standalone
The HLS project is independent of HLS-specific methods or classes, enabling standalone builds using `CMakeLists.txt`. To build and run:

```
cd <path-to-repo>/hw/hls/source
mkdir build && cd build
cmake ..
make
./main_hls_model or ./main_hls_pseudorand or ./main_hls_pso or ./main_hls_system
```

This ensures the project compiles and runs without HLS-specific dependencies.

## Test memory checks
To test with debug flags and perform memory checks:

```
g++ main_hls_system.cpp -g -Og -std=c++11 -DDEBUG_SYSTEM -o test
valgrind ./test
G_SLICE=always-malloc G_DEBUG=gc-friendly valgrind -v --tool=memcheck --leak-check=full --num-callers=40 --log-file=valgrind.log ./test_pso.out ./config/sniffbot/project_config.txt ./config/sniffbot/simulation_config_ring.txt
```
