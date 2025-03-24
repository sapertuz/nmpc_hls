## Test system and memory check with debug flags
```
$ g++ main_hls_system.cpp -g -Og -std=c++11 -DDEBUG_SYSTEM -o test
$ valgrind ./test
$ G_SLICE=always-malloc G_DEBUG=gc-friendly  valgrind -v --tool=memcheck --leak-check=full --num-callers=40 --log-file=valgrind.log ./test_pso.out ./config/sniffbot/project_config.txt ./config/sniffbot/simulation_config_ring.txt
```