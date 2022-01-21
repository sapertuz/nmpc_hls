## Test system and memory check with debug flags
```
$ g++ main_hls_system.cpp -g -Og -std=c++11 -DDEBUG_SYSTEM -o test
$ valgrind ./test
```