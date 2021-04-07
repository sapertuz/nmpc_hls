Estabas a punto de compatibilizar con hls poniendo el tamanho de los vectores en las llamadas de las funciones pero te acordaste que segun la parametrizacion N, Nc y Nu cambian... tenemos que resolver esto para pensar en su implementacion en un futuro.

# Instructions
```
cd build && cmake ..
cmake --build .
cd ..
./build/nmpc_pso config/inverted_pendulum/project_config.txt
```