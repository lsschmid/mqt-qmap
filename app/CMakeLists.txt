
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_COMPILER "/tool/gcc/install/8.1.0/bin/gcc")
set(CMAKE_CXX_COMPILER "/tool/gcc/install/8.1.0/bin/g++")
set(CMAKE_CXX_FLAGS "-std=c++11 -static-libstdc++ -O3 -fopenmp -fpermissive -DNDEBUG -DIL_STD -pg -Wno-ctor-dtor-privacy")

add_executable(NaMain main.cpp)
target_link_libraries(NaMain qmap_namap_lib stdc++fs)
