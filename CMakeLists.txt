cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)

include_directories("stk")
add_executable(dpcmc 
    main.cpp
    fft4g.c
    stk/FileRead.cpp
    stk/FileWrite.cpp
    stk/Stk.cpp
)
