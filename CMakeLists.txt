cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)

if (MSVC)
add_definitions(-D__LITTLE_ENDIAN__)
add_definitions(-D_USE_MATH_DEFINES)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
add_definitions("/wd4309")
set(CMAKE_SUPPRESS_REGENERATION true)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
endif ()

include_directories("stk")
add_executable(dpcmc 
    main.cpp
    fft4g.c
    stk/FileRead.cpp
    stk/FileWrite.cpp
    stk/Stk.cpp
)
