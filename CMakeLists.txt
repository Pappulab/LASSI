cmake_minimum_required(VERSION 3.14)
project(lassi C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Og")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og")

include_directories(.)

add_executable(lassi
        src/cluster.c
        src/cluster.h
        src/energy.c
        src/energy.h
        src/global.h
        src/initialize.c
        src/initialize.h
        src/main.c
        src/mcmove.c
        src/mcmove.h
        src/parsekey.c
        src/parsekey.h
        src/print.c
        src/print.h
        src/structure.c
        src/structure.h)
