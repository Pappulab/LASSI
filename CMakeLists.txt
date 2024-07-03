cmake_minimum_required(VERSION 3.10)
project(LASSI_SRC C)

set(CMAKE_CXX_COMPILER "gcc")
set(CMAKE_C_COMPILER "gcc")

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Og -g3 -march=native -mtune=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og -g3 -march=native -mtune=native -Wall -Wextra -pedantic")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -mtune=native -fno-omit-frame-pointer")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mtune=native -fno-omit-frame-pointer")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -march=native -mtune=native")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -mtune=native")


include_directories(./src)

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

target_link_libraries(lassi m)