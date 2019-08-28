cmake_minimum_required(VERSION 3.14)
project(LASSI_SRC C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(LASSI_SRC
        cluster.c
        cluster.h
        energy.c
        energy.h
        global.h
        initialize.c
        initialize.h
        main.c
        mcmove.c
        mcmove.h
        parsekey.c
        parsekey.h
        print.c
        print.h
        structure.c
        structure.h)
