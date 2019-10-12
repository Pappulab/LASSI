gcc -g -pg main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c -o main -lm
cp -u main /work/fdar/LatticeModels/LASSI_GIT/TESTS/GDBTest/
gcc main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c -O3 main -lm
cp -u main runs/
cp -u main runs/MYTEST/
