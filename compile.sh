gcc -Og -g -pg main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c -o main -lm
cp -u main /work/fdar/LatticeModels/LASSI_GIT/TESTS/GDBTest/
gcc -O3 main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c -o main -lm
cp -u main runs/
cp -u main runs/MYTEST/
