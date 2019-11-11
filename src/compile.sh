gcc -g -pg main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c -o lassi_gdb -lm
gcc -O3 main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c -o lassi -lm
cp  lassi ../runs/Test_Simulation/
