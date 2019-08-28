gcc -g -pg main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c -o main -lm
cp -u main runs/GDB
gcc main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c -o main -lm
cp -u main runs/
cp -u main runs/MYTEST/
