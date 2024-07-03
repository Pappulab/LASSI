gcc -march=x86-64 -O3 -std=gnu11 -lm -o lassi main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c

gcc -march=x86-64 -Og -ggdb -std=gnu11 -lm -o lassi_gdb main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c

cp -u lassi /project/fava/packages/bin/lassi
cp -u lassi_gdb /project/fava/packages/bin/lassi_gdb

cp  lassi ../runs/Test_Simulation/
