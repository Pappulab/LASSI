# Only uncomment if you use Intel or AMD
#gcc -march=x86-64 
gcc -O3 -std=gnu11 -lm -o lassi main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c

# Only uncomment if you use Intel or AMD
#gcc -march=x86-64 
gcc -Og -ggdb -std=gnu11 -lm -o lassi_gdb main.c parsekey.c print.c initialize.c structure.c energy.c cluster.c mcmove.c

# Specific to the Pappu Lab's computing infrastructure
#cp -u lassi /project/fava/packages/bin/lassi
#cp -u lassi_gdb /project/fava/packages/bin/lassi_gdb

cp lassi ../runs/Test_Simulation/
