#ifndef _PRINT_H_   // include guard
#define _PRINT_H_

#include <stdio.h>
#include <stdlib.h>

void print_log(long nGen);
void write_mcmove(char* filename, long nGen, float fMCTemp);
void write_energy(char* filename, long nGen);
void print_key(void);
void write_trajectory(char* filename, long nGen);
void write_topofile(char* filename);
void write_SysProp(char* filename);
void Print_Data(long nGen, int run_it);


#endif // _PRINT_H_
