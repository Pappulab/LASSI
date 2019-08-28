#ifndef _ENERGY_H_   // include guard
#define _ENERGY_H_

#include <stdio.h>
#include <stdlib.h>

void Calc_Tot_Energy();
float energy_SC(int beadID);
float energy_cont_and_ovlp(int beadID);
float energy_chain(int chainID);

#endif // _ENERGY_H_
