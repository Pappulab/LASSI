#ifndef _ENERGY_H_   // include guard
#define _ENERGY_H_

#include <stdio.h>
#include <stdlib.h>

void Energy_Total_System(void);

float Energy_Anisotropic(int beadID);

float Energy_Isotropic(int beadID);

float Energy_Of_Chain(int chainID);

float Energy_InitPotential(int beadID);

#endif // _ENERGY_H_
