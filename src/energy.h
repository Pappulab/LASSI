#ifndef _ENERGY_H_   // include guard
#define _ENERGY_H_

#include <stdio.h>
#include <stdlib.h>

void Energy_Total_System(void);

float Energy_Anisotropic(int beadID);

float Energy_Anisotropic_Self(int beadID);

float Energy_Anisotropic_For_Chain(int beadID);

float Energy_Anisotropic_Contiguous_Range(int beadID, int smallest_bead, int largest_bead);

float Energy_Anisotropic_With_List(const int beadID, const int *bead_list, const int list_size);

float Energy_Isotropic(int beadID);

float Energy_Isotropic_Self(int beadID);

float Energy_Isotropic_For_Chain(int beadID);

float Energy_Isotropic_Contiguous_Range(int beadID, int smallest_bead, int largest_bead);

float Energy_Isotropic_With_List(const int beadID, const int *bead_list, const int list_size);

float Energy_Of_Chain(int chainID);

float Energy_Of_Chain_Self(int chainID);

float Energy_InitPotential(int beadID);

#endif // _ENERGY_H_
