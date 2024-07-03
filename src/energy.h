#ifndef _ENERGY_H_ // include guard
#define _ENERGY_H_

#include "global.h"

void Energy_Total_System(void);

float Energy_Anisotropic(const int beadID);

float Energy_Anisotropic_Self(const int beadID);

float Energy_Anisotropic_For_Chain(const int beadID);

float Energy_Anisotropic_For_Range(const int beadID, const int smBead, const int lgBead);

float Energy_Anisotropic_Contiguous_Range(const int beadID, const int smallest_bead, const int largest_bead);

float Energy_Anisotropic_For_List(const int beadID, const int listSize, const int beadList[MAX_BONDS + 1]);

float Energy_Anisotropic_With_List(const int beadID, const int* bead_list, const int list_size);

float Energy_Iso_Ovlp(int const beadType1, int const beadType2, float const xDis);

float Energy_Iso_Cont(int const beadType1, int const beadType2, float const xDis);

float Energy_Iso_fSol(int const beadType);

float Energy_Topo_Angle(const int beadID);

float Energy_OfOvlp_wNeighList(int const beadID, const int* neighList, int const neighNum);

float Energy_OfCont_wNeighList(int const beadID, const int* neighList, int const neighNum);

float Energy_ofSol_wNeighList(const int* neighList, int const neighNum);

float Energy_ofPairs_wNeighList(int const beadID, const int* neighList, int const neighNum);

float Energy_Isotropic(const int beadID);

float Energy_Isotropic_Self(const int beadID);

float Energy_Isotropic_For_Chain(const int beadID);

float Energy_Isotropic_Contiguous_Range(const int beadID, const int smallest_bead, const int largest_bead);

float Energy_Isotropic_With_List(const int beadID, const int* bead_list, const int list_size);

float Energy_Of_Chain(const int chainID);

float Energy_Of_Chain_Self(const int chainID);

float Energy_BiasingPotential(const int beadID);

void Energy_Iso_ForLocal(const int beadID, const int resi, const int* r_pos0, long double* oldEn, long double* newEn,
                         int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs);

void Energy_Iso_ForLocalEquil(const int beadID, const int resi, const int* r_pos0, long double* oldEn,
                              long double* newEn, int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs);

float Energy_OfOvlp_wNeighList_ForChains(int const beadID, const int* neighList, int const neighNum);

float Energy_OfCont_wNeighList_ForChains(int const beadID, const int* neighList, int const neighNum);

float Energy_OfOvlp_wNeighList_ForRange(int const beadID, const int loBead, const int hiBead, const int* neighList,
                                        int const neighNum);

float Energy_OfCont_wNeighList_ForRange(int const beadID, const int loBead, const int hiBead, const int* neighList,
                                        int const neighNum);

float Energy_OfOvlp_wNeighList_ForLists(const int beadID, const int listSize, const int beadList[MAX_BONDS + 1],
                                        const int* neighList, const int neighNum);

float Energy_OfCont_wNeighList_ForLists(const int beadID, const int listSize, const int beadList[MAX_BONDS + 1],
                                        const int* neighList, const int neighNum);

void Energy_Iso_ForChains(const int beadID, long double* oldEn, long double* newEn, int* ovlp_num, int* cont_num,
                          int* ovlp_neighs, int* cont_neighs);

void Energy_Iso_ForChainsEquil(const int beadID, long double* oldEn, long double* newEn, int* ovlp_num, int* cont_num,
                               int* ovlp_neighs, int* cont_neighs);

void Energy_Iso_ForCoLocal(const int thisBeadID, const int otherBeadID, const int* r_pos0, long double* oldEn,
                           long double* newEn, int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs);

void Energy_Iso_ForLists(const int beadIdx, int const listSize, const int beadList[MAX_BONDS + 1],
                         const int beadPos[MAX_BONDS + 1][POS_MAX], long double* oldEn, long double* newEn,
                         int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs);

void Energy_Iso_ForRange(const int beadID, const int smallestBead, const int largestBead, long double* oldEn,
                         long double* newEn, int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs);

void Energy_Iso_ForRangeEquil(const int beadID, const int smallestBead, const int largestBead, long double* oldEn,
                              long double* newEn, int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs);

float Energy_Topo_Angle_ForList(const int bondNum, const int* bondList);

#endif // _ENERGY_H_
