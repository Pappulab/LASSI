#ifndef _MCMOVE_H_ // include guard
#define _MCMOVE_H_

#include "global.h"

int MC_Step(float fMCTemp);

int MC_Step_Equil(float fMCTemp);

int Move_Rot(int beadID, float MyTemp);

int Move_Local(int beadID, float MyTemp);

int Move_Snake(int chainID, float MyTemp);

int Move_Trans(int chainID, float MyTemp);

int Move_Clus_Network(float MyTemp);

int Move_SmallClus_Network(int chainID, float MyTemp);

int Move_DbPvt(int beadID, float myTemp);

int Move_CoLocal(int thisBeadID, float MyTemp);

int Move_MultiLocal(int beadID, float MyTemp);

int Move_Pivot(int chainID, float MyTemp);

int Move_BranchedRot(int chainID, float MyTemp);

int Move_SmallClus_Proximity(const int chainID, const float myTemp);

int Move_Clus_Proximity(const float myTemp);

int Move_Local_Equil(int beadID, float MyTemp);

int Move_Snake_Equil(int chainID, float MyTemp);

int Move_Trans_Equil(int chainID, float MyTemp);

int Move_MultiLocal_Equil(int beadID, float MyTemp);

int Move_Pivot_Equil(int chainID, float MyTemp);

int Move_BranchedRot_Equil(int chainID, float MyTemp);

int Check_ChainDisp(const int chainID, const int* vec_disp);

int Check_MoveBeadTo(const int* newPos);

void OP_System_DispChain(int chainID, const int* movR);

void OP_System_DispChain_ForTrans(int chainID, const int* movR);

void OP_System_RestoreChain(int chainID);

void OP_System_RestoreChain_ForTrans(int chainID);

void OP_CopyBead(int* copy_arr, const int* orig_arr);

void OP_CopyBeadsToOld(const int firstB, const int lastB);

void OP_RestoreBeadsFromOld(const int firstB, const int lastB);

void OP_RestoreChain_ForSnake(const int fB, const int lB);

void OP_System_MoveBeadTo(const int beadID, const int* newPos);

void OP_System_MoveBeadTo_Inv(int beadID);

void OP_MoveBeadTo_ForMTLocal(int beadID, const int* newPos);

void OP_SwapBeads(int bead1, int bead2);

void OP_Rotation(int PivotM, int beadID, int* tmpR);

void Rot_X_90(int beadID, const int tmpR[]);

void Rot_X_180(int beadID, const int tmpR[]);

void Rot_X_270(int beadID, const int tmpR[]);

void Rot_Y_90(int beadID, const int tmpR[]);

void Rot_Y_180(int beadID, const int tmpR[]);

void Rot_Y_270(int beadID, const int tmpR[]);

void Rot_Z_90(int beadID, const int tmpR[]);

void Rot_Z_180(int beadID, const int tmpR[]);

void Rot_Z_270(int beadID, const int tmpR[]);

void OP_ShuffleRotIndecies(void);

void OP_ShuffleArray(const int arr_size, int* dum_arr);

int Check_RotStates_wNeighList(int const beadID, int const resi, const int* neighList, int const neighNum);

int Check_RotStatesOld(int const beadID, int const resi, float const MyTemp);

int Check_RotStatesNew(int const beadID, int const resi, float const MyTemp);

void OP_NormalizeRotState(const int beadVal, const int CandNums);

int OP_PickRotState(int CandNums);

void OP_System_Snake_SlitherFwd(const int firstB, const int lastB, const int* r_posNew);

void OP_System_Snake_SlitherBck(const int firstB, const int lastB, const int* r_posNew);

void OP_Beads_CopyBeadsInListToOld(const int listSize, const int* beadList);

void OP_Beads_CopyBeadsInListFromOld(const int listSize, const int* beadList);

void OP_Beads_CopyBeadsInListToPosList(const int listSize, const int* beadList, int (*beadPos)[POS_MAX]);

void OP_Lattice_EmptySitesForListOfBeads(const int listSize, const int* beadList);

void OP_Lattice_PlaceBeadsInList(const int listSize, const int* beadList);

void OP_Lattice_EmptySitesForListOfPos(const int listSize, const int (*beadPos)[POS_MAX]);

void OP_Beads_MoveBeadsInListToPos(const int listSize, const int* beadList, const int (*newPos)[POS_MAX]);

void OP_Inv_MoveBeads_InList_ToPos(const int listSize, const int* beadList);

lLDub OP_GenMHValue(lLDub fRos, lLDub bRos, lLDub Delta_En, lLDub Cur_Temp);

lLDub MC_RosenbluthSampling_ForLocal_AtOld(const int beadID, const int resi, long double* oldEn, const int neigh_num);

lLDub MC_RosenbluthSampling_ForLocal_AtNew(const int beadID, const int resi, int* bead_part, long double* newEn,
                                           const int neigh_num);

lLDub MC_RosenbluthSampling_ForChains_AtOld(const int beadID, const int resi, long double* oldEn, const int neigh_num);

lLDub MC_RosenbluthSampling_ForChains_AtNew(const int beadID, const int resi, int* bead_part, long double* newEn,
                                            const int neigh_num);

lLDub MC_RosenbluthSampling_ForLists_AtOld(const int beadIdx, const int listSize, const int beadList[MAX_BONDS + 1],
                                           lLDub* oldEn, const int neigh_num);

lLDub MC_RosenbluthSampling_ForLists_AtNew(const int beadIdx, const int listSize, const int beadList[MAX_BONDS + 1],
                                           int* bead_part, lLDub* oldEn, const int neigh_num);

lLDub MC_RosenbluthSampling_ForRange_AtOld(const int beadID, const int resi, const int smallestBead,
                                           const int largestBead, lLDub* oldEn, const int neigh_num);

lLDub MC_RosenbluthSampling_ForRange_AtNew(const int beadID, const int resi, int* bead_part, lLDub* newEn,
                                           const int neigh_num);

void OP_System_MoveBeadsInListToPos(const int listSize, const int* beadList, const int (*newPos)[POS_MAX]);

void OP_Beads_BreakBondsInList(const int listSize, const int* beadList);

void OP_Beads_BreakBond(const int beadID);

void OP_Beads_RestoreBondsInList(const int listSize, const int* beadList);

void OP_Beads_RestoreBond(const int beadID);

#endif // _MCMOVE_H_
