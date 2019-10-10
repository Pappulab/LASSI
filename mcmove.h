#ifndef _MCMOVE_H_   // include guard
#define _MCMOVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int MC_Step(float fMCTemp);
int MC_Step_Equil(float fMCTemp);

int Move_Rot(int beadID, float MyTemp);
int Move_Local(int beadID, float MyTemp);
int Move_Snake(int chainID, float MyTemp);
int Move_Trans(int chainID, float MyTemp);
int Move_Clus(float MyTemp);
int Move_SmallClus(int chainID, float MyTemp);
int Move_DbPvt(int beadID, float MyTemp);
int Move_CoLocal(int beadID, float MyTemp);
int Move_MultiLocal(int beadID, float MyTemp);
int Move_Pivot(int chainID, float MyTemp);
int Move_BranchedRot(int chainID, float MyTemp);

int Move_Local_Equil(int beadID, float MyTemp);
int Move_Snake_Equil(int chainID, float MyTemp);
int Move_Trans_Equil(int chainID, float MyTemp);
int Move_MultiLocal_Equil(int beadID, float MyTemp);
int Move_Pivot_Equil(int chainID, float MyTemp);
int PivotMCMove_Equil2(int chainID, float MyTemp);
int Move_BranchedRot_Equil(int chainID, float MyTemp);

void disp_chain(int chainID, const int movR[]);
void trans_disp_chain(int chainID, const int movR[]);
void restore_chain(int chainID);
void trans_restore_chain(int chainID);
int check_disp_chain(int chainID, const int tR[]);
void move_bead_to(int beadID, const int newPos[]);
void undo_move_bead_to(int beadID);
void move_bead_to_shake(int beadID, const int newPos[]);
int check_move_bead_to(int newPos[]);
void swap_beads(int bead1, int bead2);

void RotOperation(int PivotM, int beadID, int tmpR[]);
void Rot_X_90(int beadID,  const int tmpR[]);
void Rot_X_180(int beadID, const int tmpR[]);
void Rot_X_270(int beadID, const int tmpR[]);
void Rot_Y_90(int beadID,  const int tmpR[]);
void Rot_Y_180(int beadID, const int tmpR[]);
void Rot_Y_270(int beadID, const int tmpR[]);
void Rot_Z_90(int beadID,  const int tmpR[]);
void Rot_Z_180(int beadID, const int tmpR[]);
void Rot_Z_270(int beadID, const int tmpR[]);

void ShuffleRotIndecies(void);
int CheckRotStatesOld(int beadID, int resi, float MyTemp);
int CheckRotStatesNew(int beadID, int resi, float MyTemp);
void NormalizeRotState(int beadVal, int CandNums);
int PickRotState(int CandNums);
#endif // _MCMOVE_H_
