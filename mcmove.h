#ifndef _MCMOVE_H_   // include guard
#define _MCMOVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int MC_Step(float fMCTemp);
int MC_Step_Equil(float fMCTemp);

int RotMCMove(int beadID, float MyTemp);
int LocalMCMove(int beadID, float MyTemp);
int SlitherMCMove(int chainID, float MyTemp);
int TransMCMove(int chainID, float MyTemp);
int ClusMCMove(float MyTemp);
int SmallClusMCMove(int chainID, float MyTemp);
int DbPvtMCMove(int beadID, float MyTemp);
int CoLocalMove(int beadID, float MyTemp);
int ShakeMove(int beadID, float MyTemp);
int PivotMCMove(int chainID, float MyTemp);
int BranchedRotMCMove(int chainID, float MyTemp);

int LocalMCMove_Equil(int beadID, float MyTemp);
int SlitherMCMove_Equil(int chainID, float MyTemp);
int TransMCMove_Equil(int chainID, float MyTemp);
int ShakeMove_Equil(int beadID, float MyTemp);
int PivotMCMove_Equil(int chainID, float MyTemp);
int PivotMCMove_Equil2(int chainID, float MyTemp);
int BranchedRotMCMove_Equil(int chainID, float MyTemp);

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
