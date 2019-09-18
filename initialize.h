#ifndef _INITIALIZE_H_   // include guard
#define _INITIALIZE_H_

void MemoryAndBookkeepingInit(void);
void Initialize_Dilute(void);
void Initialize_Dilute_SM(void);
float Temperature_Function(int mode, long nGen);
void Calculate_Rot_Bias(float CurrentTemp);
void Reset_Global_Arrays(void);
#endif // _INITIALIZE_H_
