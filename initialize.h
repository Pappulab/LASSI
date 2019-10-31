#ifndef _INITIALIZE_H_   // include guard
#define _INITIALIZE_H_

void Memory_Initialization_AtStart(void);

void Global_Array_Initialization_AtStart(void);

void Reset_Global_Arrays(void);

void Initial_Conditions_Simple(void);

void Initial_Conditions_FromFile(void);

void Initial_Conditions_BreakBonds(void);

float Temperature_Function(int mode, long nGen);

void Calculate_Rot_Bias(float CurrentTemp);

#endif // _INITIALIZE_H_
