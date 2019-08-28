#ifndef _STRUCTURE_H_   // include guard
#define _STRUCTURE_H_

#include "global.h"

int LtInd(int i, int j, int k);
int LtIndV(int xArr[POS_MAX]);
float dist(int bead1, int bead2);
float distf(float f1[], float f2[]);
float distInt(int f1[], int f2[]);
float vectorMag(const int f1[]);
int check_linker_constraint(int beadID, int tmpR[]);
int ShakeConstraint(int beadID, int tmpR[MAX_VALENCY][POS_MAX]);
int check_structure_topo(void);
void avg_rdf_split(void);
void CalcTotGyrRad(void);
void CenterMySystem(void);
#endif // _STRUCTURE_H_
