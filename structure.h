#ifndef _STRUCTURE_H_   // include guard
#define _STRUCTURE_H_

#include "global.h"

int Lat_Ind_FromCoords(int i, int j, int k);
int Lat_Ind_FromVec(int *xArr);

float Dist_BeadToBead(int bead1, int bead2);
float Dist_PointTotPoint_Float(float *f1, float *f2);
float Dist_PointToPoint(int *f1, int *f2);
float Dist_VecMag(const int *f1);

int Check_LinkerConstraint(int beadID, int *tmpR);
int Check_MTLinkerConstraint(int beadID, int **tmpR);
int Check_System_Structure(void);

void RDF_ComponentWise_Avg(void);
void GyrTensor_GyrRad_Avg(void);

#endif // _STRUCTURE_H_
