#ifndef _CLUSTER_H_   // include guard
#define _CLUSTER_H_

#include "global.h"

int Clus_ChainNetwork_General(int chainID, int *naList);
int Clus_ChainNetwork_ForTotal(int chainID, int *naList);
int Clus_LimitedCluster(int chainID, int *naList);
void Clus_Distribution_Avg(int *naList);
int Clus_SecondLargestCluster(int *naList, int **naCluster);
void Clus_TotalAnalysis(int *naList, int **naCluster);

#endif // _CLUSTER_H_
