#ifndef _CLUSTER_H_   // include guard
#define _CLUSTER_H_

#include "global.h"

int chain_network(int chainID, int naList[MAX_CHAINS]);
int chain_network_for_tot(int chainID, int *naList);
int chain_network_small(int chainID, int naList[MAX_CHAINS]);
void avg_clus_dist(int naList[MAX_CHAINS]);
int clus_network_analysis(int naList[MAX_CHAINS], int naCluster[MAX_CHAINS][MAX_CHAINS]);

#endif // _CLUSTER_H_
