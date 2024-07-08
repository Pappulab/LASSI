#ifndef _CLUSTER_H_ // include guard
#define _CLUSTER_H_

#include "global.h"

int Clus_Network_ChainCluster_General(int const chainID);

int Clus_Network_ChainCluster_ForTotal(int const chainID);

int Clus_Network_LimitedCluster(int const chainID);

void Clus_Network_Distribution_Avg(void);

void Clus_Network_Distribution_MolWise_Avg(void);

int Clus_Network_SecondLargestCluster(void);

void Clus_Network_TotalAnalysis(void);

void Clus_Network_MolWise_LargestClusters(void);

int Clus_Proximity_ChainCluster_ForTotal_All(int const chainID);

int Clus_Proximity_ChainCluster_ForTotal_IntOnly(int const chainID);

void Clus_Proximity_TotalAnalysis(void);

int Clus_Proximity_SecondLargestCluster(void);

int Clus_Proximity_LimitedCluster_IntOnly(int const chainID);

int Clus_Proximity_LimitedCluster_IntOnly_Check(int const chainID, int const* OldList);

int Clus_Proximity_LimitedCluster_All(int const chainID, int* clusList);

int Clus_Proximity_LimitedCluster_All_Check(int const chainID, int const* OldList, int* NewList);

void Clus_Proximity_Distribution_Avg(void);

void Clus_Proximity_Distribution_IntOnly_MolWise_Avg(void);

void Clus_Proximity_Distribution_All_MolWise_Avg(void);

void Clus_Proximity_IntOnly_MolWise_LargestClusters(void);

void ClusAnalysis_Perform_Analysis(const int nMode);

int Clus_Perform_ChainCluster_ForTotal(int const chainID);

void Clus_Perform_MolWise_LargestClusters(void);

void Clus_Find_LargestClusters(void);

int ClusUtil_GetOvlpNeighborBeads_ForBead(const int beadID, int* restrict neighList);

int ClusUtil_GetOvlpNeighborChains_ForBead(const int beadID, int* restrict neighList);

int ClusUtil_AddOvlpCluster_OfBead(const int beadID, char* restrict caTotClusTable, int* restrict clusList,
                                   int* restrict clusSize);

int ClusUtil_AddOvlpCluster_OfChain(const int chainID, char* restrict caTotClusTable, int* restrict clusList,
                                    int* restrict clusSize);

int ClusUtil_OvlpCluster_OfChain(const int chainID, char* restrict caTotClusTable, int* restrict naClusList);

int ClusUtil_OvlpCluster_OfChain_wMaxSize(const int chainID, char* restrict caTotClusTable, int* restrict clusList,
                                          int const maxSize);

int ClusUtil_AddOvlpCluster_OfBead_CheckForSame(const int beadID, const char* restrict const caTotClusTable);

int ClusUtil_AddOvlpCluster_OfChain_CheckForSame(const int chainID, const char* restrict const caTotClusTable);

int ClusUtil_OvlpCluster_OfChain_CheckForSame(const char* const caTotClusTable, const int* const clusList,
                                              int const clusSize);

int ClusUtil_AddAnisoCluster_OfBead(const int beadID, char* restrict caTotClusTable, int* restrict clusList,
                                    int* restrict clusSize);

int ClusUtil_AddAnisoCluster_OfChain(const int chainID, char* restrict caTotClusTable, int* restrict clusList,
                                     int* restrict clusSize);

int ClusUtil_AnisoCluster_OfChain(const int chainID, char* restrict caTotClusTable, int* restrict clusList);

int ClusUtil_AnisoCluster_OfChain_wMaxSize(const int chainID, char* restrict caTotClusTable, int* restrict clusList,
                                           int const maxSize);

int ClusUtil_NextUnvisitedChain(int const nChainID, const char* const caTotClusTable);

int ClusUtil_AnisoClusters_OfSystem(int* const naFullClusList, int* const naCumClusSizes);

int ClusUtil_OvlpClusters_OfSystem(int* const naFullClusList, int* const naCumClusSizes);

void ClusUtil_GenClusSizesFromCumulativeSizes(int* const naSizeList, const int* const naCumClusSizes,
                                              int const nClusNum);

void ClusHistUtil_AddToHist_FromCountsList(lLong* const laHist, const int* const naCounts, const int nBins);

int ClusUtil_GetCluster_FromFullClusAndCumSizes(const int nClusID, int* const naOutClusList,
                                                const int* const naFullClusList, const int* const naCumSizes,
                                                const int* const naClusSizes, const int nClusNum);

void ClusHistUtil_AddToHist_MolTypeWiseDecomp_FromCluster(lLDub* const restrict ldaMolWiseHist,
                                                          const int* const restrict naClusList, const int nClusSize);

void ClusHistUtil_AddToHist_MolTypeWiseDecomp_FromFullClusAndCumSizes(lLDub* const ldaMolWiseHist,
                                                                      const int* const naClusChainIDs,
                                                                      const int* const naCumClusSizes,
                                                                      const int* const naClusSizes, int const nClusNum);

int ClusAnalysis_Ovlp_ForSystem_MolTypeWiseDecompAndSizes(lLong* const restrict naSizeHist,
                                                          lLDub* const restrict ldaMolWiseHist);

int ClusAnalysis_Aniso_ForSystem_MolTypeWiseDecompAndSizes(lLong* const restrict naSizeHist,
                                                           lLDub* const restrict ldaMolWiseHist);

int ClusUtil_GetLargestCluster_FromFullClusAndCumSizes(int* const naOutClusList, const int* const naFullClusList,
                                                       const int* const naCumSizes, const int* const naClusSizes,
                                                       const int nClusNum);

int ClusUtil_GetSecondLargestCluster_FromFullClusAndCumSizes(int* const naOutClusList, const int* const naFullClusList,
                                                             const int* const naCumSizes, const int* const naClusSizes,
                                                             const int nClusNum);

int ClusUtil_OvlpCluster_OfSystem_SecondLargest(int* naOutClusList);

int ClusUtil_OvlpCluster_OfSystem_SecondLargest_ForMCMove(int* naOutClusList);

int ClusUtil_AnisoCluster_OfSystem_SecondLargest(int* naOutClusList);

int ClusUtil_AnisoCluster_OfSystem_SecondLargest_ForMCMove(int* naOutClusList);

int ClusUtil_OfSystem_SecondLargest(int* naOutClusList, const int nMode);

void ClusUtil_MolWise_FindLargestClusters(int* const restrict naClusIDs_out, const int* const restrict naChainTypes_in,
                                          const int* const restrict naClusSizes_in,
                                          const int* const restrict naCumClusSizes_in, const int nClusNum_in);

int ClusUtil_OfSystem_MolWise_GetLargestClusters(int* const naClusIDsList_out, int* const naFullClusList_out,
                                                 int* const naClusSizes_out, int* const naCumSizes_out);

void ClusUtil_SelectSubsetOfClusters(int* const restrict naSubsetClusChains_out,
                                     int* const restrict naSubsetClusSizes_out,
                                     int* const restrict naSubsetCumSizes_out, const int nSubsetNum_in,
                                     const int* const restrict naSubSetClusIDs_in,
                                     const int* const restrict naFullClusList_in,
                                     const int* const restrict naClusSizes_in,
                                     const int* const restrict naCumClusSizes_in);

#endif // _CLUSTER_H_
