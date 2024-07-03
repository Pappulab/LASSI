#ifndef _STRUCTURE_H_ // include guard
#define _STRUCTURE_H_

#include "global.h"

int Lat_Ind_FromCoords(const int i, const int j, const int k);

int Lat_Ind_FromVec(const int* xArr);

int Lat_Ind_OfBead(const int beadID);

float Dist_BeadToBead(const int bead1, const int bead2);

float Dist_BeadToPoint(const int bead1, const int* f1);

float Dist_PointToPoint_Float(const float* f1, const float* f2);

float Dist_PointToPoint(const int* f1, const int* f2);

float Dist_Vec3n(const int* f1);

int Dist_VecMagSq(const int* f1);

int Check_LinkerConstraint(const int beadID, const int* tmpR);

int Check_MTLinkerConstraint_OLD(int beadID, int (*tmpR)[POS_MAX]);

int Check_LinkerConstraints_ForBeadList(const int listSize, const int* beadList);

int Check_BeadID_InList(const int thisBeadID, const int listSize, const int beadList[MAX_BONDS + 1]);

void PerformRuntimeSanityChecks(long nGen, int run_cycle);

void RDF_ComponentWise_Avg(void);

int RDF_ComponentIndex(const int i, const int j);

int RDFArr_Index(const int run_cycle, const int rdf_comp, const int x_pos);

int RadDen_ComponentIndex(const int i, const int j);

int RadDenArr_Index(const int run_cycle, const int rad_comp, const int x_pos);

int MolClusArr_Index(const int run_cycle, const int chain_type, const int clus_size);

void Calc_SystemCenterOfMass(lDub* tmpR);

void Calc_CenterOfMass_OfCluster(lDub* naCOM_out, const int cluster_size, const int* const restrict naClusList_in);

void Calc_CenterOfMass_OfSystem_OfMolType(lDub* tmpR, const int thisType);

void Calc_CenterOfMass_OfSystem_WithoutMolType(lDub* tmpR, const int thisType);

void RadDen_Avg_MolTypeWise_FromSysCen(void);

void RadDenHistUtil_ForSystem_FromLargestClusterOfMolTypes(void);

void RadDenHistUtil_ForSystem_FromCenterOfMassOfMolTypes(void);

void RadialDensityAnalysis_Perform_Analysis(void);

void GyrTensor_GyrRad_Avg(void);

void Calculate_Distances_For_Radius(float* thisList, const int nRad);

int NeighborSearch_ForOvlp(int const beadID, const int* startVec, int* neighList);

int NeighborSearch_ForCont(int const beadID, const int* startVec, int* contList, int* ovlpList, int* ovlpNum);

int NeighborSearch_SolSitesAround(const int* startVec);

int NeighborSearch_AroundPoint_wRad_IgnBead(const int beadID, const int* startVec, const int nRad, int* neighList);

int NeighborSearch_AroundPoint_wRad_wDists(const int beadID, const int* startVec, const int nRad, int* neighList,
                                           float* distList);

int NeighborSearch_ForCluster_Ovlp(int const beadID, const int* startVec, int* neighList);

void LatPos_copy(int* copy_vec, const int* in_vec);

void LatPos_add_wPBC(int* outVec, const int* firVec, const int* secVec);

void BeadPos_sub_wPBC(int* outVec, const int* firVec, const int* secVec);

void LatPos_add_wPBC_ofComp(int* outVec, const int* firVec, const int* secVec, const int compNum);

void LatPos_add_noPBC(int* outVec, const int* firVec, const int* secVec);

void LatPos_gen_rand_wRad(int outVec[POS_MAX], const int nRadius);

int OP_GetTopoBonds(const int beadID, int* dum_list);

int Vec3n_DotProd(const int* vec_1, const int* vec_2);

float Vec3n_CosTheta(const int* v1, const int* v2);

float Vec3n_AngleBetweenVecs(const int* vec_1, const int* vec_2);

float BeadPos_CosThetaOfBeads(const int bead1, const int bead2);

void BeadListOP_GetChainIDs(const int beadNum, const int* beadList, int* chainList);

void BeadListOP_GetChainTypes(const int beadNum, const int* beadList, int* chainList);

int BeadListOP_Filter_wrt_SecondList(const int beadNum, int* beadList, const int* propList, const int prop_val);

int BeadListOP_InvFilter_wrt_SecondList(const int beadNum, int* beadList, const int* propList, const int prop_val);

void BeadListOP_GetIntraChainID(const int beadNum, const int* beadList, int* chainList);

int BeadListOP_Filter_DbPvtLinkerConFwd(const int beadNum, int* beadList, const int thisBead);

int BeadListOP_Filter_DbPvtLinkerConBck(const int beadNum, int* beadList, const int thisBead);

int ListOP_UniqueElementsOfSortedList_Int(const int size, int* sorted_list);

int BeadList_AppendBeads(const int old_size, int* old_list, const int* app_list, const int app_size);

int UtilFunc_CompareInts(const void* a, const void* b);

int BeadList_CanTopoAngle(const int size, int* beadList);

int ChainListOP_AddUniqueChains_wSize(const int clusSize, int* clusList, const int newChainNums, const int* chainList,
                                      const int maxClus);

int ClusListOP_AddChainID(const int clusSize, int* clusList, const int chainID);

int ClusListOP_AddIfUniqueChainID(const int clusSize, int* clusList, const int chainID);

int ChainListOP_AddUniqueChains_wSize_Check(const int clusSize, int* clusList, const int newChainNums,
                                            const int* chainList, const int maxClus, const int* oldList);

void ChainListOP_GetChainTypes(int* const naTypesList, const int* const naChainIDs, const int nListSize);

void LatticeUtil_GetOvlpLattIndecies(const int* r_pos0, int* numList);

void LatticeUtil_GetLattVals_FromList(const int* nIndexList, int* nLatValsList, const int listSize);

int LatticeUtil_GetBeadIDs_FromList(const int* nLatValsList, int* nBeadsList, const int listSize);

int LatticeUtil_GetNeighBeads_AtPos(const int* restrict r_pos0, int* restrict nBeadsList);

int ListOP_GetMaxVal_Int(const int nListSize, const int* const naList);

int ListOP_IndicesForVal_Int(const int nVal, int* const naOutList, const int nListSize, const int* const naInputList);

int ListOP_GetRandomIndexForVal_Int(const int nVal, const int nListSize, const int* const naList);

int ListOP_Get2ndLargestVal_Int(const int nListSize, const int* const naInputList);

void ListOP_GenHistFromCounts_Int(int* const restrict naOutFreqHist, const int* const restrict naInputList,
                                  const int nSize);

#endif // _STRUCTURE_H_
