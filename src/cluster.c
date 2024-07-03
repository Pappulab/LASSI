#include "cluster.h"
#include "global.h"
#include "structure.h"

/// Clus_ChainNetwork_General - calculates the cluster chainID is a part of.
/// In summary generate a tree diagram of the bonding chains, and then keep
/// going down the branches to generate more sub-branches iteratively. Gets the
/// exhaustive list of the total network chainID is part of.
/// \param chainID
/// \return ClusSize - the size of this cluster+1 (the +1 is for looping)
/// This version isn't used yet, but was before when the Cluster move moved a
/// cluster, rather than the two new moves.
int Clus_Network_ChainCluster_General(int const chainID)
{
    // Updates naList_glb to have all proteins bound to  chainID and it's cluster
    // The idea is to check every bead and see if there is a unique bonded
    // chain, and to add it to naList_glb

    int i, j; // Loop iterators
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naList_glb[i]           = -1; // Initialize the list where -1 means empty.
            naChainCheckList_glb[i] = -1; // Initialize the list where no chain has been visited.
        }

    int list_it  = 0; // Iterator for naList_glb
    int clusSize = 0; // Index to track the cluster size.
    int curID;        // Index to track the current chain being looked at.
    curID                       = chainID;
    naList_glb[clusSize++]      = curID; // The cluster contains chainID by definition, and ClusSize = 1
    naChainCheckList_glb[curID] = 1;     // This means that curID has been checked!
    int fB, lB;                          // Indecies to track the first and last bead of chains.
    int chainPart;

    while (curID != -1)
        {                                                  // Keep going through naList_glb till it is exhausted.
            fB = chain_info_glb[curID][CHAIN_START];       // First bead of this chain.
            lB = fB + chain_info_glb[curID][CHAIN_LENGTH]; // Last bead+1 of this chain. Makes
                                                           // for-loops easier this way.
            for (i = fB; i < lB; i++)
                { // Loop over all the beads in this chain
                  // and see if there is a physical bond.
                    if (bead_info_glb[i][BEAD_FACE] != -1)
                        { // This means we have a bonding partner.
                            chainPart = bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_CHAINID];
                            if (naChainCheckList_glb[chainPart] == -1)
                                { // This is a unique chain.
                                    naList_glb[clusSize++]          = chainPart;
                                    naChainCheckList_glb[chainPart] = 1; // Recording that this chain has been looked at
                                }
                        }
                    // Moving on to the next bead in this chain
                }
            // Done with this chain, so let's move to the next chain, if it exists.
            list_it++; // Going one forward in naList_glb to completely exhaust the
                       // tree.
            curID = naList_glb[list_it];
        }
    return clusSize;
}

/// Clus_ChainNetwork_ForTotal - cluster colculation, specifically for
/// Clus_SecondLargestCluster and Clus_TotalAnalysis. Note that this variant of
/// the clustering does not reset naList_glb and is meant to only be used with total
/// clustering analyses. For systems with a lot of molecules, just the
/// initialization of naList_glb can take too long, so it is only done once before
/// this function is used repeatedly.
/// \param chainID
/// \return
int Clus_Network_ChainCluster_ForTotal(int const chainID)
{
    // Updates naList_glb to have all proteins bound to  chainID and it's cluster
    // The idea is to check every bead and see if there is a unique bonded
    // chain, and to add it to naList_glb This one is specifically made to be used
    // in total system network analyses, so naList_glb and naChainCheckList_glb aren't
    // reinitialized to -1. Instead, we just use naChainCheckList_glb fully.
    int i;            // Loop iterators
    int list_it  = 0; // Iterator for naList_glb
    int clusSize = 0; // Index to track the cluster size.
    int curID;        // Index to track the current chain being looked at.
    curID                       = chainID;
    naList_glb[clusSize++]      = curID; // The cluster contains chainID by definition, and ClusSize = 1
    naChainCheckList_glb[curID] = 1;     // This means that curID has been checked!
    int fB, lB;                          // Indecies to track the first and last bead of chains.
    int chainPart;

    while (curID != -1)
        {                                                  // Keep going through naList_glb till it is exhausted.
            fB = chain_info_glb[curID][CHAIN_START];       // First bead of this chain.
            lB = fB + chain_info_glb[curID][CHAIN_LENGTH]; // Last bead+1 of this chain. Makes
                                                           // for-loops easier this way.
            for (i = fB; i < lB; i++)
                { // Loop over all the beads in this chain
                  // and see if there is a physical bond.
                    if (bead_info_glb[i][BEAD_FACE] != -1)
                        { // This means we have a bonding partner.
                            chainPart = bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_CHAINID];
                            if (naChainCheckList_glb[chainPart] == -1)
                                { // This is a unique chain.
                                    naList_glb[clusSize++]          = chainPart;
                                    naChainCheckList_glb[chainPart] = 1; // Recording that this chain has been looked at
                                }
                        }
                    // Moving on to the next bead in this chain
                }
            // Done with this chain, so let's move to the next chain in the cluster,
            // if it exists.
            list_it++; // Going one forward in naList_glb to completely exhaust the
                       // tree.
            curID = naList_glb[list_it];
        }
    return clusSize;
}

/// Clus_TotalAnalysis - calculates the total networking/clustering of the
/// system and stores cluster information in naClusterMatrix_glb[][]
void Clus_Network_TotalAnalysis(void)
{
    int curID, Cluster_length, currentLargest, i, j;
    int ClusNum  = 0;
    int IsUnique = 1;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naChainCheckList_glb[i] = -1;
            naList_glb[i]           = -1;
        }
    curID          = 0; // Start with the 0th chain
    currentLargest = 0;
    j              = 0;
    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Network_ChainCluster_ForTotal(curID); // This is the length of curID cluster
            if (Cluster_length > currentLargest)
                {
                    currentLargest = Cluster_length;
                    j              = ClusNum;
                }

            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    naClusterMatrix_glb[ClusNum][i + 1] = naList_glb[i];
                    naList_glb[i]                       = -1;
                }
            naClusterMatrix_glb[ClusNum++][0] = Cluster_length;
            IsUnique                          = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }
}

/// Clus_SecondLargestCluster - does what Clus_TotalAnalysis does, and then
/// finds the second largest cluster. Note that naList_glb[] now has the chainIDs of
/// the second largest cluster! Furthermore, in the case where we have only one
/// cluster, the function returns -1, which causes SmallClusMCMove to fail, and
/// if we have multiple smallest clusters, randomly pick one.
/// \return naClusterMatrix_glb[clusID][0] - the size of the second largest cluster
int Clus_Network_SecondLargestCluster(void)
{
    /*
    Calculates the complete networking for the system and returns naList_glb which
    contains the chainIDs for the second largest cluster, if it exists.
    */
    // printf("\nStarting clus analysis\n");
    int curID, Cluster_length, currentLargest, i, j;
    int ClusNum  = 0;
    int IsUnique = 1;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naChainCheckList_glb[i] = -1;
            naList_glb[i]           = -1;
        }
    curID          = 0; // Start with the 0th chain
    currentLargest = 0;
    j              = 0;
    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Network_ChainCluster_ForTotal(curID); // This is the length of curID cluster
            if (Cluster_length > currentLargest)
                {
                    currentLargest = Cluster_length;
                    j              = ClusNum;
                }

            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    naClusterMatrix_glb[ClusNum][i + 1] = naList_glb[i];
                    naList_glb[i]                       = -1;
                }
            naClusterMatrix_glb[ClusNum++][0] = Cluster_length;
            IsUnique                          = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }
    // printf("ClusNum is %d\n", ClusNum);
    if (currentLargest == 1)
        { // Just single chains
            // printf("Only single chains\t%d\n", ClusNum);
            curID = rand() % ClusNum;
        }
    else
        { // Find second largest
            curID = 0;
            for (i = 0; i < ClusNum; i++)
                {
                    if (naClusterMatrix_glb[i][0] < currentLargest &&
                        naClusterMatrix_glb[i][0] >= naClusterMatrix_glb[curID][0])
                        {
                            curID = i;
                        }
                }
            // printf("%d:\t\t%d\t%d\n", ClusNum, currentLargest,
            // naClusterMatrix_glb[curID][0]);
            if (curID == j)
                { // Reject if only one cluster
                    // printf("\t\tTOO BIG\n");
                    return -1;
                }
        }
    for (i = 0; i < naClusterMatrix_glb[curID][0]; i++)
        { // Reupdate naList_glb to have IDs for the second largest cluster
            naList_glb[i] = naClusterMatrix_glb[curID][i + 1];
        }
    return naClusterMatrix_glb[curID][0];
}

/// Clus_LimitedCluster - performs clustering analysis on chainID and
/// immediately ends if the cluster is larger than 4. In this version of the
/// clustering, since the cluster size is going to be at most 5, it is faster to
/// recursively check through naList_glb[] to see if the newly proposed molecule
/// should be added or not, rather than use naChainCheckList_glb which acts as sort
/// of a hash table! \param chainID \return clusSize - the size of the cluster;
/// also naList_glb[] now has the chainIDs. -1 if clusSize >= 5
int Clus_Network_LimitedCluster(int const chainID)
{
    // Updates naList_glb to have all proteins bound to  chainID and it's cluster
    // The idea is to check every bead and see if there is a unique bonded
    // chain, and to add it to naList_glb If Cluster becomes larger than 5, exit and
    // return -1
    const int ClusterLimit = nLimitedClusterSize_glb;
    int i, j; // Loop iterators
    for (i = 0; i < ClusterLimit; i++)
        {
            naList_glb[i] = -1;
        }

    int list_it  = 0; // Iterator for naList_glb
    int clusSize = 0; // Index to track the cluster size.
    int curID;        // Index to track the current chain being looked at.
    curID                  = chainID;
    naList_glb[clusSize++] = curID; // The cluster contains chainID by definition, and ClusSize = 1
    int fB, lB;                     // Indecies to track the first and last bead of chains.
    int chainPart;
    int IsUnique = 1; // Tracks if a chain is unique or not. 0 is non-unique,
                      // and 1 is unique.

    while (curID != -1)
        {                                                  // Keep going through naList_glb till it is exhausted.
            fB = chain_info_glb[curID][CHAIN_START];       // First bead of this chain.
            lB = fB + chain_info_glb[curID][CHAIN_LENGTH]; // Last bead+1 of this chain. Makes
                                                           // for-loops easier this way.

            for (i = fB; i < lB; i++)
                { // Loop over all the beads in this chain
                  // and see if there is a physical bond.
                    if (bead_info_glb[i][BEAD_FACE] != -1)
                        { // This means we have a bonding partner.
                            chainPart = bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_CHAINID];
                            // Checking if this chain is unique
                            IsUnique = 1;
                            for (j = 0; j < clusSize; j++)
                                {
                                    if (chainPart == naList_glb[j])
                                        {
                                            IsUnique = 0;
                                            break;
                                        }
                                }
                            if (IsUnique == 1)
                                {
                                    naList_glb[clusSize++] = chainPart;
                                }
                            if (clusSize >= ClusterLimit)
                                {
                                    return -1;
                                }
                        }
                    // Moving on to the next bead in this chain
                }
            // Done with this chain, so let's move to the next chain, if it exists.
            list_it++; // Going one forward in naList_glb to completely exhaust the
                       // tree.
            curID = naList_glb[list_it];
        }
    // printf("%d\n", list_it);
    return clusSize;
}

/// Clus_Distribution_Avg - do what Clus_TotalAnalysis does but instead of
/// remembering clusters in naClusterMatrix_glb, update laClusHistList_glb[] and keep making
/// the total cluster histogram for the system.
void Clus_Network_Distribution_Avg(void)
{
    /*
    Calculates the cluster distribution using total_network_analysis framework,
    but keeps adding to laClusHistList_glb for total averaging at the end. Read
    total_network_analysis() for what's happening here LOL
    */
    int curID, Cluster_length, currentLargest, i;
    int IsUnique = 1;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naList_glb[i]           = -1;
            naChainCheckList_glb[i] = -1;
        }
    curID          = 0; // Start with the 0th chain
    currentLargest = 0;
    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Network_ChainCluster_ForTotal(curID); // This is the length of curID cluster
            // printf("Clus Len: %d\n", Cluster_length);
            laClusHistList_glb[Cluster_length]++; // Adding to that cluster-size bin
            if (Cluster_length > currentLargest)
                {
                    currentLargest = Cluster_length;
                }
            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    naList_glb[i] = -1;
                }
            IsUnique = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }
    nTotClusCounter_glb++;
    nLargestClusterRightNow_glb += currentLargest;
}

/// Clus_DistributionMolWise_Avg - do what Clus_TotalAnalysis does but instead
/// of remembering clusters in naClusterMatrix_glb, update laClusHistList_glb[] and keep
/// making the total cluster histogram for the system.
void Clus_Network_Distribution_MolWise_Avg(void)
{
    /*
    Calculates the cluster distribution using total_network_analysis framework,
    but keeps adding to laClusHistList_glb for total averaging at the end. Read
    total_network_analysis() for what's happening here LOL
    */
    int curID, Cluster_length, currentLargest, i;
    int IsUnique = 1;
    int thisType = 0;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naList_glb[i]           = -1;
            naChainCheckList_glb[i] = -1;
        }
    curID          = 0; // Start with the 0th chain
    currentLargest = 0;
    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Network_ChainCluster_ForTotal(curID); // This is the length of curID cluster
            // printf("Clus Len: %d\n", Cluster_length);
            laClusHistList_glb[Cluster_length]++; // Adding to that cluster-size bin
            if (Cluster_length > currentLargest)
                {
                    currentLargest = Cluster_length;
                }
            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    thisType = chain_info_glb[naList_glb[i]][CHAIN_TYPE];
                    ldaMOLCLUS_Arr_glb[MolClusArr_Index(0, thisType, Cluster_length - 1)]++;
                    naList_glb[i] = -1;
                    // printf("%d\n", thisType);
                }
            IsUnique = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }
    nLargestClusterRightNow_glb += currentLargest;
}

/*
 * TODO: I need to update all the clustering functions with the newer neighbor-search-like functions.
 */

/// Clus_MolWiseLargestCluster - performs a total clustering analysis of the
/// system. Then finds out the largest cluster for each MolType, where
/// redundancy is allowed. naList_glb contains the cluster IDs of the chain types.
/// Note that naList_glb[0] contains the system's overall largest cluster, naList_glb[1]
/// is for molType=0 and so on.
void Clus_Network_MolWise_LargestClusters(void)
{
    int curID, Cluster_length, i;
    int ClusNum  = 0;
    int IsUnique = 1;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naChainCheckList_glb[i] = -1;
            naList_glb[i]           = -1;
        }
    curID = 0; // Start with the 0th chain

    int molType = 0;
    int* clus_numof_MolType;  // Tracks the number breakdown in each cluster
    int* clusID_of_type;      // Tracks the clusterID of largest of each type
    int* largestClus_of_type; // Tracks the largest number seen for each type

    clus_numof_MolType  = malloc((tot_chain_types_glb + 1) * sizeof(lInt));
    clusID_of_type      = malloc((tot_chain_types_glb + 1) * sizeof(lInt));
    largestClus_of_type = malloc((tot_chain_types_glb + 1) * sizeof(lInt));

    for (i = 0; i <= tot_chain_types_glb; i++)
        {
            clus_numof_MolType[i]  = 0;
            clusID_of_type[i]      = 0;
            largestClus_of_type[i] = 0;
        }

    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Network_ChainCluster_ForTotal(curID); // This is the length of curID cluster

            for (i = 0; i < tot_chain_types_glb; i++)
                {
                    clus_numof_MolType[i] = 0;
                }

            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    naClusterMatrix_glb[ClusNum][i + 1] = naList_glb[i];
                    molType                             = chain_info_glb[naList_glb[i]][CHAIN_TYPE];
                    clus_numof_MolType[molType]++;
                    naList_glb[i] = -1;
                }

            for (i = 0; i < tot_chain_types_glb; i++)
                {
                    if (clus_numof_MolType[i] > largestClus_of_type[i])
                        {
                            largestClus_of_type[i] = clus_numof_MolType[i];
                            clusID_of_type[i]      = ClusNum;
                        }
                    clus_numof_MolType[i] = 0;
                }

            naClusterMatrix_glb[ClusNum++][0] = Cluster_length;
            IsUnique                          = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }

    for (i = 0; i < tot_chain_types_glb; i++)
        {
            // printf("(%d %d)\t", largestClus_of_type[i],clusID_of_type[i]);
            naList_glb[i] = clusID_of_type[i]; // Now naList_glb contains the cluster ID's.
        }
    // printf("\n");
    free(clus_numof_MolType);
    free(clusID_of_type);
    free(largestClus_of_type);
}

int Clus_Proximity_ChainCluster_ForTotal_IntOnly(int const chainID)
{
    // Updates naList_glb to have all proteins close to chainID
    // The idea is to check every bead and see if there is a unique chain in the
    // (+-1,+-1) cube, and to add it to naList_glb This one is specifically made to
    // be used in total system proximity cluster analyses, so naList_glb and
    // naChainCheckList_glb aren't reinitialized to -1. Instead, we just use
    // naChainCheckList_glb fully. This version defines a bond only when the
    // neighboring beads are interacting with each other favorably: E_{ij} < 0.
    int i, j, k;      // Loop iterators
    int list_it  = 0; // Iterator for naList_glb
    int clusSize = 0; // Index to track the cluster size.
    int curID;        // Index to track the current chain being looked at.
    curID                       = chainID;
    naList_glb[clusSize++]      = curID; // The cluster contains chainID by definition, and ClusSize = 1
    naChainCheckList_glb[curID] = 1;     // This means that curID has been checked!
    int fB, lB;                          // Indecies to track the first and last bead of chains.
    int chainPart;
    int tmpBead = 0;
    int resi, resj; // Tracking the type of the bead to check if they are
                    // interacting via E_OVLP
    int tmpR[POS_MAX]  = {0};
    int tmpR2[POS_MAX] = {0};
    while (curID != -1)
        {                                                  // Keep going through naList_glb till it is exhausted.
            fB = chain_info_glb[curID][CHAIN_START];       // First bead of this chain.
            lB = fB + chain_info_glb[curID][CHAIN_LENGTH]; // Last bead+1 of this chain. Makes
                                                           // for-loops easier this way.
            for (i = fB; i < lB; i++)
                { // Loop over all the beads in this chain and see if there
                  // is another bead around.
                    for (j = 0; j < POS_MAX; j++)
                        { // This is where I am right now
                            tmpR[j] = bead_info_glb[i][j];
                        }
                    resi = bead_info_glb[i][BEAD_TYPE];
                    if (bead_info_glb[i][BEAD_FACE] != -1)
                        { // This means we have a bonding partner.
                            chainPart = bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_CHAINID];
                            // Checking if this chain is unique
                            if (naChainCheckList_glb[chainPart] == -1)
                                { // This is a unique chain.
                                    naList_glb[clusSize++]          = chainPart;
                                    naChainCheckList_glb[chainPart] = 1; // Recording that this chain has been looked at
                                }
                        }
                    if (nBeadTypeCanOvlp_glb[resi] == 0)
                        {
                            continue;
                        }
                    for (k = 0; k < MAX_ROTSTATES - 1; k++)
                        {
                            tmpBead = naRot_IndArr_glb[k];
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    tmpR2[j] =
                                        (tmpR[j] + naLocalArr_glb[tmpBead][j] + naBoxSize_glb[j]) % naBoxSize_glb[j];
                                }
                            tmpBead = Lat_Ind_FromVec(tmpR2);
                            tmpBead = naTotLattice_glb[tmpBead];
                            if (tmpBead != -1)
                                {
                                    resj = bead_info_glb[tmpBead][BEAD_TYPE];
                                    if (faEnergy_glb[resi][resj][E_OVLP] >= 0.)
                                        { // If not interacting, or repelling, it is not a
                                          // bond
                                            continue;
                                        }
                                    chainPart = bead_info_glb[tmpBead][BEAD_CHAINID];
                                    if (naChainCheckList_glb[chainPart] == -1)
                                        { // This is a unique chain.
                                            naList_glb[clusSize++] = chainPart;
                                            naChainCheckList_glb[chainPart] =
                                                1; // Recording that this chain has been looked at
                                        }
                                }
                        }
                    // Moving on to the next bead in this chain
                }
            // Done with this chain, so let's move to the next chain in the cluster,
            // if it exists.
            list_it++; // Going one forward in naList_glb to completely exhaust the
                       // tree.
            curID = naList_glb[list_it];
        }
    return clusSize;
}

int Clus_Proximity_ChainCluster_ForTotal_All(int const chainID)
{
    // Updates naList_glb to have all proteins close to chainID
    // The idea is to check every bead and see if there is a unique chain in the
    // (+-1,+-1) cube, and to add it to naList_glb This one is specifically made to
    // be used in total system proximity cluster analyses, so naList_glb and
    // naChainCheckList_glb aren't reinitialized to -1. Instead, we just use
    // naChainCheckList_glb fully. As opposed to the IntOnly version, this version
    // defines a bond as being within the (+-1,+-1) only, even if there is no
    // interaction energy.
    int i, j, k;      // Loop iterators
    int list_it  = 0; // Iterator for naList_glb
    int clusSize = 0; // Index to track the cluster size.
    int curID;        // Index to track the current chain being looked at.
    curID                       = chainID;
    naList_glb[clusSize++]      = curID; // The cluster contains chainID by definition, and ClusSize = 1
    naChainCheckList_glb[curID] = 1;     // This means that curID has been checked!
    int fB, lB;                          // Indecies to track the first and last bead of chains.
    int chainPart;
    int tmpBead        = 0;
    int tmpR[POS_MAX]  = {0};
    int tmpR2[POS_MAX] = {0};
    while (curID != -1)
        {                                                  // Keep going through naList_glb till it is exhausted.
            fB = chain_info_glb[curID][CHAIN_START];       // First bead of this chain.
            lB = fB + chain_info_glb[curID][CHAIN_LENGTH]; // Last bead+1 of this chain. Makes
                                                           // for-loops easier this way.
            for (i = fB; i < lB; i++)
                { // Loop over all the beads in this chain and see if there
                  // is another bead around.
                    for (j = 0; j < POS_MAX; j++)
                        { // This is where I am right now
                            tmpR[j] = bead_info_glb[i][j];
                        }
                    for (k = 0; k < MAX_ROTSTATES - 1; k++)
                        {
                            tmpBead = naRot_IndArr_glb[k];
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    tmpR2[j] =
                                        (tmpR[j] + naLocalArr_glb[tmpBead][j] + naBoxSize_glb[j]) % naBoxSize_glb[j];
                                }
                            tmpBead = Lat_Ind_FromVec(tmpR2);
                            tmpBead = naTotLattice_glb[tmpBead];
                            if (tmpBead != -1)
                                {
                                    chainPart = bead_info_glb[tmpBead][BEAD_CHAINID];
                                    if (naChainCheckList_glb[chainPart] == -1)
                                        { // This is a unique chain.
                                            naList_glb[clusSize++] = chainPart;
                                            naChainCheckList_glb[chainPart] =
                                                1; // Recording that this chain has been looked at
                                        }
                                }
                        }
                    // Moving on to the next bead in this chain
                }
            // Done with this chain, so let's move to the next chain in the cluster,
            // if it exists.
            list_it++; // Going one forward in naList_glb to completely exhaust the
                       // tree.
            curID = naList_glb[list_it];
        }
    return clusSize;
}

void Clus_Proximity_TotalAnalysis(void)
{
    int curID, Cluster_length, currentLargest, i, j;
    int ClusNum  = 0;
    int IsUnique = 1;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naChainCheckList_glb[i] = -1;
            naList_glb[i]           = -1;
        }
    curID          = 0; // Start with the 0th chain
    currentLargest = 0;
    j              = 0;
    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Proximity_ChainCluster_ForTotal_IntOnly(curID); // This is the length of curID cluster
            if (Cluster_length > currentLargest)
                {
                    currentLargest = Cluster_length;
                    j              = ClusNum;
                }

            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    naClusterMatrix_glb[ClusNum][i + 1] = naList_glb[i];
                    naList_glb[i]                       = -1;
                }
            naClusterMatrix_glb[ClusNum++][0] = Cluster_length;
            IsUnique                          = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }
}

int Clus_Proximity_SecondLargestCluster(void)
{
    /*
    Calculates the complete networking for the system and returns naList_glb which
    contains the chainIDs for the second largest cluster, if it exists.
    */
    // printf("\nStarting clus analysis\n");
    int curID, Cluster_length, currentLargest, i, j;
    int ClusNum  = 0;
    int IsUnique = 1;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naChainCheckList_glb[i] = -1;
            naList_glb[i]           = -1;
        }
    curID          = 0; // Start with the 0th chain
    currentLargest = 0;
    j              = 0;
    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Proximity_ChainCluster_ForTotal_IntOnly(curID); // This is the length of curID cluster
            if (Cluster_length > currentLargest)
                {
                    currentLargest = Cluster_length;
                    j              = ClusNum;
                }

            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    naClusterMatrix_glb[ClusNum][i + 1] = naList_glb[i];
                    naList_glb[i]                       = -1;
                }
            naClusterMatrix_glb[ClusNum++][0] = Cluster_length;
            IsUnique                          = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }
    // printf("ClusNum is %d\n", ClusNum);
    if (currentLargest == 1)
        { // Just single chains
            // printf("Only single chains\t%d\n", ClusNum);
            curID = rand() % ClusNum;
        }
    else
        { // Find second largest
            curID = 0;
            for (i = 0; i < ClusNum; i++)
                {
                    if (naClusterMatrix_glb[i][0] < currentLargest &&
                        naClusterMatrix_glb[i][0] >= naClusterMatrix_glb[curID][0])
                        {
                            curID = i;
                        }
                }
            // printf("%d:\t\t%d\t%d\n", ClusNum, currentLargest,
            // naClusterMatrix_glb[curID][0]);
            if (curID == j)
                { // Reject if only one cluster
                    // printf("\t\tTOO BIG\n");
                    return -1;
                }
        }
    for (i = 0; i < naClusterMatrix_glb[curID][0]; i++)
        { // Reupdate naList_glb to have IDs for the second largest cluster
            naList_glb[i] = naClusterMatrix_glb[curID][i + 1];
        }
    return naClusterMatrix_glb[curID][0];
}

int Clus_Proximity_LimitedCluster_IntOnly(int const chainID)
{
    // Updates naList_glb to have all proteins bound to  chainID and it's cluster
    // The idea is to check every bead and see if there is a unique bonded
    // chain, and to add it to naList_glb If Cluster becomes larger than 15, exit
    // and return -1
    int ClusterLimit = 15;
    int i, j, k; // Loop iterators
    for (i = 0; i < ClusterLimit; i++)
        {
            naList_glb[i] = -1;
        }

    int list_it  = 0; // Iterator for naList_glb
    int clusSize = 0; // Index to track the cluster size.
    int curID;        // Index to track the current chain being looked at.
    curID                  = chainID;
    naList_glb[clusSize++] = curID; // The cluster contains chainID by definition, and ClusSize = 1
    int fB, lB;                     // Indecies to track the first and last bead of chains.
    int chainPart;
    int tmpBead = 0;
    int resi, resj; // Tracking the type of the bead to check if they are
                    // interacting via E_OVLP
    int tmpR[POS_MAX]  = {0};
    int tmpR2[POS_MAX] = {0};
    int IsUnique       = 1; // Tracks if a chain is unique or not. 0 is non-unique,
                            // and 1 is unique.

    while (curID != -1)
        {                                                  // Keep going through naList_glb till it is exhausted.
            fB = chain_info_glb[curID][CHAIN_START];       // First bead of this chain.
            lB = fB + chain_info_glb[curID][CHAIN_LENGTH]; // Last bead+1 of this chain. Makes
                                                           // for-loops easier this way.

            for (i = fB; i < lB; i++)
                { // Loop over all the beads in this chain and see if there
                  // is another bead around.
                    for (j = 0; j < POS_MAX; j++)
                        { // This is where I am right now
                            tmpR[j] = bead_info_glb[i][j];
                        }
                    if (bead_info_glb[i][BEAD_FACE] != -1)
                        { // This means we have a bonding partner.
                            chainPart = bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_CHAINID];
                            // Checking if this chain is unique
                            IsUnique = 1;
                            for (j = 0; j < clusSize; j++)
                                {
                                    if (chainPart == naList_glb[j])
                                        {
                                            IsUnique = 0;
                                            break;
                                        }
                                }
                            if (IsUnique == 1)
                                {
                                    naList_glb[clusSize++] = chainPart;
                                }
                            if (clusSize >= ClusterLimit)
                                {
                                    return -1;
                                }
                        }
                    resi = bead_info_glb[i][BEAD_TYPE];
                    if (nBeadTypeCanOvlp_glb[resi] == 0)
                        {
                            continue;
                        }
                    for (k = 0; k < MAX_ROTSTATES - 1; k++)
                        {
                            tmpBead = naRot_IndArr_glb[k];
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    tmpR2[j] =
                                        (tmpR[j] + naLocalArr_glb[tmpBead][j] + naBoxSize_glb[j]) % naBoxSize_glb[j];
                                }
                            tmpBead = Lat_Ind_FromVec(tmpR2);
                            tmpBead = naTotLattice_glb[tmpBead];
                            if (tmpBead != -1)
                                {
                                    resj = bead_info_glb[tmpBead][BEAD_TYPE];
                                    if (faEnergy_glb[resi][resj][E_OVLP] >= 0.)
                                        { // If not interacting, or repelling, it is not a
                                          // bond
                                            continue;
                                        }
                                    chainPart = bead_info_glb[tmpBead][BEAD_CHAINID];
                                    // Checking if this chain is unique
                                    IsUnique = 1;
                                    for (j = 0; j < clusSize; j++)
                                        {
                                            if (chainPart == naList_glb[j])
                                                {
                                                    IsUnique = 0;
                                                    break;
                                                }
                                        }
                                    if (IsUnique == 1)
                                        {
                                            naList_glb[clusSize++] = chainPart;
                                        }
                                    if (clusSize >= ClusterLimit)
                                        {
                                            return -1;
                                        }
                                }
                        }
                    // Moving on to the next bead in this chain
                }
            // Done with this chain, so let's move to the next chain, if it exists.
            list_it++; // Going one forward in naList_glb to completely exhaust the
                       // tree.
            curID = naList_glb[list_it];
        }
    // printf("%d\n", list_it);
    return clusSize;
}

int Clus_Proximity_LimitedCluster_IntOnly_Check(int const chainID, int const* OldList)
{
    // Updates naList_glb to have all proteins bound to  chainID and it's cluster
    // The idea is to check every bead and see if there is a unique bonded
    // chain, and to add it to naList_glb If Cluster becomes larger than 15, exit
    // and return -1 Furthermore, it checks while it goes to see if naList_glb[i] ==
    // OldList[i]. If not, return -1.
    int ClusterLimit = 15;
    int i, j, k; // Loop iterators
    for (i = 0; i < ClusterLimit; i++)
        {
            naList_glb[i] = -1;
        }

    int list_it  = 0; // Iterator for naList_glb
    int clusSize = 0; // Index to track the cluster size.
    int curID;        // Index to track the current chain being looked at.
    curID                  = chainID;
    naList_glb[clusSize++] = curID; // The cluster contains chainID by definition, and ClusSize = 1
    int fB, lB;                     // Indecies to track the first and last bead of chains.
    int chainPart;
    int resi, resj; // Tracking the type of the bead to check if they are
                    // interacting via E_OVLP
    int tmpBead        = 0;
    int tmpR[POS_MAX]  = {0};
    int tmpR2[POS_MAX] = {0};
    int IsUnique       = 1; // Tracks if a chain is unique or not. 0 is non-unique,
                            // and 1 is unique.

    while (curID != -1)
        {                                                  // Keep going through naList_glb till it is exhausted.
            fB = chain_info_glb[curID][CHAIN_START];       // First bead of this chain.
            lB = fB + chain_info_glb[curID][CHAIN_LENGTH]; // Last bead+1 of this chain. Makes
                                                           // for-loops easier this way.

            for (i = fB; i < lB; i++)
                { // Loop over all the beads in this chain and see if there
                  // is another bead around.
                    for (j = 0; j < POS_MAX; j++)
                        { // This is where I am right now
                            tmpR[j] = bead_info_glb[i][j];
                        }
                    if (bead_info_glb[i][BEAD_FACE] != -1)
                        { // This means we have a bonding partner.
                            chainPart = bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_CHAINID];
                            // Checking if this chain is unique
                            IsUnique = 1;
                            for (j = 0; j < clusSize; j++)
                                {
                                    if (chainPart == naList_glb[j])
                                        {
                                            IsUnique = 0;
                                            break;
                                        }
                                }
                            if (IsUnique == 1)
                                {
                                    naList_glb[clusSize++] = chainPart;
                                }
                            if (naList_glb[clusSize - 1] != OldList[clusSize - 1])
                                {
                                    return -1;
                                }
                            if (clusSize >= ClusterLimit)
                                {
                                    return -1;
                                }
                        }
                    resi = bead_info_glb[i][BEAD_TYPE];
                    if (nBeadTypeCanOvlp_glb[resi] == 0)
                        {
                            continue;
                        }
                    for (k = 0; k < MAX_ROTSTATES - 1; k++)
                        {
                            tmpBead = naRot_IndArr_glb[k];
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    tmpR2[j] =
                                        (tmpR[j] + naLocalArr_glb[tmpBead][j] + naBoxSize_glb[j]) % naBoxSize_glb[j];
                                }
                            tmpBead = Lat_Ind_FromVec(tmpR2);
                            tmpBead = naTotLattice_glb[tmpBead];
                            if (tmpBead != -1)
                                {
                                    resj = bead_info_glb[tmpBead][BEAD_TYPE];
                                    if (faEnergy_glb[resi][resj][E_OVLP] >= 0.)
                                        { // If not interacting, or repelling, it is not a
                                          // bond
                                            continue;
                                        }
                                    chainPart = bead_info_glb[tmpBead][BEAD_CHAINID];
                                    // Checking if this chain is unique
                                    IsUnique = 1;
                                    for (j = 0; j < clusSize; j++)
                                        {
                                            if (chainPart == naList_glb[j])
                                                {
                                                    IsUnique = 0;
                                                    break;
                                                }
                                        }
                                    if (IsUnique == 1)
                                        {
                                            naList_glb[clusSize++] = chainPart;
                                        }
                                    if (naList_glb[clusSize - 1] != OldList[clusSize - 1])
                                        {
                                            // printf("OOOOOOP %d\n", clusSize);
                                            return -1;
                                        }
                                    if (clusSize >= ClusterLimit)
                                        {
                                            return -1;
                                        }
                                }
                        }
                    // Moving on to the next bead in this chain
                }
            // Done with this chain, so let's move to the next chain, if it exists.
            list_it++; // Going one forward in naList_glb to completely exhaust the
                       // tree.
            curID = naList_glb[list_it];
        }
    // printf("%d\n", list_it);
    return clusSize;
}

///
/// \param chainID
/// \return
int Clus_Proximity_LimitedCluster_All(int const chainID, int* clusList)
{
    // Updates naList_glb to have all proteins bound to  chainID and it's cluster
    // The idea is to check every bead and see if there is a unique bonded
    // chain, and to add it to naList_glb If Cluster becomes larger than 15, exit
    // and return -1 If the number of chains is lower than 15, we use
    // tot_chains_glb/2.
    const int ClusterLimit = nLimitedClusterSize_glb;
    int i; // Loop iterators

    int list_it  = 0; // Iterator for naList_glb
    int clusSize = 0; // Index to track the cluster size.
    int fB, lB;       // Indecies to track the first and last bead of chains.

    int r_pos_curBead[POS_MAX] = {0};

    int tmpNeighBeadList[MAX_ROTSTATES + 1];
    int curBeadNeighNum;
    int curChainNeighNum;
    int tmpNeighChainList[MAX_ROTSTATES + 1];

    int curID            = chainID;
    clusList[clusSize++] = curID;

    while (list_it < clusSize)
        {
            curID = clusList[list_it];
            list_it++;

            fB = chain_info_glb[curID][CHAIN_START];       // First bead of this chain.
            lB = fB + chain_info_glb[curID][CHAIN_LENGTH]; // Last bead+1 of this chain. Makes
                                                           // for-loops easier this way.

            for (i = fB; i < lB; i++)
                { // Loop over all the beads in this chain. Check for neighbor chains
                    LatPos_copy(r_pos_curBead, bead_info_glb[i]);
                    curBeadNeighNum = NeighborSearch_ForCluster_Ovlp(i, r_pos_curBead, tmpNeighBeadList);
                    BeadListOP_GetChainIDs(curBeadNeighNum, tmpNeighBeadList, tmpNeighChainList);
                    qsort(tmpNeighChainList, curBeadNeighNum, sizeof(int), UtilFunc_CompareInts);
                    curChainNeighNum = ListOP_UniqueElementsOfSortedList_Int(curBeadNeighNum, tmpNeighChainList);
                    clusSize         = ChainListOP_AddUniqueChains_wSize(clusSize, clusList, curChainNeighNum,
                                                                 tmpNeighChainList, ClusterLimit);
                    if (clusSize >= ClusterLimit)
                        {
                            return -1;
                        }
                    // Moving on to the next bead in this chain
                }
            // Done with this chain, so let's move to the next chain, if it exists.
        }
    // printf("%d\n", list_it);
    clusList[clusSize] = -1;
    return clusSize;
}

int Clus_Proximity_LimitedCluster_All_Check(int const chainID, int const* OldList, int* NewList)
{
    // Updates naList_glb to have all proteins bound to  chainID and it's cluster
    // The idea is to check every bead and see if there is a unique bonded
    // chain, and to add it to naList_glb If Cluster becomes larger than 15, exit
    // and return -1 Furthermore, it checks while it goes to see if naList_glb[i] ==
    // OldList[i]. If not, return -1.
    const int ClusterLimit = nLimitedClusterSize_glb;
    int i; // Loop iterators

    int list_it  = 0; // Iterator for naList_glb
    int clusSize = 0; // Index to track the cluster size.
    int fB, lB;       // Indecies to track the first and last bead of chains.

    int r_pos_curBead[POS_MAX] = {0};

    int tmpNeighBeadList[MAX_ROTSTATES + 1];
    int curBeadNeighNum;
    int curChainNeighNum;
    int tmpNeighChainList[MAX_ROTSTATES + 1];

    int curID           = chainID;
    NewList[clusSize++] = curID;

    while (list_it < clusSize)
        {
            curID = NewList[list_it];
            list_it++;

            fB = chain_info_glb[curID][CHAIN_START];       // First bead of this chain.
            lB = fB + chain_info_glb[curID][CHAIN_LENGTH]; // Last bead+1 of this chain. Makes
                                                           // for-loops easier this way.

            for (i = fB; i < lB; i++)
                { // Loop over all the beads in this chain and see if there
                  // is another bead around.
                    LatPos_copy(r_pos_curBead, bead_info_glb[i]);
                    curBeadNeighNum = NeighborSearch_ForCluster_Ovlp(i, r_pos_curBead, tmpNeighBeadList);
                    BeadListOP_GetChainIDs(curBeadNeighNum, tmpNeighBeadList, tmpNeighChainList);
                    qsort(tmpNeighChainList, curBeadNeighNum, sizeof(int), UtilFunc_CompareInts);
                    curChainNeighNum = ListOP_UniqueElementsOfSortedList_Int(curBeadNeighNum, tmpNeighChainList);
                    clusSize         = ChainListOP_AddUniqueChains_wSize_Check(clusSize, NewList, curChainNeighNum,
                                                                       tmpNeighChainList, ClusterLimit, OldList);
                    if (clusSize >= ClusterLimit)
                        {
                            return -1;
                        }
                    // Moving on to the next bead in this chain
                }
            // Done with this chain, so let's move to the next chain, if it exists.
        }
    NewList[clusSize] = -1;
    // printf("%d\n", list_it);
    return clusSize;
}

void Clus_Proximity_Distribution_Avg(void)
{
    /*
    Calculates the cluster distribution using total_network_analysis framework,
    but keeps adding to laClusHistList_glb for total averaging at the end. Read
    total_network_analysis() for what's happening here LOL
    */
    int curID, Cluster_length, currentLargest, i;
    int IsUnique = 1;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naList_glb[i]           = -1;
            naChainCheckList_glb[i] = -1;
        }
    curID          = 0; // Start with the 0th chain
    currentLargest = 0;
    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Proximity_ChainCluster_ForTotal_IntOnly(curID); // This is the length of curID cluster
            // printf("Clus Len: %d\n", Cluster_length);
            laClusHistList_glb[Cluster_length]++; // Adding to that cluster-size bin
            if (Cluster_length > currentLargest)
                {
                    currentLargest = Cluster_length;
                }
            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    naList_glb[i] = -1;
                }
            IsUnique = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }
    nTotClusCounter_glb++;
    nLargestClusterRightNow_glb += currentLargest;
}

void Clus_Proximity_Distribution_IntOnly_MolWise_Avg(void)
{
    /*
    Calculates the cluster distribution using total_network_analysis framework,
    but keeps adding to laClusHistList_glb for total averaging at the end. Read
    total_network_analysis() for what's happening here LOL
    */
    int curID, Cluster_length, currentLargest, i;
    int IsUnique = 1;
    int thisType = 0;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naList_glb[i]           = -1;
            naChainCheckList_glb[i] = -1;
        }
    curID          = 0; // Start with the 0th chain
    currentLargest = 0;
    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Proximity_ChainCluster_ForTotal_IntOnly(curID); // This is the length of curID cluster
            // printf("Clus Len: %d\n", Cluster_length);
            laClusHistList_glb[Cluster_length]++; // Adding to that cluster-size bin
            if (Cluster_length > currentLargest)
                {
                    currentLargest = Cluster_length;
                }
            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    thisType = chain_info_glb[naList_glb[i]][CHAIN_TYPE];
                    ldaMOLCLUS_Arr_glb[MolClusArr_Index(0, thisType, Cluster_length - 1)]++;
                    naList_glb[i] = -1;
                    // printf("%d\n", thisType);
                }
            IsUnique = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }
    nLargestClusterRightNow_glb += currentLargest;
}

void Clus_Proximity_Distribution_All_MolWise_Avg(void)
{
    /*
    Calculates the cluster distribution using total_network_analysis framework,
    but keeps adding to laClusHistList_glb for total averaging at the end. Read
    total_network_analysis() for what's happening here LOL
    */
    int curID, Cluster_length, currentLargest, i;
    int IsUnique = 1;
    int thisType = 0;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naList_glb[i]           = -1;
            naChainCheckList_glb[i] = -1;
        }
    curID          = 0; // Start with the 0th chain
    currentLargest = 0;
    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Proximity_ChainCluster_ForTotal_All(curID); // This is the length of curID cluster
            // printf("Clus Len: %d\n", Cluster_length);
            laClusHistList_glb[Cluster_length]++; // Adding to that cluster-size bin
            if (Cluster_length > currentLargest)
                {
                    currentLargest = Cluster_length;
                }
            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    thisType = chain_info_glb[naList_glb[i]][CHAIN_TYPE];
                    ldaMOLCLUS_Arr_glb[MolClusArr_Index(0, thisType, Cluster_length - 1)]++;
                    naList_glb[i] = -1;
                    // printf("%d\n", thisType);
                }
            IsUnique = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }
    nLargestClusterRightNow_glb += currentLargest;
}

void Clus_Proximity_IntOnly_MolWise_LargestClusters(void)
{
    int curID, Cluster_length, i;
    int ClusNum  = 0;
    int IsUnique = 1;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naChainCheckList_glb[i] = -1;
            naList_glb[i]           = -1;
        }
    curID = 0; // Start with the 0th chain

    int molType = 0;
    int* clus_numof_MolType;  // Tracks the number breakdown in each cluster
    int* clusID_of_type;      // Tracks the clusterID of largest of each type
    int* largestClus_of_type; // Tracks the largest number seen for each type

    clus_numof_MolType  = malloc(tot_chain_types_glb * sizeof(lInt));
    clusID_of_type      = malloc(tot_chain_types_glb * sizeof(lInt));
    largestClus_of_type = malloc(tot_chain_types_glb * sizeof(lInt));

    for (i = 0; i < tot_chain_types_glb; i++)
        {
            clus_numof_MolType[i]  = 0;
            clusID_of_type[i]      = 0;
            largestClus_of_type[i] = 0;
        }

    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Proximity_ChainCluster_ForTotal_IntOnly(curID); // This is the length of curID cluster

            for (i = 0; i < tot_chain_types_glb; i++)
                {
                    clus_numof_MolType[i] = 0;
                }

            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    naClusterMatrix_glb[ClusNum][i + 1] = naList_glb[i];
                    molType                             = chain_info_glb[naList_glb[i]][CHAIN_TYPE];
                    clus_numof_MolType[molType]++;
                    naList_glb[i] = -1;
                }

            for (i = 0; i < tot_chain_types_glb; i++)
                {
                    if (clus_numof_MolType[i] > largestClus_of_type[i])
                        {
                            largestClus_of_type[i] = clus_numof_MolType[i];
                            clusID_of_type[i]      = ClusNum;
                        }
                    clus_numof_MolType[i] = 0;
                }

            naClusterMatrix_glb[ClusNum++][0] = Cluster_length;
            IsUnique                          = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }

    for (i = 0; i < tot_chain_types_glb; i++)
        {
            // printf("(%d %d)\t", largestClus_of_type[i],clusID_of_type[i]);
            naList_glb[i] = clusID_of_type[i]; // Now naList_glb contains the cluster ID's.
        }
    // printf("\n");
    free(clus_numof_MolType);
    free(clusID_of_type);
    free(largestClus_of_type);
}

void Clus_Proximity_All_MolWise_LargestClusters(void)
{
    int curID, Cluster_length, i;
    int ClusNum  = 0;
    int IsUnique = 1;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naChainCheckList_glb[i] = -1;
            naList_glb[i]           = -1;
        }
    curID = 0; // Start with the 0th chain

    int molType = 0;
    int* clus_numof_MolType;  // Tracks the number breakdown in each cluster
    int* clusID_of_type;      // Tracks the clusterID of largest of each type
    int* largestClus_of_type; // Tracks the largest number seen for each type

    clus_numof_MolType  = malloc(tot_chain_types_glb * sizeof(lInt));
    clusID_of_type      = malloc(tot_chain_types_glb * sizeof(lInt));
    largestClus_of_type = malloc(tot_chain_types_glb * sizeof(lInt));

    for (i = 0; i < tot_chain_types_glb; i++)
        {
            clus_numof_MolType[i]  = 0;
            clusID_of_type[i]      = 0;
            largestClus_of_type[i] = 0;
        }

    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Proximity_ChainCluster_ForTotal_All(curID); // This is the length of curID cluster

            for (i = 0; i < tot_chain_types_glb; i++)
                {
                    clus_numof_MolType[i] = 0;
                }

            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    naClusterMatrix_glb[ClusNum][i + 1] = naList_glb[i];
                    molType                             = chain_info_glb[naList_glb[i]][CHAIN_TYPE];
                    clus_numof_MolType[molType]++;
                    naList_glb[i] = -1;
                }

            for (i = 0; i < tot_chain_types_glb; i++)
                {
                    if (clus_numof_MolType[i] > largestClus_of_type[i])
                        {
                            largestClus_of_type[i] = clus_numof_MolType[i];
                            clusID_of_type[i]      = ClusNum;
                        }
                    clus_numof_MolType[i] = 0;
                }

            naClusterMatrix_glb[ClusNum++][0] = Cluster_length;
            IsUnique                          = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }

    for (i = 0; i < tot_chain_types_glb; i++)
        {
            // printf("(%d %d)\t", largestClus_of_type[i],clusID_of_type[i]);
            naList_glb[i] = clusID_of_type[i]; // Now naList_glb contains the cluster ID's.
        }
    // printf("\n");
    free(clus_numof_MolType);
    free(clusID_of_type);
    free(largestClus_of_type);
}

///
/// \param nMode
void ClusAnalysis_Perform_Analysis(const int nMode)
{
    // Just a function that picks the right analysis routine given the mode

    int largestNow;
    switch (nMode)
        {
            case 1:
                largestNow =
                    ClusAnalysis_Ovlp_ForSystem_MolTypeWiseDecompAndSizes(laClusHistList_glb, ldaMOLCLUS_Arr_glb);
                break;
            case 2:
                largestNow =
                    ClusAnalysis_Ovlp_ForSystem_MolTypeWiseDecompAndSizes(laClusHistList_glb, ldaMOLCLUS_Arr_glb);
                break;
            default:
                largestNow =
                    ClusAnalysis_Aniso_ForSystem_MolTypeWiseDecompAndSizes(laClusHistList_glb, ldaMOLCLUS_Arr_glb);
                break;
        }

    nTotClusCounter_glb++;
    nLargestClusterRightNow_glb += largestNow;
}

int Clus_Perform_ChainCluster_ForTotal(int const chainID)
{
    int clus_size;

    switch (nClusteringMode_glb)
        {
            case 1:
                clus_size = Clus_Proximity_ChainCluster_ForTotal_IntOnly(chainID);
                break;
            case 2:
                clus_size = Clus_Proximity_ChainCluster_ForTotal_All(chainID);
                break;
            default:
                clus_size = Clus_Network_ChainCluster_ForTotal(chainID);
                break;
        }
    return clus_size;
}

/// Clus_Perform_MolWise_LargestClusters - performs a total clustering analysis
/// of the system. Then finds out the largest cluster for each MolType, where
/// redundancy is allowed. naList_glb contains the cluster IDs of the chain types.
/// Note that naList_glb[0] contains the system's overall largest cluster, naList_glb[1]
/// is for molType=0 and so on. Furthermore, naClusterMatrix_glb[clusID][0] = Cluster
/// size, and naClusterMatrix_glb[1:ClusterSize] are all the chainID's that comprise the
/// particular cluster.
void Clus_Perform_MolWise_LargestClusters(void)
{
    int curID, Cluster_length, i;
    int ClusNum  = 0;
    int IsUnique = 1;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naChainCheckList_glb[i] = -1;
            naList_glb[i]           = -1;
        }
    curID = 0; // Start with the 0th chain

    int molType = 0;
    int* clus_numof_MolType;  // Tracks the number breakdown in each cluster
    int* clusID_of_type;      // Tracks the clusterID of largest of each type
    int* largestClus_of_type; // Tracks the largest number seen for each type
    clus_numof_MolType  = (int*) calloc((tot_chain_types_glb + 1), sizeof(int));
    clusID_of_type      = (int*) calloc((tot_chain_types_glb + 1), sizeof(int));
    largestClus_of_type = (int*) calloc((tot_chain_types_glb + 1), sizeof(int));

    for (i = 0; i <= tot_chain_types_glb; i++)
        {
            clus_numof_MolType[i]  = 0;
            clusID_of_type[i]      = 0;
            largestClus_of_type[i] = 0;
        }

    while (curID < tot_chains_glb && IsUnique == 1)
        {
            Cluster_length = Clus_Perform_ChainCluster_ForTotal(curID); // This is the length of curID cluster

            for (i = 0; i <= tot_chain_types_glb; i++)
                {
                    clus_numof_MolType[i] = 0;
                }

            clus_numof_MolType[0] = Cluster_length;
            for (i = 0; i < Cluster_length; i++)
                { // Recording the chains in this cluster
                    naClusterMatrix_glb[ClusNum][i + 1] = naList_glb[i];
                    molType                             = chain_info_glb[naList_glb[i]][CHAIN_TYPE];
                    clus_numof_MolType[molType + 1]++;
                    naList_glb[i] = -1;
                }

            for (i = 0; i <= tot_chain_types_glb; i++)
                {
                    if (clus_numof_MolType[i] > largestClus_of_type[i])
                        {
                            largestClus_of_type[i] = clus_numof_MolType[i];
                            clusID_of_type[i]      = ClusNum;
                        }
                    clus_numof_MolType[i] = 0;
                }

            naClusterMatrix_glb[ClusNum++][0] = Cluster_length;
            IsUnique                          = 0; // Assume not unique -- just got analyzed.
            while (curID < tot_chains_glb && IsUnique == 0)
                { // Finding the next chainID that hasn't been analyzed.
                    curID++;
                    IsUnique = 0; // Assume not unique.
                    if (naChainCheckList_glb[curID] == -1)
                        {
                            IsUnique = 1;
                        }
                }
        }

    for (i = 0; i <= tot_chain_types_glb; i++)
        {
            naList_glb[i] = clusID_of_type[i]; // Now naList_glb contains the cluster ID's.
        }

    free(clus_numof_MolType);
    free(clusID_of_type);
    free(largestClus_of_type);
}

void Clus_Find_LargestClusters(void)
{
    // Just a function that picks the right analysis routine given the
    // clustering mode

    switch (nClusteringMode_glb)
        {
            case 1:
                Clus_Network_MolWise_LargestClusters();
                break;
            case 2:
                Clus_Proximity_IntOnly_MolWise_LargestClusters();
                break;
            default:
                Clus_Proximity_All_MolWise_LargestClusters();
                break;
        }
}

/// ClusUtil_GetOvlpNeighborBeads_ForBead: Given the beadID, we get the chainIDs of the neighboring beads and store
/// the values in neighList. Note that the list have beadID in it as well.
/// \param beadID
/// \param neighList
/// \return The number of neighbors for this bead + 1 since we always have beadID.
int ClusUtil_GetOvlpNeighborBeads_ForBead(const int beadID, int* restrict neighList)
{
    const int r_pos0[POS_MAX] = {bead_info_glb[beadID][POS_X], bead_info_glb[beadID][POS_Y],
                                 bead_info_glb[beadID][POS_Z]};

    return LatticeUtil_GetNeighBeads_AtPos(r_pos0, neighList);
}

/// ClusUtil_GetOvlpNeighborChains_ForBead - Given this beadID, we get all the chainIDs of it's neighbors,
/// including it's own chain. We return the number of chains.
/// \param beadID
/// \param neighList
/// \return Number of ChainIDs. This should always be >= 1.
int ClusUtil_GetOvlpNeighborChains_ForBead(const int beadID, int* restrict neighList)
{

    int tmpBeadList[CLUS_CONTACT_NEIGHS + 1];

    const int tmpNum = ClusUtil_GetOvlpNeighborBeads_ForBead(beadID, tmpBeadList);

    BeadListOP_GetChainIDs(tmpNum, tmpBeadList, neighList);

    return tmpNum;
}

/// ClusUtil_AddOvlpCluster_OfBead - Given this beadID, and a hash-table like structure caTotClusTable we add the new
/// and unique chains that are part of the cluster. caTotClusTable should be an array that is 'tot_chains_glb' long
/// where caTotClusTable[chainID] 0 means a unique chain, and 1 means that the chain has already been added.
/// Furthermore, also adds the chains into clusList, while increasing the total clusSize.
/// Lastly, the function returns the number of unique chains found.
/// \param beadID
/// \param caTotClusTable
/// \param clusList
/// \param clusSize
/// \return
int ClusUtil_AddOvlpCluster_OfBead(const int beadID, char* restrict caTotClusTable, int* restrict clusList,
                                   int* restrict clusSize)
{
    int tmpChainList[CLUS_CONTACT_NEIGHS];
    const int dumChainNum = ClusUtil_GetOvlpNeighborChains_ForBead(beadID, tmpChainList);

    int i;
    const int startSize = *clusSize;
    int tmpChain;
    for (i = 0; i < dumChainNum; ++i)
        {
            tmpChain = tmpChainList[i];
            if (! caTotClusTable[tmpChain])
                {
                    caTotClusTable[tmpChain] = 1;
                    clusList[*clusSize]      = tmpChain;
                    (*clusSize)++;
                }
        }

    const int endSize = *clusSize - startSize;

    return endSize;
}

/// ClusUtil_AddOvlpCluster_OfBead_CheckForSame -
/// \param beadID
/// \param caTotClusTable
/// \return
int ClusUtil_AddOvlpCluster_OfBead_CheckForSame(const int beadID, const char* const caTotClusTable)
{
    int tmpChainList[CLUS_CONTACT_NEIGHS];
    const int dumChainNum = ClusUtil_GetOvlpNeighborChains_ForBead(beadID, tmpChainList);

    int i;
    int tmpChain;
    for (i = 0; i < dumChainNum; ++i)
        {
            tmpChain = tmpChainList[i];
            if (! caTotClusTable[tmpChain])
                {
                    return -1;
                }
        }

    return 1;
}

/// ClusUtil_AddOvlpCluster_OfChain - Adds all the unique chains within +-1 of each bead of this chain.
/// caTotClusTable is a sort of hash-table used to track if a chainID has been added or not.
/// Returns the number of chains that are _unique_ given caTotClusTable.
/// \param chainID
/// \param caTotClusTable
/// \param clusList
/// \param clusSize
/// \return
int ClusUtil_AddOvlpCluster_OfChain(const int chainID, char* restrict caTotClusTable, int* restrict clusList,
                                    int* restrict clusSize)
{
    const int startSize = *clusSize;

    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chain_info_glb[chainID][CHAIN_LENGTH];

    int i;
    for (i = firstB; i < lastB; ++i)
        {
            ClusUtil_AddOvlpCluster_OfBead(i, caTotClusTable, clusList, clusSize);
        }

    const int endSize = *clusSize - startSize;

    return endSize;
}

///
/// \param chainID
/// \param caTotClusTable
/// \return
int ClusUtil_AddOvlpCluster_OfChain_CheckForSame(const int chainID, const char* restrict const caTotClusTable)
{
    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chain_info_glb[chainID][CHAIN_LENGTH];

    int i;
    int tmpVal;
    for (i = firstB; i < lastB; ++i)
        {
            tmpVal = ClusUtil_AddOvlpCluster_OfBead_CheckForSame(i, caTotClusTable);
            if (tmpVal == -1)
                {
                    return -1;
                }
        }

    return 1;
}

/// ClusUtil_OvlpCluster_OfChain Calculate the Ovlp-bases cluster of molecule chainID. Returns the size of the cluster.
/// The cluster will be written on naClusList, and caTotClusTable will contain the cluster for chainID. \param chainID
/// \param caTotClusTable Boolean-like array that contains which chains have been already looked at in the system.
/// \param naClusList Array that will contain the chain indices of the chains in the cluster.
/// \return Total size of the cluster chainID is a part of.
int ClusUtil_OvlpCluster_OfChain(const int chainID, char* restrict caTotClusTable, int* restrict naClusList)
{
    naClusList[0]           = chainID;
    caTotClusTable[chainID] = 1;

    int clusSize = 1;

    int tmpChain;
    int i;
    for (i = 0; i < clusSize; ++i)
        {
            tmpChain = naClusList[i];
            ClusUtil_AddOvlpCluster_OfChain(tmpChain, caTotClusTable, naClusList, &clusSize);
        }

    return clusSize;
}

///
/// \param chainID
/// \param caTotClusTable
/// \param clusList
/// \param maxSize
/// \return
int ClusUtil_OvlpCluster_OfChain_wMaxSize(const int chainID, char* restrict caTotClusTable, int* restrict clusList,
                                          int const maxSize)
{
    clusList[0]             = chainID;
    caTotClusTable[chainID] = 1;

    int clusSize = 1;

    int tmpChain;
    int i;
    for (i = 0; i < clusSize; ++i)
        {
            tmpChain = clusList[i];
            ClusUtil_AddOvlpCluster_OfChain(tmpChain, caTotClusTable, clusList, &clusSize);
            if (clusSize >= maxSize)
                {
                    return 0;
                }
        }

    return clusSize;
}

///
/// \param caTotClusTable
/// \param clusList
/// \param clusSize
/// \return
int ClusUtil_OvlpCluster_OfChain_CheckForSame(const char* const caTotClusTable, const int* const clusList,
                                              int const clusSize)
{

    int tmpChain;
    int tmpVal;
    int i;
    for (i = 0; i < clusSize; ++i)
        {
            tmpChain = clusList[i];
            tmpVal   = ClusUtil_AddOvlpCluster_OfChain_CheckForSame(tmpChain, caTotClusTable);
            if (tmpVal == -1)
                {
                    return -1;
                }
        }

    return 1;
}

///
/// \param beadID
/// \param caTotClusTable
/// \param clusList
/// \param clusSize
/// \return
int ClusUtil_AddAnisoCluster_OfBead(const int beadID, char* restrict caTotClusTable, int* restrict clusList,
                                    int* restrict clusSize)
{
    int const beadPartner = bead_info_glb[beadID][BEAD_FACE];

    if (beadPartner != -1)
        {
            int const tmpChain = bead_info_glb[beadPartner][BEAD_CHAINID];
            if (! caTotClusTable[tmpChain])
                {
                    caTotClusTable[tmpChain] = 1;
                    clusList[*clusSize]      = tmpChain;
                    (*clusSize)++;
                    return 1;
                }
        }
    return 0;
}

///
/// \param chainID
/// \param caTotClusTable
/// \param clusList
/// \param clusSize
/// \return
int ClusUtil_AddAnisoCluster_OfChain(const int chainID, char* restrict caTotClusTable, int* restrict clusList,
                                     int* restrict clusSize)
{
    const int startSize = *clusSize;

    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chain_info_glb[chainID][CHAIN_LENGTH];

    int i;
    for (i = firstB; i < lastB; ++i)
        {
            ClusUtil_AddAnisoCluster_OfBead(i, caTotClusTable, clusList, clusSize);
        }

    const int endSize = *clusSize - startSize;

    return endSize;
}

///
/// \param chainID
/// \param caTotClusTable
/// \param clusList
/// \return
int ClusUtil_AnisoCluster_OfChain(const int chainID, char* restrict caTotClusTable, int* restrict clusList)
{
    clusList[0]             = chainID;
    caTotClusTable[chainID] = 1;

    int clusSize = 1;

    int tmpChain;
    int i;
    for (i = 0; i < clusSize; ++i)
        {
            tmpChain = clusList[i];
            ClusUtil_AddAnisoCluster_OfChain(tmpChain, caTotClusTable, clusList, &clusSize);
        }

    return clusSize;
}

///
/// \param chainID
/// \param caTotClusTable
/// \param clusList
/// \param maxSize
/// \return
int ClusUtil_AnisoCluster_OfChain_wMaxSize(const int chainID, char* restrict caTotClusTable, int* restrict clusList,
                                           int const maxSize)
{
    clusList[0]             = chainID;
    caTotClusTable[chainID] = 1;

    int clusSize = 1;

    int tmpChain;
    int i;
    for (i = 0; i < clusSize; ++i)
        {
            tmpChain = clusList[i];
            ClusUtil_AddAnisoCluster_OfChain(tmpChain, caTotClusTable, clusList, &clusSize);
            if (clusSize >= maxSize)
                {
                    return 0;
                }
        }

    return clusSize;
}

/// ClusUtil_NextUnvisitedChain - Given that we are on chain nChainID, we use caTotClusTable to find the index of the
/// next unvisited chain.
/// \param nChainID Starting chain.
/// \param caTotClusTable Boolean-like array that stores if a chain has been visited before.
/// \return Index of an unvisited chain.
int ClusUtil_NextUnvisitedChain(int const nChainID, const char* const caTotClusTable)
{
    int nNewChainID = nChainID;

    while (caTotClusTable[nNewChainID] == 1)
        {
            nNewChainID++;
        }

    return nNewChainID;
}

/// ClusUtil_AnisoClusters_OfSystem - Calculates the Aniso-based cluster of the whole system. Writes over the provided
/// arrays to save the cluster-list and the cumulative sizes. \param naFullClusList Array that will contain all the
/// clusters. \param naCumClusSizes Array that will contain the cumulative sizes of the clusters. Allows for cluster
/// extraction later by doing some index-based arithmetic. \return The total number of clusters in the system. 1 would
/// mean every molecule is part of a single cluster. While tot_chains_glb would mean that every molecule is only a part
/// of its own cluster.
int ClusUtil_AnisoClusters_OfSystem(int* const restrict naFullClusList, int* const restrict naCumClusSizes)
{
    int nCumulativeSize = 0;
    int nClusNum        = 0;

    char* caTmpClusCheckList = (char*) calloc(tot_chains_glb + 1, sizeof(char));

    int nThisChainID = 0;
    int nThisClusSize;

    while (nThisChainID < tot_chains_glb)
        {
            nThisClusSize =
                ClusUtil_AnisoCluster_OfChain(nThisChainID, caTmpClusCheckList, naFullClusList + nCumulativeSize);

            nCumulativeSize += nThisClusSize;
            nClusNum++;
            naCumClusSizes[nClusNum] = nCumulativeSize;

            nThisChainID = ClusUtil_NextUnvisitedChain(nThisChainID, caTmpClusCheckList);
        }

    free(caTmpClusCheckList);

    return nClusNum;
}

/// ClusUtil_OvlpClusters_OfSystem - Calculates the Ovlp-based cluster of the whole system. Writes over the provided
/// arrays to save the cluster-list and the cumulative sizes. \param naFullClusList Array that will contain all the
/// clusters. \param naCumClusSizes Array that will contain the cumulative sizes of the clusters. Allows for cluster
/// extraction later by doing some index-based arithmetic. \return The total number of clusters in the system. 1 would
/// mean every molecule is part of a single cluster. While tot_chains_glb would mean that every molecule is only a part
/// of its own cluster.
int ClusUtil_OvlpClusters_OfSystem(int* const naFullClusList, int* const naCumClusSizes)
{
    int nCumulativeSize = 0;
    int nClusNum        = 0;

    char* caTmpClusCheckList = (char*) calloc(tot_chains_glb + 1, sizeof(char));

    int nThisChainID = 0;
    int nThisClusSize;

    while (nThisChainID < tot_chains_glb)
        {
            nThisClusSize =
                ClusUtil_OvlpCluster_OfChain(nThisChainID, caTmpClusCheckList, naFullClusList + nCumulativeSize);

            nCumulativeSize += nThisClusSize;
            nClusNum++;
            naCumClusSizes[nClusNum] = nCumulativeSize;

            nThisChainID = ClusUtil_NextUnvisitedChain(nThisChainID, caTmpClusCheckList);
        }

    free(caTmpClusCheckList);

    return nClusNum;
}

/// ClusUtil_GenClusSizesFromCumulativeSizes - Given the total number of clusters nClusNum, and the cumulative-sizes of
/// the clusters naCumClusSizes, we update the naSizeList to contain all the cluster sizes. Remember that the sizes
/// will be unordered. size[i] = cum_size[i+1]-cum_size[i]
/// \param naSizeList Will contain sizes of clusters.
/// \param naCumClusSizes Contains cumulative sizes.
/// \param nClusNum Total number of clusters.
void ClusUtil_GenClusSizesFromCumulativeSizes(int* const naSizeList, const int* const naCumClusSizes,
                                              int const nClusNum)
{
    int i;

    for (i = 0; i < nClusNum; i++)
        {
            naSizeList[i] = naCumClusSizes[i + 1] - naCumClusSizes[i];
        }
}

/// ClusHistUtil_AddToHist_FromCountsList - Using the occurences in naCounts, we add to the histogram laHist, which
/// is a long-double. naCounts[i] represents the index of the histogram.
/// \param laHist Histogram to be added on to.
/// \param naCounts Array that has the occurences.
/// \param nBins Total number of bins to look at.
void ClusHistUtil_AddToHist_FromCountsList(lLong* const laHist, const int* const naCounts, const int nBins)
{
    int i;
    for (i = 0; i < nBins; i++)
        {
            laHist[naCounts[i]]++;
        }
}

///
/// \param naClusSizeHist
/// \param naCumClusSizes
/// \param nClusNum
void ClusUtil_GenHistFromCumulativeSizes(lLong* const restrict naClusSizeHist, const int* const restrict naCumClusSizes,
                                         int const nClusNum)
{
    int* naTmpSizeList = (int*) calloc((nClusNum), sizeof(int));

    ClusUtil_GenClusSizesFromCumulativeSizes(naTmpSizeList, naCumClusSizes, nClusNum);

    ClusHistUtil_AddToHist_FromCountsList(naClusSizeHist, naTmpSizeList, nClusNum);

    free(naTmpSizeList);
}

/// ClusUtil_GetCluster_FromFullClusAndCumSizes - Given the cluster-list, cumulative sizes, cluster sizes and total
/// number of clusters, we get the cluster specific to nClusID and also return the size.
/// \param nClusID Index of cluster in naFullClusList
/// \param naOutClusList Array that will contain the chainIDs of the cluster
/// \param naFullClusList Array that contains the clusters
/// \param naCumSizes Array containing the cumulative sizes of the clusters.
/// \param naClusSizes Array containing the cluster sizes.
/// \param nClusNum Total number of clusters. Mostly used for debugging.
/// \return The size of this particular cluster.
int ClusUtil_GetCluster_FromFullClusAndCumSizes(const int nClusID, int* const naOutClusList,
                                                const int* const naFullClusList, const int* const naCumSizes,
                                                const int* const naClusSizes, const int nClusNum)
{
#if DEBUG_BUILD
    if (nClusID > nClusNum)
        {
            fprintf(stderr, "Incorrect cluster-ID (%d) is greater than number of clusters (%d)\n", nClusID, nClusNum);
            fprintf(stderr, "Crashing.");
            exit(1);
        }
    if (nClusID < 0)
        {
            fprintf(stderr, "Negative cluster-ID (%d)\n", nClusID);
            fprintf(stderr, "Crashing.");
            exit(1);
        }
#endif

    int i;
    const int thisSize = naClusSizes[nClusID];
    const int nStart   = naCumSizes[nClusID];

    for (i = 0; i < thisSize; i++)
        {
            naOutClusList[i] = naFullClusList[i + nStart];
        }

    return thisSize;
}

/// ClusHistUtil_AddToHist_MolTypeWiseDecomp_FromCluster - For the cluster in naClusList (of size nClusSize) we
/// calculate the cluster composition, and add to the ldaMolWiseHist array.
/// \param ldaMolWiseHist Histogram array. Bins are calculated using MolClusArr_Index()
/// \param naClusList The chainIDs of this cluster.
/// \param nClusSize The size of this cluster.
void ClusHistUtil_AddToHist_MolTypeWiseDecomp_FromCluster(lLDub* const restrict ldaMolWiseHist,
                                                          const int* const restrict naClusList, const int nClusSize)
{
#ifdef DEBUG_BUILD
    if (nClusSize < 1)
        {
            fprintf(stderr, "Cluster size is %d\n", nClusSize);
            fputs("Cluster sizes have to be positive integers\n", stderr);
            fputs("Crashing!\n", stderr);
            exit(1);
        }
#endif
    int i;
    const int clusSizeBin = nClusSize - 1;
    int chainID, chainType, histBin;
    for (i = 0; i < nClusSize; i++)
        {
            chainID   = naClusList[i];
            chainType = chain_info_glb[chainID][CHAIN_TYPE];
            histBin   = MolClusArr_Index(0, chainType, clusSizeBin);
            ldaMolWiseHist[histBin]++;
        }
}

/// ClusHistUtil_AddToHist_MolTypeWiseDecomp_FromFullClusAndCumSizes - Goes over every cluster and adds the cluster
/// composition to ldaMolWiseHist array.
/// \param ldaMolWiseHist
/// \param naClusChainIDs Contains the clusters (chainIDs within every cluster)
/// \param naCumClusSizes  Contains the cumulative sizes given the clusters.
/// \param naClusSizes Contains the cluster sizes.
/// \param nClusNum Contains the total number of clusters.
void ClusHistUtil_AddToHist_MolTypeWiseDecomp_FromFullClusAndCumSizes(lLDub* const ldaMolWiseHist,
                                                                      const int* const naClusChainIDs,
                                                                      const int* const naCumClusSizes,
                                                                      const int* const naClusSizes, int const nClusNum)
{
    int* thisCluster = calloc(tot_chains_glb + 1, sizeof(int));

    int i;
    int thisSize;
    for (i = 0; i < nClusNum; i++)
        {
            thisSize = ClusUtil_GetCluster_FromFullClusAndCumSizes(i, thisCluster, naClusChainIDs, naCumClusSizes,
                                                                   naClusSizes, nClusNum);
            ClusHistUtil_AddToHist_MolTypeWiseDecomp_FromCluster(ldaMolWiseHist, thisCluster, thisSize);
        }

    free(thisCluster);
}

/// ClusAnalysis_Ovlp_ForSystem_MolTypeWiseDecompAndSizes - Performs a total clustering analysis of the system, and also
/// gets the composition of each cluster. naSizeHist is the histogram for the cluster sizes. ldaMolWiseHist is the
/// histogram for the composition of each cluster. The composition corresponds to the number of different chain types
/// within each cluster (or cluster-size)
/// \param naSizeHist
/// \param ldaMolWiseHist
int ClusAnalysis_Ovlp_ForSystem_MolTypeWiseDecompAndSizes(lLong* const restrict naSizeHist,
                                                          lLDub* const restrict ldaMolWiseHist)
{
    int* naFullClusList = (int*) calloc((tot_chains_glb + 1), sizeof(int));
    int* naCumClusSizes = (int*) calloc((tot_chains_glb + 1), sizeof(int));

    const int nClusNum = ClusUtil_OvlpClusters_OfSystem(naFullClusList, naCumClusSizes);

    int* naClusSizes = (int*) calloc((nClusNum + 1), sizeof(int));

    ClusUtil_GenClusSizesFromCumulativeSizes(naClusSizes, naCumClusSizes, nClusNum);

    ClusHistUtil_AddToHist_FromCountsList(naSizeHist, naClusSizes, nClusNum);

    ClusHistUtil_AddToHist_MolTypeWiseDecomp_FromFullClusAndCumSizes(ldaMolWiseHist, naFullClusList, naCumClusSizes,
                                                                     naClusSizes, nClusNum);

    const int currentLargest = ListOP_GetMaxVal_Int(nClusNum, naClusSizes);

    free(naFullClusList);
    free(naCumClusSizes);
    free(naClusSizes);

    return currentLargest;
}

/// ClusAnalysis_Aniso_ForSystem_MolTypeWiseDecompAndSizes - Performs a total clustering analysis of the system, and
/// also gets the composition of each cluster. naSizeHist is the histogram for the cluster sizes. ldaMolWiseHist is the
/// histogram for the composition of each cluster. The composition corresponds to the number of different chain types
/// within each cluster (or cluster-size)
/// \param naSizeHist
/// \param ldaMolWiseHist
int ClusAnalysis_Aniso_ForSystem_MolTypeWiseDecompAndSizes(lLong* const restrict naSizeHist,
                                                           lLDub* const restrict ldaMolWiseHist)
{
    int* naFullClusList = (int*) calloc((tot_chains_glb + 1), sizeof(int));
    int* naCumClusSizes = (int*) calloc((tot_chains_glb + 1), sizeof(int));

    const int nClusNum = ClusUtil_AnisoClusters_OfSystem(naFullClusList, naCumClusSizes);

    int* naClusSizes = (int*) calloc((nClusNum + 1), sizeof(int));

    ClusUtil_GenClusSizesFromCumulativeSizes(naClusSizes, naCumClusSizes, nClusNum);

    ClusHistUtil_AddToHist_FromCountsList(naSizeHist, naClusSizes, nClusNum);

    ClusHistUtil_AddToHist_MolTypeWiseDecomp_FromFullClusAndCumSizes(ldaMolWiseHist, naFullClusList, naCumClusSizes,
                                                                     naClusSizes, nClusNum);

    const int currentLargest = ListOP_GetMaxVal_Int(nClusNum, naClusSizes);

    free(naFullClusList);
    free(naCumClusSizes);
    free(naClusSizes);

    return currentLargest;
}

///
/// \param naOutClusList
/// \param naFullClusList
/// \param naCumSizes
/// \param naClusSizes
/// \param nClusNum
/// \return
int ClusUtil_GetLargestCluster_FromFullClusAndCumSizes(int* const restrict naOutClusList,
                                                       const int* const restrict naFullClusList,
                                                       const int* const restrict naCumSizes,
                                                       const int* const restrict naClusSizes, const int nClusNum)
{
    int i;
    const int thisSize = ListOP_GetMaxVal_Int(nClusNum, naClusSizes);
    const int nClusID  = ListOP_GetRandomIndexForVal_Int(thisSize, nClusNum, naClusSizes);
#if DEBUG_BUILD
    if (nClusID > nClusNum)
        {
            fprintf(stderr, "Incorrect cluster-ID (%d) is greater than number of clusters (%d)\n", nClusID, nClusNum);
            fprintf(stderr, "Crashing.");
            exit(1);
        }
    if (nClusID < 0)
        {
            fprintf(stderr, "Negative cluster-ID (%d)\n", nClusID);
            fprintf(stderr, "Crashing.");
            exit(1);
        }
#endif

    const int nStart = naCumSizes[nClusID];

    for (i = 0; i < thisSize; i++)
        {
            naOutClusList[i] = naFullClusList[i + nStart];
        }

    return thisSize;
}

///
/// \param naOutClusList
/// \param naFullClusList
/// \param naCumSizes
/// \param naClusSizes
/// \param nClusNum
/// \return
int ClusUtil_GetSecondLargestCluster_FromFullClusAndCumSizes(int* const restrict naOutClusList,
                                                             const int* const restrict naFullClusList,
                                                             const int* const restrict naCumSizes,
                                                             const int* const restrict naClusSizes, const int nClusNum)
{
    int i;

    const int thisSize = ListOP_Get2ndLargestVal_Int(nClusNum, naClusSizes);
    const int nClusID  = ListOP_GetRandomIndexForVal_Int(thisSize, nClusNum, naClusSizes);

#if DEBUG_BUILD
    if (nClusID > nClusNum)
        {
            fprintf(stderr, "Incorrect cluster-ID (%d) is greater than number of clusters (%d)\n", nClusID, nClusNum);
            fprintf(stderr, "Crashing.");
            exit(1);
        }
    if (nClusID < 0)
        {
            fprintf(stderr, "Negative cluster-ID (%d)\n", nClusID);
            fprintf(stderr, "Crashing.");
            exit(1);
        }
#endif
    const int nStart = naCumSizes[nClusID];

    for (i = 0; i < thisSize; i++)
        {
            naOutClusList[i] = naFullClusList[i + nStart];
        }

    return thisSize;
}

int ClusUtil_OfSystem_SecondLargest(int* const naOutClusList, const int nMode)
{
    int clusterSize;
    switch (nMode)
        {
            case 1:
                clusterSize = ClusUtil_OvlpCluster_OfSystem_SecondLargest(naOutClusList);
                break;
            case 2:
                clusterSize = ClusUtil_OvlpCluster_OfSystem_SecondLargest(naOutClusList);
                break;
            default:
                clusterSize = ClusUtil_AnisoCluster_OfSystem_SecondLargest(naOutClusList);
                break;
        }

    return clusterSize;
}

/// ClusUtil_OvlpCluster_OfSystem_SecondLargest
/// \param naOutClusList
/// \return Size of the cluster
int ClusUtil_OvlpCluster_OfSystem_SecondLargest(int* naOutClusList)
{
    int thisSize;

    int* naFullClusList = (int*) calloc((tot_chains_glb + 1), sizeof(int));
    int* naCumClusSizes = (int*) calloc((tot_chains_glb + 1), sizeof(int));

    const int nClusNum = ClusUtil_OvlpClusters_OfSystem(naFullClusList, naCumClusSizes);

    int* naClusSizes = (int*) calloc((nClusNum + 1), sizeof(int));

    ClusUtil_GenClusSizesFromCumulativeSizes(naClusSizes, naCumClusSizes, nClusNum);

    thisSize = ClusUtil_GetSecondLargestCluster_FromFullClusAndCumSizes(naOutClusList, naFullClusList, naCumClusSizes,
                                                                        naClusSizes, nClusNum);

    free(naFullClusList);
    free(naCumClusSizes);
    free(naClusSizes);

    return thisSize;
}

///
/// \param naOutClusList
/// \return Size of the cluster. If there is only 1 cluster, return -1 so the move fails.
int ClusUtil_OvlpCluster_OfSystem_SecondLargest_ForMCMove(int* naOutClusList)
{
    int thisSize;

    int* naFullClusList = (int*) calloc((tot_chains_glb + 1), sizeof(int));
    int* naCumClusSizes = (int*) calloc((tot_chains_glb + 1), sizeof(int));

    const int nClusNum = ClusUtil_OvlpClusters_OfSystem(naFullClusList, naCumClusSizes);

    int* naClusSizes = (int*) calloc((nClusNum + 1), sizeof(int));

    if (nClusNum == 1)
        {
            thisSize = -1;
        }
    else
        {
            ClusUtil_GenClusSizesFromCumulativeSizes(naClusSizes, naCumClusSizes, nClusNum);

            thisSize = ClusUtil_GetSecondLargestCluster_FromFullClusAndCumSizes(naOutClusList, naFullClusList,
                                                                                naCumClusSizes, naClusSizes, nClusNum);
        }

    free(naFullClusList);
    free(naCumClusSizes);
    free(naClusSizes);

    return thisSize;
}

///
/// \param naOutClusList
/// \return Size of the cluster
int ClusUtil_AnisoCluster_OfSystem_SecondLargest(int* naOutClusList)
{
    int thisSize;

    int* naFullClusList = (int*) calloc((tot_chains_glb + 1), sizeof(int));
    int* naCumClusSizes = (int*) calloc((tot_chains_glb + 1), sizeof(int));

    const int nClusNum = ClusUtil_AnisoClusters_OfSystem(naFullClusList, naCumClusSizes);

    int* naClusSizes = (int*) calloc((nClusNum + 1), sizeof(int));

    ClusUtil_GenClusSizesFromCumulativeSizes(naClusSizes, naCumClusSizes, nClusNum);

    thisSize = ClusUtil_GetSecondLargestCluster_FromFullClusAndCumSizes(naOutClusList, naFullClusList, naCumClusSizes,
                                                                        naClusSizes, nClusNum);

    free(naFullClusList);
    free(naCumClusSizes);
    free(naClusSizes);

    return thisSize;
}

///
/// \param naOutClusList
/// \return Size of the cluster. If there is only 1 cluster, return -1 so the move fails.
int ClusUtil_AnisoCluster_OfSystem_SecondLargest_ForMCMove(int* naOutClusList)
{
    int thisSize;

    int* naFullClusList = (int*) calloc((tot_chains_glb + 1), sizeof(int));
    int* naCumClusSizes = (int*) calloc((tot_chains_glb + 1), sizeof(int));

    const int nClusNum = ClusUtil_AnisoClusters_OfSystem(naFullClusList, naCumClusSizes);

    int* naClusSizes = (int*) calloc((nClusNum + 1), sizeof(int));

    if (nClusNum == 1)
        {
            thisSize = -1;
        }
    else
        {
            ClusUtil_GenClusSizesFromCumulativeSizes(naClusSizes, naCumClusSizes, nClusNum);

            thisSize = ClusUtil_GetSecondLargestCluster_FromFullClusAndCumSizes(naOutClusList, naFullClusList,
                                                                                naCumClusSizes, naClusSizes, nClusNum);
        }

    free(naFullClusList);
    free(naCumClusSizes);
    free(naClusSizes);

    return thisSize;
}

///
/// \param naClusIDsList_out
/// \param naFullClusList_out
/// \param naClusSizes_out
/// \param naCumClusSizes_out
/// \return
int ClusUtil_OfSystem_MolWise_GetLargestClusters(int* const restrict naClusIDsList_out,
                                                 int* const restrict naFullClusList_out,
                                                 int* const restrict naClusSizes_out,
                                                 int* const restrict naCumClusSizes_out)
{

    const int nClusNum = ClusUtil_OvlpClusters_OfSystem(naFullClusList_out, naCumClusSizes_out);

    ClusUtil_GenClusSizesFromCumulativeSizes(naClusSizes_out, naCumClusSizes_out, nClusNum);

    int* const naChainTypes = (int*) calloc((tot_chains_glb + 1), sizeof(int));

    ChainListOP_GetChainTypes(naChainTypes, naFullClusList_out, tot_chains_glb);

    ClusUtil_MolWise_FindLargestClusters(naClusIDsList_out, naChainTypes, naClusSizes_out, naCumClusSizes_out,
                                         nClusNum);

    free(naChainTypes);

    return nClusNum;
}

/// ClusUtil_SelectSubsetOfClusters - Assuming that we have a complete cluster list, cumulative sizes, and cluster
/// sizes, and a small set of cluster IDs from that complete list, we generate a new _cluster_ list that only has the
/// subset of cluster IDs supplied. The new list will have chainIDs for each of the subset of clusters, and we will also
/// write the cumulative size and cluster size of each of those. Thus, we can use
/// ClusUtil_GetCluster_FromFullClusAndCumSizes on the new cluster with the new sizes. Since this appends the clusters
/// one after the other, the newer cluster list _could_ be larger than the total number of chains in the system by
/// having redundant clusters.
///
/// \param naSubsetClusChains_out
/// \param naSubsetClusSizes_out
/// \param naSubsetCumSizes_out
/// \param nSubsetNum_in
/// \param naSubSetClusIDs_in
/// \param naFullClusList_in
/// \param naClusSizes_in
/// \param naCumClusSizes_in
void ClusUtil_SelectSubsetOfClusters(int* const restrict naSubsetClusChains_out,
                                     int* const restrict naSubsetClusSizes_out,
                                     int* const restrict naSubsetCumSizes_out, const int nSubsetNum_in,
                                     const int* const restrict naSubSetClusIDs_in,
                                     const int* const restrict naFullClusList_in,
                                     const int* const restrict naClusSizes_in,
                                     const int* const restrict naCumClusSizes_in)
{
    int i, thisSize;

    int cumSize             = 0;
    naSubsetCumSizes_out[0] = 0;
    for (i = 0; i < nSubsetNum_in; i++)
        {
            thisSize = ClusUtil_GetCluster_FromFullClusAndCumSizes(naSubSetClusIDs_in[i],
                                                                   naSubsetClusChains_out + cumSize, naFullClusList_in,
                                                                   naCumClusSizes_in, naClusSizes_in, nSubsetNum_in);
            cumSize += thisSize;
            naSubsetClusSizes_out[i]    = thisSize;
            naSubsetCumSizes_out[i + 1] = cumSize;
        }
}

///
/// \param naClusIDs_out
/// \param naChainTypes_in
/// \param naClusSizes_in
/// \param naCumClusSizes_in
/// \param nClusNum_in
void ClusUtil_MolWise_FindLargestClusters(int* const restrict naClusIDs_out, const int* const restrict naChainTypes_in,
                                          const int* const restrict naClusSizes_in,
                                          const int* const restrict naCumClusSizes_in, const int nClusNum_in)
{

    int i;
    int j;
    int thisSize = 0;

    int* const thisClus       = malloc((tot_chains_glb + 1) * sizeof(int));
    int* const naLargestSizes = calloc(tot_chain_types_glb + 1, sizeof(int));
    int* const tmpSizesList   = calloc(tot_chain_types_glb + 1, sizeof(int));

    for (i = 0; i < nClusNum_in; i++)
        {
            thisSize = ClusUtil_GetCluster_FromFullClusAndCumSizes(i, thisClus, naChainTypes_in, naCumClusSizes_in,
                                                                   naClusSizes_in, nClusNum_in);
            tmpSizesList[0] = thisSize;
            ListOP_GenHistFromCounts_Int(tmpSizesList + 1, thisClus, thisSize);

            for (j = 0; j < tot_chain_types_glb + 1; j++)
                {
                    if (tmpSizesList[j] >= naLargestSizes[j])
                        {
                            naLargestSizes[j] = tmpSizesList[j];
                            naClusIDs_out[j]  = i;
                        }
                    tmpSizesList[j] = 0;
                }
        }

    free(thisClus);
    free(naLargestSizes);
    free(tmpSizesList);
}

/// ListOP_ReplaceIfLarger_Int - naOutList will have the max(naOutList[i], naInList[i]).
/// \param naOutList
/// \param naInList
/// \param nSize
void ListOP_ReplaceWithLarger_Int(int* const restrict naOutList, const int* const restrict naInList, const int nSize)
{
    int i;
    for (i = 0; i < nSize; i++)
        {
            if (naInList[i] >= naOutList[i])
                {
                    naOutList[i] = naInList[i];
                }
        }
}
