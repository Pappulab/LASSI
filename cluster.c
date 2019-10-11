#include "global.h"
#include "cluster.h"
#include "structure.h"

int Clus_ChainNetwork_General(int chainID){
  //Updates naList to have all proteins bound to  chainID and it's cluster
  //The idea is to check every bead and see if there is a unique bonded chain, and to add it to naList

 int i, j; //Loop iterators
 for(i=0; i <= tot_chains; i++){
   naList[i] = -1;//Initialize the list where -1 means empty.
   naChainCheckList[i] = -1;//Initialize the list where no chain has been visited.
 }


 int list_it = 0;//Iterator for naList
 int clusSize = 0;//Index to track the cluster size.
 int curID;//Index to track the current chain being looked at.
 curID = chainID; naList[clusSize++] = curID;//The cluster contains chainID by definition, and ClusSize = 1
 naChainCheckList[curID] = 1;//This means that curID has been checked!
 int fB, lB;//Indecies to track the first and last bead of chains.
 int chainPart;

 while (curID != -1){//Keep going through naList till it is exhausted.
  fB = chain_info[curID][CHAIN_START];//First bead of this chain.
  lB  = fB + chain_info[curID][CHAIN_LENGTH];//Last bead+1 of this chain. Makes for-loops easier this way.
  for(i=fB; i<lB; i++){//Loop over all the beads in this chain and see if there is a physical bond.
    if(bead_info[i][BEAD_FACE] != -1){//This means we have a bonding partner.
      chainPart = bead_info[bead_info[i][BEAD_FACE]][BEAD_CHAINID];
      if(naChainCheckList[chainPart] == -1 ){//This is a unique chain.
        naList[clusSize++] = chainPart;
          naChainCheckList[chainPart] = 1;//Recording that this chain has been looked at
      }
     }
     //Moving on to the next bead in this chain
    }
    //Done with this chain, so let's move to the next chain, if it exists.
    list_it++;//Going one forward in naList to completely exhaust the tree.
    curID = naList[list_it];
  }
  return clusSize;
}

int Clus_ChainNetwork_ForTotal(int chainID){
  //Updates naList to have all proteins bound to  chainID and it's cluster
  //The idea is to check every bead and see if there is a unique bonded chain, and to add it to naList
  //This one is specifically made to be used in total system network analyses, so naList and naChainCheckList aren't
  //reinitialized to -1. Instead, we just use naChainCheckList fully.
 int i; //Loop iterators
 int list_it = 0;//Iterator for naList
 int clusSize = 0;//Index to track the cluster size.
 int curID;//Index to track the current chain being looked at.
 curID = chainID; naList[clusSize++] = curID;//The cluster contains chainID by definition, and ClusSize = 1
 naChainCheckList[curID] = 1;//This means that curID has been checked!
 int fB, lB;//Indecies to track the first and last bead of chains.
 int chainPart;

 while (curID != -1){//Keep going through naList till it is exhausted.
  fB = chain_info[curID][CHAIN_START];//First bead of this chain.
  lB  = fB + chain_info[curID][CHAIN_LENGTH];//Last bead+1 of this chain. Makes for-loops easier this way.
  for(i=fB; i<lB; i++){//Loop over all the beads in this chain and see if there is a physical bond.
    if(bead_info[i][BEAD_FACE] != -1){//This means we have a bonding partner.
      chainPart = bead_info[bead_info[i][BEAD_FACE]][BEAD_CHAINID];
      if(naChainCheckList[chainPart] == -1 ){//This is a unique chain.
        naList[clusSize++] = chainPart;
          naChainCheckList[chainPart] = 1;//Recording that this chain has been looked at
      }
     }
     //Moving on to the next bead in this chain
    }
    //Done with this chain, so let's move to the next chain in the cluster, if it exists.
    list_it++;//Going one forward in naList to completely exhaust the tree.
    curID = naList[list_it];
  }
  return clusSize;
}

void Clus_Distribution_Avg(void){
  /*
  Calculates the cluster distribution using total_network_analysis framework, but keeps adding to
  naClusHistList for total averaging at the end. Read total_network_analysis() for what's happening here LOL
  */
  int curID, Cluster_length, currentLargest, i;
  int IsUnique = 1;
  for(i=0; i <= tot_chains; i++){
    naList[i]           = -1;
    naChainCheckList[i] = -1;
  }
  curID = 0;//Start with the 0th chain
  currentLargest = 0;
  while(curID < tot_chains && IsUnique == 1){
    Cluster_length = Clus_ChainNetwork_ForTotal(curID);//This is the length of curID cluster
    //printf("Clus Len: %d\n", Cluster_length);
    naClusHistList[Cluster_length]++; //Adding to that cluster-size bin
    if(Cluster_length > currentLargest){
      currentLargest = Cluster_length;
    }
    IsUnique = 0;//Assume not unique -- just got analyzed.
    while (curID < tot_chains && IsUnique == 0){//Finding the next chainID that hasn't been analyzed.
    curID++;
    IsUnique = 0;//Assume not unique.
    if (naChainCheckList[curID] == -1){
      IsUnique = 1;
    }
    }
  }
  nTotClusCounter++;
  nLargestClusterRightNow += currentLargest;
}

int Clus_LimitedCluster(int chainID){
  //Updates naList to have all proteins bound to  chainID and it's cluster
  //The idea is to check every bead and see if there is a unique bonded chain, and to add it to naList
  //If Cluster becomes larger than 5, exit and return -1

 int i, j; //Loop iterators
 for(i=0; i<5; i++){
   naList[i] = -1;
  }

 int list_it = 0;//Iterator for naList
 int clusSize = 0;//Index to track the cluster size.
 int curID;//Index to track the current chain being looked at.
 curID = chainID; naList[clusSize++] = curID;//The cluster contains chainID by definition, and ClusSize = 1
 int fB, lB;//Indecies to track the first and last bead of chains.
 int chainPart;
 int IsUnique = 1;//Tracks if a chain is unique or not. 0 is non-unique, and 1 is unique.

 while (curID != -1){//Keep going through naList till it is exhausted.
  fB = chain_info[curID][CHAIN_START];//First bead of this chain.
  lB  = fB + chain_info[curID][CHAIN_LENGTH];//Last bead+1 of this chain. Makes for-loops easier this way.

  for(i=fB; i<lB; i++){//Loop over all the beads in this chain and see if there is a physical bond.
    if(bead_info[i][BEAD_FACE] != -1){//This means we have a bonding partner.
      chainPart = bead_info[bead_info[i][BEAD_FACE]][BEAD_CHAINID];
      //Checking if this chain is unique
      IsUnique = 1;
      for(j=0; j<clusSize; j++){
        if(chainPart == naList[j]){
          IsUnique = 0;
          break;
        }
      }
      if (IsUnique == 1){
        naList[clusSize++] = chainPart;
      }
      if (clusSize >= 5){
        return -1;
      }
     }
     //Moving on to the next bead in this chain
    }
    //Done with this chain, so let's move to the next chain, if it exists.
    list_it++;//Going one forward in naList to completely exhaust the tree.
    curID = naList[list_it];
  }
  //printf("%d\n", list_it);
  return clusSize;
}

int Clus_SecondLargestCluster(void){
    /*
    Calculates the complete networking for the system and returns naList which contains the chainIDs
    for the second largest cluster, if it exists.
    */
    //printf("\nStarting clus analysis\n");
    int curID, Cluster_length, currentLargest, i, j;
    int ClusNum = 0;
    int IsUnique = 1;
    for(i=0;i<=tot_chains;i++){
        naChainCheckList[i] = -1;
        naList[i]           = -1;
    }
    curID = 0;//Start with the 0th chain
    currentLargest = 0; j = 0;
    while(curID < tot_chains && IsUnique == 1){
        Cluster_length = Clus_ChainNetwork_ForTotal(curID);//This is the length of curID cluster
        if(Cluster_length > currentLargest) {
            currentLargest = Cluster_length;
            j = ClusNum;
        }

        for(i=0; i<Cluster_length; i++){//Recording the chains in this cluster
            naCluster[ClusNum][i+1] = naList[i];
            naList[i] = -1;
        }
        naCluster[ClusNum++][0] = Cluster_length;
        IsUnique = 0;//Assume not unique -- just got analyzed.
        while (curID<tot_chains && IsUnique == 0){//Finding the next chainID that hasn't been analyzed.
            curID++;
            IsUnique = 0;//Assume not unique.
            if (naChainCheckList[curID] == -1){
                IsUnique = 1;
            }
        }
    }
    //printf("ClusNum is %d\n", ClusNum);
    if (currentLargest == 1){//Just single chains
        //printf("Only single chains\t%d\n", ClusNum);
        curID = rand() % ClusNum;
    }
    else{//Find second largest
        curID = 0;
        for(i=0; i<ClusNum; i++){
            if(naCluster[i][0] < currentLargest && naCluster[i][0] >= naCluster[curID][0]){
                curID = i;
            }
        }
        //printf("%d:\t\t%d\t%d\n", ClusNum, currentLargest, naCluster[curID][0]);
        if (curID == j){//Reject if only one cluster
            //printf("\t\tTOO BIG\n");
            return -1;
        }
    }
    for(i=0; i<naCluster[curID][0]; i++){//Reupdate naList to have IDs for the second largest cluster
        naList[i] = naCluster[curID][i+1];
    }
    return naCluster[curID][0];
}

void Clus_TotalAnalysis(void){
    int curID, Cluster_length, currentLargest, i, j;
    int ClusNum = 0;
    int IsUnique = 1;
    for(i=0;i<=tot_chains;i++){
        naChainCheckList[i] = -1;
        naList[i]           = -1;
    }
    curID = 0;//Start with the 0th chain
    currentLargest = 0; j = 0;
    while(curID < tot_chains && IsUnique == 1){
        Cluster_length = Clus_ChainNetwork_ForTotal(curID);//This is the length of curID cluster
        if(Cluster_length > currentLargest) {
            currentLargest = Cluster_length;
            j = ClusNum;
        }

        for(i=0; i<Cluster_length; i++){//Recording the chains in this cluster
            naCluster[ClusNum][i+1] = naList[i];
            naList[i] = -1;
        }
        naCluster[ClusNum++][0] = Cluster_length;
        IsUnique = 0;//Assume not unique -- just got analyzed.
        while (curID<tot_chains && IsUnique == 0){//Finding the next chainID that hasn't been analyzed.
            curID++;
            IsUnique = 0;//Assume not unique.
            if (naChainCheckList[curID] == -1){
                IsUnique = 1;
            }
        }
    }
}
