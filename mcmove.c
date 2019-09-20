#include "global.h"
#include "mcmove.h"
#include "structure.h"
#include "cluster.h"
#include "energy.h"

int MC_Step(float fMCTemp) {
  int mode=0; //Which move to do
  int nAccept; //Used in MC steps
  float prob = (float)rand()/(float)RAND_MAX;
  int i;//Internal iterator, and int.
  for (i=0; i<MAX_MV; i++) {
    if (prob < fMCFreq[i]){//Pick this move!
      mode = i;
      break;
    }
  }
    nAccept = 0;

  //Performing the appropriate move.
  switch (mode) {
    //Translation
    case MV_TRANS:
      i = rand() % tot_chains; //Pick random chain
      nAccept = TransMCMove(i, fMCTemp);
      break;

    //Cluster translation
    case MV_CLSTR:
        nAccept = ClusMCMove(fMCTemp);
      break;

    //Cluster translation iff ClusSize <= 5
    case MV_SMCLSTR:
      i = rand() % tot_chains;//Pick a random chain
      nAccept = SmallClusMCMove(i, fMCTemp);
      break;

    //Change Rotational State
    case MV_FACEC:
      i = rand() % tot_beads;
      nAccept = RotMCMove(i, fMCTemp);
      break;

    //Local Move
    case MV_LOCAL:
      i = rand() % tot_beads;
      nAccept = LocalMCMove(i, fMCTemp);
      break;

    //Slithering Snake
    case MV_SNAKE:
      i = rand() % tot_chains;
      nAccept = SlitherMCMove(i, fMCTemp);
      break;

    //Double Pivot
    case MV_DBPVT:
      i = rand() % tot_beads;//Pick a random bead
      nAccept = DbPvtMCMove(i, fMCTemp);
      break;

    //Co-Local
    case MV_COLOCAL:
      i = rand() % tot_beads;//Pick a random bead
      nAccept = CoLocalMove(i, fMCTemp);
      break;

    //Shake Move
    case MV_SHAKE:
      i = rand() % tot_beads;
      nAccept = ShakeMove(i, fMCTemp);
      break;

    //Pivot
    case MV_PIVOT:
      i = rand() % tot_chains;
      nAccept = PivotMCMove(i, fMCTemp);
      break;

    //Branched Rotation
    case MV_BRROT:
      i = rand() % tot_chains;
      nAccept = BranchedRotMCMove(i, fMCTemp);
      break;

    default:
        nAccept = 0;
      break;
  }

  MCAccepMat[nAccept][mode]++;//Just adding which move got accepted/rejected.
  return mode*12+nAccept;
}

int MC_Step_Equil(float fMCTemp) {
  int mode=0; //Which move to do
  int nAccept; //Used in MC steps
  float prob = (float)rand()/(float)RAND_MAX;
  int i;//Internal iterator, and int.
  for (i=0; i<MAX_MV; i++) {
    if (prob < fMCFreq[i]){//Pick this move!
      mode = i;
      break;
    }
  }
  //Start by assuming failure.
  nAccept = 0;

  //Performing the appropriate move.
  switch (mode) {
    // translation of a random chain
    case MV_TRANS:
      i = rand() % tot_chains; //Pick random chain
      nAccept = TransMCMove_Equil(i, fMCTemp);
      break;

    // cluster translation: moves largest cluster to another spot.
    case MV_CLSTR:
        nAccept = ClusMCMove(fMCTemp);
      break;

    case MV_SMCLSTR:
      i = rand() % tot_chains;//Pick a random chain
      nAccept = SmallClusMCMove(i, fMCTemp);
      break;

    // face change
    case MV_FACEC:
      //  printf("FACE MOVE\n");
      nAccept = 0;
      break;

    // local moves
    case MV_LOCAL:
      i = rand() % tot_beads;
          nAccept = LocalMCMove_Equil(i, fMCTemp);
      break;

    // slithering snake
    case MV_SNAKE:
        i = rand() % tot_chains;
          nAccept = SlitherMCMove_Equil(i, fMCTemp);
        break;

    // double pivot
    case MV_DBPVT:
      i = rand() % tot_beads;//Pick a random bead
      nAccept = DbPvtMCMove(i, fMCTemp);
      break;

    // co-local
    case MV_COLOCAL:
        nAccept = 0;
      break;

    // shake
    case MV_SHAKE:
      i = rand() % tot_beads;
          nAccept = ShakeMove_Equil(i, fMCTemp);
      break;

    // pivot
    case MV_PIVOT:
    i = rand() % tot_chains;
          nAccept = PivotMCMove_Equil(i, fMCTemp);
    break;

    // branched rotate
    case MV_BRROT:
        i = rand() % tot_chains;
          nAccept = BranchedRotMCMove_Equil(i, fMCTemp);
        break;

    default:
        nAccept = 0;
      break;
  }

  MCAccepMat[nAccept][mode]++;//Just adding which move got accepted/rejected.
  return mode*12+nAccept;
}

int RotMCMove(int beadID, float MyTemp){
  //Performs a rotational MC-Move on beadID
  int bAccept; //Used in MC steps
  int resi, resj; //To track bead types.
  //Firstly -- make sure that bead i can even have rotational states.
  resi = bead_info[beadID][BEAD_TYPE];
  //printf("Beginning ROT\n");
  if(TypeCanRot[resi] == 0){//Skip beads that cannot rotate!
    bAccept = 0;
    return bAccept;
    }

  float MCProb, oldEn, newEn; //For Metropolis Hastings
  oldEn = 0.;
  newEn = 0.;
  int i; //General looping iterators
  int yTemp;
  int FWWeight;

  if (bead_info[beadID][BEAD_FACE] != -1){
    resj  = bead_info[bead_info[beadID][BEAD_FACE]][BEAD_TYPE];//This is who I'm currently bonded to
    oldEn = fEnergy[resi][resj][E_SC_SC];
    }
  ShuffleRotIndecies();
  FWWeight = CheckRotStatesOld(beadID, resi, MyTemp);
  NormalizeRotState(0, FWWeight);

  yTemp = PickRotState(FWWeight);

  if (yTemp != -1){//There is a bead here
    resj = bead_info[yTemp][BEAD_TYPE];
    newEn = fEnergy[resi][resj][E_SC_SC];
    }
  //See if we can accept this move
  MCProb = (float)rand()/(float)RAND_MAX;
  if (MCProb < expf((oldEn-newEn)/MyTemp)){//Accept this state
    if (bead_info[beadID][BEAD_FACE] != -1){//Break old bond
      bead_info[bead_info[beadID][BEAD_FACE]][BEAD_FACE] = -1;//Breaking bond with old partner
      }
    bead_info[beadID][BEAD_FACE] = yTemp;
    if(yTemp != -1){//New bond
      bead_info[yTemp][BEAD_FACE]  = beadID;
      }
    bAccept = 1;
    return bAccept;
    }
  else{
    bAccept = 0;
    return bAccept;
    }
}

int LocalMCMove(int beadID, float MyTemp){//Performs a local translation MC-move on beadID

  int bAccept = 0; //Used in MC steps
  float MCProb, oldEn, newEn; //For Metropolis Hastings
  oldEn = 0.;
  newEn = 0.;
  int i, j;//Loop iterators
  int resi, resj;
  int xTemp, yTemp, lRadUp, lRadLow;//Random numbers to store things
  int tmpR[POS_MAX], tmpR2[POS_MAX];//Vectors to stores coordinates.
  int FWWeight, BWWeight;//Used to perform orientational bias MC
  float FWRos, BWRos;//Forwards and backwards Rosenbluth Factors
  FWRos = 1; BWRos = 1;
  for(j=0; j<POS_MAX; j++){//Initializing the vectors to where this bead is.
    tmpR[j]  = bead_info[beadID][j];
  }
  //Initialize the radii for the search of next trial location
  //For now just +-2
  lRadLow = linker_len[beadID][0];
  //lRadLow = 2;
  lRadUp  = lRadLow*2 + 1;//2*2+1

  xTemp = 0; yTemp = 0;//Initialize these guys.
  while(yTemp == 0 && xTemp < nMCMaxTrials){//Attempt to find an empty lattice point.
    for(j=0; j<POS_MAX; j++){
      tmpR2[j] = (rand() % lRadUp) - lRadLow;//Generate number between -2 and 2
      tmpR2[j] = (tmpR[j] + tmpR2[j] + nBoxSize[j]) % nBoxSize[j];
    }
    yTemp = check_move_bead_to(tmpR2);
    if(yTemp == 1){//This means we found an empty lattice site. So let's check if the linkers are okay.
      yTemp = check_linker_constraint(beadID, tmpR2);
    }
    xTemp++;
  }
  if(xTemp == nMCMaxTrials || yTemp == 0){//No space -- reject.
    bAccept = 0;
    return bAccept;
  }
  else{//Have successfully found a good lattice spot. Let's perform the usual Metropolis-Hastings shenanigans.
    resi = bead_info[beadID][BEAD_TYPE];//I want to treat linker beads differently always because they have no rotational states
    oldEn = energy_cont_and_ovlp(beadID);
    if (TypeCanRot[resi] == 1){//Only non linkers can bond
      if (bead_info[beadID][BEAD_FACE] != -1){
        resj = bead_info[bead_info[beadID][BEAD_FACE]][BEAD_TYPE];//This is type of who I'm currently bonded to
        oldEn += fEnergy[resi][resj][E_SC_SC];
      }

      //ShuffleRotIndecies(); //No need to shuffle just to check.
      BWWeight = CheckRotStatesOld(beadID, resi, MyTemp);
      NormalizeRotState(0, BWWeight);
      BWRos = bolt_norm[0];

      move_bead_to(beadID, tmpR2);//Does not break bond

      ShuffleRotIndecies();//Need to shuffle because this also acts as selecting new partner
      FWWeight = CheckRotStatesNew(beadID, resi, MyTemp);
      NormalizeRotState(0, FWWeight);
      FWRos = bolt_norm[0];

      yTemp = PickRotState(FWWeight);

      if (yTemp != -1){//There is a bead at this position in the rot_trial, so let's add the energy.
        resj = bead_info[yTemp][BEAD_TYPE];
        newEn = fEnergy[resi][resj][E_SC_SC];
      }
  }
    else{//These beads have no rotational states.
      yTemp = -1;
      move_bead_to(beadID, tmpR2);
    }
    //Now let's calculate the energy of the new state. SC-SC energy is already done.
    newEn += energy_cont_and_ovlp(beadID);
    MCProb = (float)rand()/(float)RAND_MAX;
    if (MCProb < (FWRos/BWRos)*expf((oldEn-newEn)/MyTemp)){//Accept this state
      if (bead_info[beadID][BEAD_FACE] != -1){//Breaking old bond
        bead_info[bead_info[beadID][BEAD_FACE]][BEAD_FACE] = -1;
      }
      bead_info[beadID][BEAD_FACE] = yTemp;
      if(yTemp != -1){//Making new bond
        bead_info[yTemp][BEAD_FACE]  = beadID;
      }
      bAccept = 1;
      return bAccept;
    }
    else{
        //TODO: make separate undo function
      undo_move_bead_to(beadID);
      bAccept = 0;
      return bAccept;
    }
  }
}

int SlitherMCMove(int chainID, float MyTemp){//Performs a slither MC-move on chainID

  int firstB, lastB;//Track first and last+1 bead of chainID. Makes reading easier.
  int bAccept = 0; //Used in MC steps
  //Finding the bounds for looping over the molecule/chain
  firstB = chain_info[chainID][CHAIN_START];
  lastB  = firstB + chain_info[chainID][CHAIN_LENGTH];
  if(lastB - firstB == 1){//This means we have a monomer. Reject the move, because Local or Trans
    //moves should be the ones that move monomers.
    bAccept = 0;
    return bAccept;
  }
  else{
    if(TypeIsLinear[chain_info[chainID][CHAIN_TYPE]] != 1){//If chain is not linear. Reject move because slithering will not work!
      bAccept = 0;
      return bAccept;
    }
  }
  //This chain is suitable to perform a slithering-snake move on.

  float MCProb, oldEn, newEn; //For Metropolis Hastings
  oldEn = 0.; newEn = 0.;
  int i, j;//Loop iterators
  int resi, resj;
  int xTemp, yTemp, lRadUp, lRadLow;//Random numbers to store things
  int tmpR[POS_MAX], tmpR2[POS_MAX], tmpR3[POS_MAX];//Vectors to store positions.
  int FWWeight, BWWeight;//Used to perform orientational bias MC
  double FSum, BSum;

  MCProb = (float)rand()/(float)RAND_MAX;//To decide if we slither forwards or backwards
  if(MCProb < 0.5){//Forwards slither, so lastB-1 (last bead) is anchor
  lRadUp  = 2 * linker_len[lastB - 1][0] + 1;//lastB-1 will be replaced by lastB-2
  lRadLow = linker_len[lastB - 1][0];
  yTemp = 0; xTemp = 0;//Using to track trials for placing beads
  while (xTemp < nMCMaxTrials && yTemp == 0) {
    for(j=0; j<POS_MAX; j++){
        tmpR[j] = (rand() % lRadUp) - lRadLow;
        tmpR[j] = (bead_info[lastB-1][j] + tmpR[j] + nBoxSize[j]) % nBoxSize[j];
      }
    yTemp = check_move_bead_to(tmpR);// 0: there is no space, 1: there is space
    xTemp++;
      }
    }
  else{//Backwards slither, so firstB is anchor
      lRadUp  = (int)2 * linker_len[firstB][0] + 1;//firstB will be replaced by firstB+1
      lRadLow = (int)linker_len[firstB][0];
      yTemp = 0; xTemp = 0;//Using to track trials for placing beads
      while (xTemp < nMCMaxTrials && yTemp == 0) {
        for(j=0; j<POS_MAX; j++){
          tmpR[j] = (rand() % lRadUp) - lRadLow;
          tmpR[j] = (bead_info[firstB][j] + tmpR[j] + nBoxSize[j]) % nBoxSize[j];
        }
      yTemp = check_move_bead_to(tmpR);// 0: there is no space, 1: there is space
      xTemp++;
        }
    }

    if (yTemp == 0 || xTemp == nMCMaxTrials){//Couldn't find a spot, so reject the damn move
      bAccept = 0;
      return bAccept;
    }

   //Let's remember where this chain exists.
  for(i=firstB; i<lastB; i++){
    for(j=0; j<BEADINFO_MAX; j++){
      old_bead[i][j] = bead_info[i][j];
    }
  }
    //Recording rotational degeneracy has no direction so doesn't matter if forwards or backwards
    yTemp = 0;
    for (i=firstB; i<lastB; i++){
      resi = bead_info[i][BEAD_TYPE];
      oldEn += energy_cont_and_ovlp(i);
    if (TypeCanRot[resi] == 0){//Skip beads that don't interact
      continue;
    }
    if (bead_info[i][BEAD_FACE] != -1){//I am bonded to something
      resj   = bead_info[bead_info[i][BEAD_FACE]][BEAD_TYPE];//Type of bead I'm bonded to
      oldEn += fEnergy[resi][resj][E_SC_SC];
      }

    //ShuffleRotIndecies();
    BWWeight = CheckRotStatesOld(i, resi, MyTemp);
    NormalizeRotState(yTemp, BWWeight);
    yTemp++;
    }
    //Done with checking states in the old location. Take the sum.
    BSum = 1.;
    for (i = 0; i < yTemp; i++){
        BSum += log10(bolt_norm[i]);
    }

  if(MCProb < 0.5){//Slithering the chain forwards in ID-space
  for (i=firstB; i<lastB-1; i++){
    for(j=0;j<POS_MAX;j++){
      tmpR2[j] = bead_info[i][j];
      bead_info[i][j] = old_bead[i+1][j];//Hopping over by one bead
      tmpR3[j] = bead_info[i][j];
    }
    if(i==firstB){//Only the firstB's location is empty
      naTotLattice[LtIndV(tmpR2)] = -1;
    }
      naTotLattice[LtIndV(tmpR3)] = i;
  }
  //Moving the last bead, and it has to be done independently because lastB-1 -> tmpR
  i = lastB-1;
  for(j=0;j<POS_MAX;j++){
    bead_info[i][j] = tmpR[j];
  }
      naTotLattice[LtIndV(tmpR)]   = i;
  }
  else{//Slithering backwards in ID-space
    for (i=firstB+1; i<lastB; i++){
      for(j=0; j<POS_MAX; j++){
        tmpR2[j] = bead_info[i][j];
        bead_info[i][j] = old_bead[i-1][j];//Hopping back by one bead
        tmpR3[j] = bead_info[i][j];
      }
      if(i==lastB-1){//Only the lastB-1's location is empty
        naTotLattice[LtIndV(tmpR2)] = -1;
      }
        naTotLattice[LtIndV(tmpR3)] = i;
    }
    //Moving the first bead, and it has to be done independently because firstB -> tmpR
    i = firstB;
    for(j=0; j<POS_MAX; j++){
      bead_info[i][j] = tmpR[j];
    }
      naTotLattice[LtIndV(tmpR)]    = i;
  }

  //Have slithered the chain whichever way, so let's check rotational states in the new location
  for(i=firstB; i<lastB; i++){//Break the old bonds first!
    if(bead_info[i][BEAD_FACE] != -1){//Need to break this bond
      bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE] = -1;
      bead_info[i][BEAD_FACE] = -1;
    }
  }

    yTemp = 0;
     for (i=firstB; i<lastB; i++){//Counting states in the new location
         resi = bead_info[i][BEAD_TYPE];
         newEn += energy_cont_and_ovlp(i);
         if (TypeCanRot[resi] == 0){//Skip non-bonders
         continue;
         }
        ShuffleRotIndecies();
        FWWeight = CheckRotStatesNew(i, resi, MyTemp);
        NormalizeRotState(yTemp, FWWeight);
        //Note that the bonds need to be formed in this loop so that we don't overcount!
        if(bead_info[i][BEAD_FACE] == -1){//Make sure this bead is unbonded!
        //Let's assign a rotational state to this bead
            xTemp = PickRotState(FWWeight);
            if(xTemp != -1){//An appropriate partner has been selected. Form the bonds and add the energy
                resj = bead_info[xTemp][BEAD_TYPE];
                bead_info[i][BEAD_FACE] = xTemp;
                bead_info[xTemp][BEAD_FACE] = i;
                newEn += fEnergy[resi][resj][E_SC_SC];
            }
          }
        yTemp++;//This keeps track of which residue*/
    }
    FSum = 1.;
    for (i=0; i<yTemp; i++){
          FSum += log10(bolt_norm[i]);
      }
     //Doing the Metropolis-Hastings thing
     MCProb = (float)rand()/(float)RAND_MAX;
     if ( MCProb < (FSum/BSum)*expf((oldEn-newEn)/MyTemp) ){//Accept. Bonds have been handled before!
      bAccept = 1;
      return bAccept;
    }
     else{//TODO: Make a function that resets a single chain, breaking new bonds and restoring old ones.
       for(i=firstB; i<lastB; i++){//Resetting the lattice
         if(bead_info[i][BEAD_FACE] != -1){//Need to break the newly proposed bond
           bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE] = -1;
         }
         for(j=0; j<POS_MAX; j++){//This is where I am
          tmpR2[j] = bead_info[i][j];
        }
           naTotLattice[LtIndV(tmpR2)]  = -1;
         }
       for(i=firstB; i<lastB; i++){
         for(j=0;j<BEADINFO_MAX; j++){
            bead_info[i][j] = old_bead[i][j];//Restoring
            if(j<POS_MAX){//This is where I should be
              tmpR2[j] = bead_info[i][j];
            }
         }
         if(bead_info[i][BEAD_FACE] != -1){//I was bonded so restore the bond
           bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE] = i;}
           naTotLattice[LtIndV(tmpR2)]  = i;
         }
         bAccept = 0;
         return bAccept;
        }
}

int TransMCMove(int chainID, float MyTemp){//Performs a translation move with orientational bias
  int bAccept = 0; //Used in MC steps
  float MCProb, oldEn, newEn; //For Metropolis Hastings
  oldEn = 0.; newEn = 0.;
  int i, j;//Loop iterators
  int resi, resj, firstB, lastB;
  int xTemp, yTemp, lRadUp, lRadLow;//Random numbers to store things
  int tmpR[POS_MAX];//Vectors to store coordinates.
  int FWWeight, BWWeight;//Used to perform orientational bias MC
  float FSum, BSum;//Overall Rosenbluth sums


  //Finding the bounds for looping over the molecule/chain
 firstB = chain_info[chainID][CHAIN_START];
 lastB  = firstB + chain_info[chainID][CHAIN_LENGTH];
  //Radii for translation moves. All moves are L/4 radius
  lRadLow = nBoxSize[2]/2;
  lRadUp  = 2*lRadLow + 1;

  xTemp = 0; yTemp = 0;
  while(xTemp < nMCMaxTrials && yTemp == 0){
      for (j=0; j<POS_MAX; j++) {
        tmpR[j] = (rand() % lRadUp) - lRadLow;//Random vector to move all beads within r=L/4
      }
      yTemp = check_disp_chain(chainID, tmpR);//yTemp=0 means clash
      xTemp++;
    }
  if (yTemp == 0 || xTemp == nMCMaxTrials){//We have failed to find a good spot for this chain.
    bAccept = 0;
    return bAccept;
  }
  //We now have a chain which when moved does not overlap.

  yTemp = 0;
  for (i = firstB; i < lastB; i++){//Rosenbluth in old location.
    resi = bead_info[i][BEAD_TYPE];
    oldEn += energy_cont_and_ovlp(i);
    if (TypeCanRot[resi] != 1){//Skip beads that cannot bond.
      continue;
      }

    if (bead_info[i][BEAD_FACE] != -1){//I am bonded to something
      resj = bead_info[bead_info[i][BEAD_FACE]][BEAD_TYPE];//Type of bead I am bonded to
      oldEn += fEnergy[resi][resj][E_SC_SC];//Adding the energy.
      }

    //ShuffleRotIndecies();
    BWWeight = CheckRotStatesOld(i, resi, MyTemp);
    NormalizeRotState(yTemp, BWWeight);
    yTemp++;
    }

  BSum = 1.;
  for (i = 0; i < yTemp; i++){
    BSum += log10(bolt_norm[i]);
    }

  trans_disp_chain(chainID, tmpR);//Moved the chain, broke bonds, and remembered stuff

  yTemp = 0;
  for (i = firstB; i < lastB; i++){//Rosenbluth in new location
    resi = bead_info[i][BEAD_TYPE];
    newEn += energy_cont_and_ovlp(i);
    if (TypeCanRot[resi] != 1){//Because linkers don't have rotational states
      continue;
      }

    ShuffleRotIndecies();
    FWWeight = CheckRotStatesNew(i, resi, MyTemp);
    NormalizeRotState(yTemp, FWWeight);
    //Note that the bonds need to be formed in this loop so that we don't overcount!
    if(bead_info[i][BEAD_FACE] == -1){//Make sure this bead is unbonded!
    //Let's assign a rotational state to this bead
      xTemp = PickRotState(FWWeight);
      if(xTemp != -1){//An appropriate partner has been selected. Form the bonds and add the energy
        resj = bead_info[xTemp][BEAD_TYPE];
        bead_info[i][BEAD_FACE] = xTemp;
        bead_info[xTemp][BEAD_FACE] = i;
        newEn += fEnergy[resi][resj][E_SC_SC];
        }
      }
    yTemp++;//This keeps track of which residue*/
  }

  FSum = 1.;
  for (i = 0; i < yTemp; i++){
    FSum += log10(bolt_norm[i]);
    }
   MCProb = (float)rand()/(float)RAND_MAX;
   if ( MCProb < (FSum/BSum)*expf((oldEn-newEn)/MyTemp) ){//Accept the move. Remember that the bonds were assigned above!
    bAccept = 1;
    return bAccept;
   }
   else{
       trans_restore_chain(chainID);
       bAccept = 0;
       return bAccept;
   }

}

int ClusMCMove(float MyTemp){
  //Attempts to move the second largest cluster

  int bAccept = 0; //Used in MC steps, assume that move fails initially.
  int ClusSize, i, j;//Loop iterators
  int yTemp;
  int nTemp[POS_MAX];
  int lRadLow, lRadUp;//Radii bounds.
  float oldEn, newEn, MCProb;
  oldEn = 0.0; newEn = 0.0;

  ClusSize = clus_network_analysis(naList, naCluster);//Second largest cluster;

  if(ClusSize != -1){
    //Radii for translation moves. All moves are L/2 radius
    lRadLow = nBoxSize[2]/2;
    lRadUp  = 2*lRadLow + 1;
      for(j=0;j<POS_MAX;j++){
       nTemp[j] = (rand() % lRadUp) - lRadLow;   //Random vector to displace the cluster
      }
      for(i=0; i<ClusSize; i++){
       yTemp = check_disp_chain(naList[i], nTemp);//Checking for steric clash
       if (yTemp == 0){
         bAccept = 0;
         //printf("End CLUS - No space\n");
         return bAccept;
       }
      }
      //This means that there is no steric clash when the cluster is moved.
      for (i=0; i<ClusSize; i++){
       oldEn += energy_chain(naList[i]);   //Old energy
      }
      for(i=0; i<ClusSize; i++){
       disp_chain(naList[i], nTemp);//Moving the cluster properly
      }
      for (i=0; i<ClusSize; i++){
       newEn += energy_chain(naList[i]);   //New energy
      }
      MCProb = (float)rand()/(float)RAND_MAX;
      if ( MCProb < expf((oldEn-newEn)/MyTemp)){
       bAccept = 1;//Accept the move
       //printf("End CLUS - Yes\n");
      }
      else{
       bAccept = 0;   //Reject the move so I have to restore cluster back
       for(i=0; i<ClusSize; i++){
       restore_chain(naList[i]);//Placing  the cluster back properly
        }
        //printf("End CLUS - Failed.\n");
      }
  }
  //printf("Ending CLUS w/ %d clus\n", ClusSize);
  return bAccept;
}

int SmallClusMCMove(int chainID, float MyTemp){
  //Performs a cluster move where a given chain and it's cluster are moved. No new 'bonds' are made so the move is reversible....

  int bAccept = 0; //Used in MC steps, assume that move fails initially.
  int ClusSize, i, j;//Loop iterators
  int yTemp;
  int nTemp[POS_MAX];
  int lRadLow, lRadUp;//Radii bounds.
  float oldEn, newEn, MCProb;
  oldEn = 0.0; newEn = 0.0;
  //printf("Beginning CLUS\n");
  ClusSize = chain_network_small(chainID, naList);//Looking at everything that is connected to chainID
  //Remember that naList[] contains the chainID's of the network chainID is part of from 0 - ClusSize-1.
  //printf("Done with network\t %d\n", ClusSize);
  if(ClusSize >= 1){
    //Radii for translation moves. All moves are L/4 radius
    //I guess moving single chains around as well is not a bad idea
    lRadLow = nBoxSize[2]/2;
    lRadUp  = 2*lRadLow + 1;
      for(j=0;j<POS_MAX;j++){
       nTemp[j] = (rand() % lRadUp) - lRadLow;   //Random vector to displace the cluster
      }
      for(i=0; i<ClusSize; i++){
        //printf("%d\n", naList[i]);
       yTemp = check_disp_chain(naList[i], nTemp);//Checking for steric clash
       if (yTemp == 0){
         bAccept = 0;
         //printf("End CLUS - No space\n");
         return bAccept;
       }
      }
      //This means that there is no steric clash when the cluster is moved.
      for (i=0; i<ClusSize; i++){
       oldEn += energy_chain(naList[i]);   //Old energy
      }
      for(i=0; i<ClusSize; i++){
       disp_chain(naList[i], nTemp);//Moving the cluster properly
      }
      for (i=0; i<ClusSize; i++){
       newEn += energy_chain(naList[i]);   //New energy
      }
      MCProb = (float)rand()/(float)RAND_MAX;
      if ( MCProb < expf((oldEn-newEn)/MyTemp)){
       bAccept = 1;//Accept the move
       //printf("End CLUS - Yes\n");
      }
      else{
       bAccept = 0;   //Reject the move so I have to restore cluster back
       for(i=0; i<ClusSize; i++){
       restore_chain(naList[i]);//Placing  the cluster back properly
        }
        //printf("End CLUS - Failed.\n");
      }
  }
  //printf("Ending CLUS w/ %d clus\n", ClusSize);
  return bAccept;
}

int DbPvtMCMove(int beadID, float MyTemp){//Performs a double-pivot move.
  /* Molecule MUST be LINEAR
  The move requires selecting a random bead, which is beadID. Then, we'll search the lattice in +-2 sites around beadID.
  Let i be the position of beadID along it's chain. Let i' denote same position along another chain of the same type. We want
  dist(i,i'+1) < linker_len[i] && dist (i',i+1) < linker_len[i']. We'll count however many candidates there are and select one randomly.
  Then we just swap beads from i to then end of the chain. Do a Metropolis thing, and decide.
  In other words, i'+1 becomes i+1, i'+2 becomes i+2 until N, and i+1 become i'+1 and so on.
   */
    int bAccept = 0; //Move acceptance and such LEL
    int PChainID = bead_info[beadID][BEAD_CHAINID];//The proposed chainID
    int PType  = chain_info[PChainID][CHAIN_TYPE];//Type of chain.
    int PStart = chain_info[PChainID][CHAIN_START];//Start of this chain.
    int PEnd   = PStart + chain_info[PChainID][CHAIN_LENGTH];//LastBead+1 for this chain.
    if(TypeIsLinear[PType] == 0){
      //Reject the move because the chain is not linear.
      return bAccept;
    }
    //Now make sure that the proposed bead is neither the start or end of a chain.
    if(PStart == beadID || PEnd - 1 == beadID){
      return bAccept;
    }

    int thisbead, otherbead;// The beads that will be swapped!
    int thischain, mychain, thischaintype, thischainstart;//Storing chainIDs for the eventual swapping
    int i, j, k;//Loop iterators
    int SrchLen = 2;//The length of the box we'll be searching.
    int nTemp[POS_MAX] = {0};//Storing the lattice position
    int nListLen = 0;//Tracks if there are any bead candidates
    int nListLen_back = 0;//Tracks the number of candidates for the reverse move.
    int candList[MAX_ROTSTATES] = {0};
    float MCProb;

    for (i=-SrchLen; i<=SrchLen;i++){
        nTemp[POS_X] = (bead_info[beadID][POS_X] + i + nBoxSize[POS_X] ) % nBoxSize[POS_X];
    for (j=-SrchLen; j<=SrchLen;j++){
        nTemp[POS_Y] = (bead_info[beadID][POS_Y] + j + nBoxSize[POS_Y] ) % nBoxSize[POS_Y];
    for (k=-SrchLen; k<=SrchLen;k++){
        nTemp[POS_Z] = (bead_info[beadID][POS_Z] + k + nBoxSize[POS_Z] ) % nBoxSize[POS_Z];
        thisbead = naTotLattice[LtIndV(nTemp)];//Seeing what is at that location
        //thisbead is i'+1 from the above explanation
        if(thisbead != beadID && thisbead != -1 && thisbead != beadID+1){//Different bead from me
        if (PType ==
            chain_info[bead_info[thisbead][BEAD_CHAINID]][CHAIN_TYPE]
           && PChainID != bead_info[thisbead][BEAD_CHAINID]){//Making sure that they are the same chain types, and that the chains are different
        if(beadID - PStart + 1 ==
           thisbead - chain_info[bead_info[thisbead][BEAD_CHAINID]][CHAIN_START]){//Remember that we are searching for i'+1, not i'
            //NOW we can see if there is room to make bridges
            //Remember that linker_len[beadID][0] is the linker dist for beadID-1,
            //and linker_len[beadID][1] is the linker dist for beadID+1
               if(dist(beadID, thisbead) < 1.74 * linker_len[beadID][1] &&
                   dist(beadID+1, thisbead-1 ) < 1.74 * linker_len[thisbead - 1][1] ){//The linker lengths are correct to change connectivity
                       candList[nListLen] = thisbead;//Storing which bead it is
                       nListLen++;//Onto the next bead
                       if(nListLen == MAX_ROTSTATES){goto FoundMax;}//MAX_ROTSTATES is the size of candList[]
                                }
                            }
                    }
                }
            }

        }
    }
    //Should have found all the candidates
    FoundMax:

    if (nListLen == 0){
      bAccept = 0;
      return bAccept;
    }

    thisbead       = rand() % nListLen;//Randomly select a candidate
    thisbead       = candList[thisbead];//Pick the ID of the candidate
    thischain      = bead_info[thisbead][BEAD_CHAINID];//The chain ID of the candidate
    mychain        = bead_info[beadID][BEAD_CHAINID];
    thischaintype  = chain_info[thischain][CHAIN_TYPE];
    thischainstart = chain_info[thischain][CHAIN_START];

    //For detailed balance, we need to count how many candidates thisbead-1 has!

    for (i=-SrchLen; i<=SrchLen;i++){
        nTemp[POS_X] = (bead_info[thisbead-1][POS_X] + i + nBoxSize[POS_X] ) % nBoxSize[POS_X];
    for (j=-SrchLen; j<=SrchLen;j++){
        nTemp[POS_Y] = (bead_info[thisbead-1][POS_Y] + j + nBoxSize[POS_Y] ) % nBoxSize[POS_Y];
    for (k=-SrchLen; k<=SrchLen;k++){
        nTemp[POS_Z] = (bead_info[thisbead-1][POS_Z] + k + nBoxSize[POS_Z] ) % nBoxSize[POS_Z];
        otherbead = naTotLattice[LtIndV(nTemp)];//Seeing what is at that location
        //otherbead is i'+1 from the above explanation
        if(otherbead != thisbead-1 && otherbead != -1 && otherbead != thisbead){//Different bead from me
        if (thischaintype ==
            chain_info[bead_info[otherbead][BEAD_CHAINID]][CHAIN_TYPE]
           && thischain != bead_info[otherbead][BEAD_CHAINID]){//Making sure that they are the same chain types, and that the chains are different
        if(thisbead - PStart ==
           otherbead - chain_info[bead_info[thisbead][BEAD_CHAINID]][CHAIN_START]){//Remember that we are searching for i'+1, not i'
            //NOW we can see if there is room to make bridges
            //Remember that linker_len[thisbead][0] is the linker dist for thisbead-1,
            //and linker_len[thisbead][1] is the linker dist for thisbead+1
               if(dist(thisbead-1, otherbead) < 1.74 * linker_len[thisbead - 1][1] &&
                   dist(thisbead, otherbead-1 ) < 1.74 * linker_len[otherbead - 1][1] ){//The linker lengths are correct to change connectivity
                       //candList[nListLen] = thisbead;//Storing which bead it is
                       nListLen_back++;//Onto the next guy
                       if(nListLen_back == MAX_ROTSTATES){goto FoundMax_back;}//MAX_ROTSTATES is the size of candList[]
                                }
                            }
                    }
                }
            }

        }
    }
    //Should have found all the candidates
    FoundMax_back:
    MCProb = (float)rand()/(float)RAND_MAX;
    //*/
    if(MCProb < (float)(nListLen+1.)/(float)(nListLen_back+1.)){
    //We can perform the swap
    j=1;//Tracks the beads
        for ( i = beadID+1; i < PEnd; i++){//Swapping from beadID+1 onwards
            swap_beads(beadID+j, thisbead+(j-1));//This is pretty dumb, but easier to read UGH
            j++;
        }
        bAccept  = 1;
        return bAccept;
      }
    else{
      bAccept = 0;
      return bAccept;
    }//*/
}

int CoLocalMove(int beadID, float MyTemp){
  /*
  Translate a bead and its partner in tandem. If no partner, reject move.
  */

  int bAccept = 0; //Used in MC steps

  if (bead_info[beadID][BEAD_FACE] == -1){
    bAccept = 0;
    //printf("Bead has no partner!\n");
    return bAccept;
  }
  int tmpR[POS_MAX];//Random translation vector.
  int tmpR1[POS_MAX], tmpR2[POS_MAX];
  int beadPart = bead_info[beadID][BEAD_FACE];
  //printf("Partner is (%d)?!\n", beadPart);
  float MCProb, oldEn, newEn; //For Metropolis Hastings.
  oldEn = 0.;
  newEn = 0.;
  int i, j;//Loop iterators
  int xTemp, yTemp;
  int lRadLow, lRadUp;
  lRadLow = 2;
  lRadUp  = lRadLow*2 + 1;//2*2+1

  xTemp = 0; yTemp = 0;//Initialize these guys.
  while(yTemp == 0 && xTemp < nMCMaxTrials){
    for(j=0; j<POS_MAX; j++){
      tmpR[j]  = (rand() % lRadUp) - lRadLow;
      tmpR1[j] = (bead_info[beadID][j]   + tmpR[j] + nBoxSize[j]) % nBoxSize[j];
      tmpR2[j] = (bead_info[beadPart][j] + tmpR[j] + nBoxSize[j]) % nBoxSize[j];
    }

    yTemp = check_move_bead_to(tmpR1)*check_move_bead_to(tmpR2);
    if(yTemp == 1){//This means we found an empty lattice site. So let's check if the linkers are okay.
      yTemp = check_linker_constraint(beadID, tmpR1)*check_linker_constraint(beadPart, tmpR2);
    }
    xTemp++;
  }
  if(xTemp == nMCMaxTrials || yTemp == 0){//This means that we have failed to find an appropriate spot for this bead to be moved to. Therefore, the move is rejected!
    bAccept = 0;
    return bAccept;
  }
  oldEn = energy_cont_and_ovlp(beadID)+energy_cont_and_ovlp(beadPart);
  move_bead_to(beadID, tmpR1); move_bead_to(beadPart, tmpR2);
  newEn = energy_cont_and_ovlp(beadID)+energy_cont_and_ovlp(beadPart);
  MCProb = (float)rand()/(float)RAND_MAX;
  if (MCProb < expf((oldEn-newEn)/MyTemp)){
    bAccept = 1;
    //printf("Accepted!\n");
    return bAccept;
  }
  else{
    undo_move_bead_to(beadID);
    undo_move_bead_to(beadPart);
    bAccept = 0;
    //printf("Rejected move; undoing!\n");
    return bAccept;
  }
}

int ShakeMove(int beadID, float MyTemp){
    int topIt; //Iterator for topo_info
    int i,j; //Loop iterators
    int curID; //current bead being looked at
    int tmpR[MAX_VALENCY][POS_MAX]; //Storing temporary locations
    int  bAccept;
    int lRadLow, lRadUp;
    lRadLow = 2;
    lRadUp  = lRadLow*2 + 1;//2*2+1
    int xTemp = 0;
    int yTemp = 0;

    curID = beadID; topIt = 0;
    while (curID != -1) {
        for (j=0; j<BEADINFO_MAX; j++) {
            old_bead[curID][j] = bead_info[curID][j];//Remembering
            if(j<POS_MAX){
                tmpR[topIt][j] = bead_info[curID][j];
            }
        }
        curID = topo_info[beadID][topIt++];
    }

    while(yTemp == 0 && xTemp < nMCMaxTrials){
        curID = beadID; topIt = 0;
        while (curID != -1){
            naTotLattice[LtIndV(bead_info[curID])] = -1;
            curID = topo_info[beadID][topIt++];
        }
        curID = beadID; topIt = 0;
        while (curID != -1){
            yTemp = 1;
            for(j = 0; j < POS_MAX; j++){
                tmpR[topIt][j] = (rand() % lRadUp) - lRadLow;
                tmpR[topIt][j] = (bead_info[curID][j] + tmpR[topIt][j] + nBoxSize[j]) % nBoxSize[j];
            }
            if (naTotLattice[LtIndV(tmpR[topIt])] != -1 ){
                yTemp = 0;
                break;
            }
            naTotLattice[LtIndV(tmpR[topIt])] = curID;
            curID = topo_info[beadID][topIt++];
        }
        if (yTemp == 1){//No steric clash so check for topology constraint
            yTemp = ShakeConstraint(beadID, tmpR);
        }
        for(i = 0; i < topIt; i++){
            naTotLattice[LtIndV(tmpR[i])] = -1;
        }
        xTemp++;
    }

    if(xTemp == nMCMaxTrials || yTemp == 0){//Linker or steric clash didn't work out
        curID = beadID; topIt = 0;
        while (curID != -1) {
            naTotLattice[LtIndV(bead_info[curID])] = curID;
            curID = topo_info[beadID][topIt++];
        }
        //printf("No space!\n");
        bAccept = 0;
        return bAccept;
    }

    int FWWeight, BWWeight;//Used to perform orientational bias MC
    float FSum, BSum;
    int resi, resj;
    //int tmpR2[POS_MAX], tmpR3[POS_MAX];
    float oldEn = 0.;
    float newEn = 0.;
    float MCProb;


    curID = beadID; topIt = 0; yTemp = 0;
    while (curID != -1){
        resi = bead_info[curID][BEAD_TYPE];
        oldEn += energy_cont_and_ovlp(curID);
    if (TypeCanRot[resi] != 1){//Skip non-interactors
        curID = topo_info[beadID][topIt++];
        continue;
    }
    if (bead_info[curID][BEAD_FACE] != -1){//I am bonded to something
        resj = bead_info[bead_info[curID][BEAD_FACE]][BEAD_TYPE];//Type of bead I'm bonded to
        oldEn += fEnergy[resi][resj][E_SC_SC];
    }

    ShuffleRotIndecies();
    BWWeight = CheckRotStatesOld(curID, resi, MyTemp);
    NormalizeRotState(yTemp, BWWeight);
    curID = topo_info[beadID][topIt++];
    yTemp++;
    }

    BSum = 1.;
    for(i=0; i<yTemp; i++){
        BSum += log10(bolt_norm[i]);
    }

    curID = beadID; topIt = 0;
    while (curID != -1){
        move_bead_to_shake(curID, tmpR[topIt]);
        curID = topo_info[beadID][topIt++];
    }

    curID = beadID; topIt = 0;
    while (curID != -1){//Breaking old bonds
        i = bead_info[curID][BEAD_FACE];
        if(i != -1){
            bead_info[i][BEAD_FACE] = -1;
            bead_info[curID][BEAD_FACE] = -1;
            }
        curID = topo_info[beadID][topIt++];
    }

    curID = beadID; topIt = 0; yTemp = 0;
    while (curID != -1){
        resi   = bead_info[curID][BEAD_TYPE];
        newEn += energy_cont_and_ovlp(curID);
        if (TypeCanRot[resi] != 1){//Skip non-interactors
            curID = topo_info[beadID][topIt++];
            continue;
        }
        ShuffleRotIndecies();
        FWWeight = CheckRotStatesNew(curID, resi, MyTemp);
        NormalizeRotState(yTemp, FWWeight);
        if(bead_info[curID][BEAD_FACE] == -1){//Make sure this bead is unbonded!
            xTemp = PickRotState(FWWeight);
                if(xTemp != -1){//An appropriate partner has been selected. Form the bonds and add the energy
                    resj = bead_info[xTemp][BEAD_TYPE];
                    bead_info[curID][BEAD_FACE] = xTemp;
                    bead_info[xTemp][BEAD_FACE] = curID;
                    newEn += fEnergy[resi][resj][E_SC_SC];
                }
        }
  curID = topo_info[beadID][topIt++];
  yTemp++;
  }

  FSum = 1.;
  for (i = 0; i < yTemp; i++){
      FSum += log10(bolt_norm[i]);
  }

  MCProb = (float)rand()/(float)RAND_MAX;
  if ( MCProb < (FSum/BSum)*expf((oldEn-newEn)/MyTemp) ){//Accept the move. Remember that the bonds were assigned above
      bAccept = 1;
      return bAccept;
  }
  else{
      curID = beadID; topIt = 0;
      while (curID != -1){
          i = bead_info[curID][BEAD_FACE];
          if(i != -1){
              bead_info[i][BEAD_FACE] = -1;
              bead_info[curID][BEAD_FACE] = -1;
          }
          naTotLattice[LtIndV(bead_info[curID])] = -1;
          for(j = 0; j < BEADINFO_MAX; j++){
              bead_info[curID][j] = old_bead[curID][j];
          }
          curID = topo_info[beadID][topIt++];
      }
      curID = beadID; topIt = 0;
      while (curID != -1){
          naTotLattice[LtIndV(bead_info[curID])] = curID;
          i = bead_info[curID][BEAD_FACE];
          if(i != -1){
              bead_info[i][BEAD_FACE] = curID;
          }
          curID = topo_info[beadID][topIt++];
      }
        bAccept = 0;
        return bAccept;
    }
}

int PivotMCMove(int chainID, float MyTemp){
  //Performs a pivot move on chainID
  /*
  Randomly pick a bead in this chain (anchorBead), and perform a symmetry operation on the beads after the anchorBead. Note that if anchorBead is the second-to-last, or last, bead, the move is rejected. Local and slither moves would be better for those.
  The protein must be linear and the move is rejected if the chain is branched.
  It is assumed that a linear chain has been passed!
  */
  int bAccept = 0;
  //Check if the chain is longer than 3 or if it is linear
  if(chain_info[chainID][CHAIN_LENGTH] <=3 || TypeIsLinear[chain_info[chainID][CHAIN_TYPE]] != 1){
    bAccept = 0;
    return bAccept;
  }
  int firstB, lastB;
  //Get the beadID's for the first and last bead
  firstB = chain_info[chainID][CHAIN_START];
  lastB  = firstB + chain_info[chainID][CHAIN_LENGTH];
  //This chain is long enough
  int chainLength = chain_info[chainID][CHAIN_LENGTH];
  //Randomly select a bead
  int anchorBead  = chainLength-2;
  anchorBead      = 2 + (rand() % anchorBead);

  int PivotDir = 0; //-1 Means backwards, +1 means forwards. Always Pivot the smaller portion
  PivotDir   = anchorBead > chainLength/2 ? 1 : -1;
  anchorBead = firstB + anchorBead;
  //printf("Bead: %d\n", anchorBead);
  //Checking if bead is first or last. If so, reject move.
  if(anchorBead == firstB || anchorBead > lastB-2){
    bAccept = 0;
    //printf("P Bad Bead\n");
    return bAccept;
  }

  int PivotM;
  PivotM = rand() % 10;
  //printf("Move: %d\n", PivotM);
  if(PivotM == 0){
    bAccept = 1;
    return bAccept;
  }

  int i,j;
  int xTemp, yTemp;
  int anchorPos[POS_MAX];
  int tmpList[MAX_CHAINLEN];
  int listLen = 0;
  for(j=0; j<POS_MAX; j++){
    anchorPos[j] = bead_info[anchorBead][j];
  }
  if(PivotDir == 1){
      for (i = anchorBead + 1; i < lastB; i++) {
          tmpList[listLen++] = i;
      }
  }
  else{
      for (i = firstB; i < anchorBead; i++) {
          tmpList[listLen++] = i;
      }
  }

  xTemp = 0; yTemp = 0;
  while (xTemp < nMCMaxTrials && yTemp == 0) {
      //for (i = anchorBead + 1; i < lastB; i++) {
      for (j = 0; j < listLen; j++){
          i  = tmpList[j];
          RotOperation(PivotM, i, anchorPos);
          yTemp = check_move_bead_to(naTempR);
          if (yTemp == 0) {
              xTemp++;
              break;
          }
      }
      xTemp++;
  }

  if(xTemp == nMCMaxTrials || yTemp == 0){
    bAccept = 0;
    //printf("P Clashout\n");
    return bAccept;
  }

  int FWWeight, BWWeight;//Used to perform orientational bias MC
  float FSum, BSum;
  int resi, resj;
  float oldEn = 0.;
  float newEn = 0.;
  float MCProb;

  yTemp = 0;
    //for (i = anchorBead + 1; i < lastB; i++) {
    for (j = 0; j < listLen; j++){
      i  = tmpList[j];
      resi = bead_info[i][BEAD_TYPE];
      oldEn += energy_cont_and_ovlp(i);
      if (TypeCanRot[resi] != 1){//Skip beads that cannot bond.
        continue;
      }
  if (bead_info[i][BEAD_FACE] != -1){//I am bonded to something
      resj = bead_info[bead_info[i][BEAD_FACE]][BEAD_TYPE];//Type of bead I am bonded to
      oldEn += fEnergy[resi][resj][E_SC_SC];//Adding the energy.
    }
    ShuffleRotIndecies();
    BWWeight = CheckRotStatesOld(i, resi, MyTemp);
    NormalizeRotState(yTemp, BWWeight);
    yTemp++;
  }
  BSum = 1.;
  for (i = 0; i < yTemp; i++){
      BSum += log10(bolt_norm[i]);
  }

    //for (i = anchorBead + 1; i < lastB; i++) {
    for (j = 0; j < listLen; j++){
      i  = tmpList[j];
      RotOperation(PivotM, i, anchorPos);
      move_bead_to(i, naTempR);
  }

    //for (i = anchorBead + 1; i < lastB; i++) {
    for (j = 0; j < listLen; j++){
      i  = tmpList[j];
      if(bead_info[i][BEAD_FACE] != -1){
          bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE] = -1;
            bead_info[i][BEAD_FACE] = -1;
      }
  }

  yTemp = 0;
    //for (i = anchorBead + 1; i < lastB; i++) {
    for (j = 0; j < listLen; j++){
      i  = tmpList[j];
      resi = bead_info[i][BEAD_TYPE];
      newEn += energy_cont_and_ovlp(i);
      if (TypeCanRot[resi] != 1){//Because linkers don't have rotational states
          continue;
      }
      ShuffleRotIndecies();
      FWWeight = CheckRotStatesNew(i, resi, MyTemp);
      NormalizeRotState(yTemp, FWWeight);
      //Note that the bonds need to be formed in this loop so that we don't overcount!
      if(bead_info[i][BEAD_FACE] == -1){//Make sure this bead is unbonded!
      //Let's assign a rotational state to this bead
        xTemp = PickRotState(FWWeight);
        if(xTemp != -1){//An appropriate partner has been selected. Form the bonds and add the energy
          resj = bead_info[xTemp][BEAD_TYPE];
          bead_info[i][BEAD_FACE] = xTemp;
          bead_info[xTemp][BEAD_FACE] = i;
          newEn += fEnergy[resi][resj][E_SC_SC];
          }
        }
      yTemp++;
  }

  FSum = 1.;
  for (i = 0; i < yTemp; i++){
      FSum += log10(bolt_norm[i]);
  }

//printf("%f\n", (FSum/BSum)*expf((oldEn-newEn)/MyTemp));
  MCProb = (float)rand()/(float)RAND_MAX;
  if ( MCProb < (FSum/BSum)*expf((oldEn-newEn)/MyTemp)){//Accept the move. Remember that the bonds were assigned above!
   bAccept = 1;
   return bAccept;
  }
  else{//Rejecting move
      //for (i = anchorBead + 1; i < lastB; i++) {
      for (j = 0; j < listLen; j++){
        i  = tmpList[j];
      if(bead_info[i][BEAD_FACE] != -1){
        bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE] = -1;
        bead_info[i][BEAD_FACE] = -1;
      }
    }
      //for (i = anchorBead + 1; i < lastB; i++) {
      for (j = 0; j < listLen; j++){
        i  = tmpList[j];
      undo_move_bead_to(i);
    }
      bAccept = 0;
      return bAccept;
  }
}

int BranchedRotMCMove(int chainID, float MyTemp){
      //Rotates a branched molecule about the branching, which is assumed to be the firstB of chainID
      //Performs a PivotMCMove() on molecule where the rotation occurs around firstB
      /*
      Set the first bead as anchorBead, and perform a symmetry operation on the beads after the anchorBead (the whole molecule). Note that if the molecule is linear, the move is outright rejected. Again, it is assumed that the first bead in that molecule is the 'node'.
      */
      int bAccept = 0;
      //Reject if the molecule is linear
      if(TypeIsLinear[chain_info[chainID][CHAIN_TYPE]] == 1){
        bAccept = 0;
        return bAccept;
      }
      int firstB, lastB;
      //Get the beadID's for the first and last bead
      firstB = chain_info[chainID][CHAIN_START];
      lastB  = firstB + chain_info[chainID][CHAIN_LENGTH];
      //Randomly select a bead
      int anchorBead  = firstB;
      anchorBead = firstB + anchorBead;

      //Randomly selecting a symmetry operation
      int PivotM;
      PivotM = rand() % 10;
      //printf("Move: %d\n", PivotM);
      if(PivotM == 0){
        bAccept = 1;
        return bAccept;
      }

      int i,j;
      int xTemp, yTemp;
      int anchorPos[POS_MAX];
      for(j=0; j<POS_MAX; j++){
        anchorPos[j] = bead_info[anchorBead][j];
      }

      xTemp = 0; yTemp = 0;
      while (xTemp < nMCMaxTrials && yTemp == 0){
        for(i=anchorBead+1; i<lastB; i++){
          RotOperation(PivotM, i, anchorPos);
          yTemp = check_move_bead_to(naTempR);
          if(yTemp == 0){
            xTemp++;
            break;
          }
        }
        xTemp++;
      }
      if(xTemp == nMCMaxTrials || yTemp == 0){
        bAccept = 0;
        return bAccept;
      }

      int FWWeight, BWWeight;//Used to perform orientational bias MC
      double FSum, BSum;
      int resi, resj;
      float oldEn = 0.;
      float newEn = 0.;
      float MCProb;

      yTemp = 0;
      for (i=anchorBead+1; i<lastB; i++){
          resi = bead_info[i][BEAD_TYPE];
          oldEn += energy_cont_and_ovlp(i);
          if (TypeCanRot[resi] != 1){//Skip beads that cannot bond.
            continue;
          }
      if (bead_info[i][BEAD_FACE] != -1){//I am bonded to something
          resj = bead_info[bead_info[i][BEAD_FACE]][BEAD_TYPE];//Type of bead I am bonded to
          oldEn += fEnergy[resi][resj][E_SC_SC];//Adding the energy.
        }
        ShuffleRotIndecies();
        BWWeight = CheckRotStatesOld(i, resi, MyTemp);
        NormalizeRotState(yTemp, BWWeight);
        yTemp++;
      }
        BSum = 1.;
        for (i = 0; i < yTemp; i++){
          BSum += log10(bolt_norm[i]);
        }

        for(i=anchorBead+1; i<lastB; i++){
          RotOperation(PivotM, i, anchorPos);
          move_bead_to(i, naTempR);
        }

        for(i=anchorBead+1; i<lastB; i++){
          if(bead_info[i][BEAD_FACE] != -1){
            bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE] = -1;
            bead_info[i][BEAD_FACE] = -1;
          }
        }

        yTemp = 0;
         for (i=anchorBead+1; i < lastB; i++){//Counting states in the new location
             resi = bead_info[i][BEAD_TYPE];
             newEn += energy_cont_and_ovlp(i);
             if (TypeCanRot[resi] != 1){//Because linkers don't have rotational states
              continue;
             }
          ShuffleRotIndecies();
          FWWeight = CheckRotStatesNew(i, resi, MyTemp);
          NormalizeRotState(yTemp, FWWeight);
          //Note that the bonds need to be formed in this loop so that we don't overcount!
          if(bead_info[i][BEAD_FACE] == -1){//Make sure this bead is unbonded!
          //Let's assign a rotational state to this bead
            xTemp = PickRotState(FWWeight);
            if(xTemp != -1){//An appropriate partner has been selected. Form the bonds and add the energy
              resj = bead_info[xTemp][BEAD_TYPE];
              bead_info[i][BEAD_FACE] = xTemp;
              bead_info[xTemp][BEAD_FACE] = i;
              newEn += fEnergy[resi][resj][E_SC_SC];
              }
            }
          yTemp++;
        }
        FSum = 1.;

        for (i = 0; i < yTemp; i++){
              FSum += log10(bolt_norm[i]);
          }

      MCProb = (float)rand()/(float)RAND_MAX;
      if ( MCProb < (FSum/BSum)*expf((oldEn-newEn)/MyTemp)){//Accept the move. Remember that the bonds were assigned above!
       bAccept = 1;
       return bAccept;
      }
      else{//Rejecting move
        for(i=anchorBead+1; i<lastB; i++){
          if(bead_info[i][BEAD_FACE] != -1){
            bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE] = -1;
            bead_info[i][BEAD_FACE] = -1;
          }
        }
        for(i=anchorBead+1; i<lastB; i++){
          undo_move_bead_to(i);
        }
          bAccept = 0;
          return bAccept;
      }
}

int LocalMCMove_Equil(int beadID, float MyTemp){//Performs a local translation MC-move on beadID

  int bAccept = 0; //Used in MC steps
  float MCProb, oldEn, newEn; //For Metropolis Hastings
  oldEn = 0.;
  newEn = 0.;
  int i, j;//Loop iterators
  int xTemp, yTemp, lRadUp, lRadLow;//Random numbers to store things
  int tmpR[POS_MAX], tmpR2[POS_MAX];//Vectors to stores coordinates.
  int FWWeight, BWWeight;//Used to perform orientational bias MC
  //printf("Beginning LOCAL\n");
  for(j=0; j<POS_MAX; j++){//Initializing the vectors to where this bead is.
    tmpR[j]  = bead_info[beadID][j];
  }
  //Initialize the radii for the search of next trial location
  lRadLow = linker_len[beadID][0];
  lRadUp  = lRadLow*2 + 1;//2*2+1

  xTemp = 0; yTemp = 0;//Initialize these guys.
  while(yTemp == 0 && xTemp < nMCMaxTrials){//Attempt to find an empty lattice point.
    for(j=0; j<POS_MAX; j++){
      tmpR2[j] = (rand() % lRadUp) - lRadLow;//Generate number between -2 and 2
      tmpR2[j] = (tmpR[j] + tmpR2[j] + nBoxSize[j]) % nBoxSize[j];
    }
    yTemp = check_move_bead_to(tmpR2);
    if(yTemp == 1){//This means we found an empty lattice site. So let's check if the linkers are okay.
      yTemp = check_linker_constraint(beadID, tmpR2);
    }
    xTemp++;
  }
  if(xTemp == nMCMaxTrials || yTemp == 0){//This means that we have failed to find an appropriate spot for this bead to be moved to. Therefore, the move is rejected!
    bAccept = 0;
    //printf("End LOCAL - No space\n");
    return bAccept;
  }
  //Have successfully found a good lattice spot.
    FWWeight = 1; BWWeight = 1;
    oldEn = energy_cont_and_ovlp(beadID);
    move_bead_to(beadID, tmpR2);
    //Now let's calculate the energy of the new state. SC-SC energy is already done.
    newEn += energy_cont_and_ovlp(beadID);
    MCProb = (float)rand()/(float)RAND_MAX;
    if (MCProb < ((float)FWWeight/(float)BWWeight)*expf((oldEn-newEn)/MyTemp)){//Accept this state
      bAccept = 1;//Accepting!
      //printf("End LOCAL - Accept\n");
      return bAccept;
    }
    else{
      undo_move_bead_to(beadID);
      bAccept = 0;
      //printf("End LOCAL - Fail\n");
      return bAccept;
    }
}

int SlitherMCMove_Equil(int chainID, float MyTemp){//Performs a slither MC-move on chainID

  int firstB, lastB;//Track first and last+1 bead of chainID. Makes reading easier.
  int bAccept = 0; //Used in MC steps
  //Finding the bounds for looping over the molecule/chain
  firstB = chain_info[chainID][CHAIN_START];
  lastB  = firstB + chain_info[chainID][CHAIN_LENGTH];
  if(lastB - firstB == 1){//This means we have a monomer. Reject the move, because Local or Trans
    //moves should be the ones that move monomers.
    bAccept = 0;
    return bAccept;
  }
  else{
    if(TypeIsLinear[chain_info[chainID][CHAIN_TYPE]] != 1){//If chain is not linear. Reject move because slithering will not werk!
      bAccept = 0;
      return bAccept;
    }
  }
  //This chain is suitable to perform a slithering-snake move on.

  float MCProb, oldEn, newEn; //For Metropolis Hastings
  oldEn = 0.; newEn = 0.;
  int i, j, k;//Loop iterators
  int xTemp, yTemp, lRadUp, lRadLow;//Random numbers to store things
  int tmpR[POS_MAX], tmpR2[POS_MAX], tmpR3[POS_MAX];//Vectors to store positions.

  MCProb = (float)rand()/(float)RAND_MAX;//To decide if we slither forwards or backwards
  if(MCProb < 0.5){//Forwards slither, so lastB-1 (last bead) is anchor
  lRadUp  = (int)2 * linker_len[lastB - 1][0] + 1;//lastB-1 will be replaced by lastB-2
  lRadLow = (int)linker_len[lastB - 1][0];
  yTemp = 0; xTemp = 0;//Using to track trials for placing beads
  while (xTemp < nMCMaxTrials && yTemp == 0) {
    for(j=0; j<POS_MAX; j++){
        tmpR[j] = (rand() % lRadUp) - lRadLow;
        tmpR[j] = (bead_info[lastB-1][j] + tmpR[j] + nBoxSize[j]) % nBoxSize[j];
      }
    yTemp = check_move_bead_to(tmpR);// 0: there is no space, 1: there is space
    xTemp++;
      }
    }
  else{//Backwards slither, so firstB is anchor
      lRadUp  = (int)2 * linker_len[firstB][0] + 1;//firstB will be replaced by firstB+1
      lRadLow = (int)linker_len[firstB][0];
      yTemp = 0; xTemp = 0;//Using to track trials for placing beads
      while (xTemp < nMCMaxTrials && yTemp == 0) {
        for(j=0; j<POS_MAX; j++){
          tmpR[j] = (rand() % lRadUp) - lRadLow;
          tmpR[j] = (bead_info[firstB][j] + tmpR[j] + nBoxSize[j]) % nBoxSize[j];
        }
      yTemp = check_move_bead_to(tmpR);// 0: there is no space, 1: there is space
      xTemp++;
        }
    }
    if (yTemp == 0 || xTemp == nMCMaxTrials){//Couldn't find a spot, so reject the damn move
      bAccept = 0;
      return bAccept;
    }
    //We should have a spot to move to! tmpR has the location

   //Let's remember where this chain exists.
  for(i=firstB; i<lastB; i++){
    for(j=0; j<BEADINFO_MAX; j++){
      old_bead[i][j] = bead_info[i][j];
    }

  }

    for (i=firstB; i<lastB; i++){//Counting states in the previous location
      oldEn += energy_cont_and_ovlp(i);
    }

  if(MCProb < 0.5){
    //Slithering the chain forwards in ID-space
  for (i=firstB; i<lastB-1; i++){
    for(j=0;j<POS_MAX;j++){
      tmpR2[j] = bead_info[i][j];
      bead_info[i][j] = old_bead[i+1][j];//Hopping over by one bead
      tmpR3[j] = bead_info[i][j];
    }
    if(i==firstB){//Only the firstB's location is empty
      naTotLattice[LtIndV(tmpR2)] = -1;
    }
      naTotLattice[LtIndV(tmpR3)] = i;
  }
  //Moving the last bead, and it has to be done independently because lastB-1 -> tmpR
  i = lastB-1;
  for(j=0;j<POS_MAX;j++){
    bead_info[i][j] = tmpR[j];
  }
      naTotLattice[LtIndV(tmpR)]    = i;
  }
  else{//Slithering backwards in ID-space
    for (i=firstB+1; i<lastB; i++){
      for(j=0; j<POS_MAX; j++){
        tmpR2[j] = bead_info[i][j];
        bead_info[i][j] = old_bead[i-1][j];//Hopping back by one bead
        tmpR3[j] = bead_info[i][j];
      }
      if(i==lastB-1){//Only the lastB-1's location is empty
        naTotLattice[LtIndV(tmpR2)] = -1;
      }
        naTotLattice[LtIndV(tmpR3)] = i;
    }
    //Moving the first bead, and it has to be done independently because firstB -> tmpR
    i = firstB;
    for(j=0; j<POS_MAX; j++){
      bead_info[i][j] = tmpR[j];
    }
      naTotLattice[LtIndV(tmpR)]    = i;
  }

     for (i=firstB; i<lastB; i++){//Counting states in the new location
         newEn += energy_cont_and_ovlp(i);
       }
     //Doing the Metropolis-Hastings thing
     MCProb = (float)rand()/(float)RAND_MAX;
     if (MCProb < expf((oldEn-newEn)/MyTemp)){//Accept. Bonds have been handled before!
      bAccept = 1;
      return bAccept;
    }
     else{
       for(i=firstB; i<lastB; i++){//Resetting the lattice
         for(j=0;j<POS_MAX; j++){//This is where I am
          tmpR2[j] = bead_info[i][j];
        }
           naTotLattice[LtIndV(tmpR2)]  = -1;
         }
       for(i=firstB; i<lastB; i++){
         for(j=0;j<BEADINFO_MAX; j++){
            bead_info[i][j] = old_bead[i][j];//Restoring
            if(j<POS_MAX){//This is where I should be
              tmpR2[j] = bead_info[i][j];
            }
         }
           naTotLattice[LtIndV(tmpR2)]  = i;
         }
         bAccept = 0;
         return bAccept;
      }
}

int TransMCMove_Equil(int chainID, float MyTemp){//Performs a translation move with orientational bias
  int bAccept = 0; //Used in MC steps
  float MCProb, oldEn, newEn; //For Metropolis Hastings
  oldEn = 0.; newEn = 0.;
  int i, j;//Loop iterators
  int firstB, lastB;
  int xTemp, yTemp, lRadUp, lRadLow;//Random numbers to store things
  int tmpR[POS_MAX];//Vectors to store coordinates.
  //Finding the bounds for looping over the molecule/chain
 firstB = chain_info[chainID][CHAIN_START];
 lastB  = firstB + chain_info[chainID][CHAIN_LENGTH];
  //Radii for translation moves. All moves are L/4 radius
  lRadLow = nBoxSize[2]/2;
  lRadUp  = 2*lRadLow + 1;
  //Initialize these iterators.
  //printf("Beginning TRANS\n");
  xTemp = 0; yTemp = 0;
  while(xTemp < nMCMaxTrials && yTemp == 0){
      for (j=0; j<POS_MAX; j++) {
        tmpR[j] = (rand() % lRadUp) - lRadLow;//Random vector to move all beads within r=L/4
      }
      yTemp = check_disp_chain(chainID, tmpR);//yTemp=0 means clash
      xTemp++;
    }
  if (yTemp == 0 || xTemp == nMCMaxTrials){//We have failed to find a good spot for this chain.
    bAccept = 0;
    //printf("Ending TRANS no space\n");
    return bAccept;
  }
  //We now have a chain which when moved does not overlap.
  //Initialize the orientational-bias sums and vectors.
  //Counting states in the previous location
  for (i = firstB; i < lastB; i++){
      oldEn += energy_cont_and_ovlp(i);
  }
   trans_disp_chain(chainID, tmpR);//Moved the chain, broke bonds, and remembered stuff
   for (i = firstB; i < lastB; i++){//Counting states in the new location
       newEn += energy_cont_and_ovlp(i);
  }
   MCProb = (float)rand()/(float)RAND_MAX;
   if ( MCProb < expf((oldEn-newEn)/MyTemp)){//Accept the move. Remember that the bonds were assigned above!
    bAccept = 1;
    return bAccept;
   }
   else{
       trans_restore_chain(chainID);
       bAccept = 0;
       return bAccept;
   }

}

int ShakeMove_Equil(int beadID, float MyTemp){

  int topIt; //Iterator for topo_info
  int i,j; //Loop iterators
  int curID; //current bead being looked at
  int tmpR[MAX_VALENCY][POS_MAX]; //Storing temporary locations
  int  bAccept;
  int lRadLow, lRadUp;
  lRadLow = 2;
  lRadUp  = lRadLow*2 + 1;//2*2+1
  int xTemp = 0;
  int yTemp = 0;
    curID = beadID; topIt = 0;
    while (curID != -1) {
        for (j=0; j<BEADINFO_MAX; j++) {
            old_bead[curID][j] = bead_info[curID][j];//Remembering
            if(j<POS_MAX){
                tmpR[topIt][j] = bead_info[curID][j];
            }
        }
        curID = topo_info[beadID][topIt++];
    }

    while(yTemp == 0 && xTemp < nMCMaxTrials){
        curID = beadID; topIt = 0;
        while (curID != -1){
            naTotLattice[LtIndV(bead_info[curID])] = -1;
            curID = topo_info[beadID][topIt++];
        }
        curID = beadID; topIt = 0;
        while (curID != -1){
            yTemp = 1;
            for(j = 0; j < POS_MAX; j++){
                tmpR[topIt][j] = (rand() % lRadUp) - lRadLow;
                tmpR[topIt][j] = (bead_info[curID][j] + tmpR[topIt][j] + nBoxSize[j]) % nBoxSize[j];
            }
            if (naTotLattice[LtIndV(tmpR[topIt])] != -1 ){
                yTemp = 0;
                break;
            }
            naTotLattice[LtIndV(tmpR[topIt])] = curID;
            curID = topo_info[beadID][topIt++];
        }
        if (yTemp == 1){//No steric clash so check for topology constraint
            yTemp = ShakeConstraint(beadID, tmpR);
        }
        for(i = 0; i < topIt; i++){
            naTotLattice[LtIndV(tmpR[i])] = -1;
        }
        xTemp++;
    }

    if(xTemp == nMCMaxTrials || yTemp == 0){//Linker or steric clash didn't work out
        curID = beadID; topIt = 0;
        while (curID != -1) {
            naTotLattice[LtIndV(bead_info[curID])] = curID;
            curID = topo_info[beadID][topIt++];
        }
        //printf("No space!\n");
        bAccept = 0;
        return bAccept;
    }

  float oldEn = 0.;
  float newEn = 0.;
  float MCProb;

  //printf("Starting BW\n");
  curID = beadID; topIt = 0;
  while (curID != -1){
    oldEn += energy_cont_and_ovlp(curID);
    curID = topo_info[beadID][topIt++];
  }

  curID = beadID; topIt = 0;
  while (curID != -1){
    move_bead_to_shake(curID, tmpR[topIt]);
    curID = topo_info[beadID][topIt++];
  }

  curID = beadID; topIt = 0;
  while (curID != -1){
    newEn += energy_cont_and_ovlp(curID);
    curID = topo_info[beadID][topIt++];
  }

    MCProb = (float)rand()/(float)RAND_MAX;
    if ( MCProb < expf((oldEn-newEn)/MyTemp)){//Accept the move. Remember that the bonds were assigned above!
     bAccept = 1;
     //printf("Win\n" );
     return bAccept;
    }
    else{
        curID = beadID; topIt = 0;
        while (curID != -1){
            naTotLattice[LtIndV(bead_info[curID])] = -1;
            curID = topo_info[beadID][topIt++];
        }
        curID = beadID; topIt = 0;
        while (curID != -1){
            for(i = 0; i < BEADINFO_MAX; i++){
                bead_info[curID][i] = old_bead[curID][i];
            }
              naTotLattice[LtIndV(bead_info[curID])] = curID;
              curID = topo_info[beadID][topIt++];
        }
        //printf("fail\n");
        bAccept = 0;
        return bAccept;
    }


}

int PivotMCMove_Equil(int chainID, float MyTemp){
  //Performs a pivot move on chainID
  /*
  Randomly pick a bead in this chain (anchorBead), and perform a symmetry operation on the beads after the anchorBead. Note that if anchorBead is the second-to-last, or last, bead, the move is rejected. Local and slither moves would be better for those.
  Obviously the protein must be linear and the move is rejected if the chain is branched.
  It is assumed that a linear chain has been passed!
  */
  int bAccept = 0;
  //Check if the chain is longer than 3 or if it is linear
  if( chain_info[chainID][CHAIN_LENGTH] <=3 || TypeIsLinear[chain_info[chainID][CHAIN_TYPE]] != 1){
    bAccept = 0;
    return bAccept;
  }
  int firstB, lastB;
  //Get the beadID's for the first and last bead
  firstB = chain_info[chainID][CHAIN_START];
  lastB  = firstB + chain_info[chainID][CHAIN_LENGTH];
  //This chain is long enough
  int chainLength = chain_info[chainID][CHAIN_LENGTH];
  //Randomly select a bead
  int anchorBead  = rand() % chainLength;
  anchorBead = firstB + anchorBead;
  //printf("Bead: %d\n", anchorBead);
  //Checking if bead is first or last. If so, reject move.
  if(anchorBead == firstB || anchorBead > lastB-2){
    bAccept = 0;
    //printf("P Bad Bead\n");
    return bAccept;
  }

  int PivotM;
  PivotM = rand() % 10;
  //printf("Move: %d\n", PivotM);
  if(PivotM == 0){
    //printf("P Accept\n");
    bAccept = 1;
    return bAccept;
  }

  int i,j;
  int xTemp, yTemp;
  int anchorPos[POS_MAX];
  for(j=0; j<POS_MAX; j++){
    anchorPos[j] = bead_info[anchorBead][j];
  }

  xTemp = 0; yTemp = 0;
  while (xTemp < nMCMaxTrials && yTemp == 0){
    for(i=anchorBead+1; i<lastB; i++){
      RotOperation(PivotM, i, anchorPos);
      yTemp = check_move_bead_to(naTempR);
      if(yTemp == 0){
        xTemp++;
        break;
      }
    }
    xTemp++;
  }
  if(xTemp == nMCMaxTrials || yTemp == 0){
    bAccept = 0;
    //printf("P Clashout\n");
    return bAccept;
  }

  float oldEn = 0.;
  float newEn = 0.;
  float MCProb;

  yTemp = 0;
  for (i = anchorBead+1; i < lastB; i++){
      oldEn += energy_cont_and_ovlp(i);
  }
    for(i=anchorBead+1; i<lastB; i++){
      RotOperation(PivotM, i, anchorPos);
      move_bead_to(i, naTempR);
    }

    yTemp = 0;
     for (i = anchorBead+1; i < lastB; i++){//Counting states in the new location
         newEn += energy_cont_and_ovlp(i);
    }

  MCProb = (float)rand()/(float)RAND_MAX;
  if ( MCProb < expf((oldEn-newEn)/MyTemp)){//Accept the move. Remember that the bonds were assigned above!
   bAccept = 1;
   //printf("P Win\n");
   return bAccept;
  }
  else{
    for(i=anchorBead+1; i<lastB; i++){
      undo_move_bead_to(i);
    }
    //printf("P Fail\n");
      bAccept = 0;
      return bAccept;
  }


}

int BranchedRotMCMove_Equil(int chainID, float MyTemp){
      //Rotates a branched molecule about the branching, which is assumed to be the firstB of chainID
      //Performs a PivotMCMove() on molecule where the rotation occurs around firstB
      /*
      Set the first bead as anchorBead, and perform a symmetry operation on the beads after the anchorBead (the whole molecule). Note that if the molecule is linear, the move is outright rejected. Again, it is assumed that the first bead in that molecule is the 'node'.
      */
      int bAccept = 0;
      //Reject if the molecule is linear
      if(TypeIsLinear[chain_info[chainID][CHAIN_TYPE]] == 1){
        bAccept = 0;
        return bAccept;
      }
      int firstB, lastB;
      //Get the beadID's for the first and last bead
      firstB = chain_info[chainID][CHAIN_START];
      lastB  = firstB + chain_info[chainID][CHAIN_LENGTH];

      //Randomly select a bead
      int anchorBead  = firstB;
      anchorBead = firstB + anchorBead;

      //Randomly selecting a symmetry operation
      int PivotM;
      PivotM = rand() % 10;
      //printf("Move: %d\n", PivotM);
      if(PivotM == 0){
        bAccept = 1;
        return bAccept;
      }

      int i,j;
      int xTemp, yTemp;
      int anchorPos[POS_MAX];
      for(j=0; j<POS_MAX; j++){
        anchorPos[j] = bead_info[anchorBead][j];
      }

      xTemp = 0; yTemp = 0;
      while (xTemp < nMCMaxTrials && yTemp == 0){
        for(i=anchorBead+1; i<lastB; i++){
          RotOperation(PivotM, i, anchorPos);
          yTemp = check_move_bead_to(naTempR);
          if(yTemp == 0){
            xTemp++;
            break;
          }
        }
        xTemp++;
      }
      if(xTemp == nMCMaxTrials || yTemp == 0){
        bAccept = 0;
        return bAccept;
      }

      float oldEn = 0.;
      float newEn = 0.;
      float MCProb;

      for (i=anchorBead+1; i<lastB; i++){
          oldEn += energy_cont_and_ovlp(i);
      }
        for(i=anchorBead+1; i<lastB; i++){
          RotOperation(PivotM, i, anchorPos);
          move_bead_to(i, naTempR);
        }
         for (i=anchorBead+1; i < lastB; i++){//Counting states in the new location
             newEn += energy_cont_and_ovlp(i);
        }

      MCProb = (float)rand()/(float)RAND_MAX;
      if ( MCProb < expf((oldEn-newEn)/MyTemp)){//Accept the move. Remember that the bonds were assigned above!
       bAccept = 1;
       return bAccept;
      }
      else{//Rejecting move
        for(i=anchorBead+1; i<lastB; i++){
          undo_move_bead_to(i);
        }
          bAccept = 0;
          return bAccept;
      }
}

int check_disp_chain(int chainID, const int tR[]){//Checks if chain can be displaced by tR
  int i,j;
  int canI = 1;
  int tmpR[POS_MAX];

  for (i=chain_info[chainID][CHAIN_START];i<chain_info[chainID][CHAIN_START]+chain_info[chainID][CHAIN_LENGTH]; i++){
      for(j=0;j<POS_MAX;j++){
      tmpR[j] = (bead_info[i][j] + tR[j] + nBoxSize[j]) % nBoxSize[j];
    }
    if (naTotLattice[LtIndV(tmpR)] != -1){
     canI = 0;// 0 means steric clash
     return canI;
    }
  }

  return canI;
}

void disp_chain(int chainID, const int movR[]){
  //Displaces current chain by movR and handles the lattice
  //Also remembers where everything was moved and saves into old_bead

  int i, l;
  int fB = chain_info[chainID][CHAIN_START];
  int lB = fB + chain_info[chainID][CHAIN_LENGTH];
  int tmpR[POS_MAX], tmpR2[POS_MAX];
    for (i=fB; i<lB; i++){
      for (l=0; l<BEADINFO_MAX; l++){
      old_bead[i][l] = bead_info[i][l];
      if(l<POS_MAX){
      tmpR[l] = old_bead[i][l];//Where we were
      bead_info[i][l] = (tmpR[l] + movR[l] + nBoxSize[l]) % nBoxSize[l];
      tmpR2[l] = bead_info[i][l];//Where we are now
    }
    }
    naTotLattice[LtIndV(tmpR)] = -1;//Removing from old place
    naTotLattice[LtIndV(tmpR2)] = i;
  }
}

void trans_disp_chain(int chainID, const int movR[]){
  //Displaces current chain by movR and handles lattice
  //Specific for TransMCMove because it breaks old bonds!
  //Also remembers where everything was moved and saves into old_bead

  int i, l;
  int fB = chain_info[chainID][CHAIN_START];
  int lB = fB + chain_info[chainID][CHAIN_LENGTH];
  int tmpR[POS_MAX], tmpR2[POS_MAX];
    for (i=fB; i<lB; i++){
      for (l=0; l<BEADINFO_MAX; l++){
      old_bead[i][l] = bead_info[i][l];
      if(l<POS_MAX){
      tmpR[l] = old_bead[i][l];//Where we were
      bead_info[i][l] = (tmpR[l] + movR[l] + nBoxSize[l]) % nBoxSize[l];
      tmpR2[l] = bead_info[i][l];//Where we are now
    }
    }
      naTotLattice[LtIndV(tmpR)]  = -1;//Removing from old place
      naTotLattice[LtIndV(tmpR2)] = i;//Adding to new place
  }
  for(i=fB; i<lB; i++){//Delete bonds AFTER old_bead has remembered things
    if(old_bead[i][BEAD_FACE] != -1){
      bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE] = -1;//Break old bond
      bead_info[i][BEAD_FACE] = -1;//Breaking this bond
    }
  }
}

void restore_chain(int chainID){//Uses old_bead to undo what disp_chain does.
  int i, l;
  int fB = chain_info[chainID][CHAIN_START];
  int lB = fB + chain_info[chainID][CHAIN_LENGTH];
  int tmpR[POS_MAX], tmpR2[POS_MAX];
  for (i=fB; i<lB; i++){
      for (l=0; l<BEADINFO_MAX; l++){
        if(l<POS_MAX){
      tmpR[l] = bead_info[i][l];//Where we now are and must be removed from
      tmpR2[l] = old_bead[i][l];//Where we were and will be moved back to
    }
      bead_info[i][l] = old_bead[i][l];//Moving back
    }
      naTotLattice[LtIndV(tmpR)] = -1;//Removing from old place
      naTotLattice[LtIndV(tmpR2)] = i;
  }
}

void trans_restore_chain(int chainID){//Uses old_bead to undo what trans_disp_chain does.
  //Note that removing and placing separately makes sure that no incorrect bonds are formed, or unformed.
  int i, l;
  int fB = chain_info[chainID][CHAIN_START];
  int lB = fB + chain_info[chainID][CHAIN_LENGTH];
  int tmpR[POS_MAX], tmpR2[POS_MAX];
  for (i=fB; i<lB; i++){//Let's remove from the new state
    if(bead_info[i][BEAD_FACE] != -1){//Destroying the newly proposed bond
      bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE] = -1;
      bead_info[i][BEAD_FACE] = -1;
    }
      for (l=0; l<POS_MAX; l++){
      tmpR[l]  = bead_info[i][l];//Where we now are and must be removed from
      tmpR2[l] = old_bead[i][l];//This is where we should be
  }
      naTotLattice[LtIndV(tmpR)] = -1;//Removing from old place
      naTotLattice[LtIndV(tmpR2)] = i;//Removing from old place
  }
  for (i=fB; i<lB; i++){
    for(l=0; l<BEADINFO_MAX; l++){
    bead_info[i][l] = old_bead[i][l];//Moving back
  }
  }
  for (i=fB; i<lB; i++){
    if(old_bead[i][BEAD_FACE] != -1){//Restoring the old bond
      bead_info[i][BEAD_FACE] = old_bead[i][BEAD_FACE];
      bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE] = i;
    }
  }
}

inline int check_move_bead_to(int newPos[]){//Checks if I can move here

  if (naTotLattice[LtIndV(newPos)] != -1){
     return 0;//There's something here already, brah
  }
  return 1;
}

void move_bead_to(int beadID, const int newPos[]){//Updates position to new one and handles lattice
  int i;
  int tmpR[POS_MAX], tmpR2[POS_MAX];
  for (i=0;i<BEADINFO_MAX;i++){
    old_bead[beadID][i]  = bead_info[beadID][i];
    if(i<POS_MAX){
        tmpR[i] = old_bead[beadID][i];
        bead_info[beadID][i] = newPos[i];
        tmpR2[i] = bead_info[beadID][i];
    }
  }
    naTotLattice[LtIndV(tmpR)] = -1;//Removing from old place
    naTotLattice[LtIndV(tmpR2)] = beadID;
}

void undo_move_bead_to(int beadID){//Undoes what the above function does
  int i;
  int tmpR[POS_MAX], tmpR2[POS_MAX];
  for (i=0;i<BEADINFO_MAX;i++){
    if(i<POS_MAX){
    tmpR[i]  = bead_info[beadID][i];
    tmpR2[i] = old_bead[beadID][i];
  }
    bead_info[beadID][i] = old_bead[beadID][i];
  }

    naTotLattice[LtIndV(tmpR)] = -1;//Removing from old place
    naTotLattice[LtIndV(tmpR2)] = beadID;//Restoring

  if(bead_info[beadID][BEAD_FACE] != -1){
    bead_info[bead_info[beadID][BEAD_FACE]][BEAD_FACE] = beadID;
  }
}

void move_bead_to_shake(int beadID, const int newPos[]){//Updates position to new one and handles lattice specifically for shake move
    int i;
    int tmpR2[POS_MAX];
    for (i=0;i<POS_MAX;i++){
        bead_info[beadID][i] = newPos[i];
        tmpR2[i] = bead_info[beadID][i];
    }
    naTotLattice[LtIndV(tmpR2)] = beadID;
    i = bead_info[beadID][BEAD_FACE];
    if (i != -1){
        bead_info[i][BEAD_FACE] = -1;
        bead_info[beadID][BEAD_FACE] = -1;
    }
}

inline int check_disp_bead(int beadID, const int movR[]){//Checks if bead can be displaced by tmpR
  int l;
  int canI = 1;
  int tmpR[POS_MAX];
  for (l=0; l<POS_MAX; l++){
  tmpR[l] = (bead_info[beadID][l] + movR[l] + nBoxSize[l]) % nBoxSize[l];//Where we want to go
  }
  if (naTotLattice[LtIndV(tmpR)] != -1){
    canI = 0;//There's something here already, brah
  }
  return canI;
}

void swap_beads(int bead1, int bead2){
  //Swaps ALL properties of the beads, and exchanged bonding partners
  int MyF1, MyF2;
  int i;
  int tmpR[POS_MAX], tmpR2[POS_MAX];
  //First the coordinates
    for(i=0; i<POS_MAX; i++){
      tmpR[i]  = bead_info[bead1][i];
      tmpR2[i] = bead_info[bead2][i];
      bead_info[bead1][i] = tmpR2[i];
      bead_info[bead2][i] = tmpR[i];
    }

    //Onto bonds and partners
    MyF1 = bead_info[bead1][BEAD_FACE];//Bead 1's bonding partner
    MyF2 = bead_info[bead2][BEAD_FACE];//Bead 2's bonding partner.
    //If they are bonded to each other, don't do anything
    if(MyF1 != bead2){
    bead_info[bead1][BEAD_FACE] = MyF2;
    bead_info[bead2][BEAD_FACE] = MyF1;
    if(MyF1 != -1){//Need to swap partners -- very 2019
        bead_info[MyF1][BEAD_FACE] = bead2;//It's bonded to bead2 now
    }
    if(MyF2 != -1){
        bead_info[MyF2][BEAD_FACE] = bead1;
    }
  }
  //Swap them on the lattice
    naTotLattice[LtIndV(tmpR)]    = bead2;
    naTotLattice[LtIndV(tmpR2)]   = bead1;
}

void Rot_X_90(int beadID, const int tmpR[POS_MAX]){

  int j;
  //Performs a rotation by 90 degress around x-axis. tmpR is the new origin.
  for(j=0; j<POS_MAX; j++){
    naTempR[j] = bead_info[beadID][j];
  }
  naTempR[POS_Y] = tmpR[POS_Y] + tmpR[POS_Z] - bead_info[beadID][POS_Z];
  naTempR[POS_Z] = bead_info[beadID][POS_Y]  - tmpR[POS_Y] + tmpR[POS_Z];
  //Adjusting for periodic boundaries
  naTempR[POS_Y] = (naTempR[POS_Y] + nBoxSize[POS_Y]) % nBoxSize[POS_Y];
  naTempR[POS_Z] = (naTempR[POS_Z] + nBoxSize[POS_Z]) % nBoxSize[POS_Z];
}

void Rot_X_180(int beadID, const int tmpR[POS_MAX]){

  int j;
  //Performs a rotation by 180 degress around x-axis. tmpR is the new origin for the rotation.
  for(j=0; j<POS_MAX; j++){
    naTempR[j] = bead_info[beadID][j];
  }
  naTempR[POS_Y] = 2*tmpR[POS_Y] - bead_info[beadID][POS_Y];
  naTempR[POS_Z] = 2*tmpR[POS_Z] - bead_info[beadID][POS_Z];
  //Adjusting for periodic boundaries
  naTempR[POS_Y] = (naTempR[POS_Y] + nBoxSize[POS_Y]) % nBoxSize[POS_Y];
  naTempR[POS_Z] = (naTempR[POS_Z] + nBoxSize[POS_Z]) % nBoxSize[POS_Z];
}

void Rot_X_270(int beadID, const int tmpR[POS_MAX]){

  int j;
  //Performs a rotation by 270 degress around x-axis. tmpR is the new origin.
  for(j=0; j<POS_MAX; j++){
    naTempR[j] = bead_info[beadID][j];
  }
  naTempR[POS_Y] =  bead_info[beadID][POS_Z] + tmpR[POS_Y] - tmpR[POS_Z];
  naTempR[POS_Z] =  tmpR[POS_Y] + tmpR[POS_Z] - bead_info[beadID][POS_Y];
  //Adjusting for periodic boundaries
  naTempR[POS_Y] = (naTempR[POS_Y] + nBoxSize[POS_Y]) % nBoxSize[POS_Y];
  naTempR[POS_Z] = (naTempR[POS_Z] + nBoxSize[POS_Z]) % nBoxSize[POS_Z];
}

void Rot_Y_90(int beadID, const int tmpR[POS_MAX]){
  int j;
  //Performs a rotation by 90 degress around x-axis. tmpR is the new origin.
  for(j=0; j<POS_MAX; j++){
    naTempR[j] = bead_info[beadID][j];
  }
  naTempR[POS_X] = (bead_info[beadID][POS_Z] - tmpR[POS_Z] + tmpR[POS_X]);
  naTempR[POS_Z] = (tmpR[POS_X] + tmpR[POS_Z] - bead_info[beadID][POS_X]);
  //Adjusting for periodic boundaries
  naTempR[POS_X] = (naTempR[POS_X] + nBoxSize[POS_X]) % nBoxSize[POS_X];
  naTempR[POS_Z] = (naTempR[POS_Z] + nBoxSize[POS_Z]) % nBoxSize[POS_Z];
}

void Rot_Y_180(int beadID, const int tmpR[POS_MAX]){
  int j;
  //Performs a rotation by 270 degress around x-axis. tmpR is the new origin.
  for(j=0; j<POS_MAX; j++){
    naTempR[j] = bead_info[beadID][j];
  }
  naTempR[POS_X] = (2*tmpR[POS_X] - bead_info[beadID][POS_X]);
  naTempR[POS_Z] = (2*tmpR[POS_Z] - bead_info[beadID][POS_Z]);
  //Adjusting for periodic boundaries
  naTempR[POS_X] = (naTempR[POS_X] + nBoxSize[POS_X]) % nBoxSize[POS_X];
  naTempR[POS_Z] = (naTempR[POS_Z] + nBoxSize[POS_Z]) % nBoxSize[POS_Z];
}

void Rot_Y_270(int beadID, const int tmpR[POS_MAX]){
  int j;
  //Performs a rotation by 270 degress around x-axis. tmpR is the new origin.
  for(j=0; j<POS_MAX; j++){
    naTempR[j] = bead_info[beadID][j];
  }
  naTempR[POS_X] = (tmpR[POS_X] + tmpR[POS_Z] - bead_info[beadID][POS_Z]);
  naTempR[POS_Z] = (bead_info[beadID][POS_X]- tmpR[POS_X] + tmpR[POS_Z]);
  //Adjusting for periodic boundaries
  naTempR[POS_X] = (naTempR[POS_X] + nBoxSize[POS_X]) % nBoxSize[POS_X];
  naTempR[POS_Z] = (naTempR[POS_Z] + nBoxSize[POS_Z]) % nBoxSize[POS_Z];
}

void Rot_Z_90(int beadID, const int tmpR[POS_MAX]){
  int j;
  //Performs a rotation by 270 degress around x-axis. tmpR is the new origin.
  for(j=0; j<POS_MAX; j++){
    naTempR[j] = bead_info[beadID][j];
  }
  naTempR[POS_X] = (tmpR[POS_X] + tmpR[POS_Y] - bead_info[beadID][POS_Y]);
  naTempR[POS_Y] = (bead_info[beadID][POS_X] - tmpR[POS_X] + tmpR[POS_Y]);
  //Adjusting for periodic boundaries
  naTempR[POS_X] = (bead_info[beadID][POS_X] + nBoxSize[POS_X]) % nBoxSize[POS_X];
  naTempR[POS_Y] = (bead_info[beadID][POS_Y] + nBoxSize[POS_Y]) % nBoxSize[POS_Y];
}

void Rot_Z_180(int beadID, const int tmpR[POS_MAX]){
  int j;
  //Performs a rotation by 270 degress around x-axis. tmpR is the new origin.
  for(j=0; j<POS_MAX; j++){
    naTempR[j] = bead_info[beadID][j];
  }
  naTempR[POS_X] = (2*tmpR[POS_X] - bead_info[beadID][POS_X]);
  naTempR[POS_Y] = (2*tmpR[POS_Y] - bead_info[beadID][POS_Y]);
  //Adjusting for periodic boundaries
  naTempR[POS_X] = (bead_info[beadID][POS_X] + nBoxSize[POS_X]) % nBoxSize[POS_X];
  naTempR[POS_Y] = (bead_info[beadID][POS_Y] + nBoxSize[POS_Y]) % nBoxSize[POS_Y];
}

void Rot_Z_270(int beadID, const int tmpR[POS_MAX]){
  int j;
  //Performs a rotation by 270 degress around x-axis. tmpR is the new origin.
  for(j=0; j<POS_MAX; j++){
    naTempR[j] = bead_info[beadID][j];
  }
  naTempR[POS_X] = (bead_info[beadID][POS_Y]  + tmpR[POS_X] - tmpR[POS_Y]);
  naTempR[POS_Y] = (tmpR[POS_X] + tmpR[POS_Y] - bead_info[beadID][POS_X]);
  //Adjusting for periodic boundaries
  naTempR[POS_X] = (bead_info[beadID][POS_X] + nBoxSize[POS_X]) % nBoxSize[POS_X];
  naTempR[POS_Y] = (bead_info[beadID][POS_Y] + nBoxSize[POS_Y]) % nBoxSize[POS_Y];
}

void RotOperation(int PivotM, int beadID, int tmpR[POS_MAX]){
  int j;
  switch(PivotM){
    case 1:
    Rot_X_90(beadID, tmpR);
    break;

    case 2:
    Rot_X_180(beadID, tmpR);
    break;

    case 3:
    Rot_X_270(beadID, tmpR);
    break;

    case 4:
    Rot_Y_90(beadID, tmpR);
    break;

    case 5:
    Rot_Y_180(beadID, tmpR);
    break;

    case 6:
    Rot_Y_270(beadID, tmpR);
    break;

    case 7:
    Rot_Z_90(beadID, tmpR);
    break;

    case 8:
    Rot_Z_180(beadID, tmpR);
    break;

    case 9:
    Rot_Z_270(beadID, tmpR);
    break;

    default:
    //printf("How\n");
    for(j=0; j<POS_MAX; j++){
      naTempR[j] = bead_info[beadID][j]; //Does nothing; should make the move fail.
    }
    break;
  }



}

void ShuffleRotIndecies(void){
  //Algorithm taken from
  //https://www.w3resource.com/c-programming-exercises/array/c-array-exercise-77.php
  //Takes the array and shuffles the numbers in it. Used to randomly, but fully, sample the
  //rotational states.

  int i, j;
  int i_val, j_val;

    for (i=MAX_ROTSTATES-2; i>0; i--)
    {
        j = rand() % (i+1);
        i_val = Rot_IndArr[i];
        j_val = Rot_IndArr[j];
        Rot_IndArr[i] = j_val;
        Rot_IndArr[j] = i_val;
    }
}

/*int CheckRotStatesOld(int beadID, int resi, float MyTemp){

  int i, j, k, tmpBead;
  int tmpR[POS_MAX], tmpR2[POS_MAX];
  for(j=0; j<POS_MAX; j++){
    tmpR[j] = bead_info[beadID][j];
    }
  for(k=0; k<MAX_ROTSTATES-1; k++){
    i = Rot_IndArr[k];
    for(j=0; j<POS_MAX; j++){
      tmpR2[j] = (tmpR[j] + LocalArr[i][j] + nBoxSize[j]) % nBoxSize[j];
      }
      tmpBead = LtIndV(tmpR2);
      tmpBead = naTotLattice[tmpBead];
      if (tmpBead != -1){
        j = bead_info[tmpBead][BEAD_TYPE];
        if (fEnergy[resi][j][E_SC_SC] != 0 && (bead_info[tmpBead][BEAD_FACE] == -1 || bead_info[tmpBead][BEAD_FACE] == beadID)){
            //bolt_fac[k] = exp(-fEnergy[resi][j][E_SC_SC] / MyTemp);
            bolt_fac[k] = dbias_bolt_fac[resi][j];
            rot_trial[0][k] = tmpBead;
          }
        else{
            bolt_fac[k] = fRot_Bias;
            rot_trial[0][k] = -1;
          }
        }
      else{
          bolt_fac[k] = fRot_Bias;
          rot_trial[0][k] = -1;
        }
    }

  return k;

}*/

int CheckRotStatesOld(int beadID, int resi, float MyTemp){

    int i, j, k, tmpBead;
    int CandNums = 0;
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    for(j=0; j<POS_MAX; j++){
        tmpR[j] = bead_info[beadID][j];
    }
    for(k=0; k<MAX_ROTSTATES-1; k++){
        i = Rot_IndArr[k];
        for(j=0; j<POS_MAX; j++){
            tmpR2[j] = (tmpR[j] + LocalArr[i][j] + nBoxSize[j]) % nBoxSize[j];
        }
        tmpBead = LtIndV(tmpR2);
        tmpBead = naTotLattice[tmpBead];
        if (tmpBead != -1){
            j = bead_info[tmpBead][BEAD_TYPE];
            if (fEnergy[resi][j][E_SC_SC] != 0 && (bead_info[tmpBead][BEAD_FACE] == -1 || bead_info[tmpBead][BEAD_FACE] == beadID)){
                bolt_fac[CandNums] = dbias_bolt_fac[resi][j];
                rot_trial[0][CandNums] = tmpBead;
                CandNums++;
            }
        }
    }

    return CandNums;

}

/*int CheckRotStatesNew(int beadID, int resi, float MyTemp){

  int i, j, k, tmpBead;
  int tmpR[POS_MAX], tmpR2[POS_MAX];
  for(j=0; j<POS_MAX; j++){
    tmpR[j] = bead_info[beadID][j];
    }
  for(k=0; k<MAX_ROTSTATES-1; k++){
    i = Rot_IndArr[k];
    for(j=0; j<POS_MAX; j++){
      tmpR2[j] = (tmpR[j] + LocalArr[i][j] + nBoxSize[j]) % nBoxSize[j];
      }
      tmpBead = LtIndV(tmpR2);
      tmpBead = naTotLattice[tmpBead];
      if (tmpBead != -1){
        j = bead_info[tmpBead][BEAD_TYPE];
        if (fEnergy[resi][j][E_SC_SC] != 0 && (bead_info[tmpBead][BEAD_FACE] == -1 || bead_info[tmpBead][BEAD_FACE] == beadID)){
            //bolt_fac[k] = exp(-fEnergy[resi][j][E_SC_SC] / MyTemp);
            bolt_fac[k] = dbias_bolt_fac[resi][j];
            rot_trial[0][k] = tmpBead;
          }
        else{
            bolt_fac[k] = fRot_Bias;
            rot_trial[0][k] = -1;
          }
        }
      else{
          bolt_fac[k] = fRot_Bias;
          rot_trial[0][k] = -1;
        }
    }

  return k;

}*/

int CheckRotStatesNew(int beadID, int resi, float MyTemp){

    int i, j, k, tmpBead;
    int CandNums = 0;
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    for(j=0; j<POS_MAX; j++){
        tmpR[j] = bead_info[beadID][j];
    }
    for(k=0; k<MAX_ROTSTATES-1; k++){
        i = Rot_IndArr[k];
        for(j=0; j<POS_MAX; j++){
            tmpR2[j] = (tmpR[j] + LocalArr[i][j] + nBoxSize[j]) % nBoxSize[j];
        }
        tmpBead = LtIndV(tmpR2);
        tmpBead = naTotLattice[tmpBead];
        if (tmpBead != -1){
            j = bead_info[tmpBead][BEAD_TYPE];
            if (fEnergy[resi][j][E_SC_SC] != 0 && (bead_info[tmpBead][BEAD_FACE] == -1 || bead_info[tmpBead][BEAD_FACE] == beadID)){
                bolt_fac[CandNums] = dbias_bolt_fac[resi][j];
                rot_trial[0][CandNums] = tmpBead;
                CandNums++;
            }
        }
    }

    return CandNums;

}

void NormalizeRotState(int beadVal, int CandNums){
  int i;
  bolt_norm[beadVal] = 0.;
  if (CandNums > 0){//There is a possible candidate, so normalize bolt_fac
      for(i=0; i<CandNums; i++){
          bolt_norm[beadVal] += bolt_fac[i];
      }
      for(i=0; i<CandNums; i++){
          bolt_fac[i] /= bolt_norm[beadVal];
      }
      for(i=1; i<CandNums; i++){
          bolt_fac[i] += bolt_fac[i - 1];
      }
  }
    bolt_norm[beadVal] += 10.;
}

int PickRotState(int CandNums){
    int newRot = -1;
    int i;
    float fProb;
    int nCheck = CandNums + 1;
    nCheck = rand() % nCheck;
    if (nCheck == 0){
        newRot = -1;
    }
    else{
        fProb = (float)rand()/(float)RAND_MAX;
        for(i=0; i<CandNums; i++){
            if (fProb < bolt_fac[i]){
                break;
            }
        }

        newRot = rot_trial[0][i];
    }
    return newRot;
}