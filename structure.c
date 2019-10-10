#include "global.h"
#include "structure.h"

int LtInd(int i, int j, int k){//Lattice index from 3D to 1D array
  return i + nBoxSize[0]*(j + nBoxSize[1]*k);
}

int LtIndV(int xArr[POS_MAX]){//Just the vector form of the above function. Easier to read sometimes.
      return xArr[POS_X] + nBoxSize[POS_X]*(xArr[POS_Y] + nBoxSize[POS_Y]*xArr[POS_Z]);
}

float distf(float f1[POS_MAX], float f2[POS_MAX]) {
  float d[POS_MAX];
  int i;
  for (i=0; i<POS_MAX; i++) {
    d[i] = fabsf(f1[i] - f2[i]);
    d[i] = d[i] > (float)nBoxSize[i]/2. ? (float)nBoxSize[i] - d[i] : d[i];
    //f[i] = (fabsf(f1[i] - f2[i]) < nBoxSize[i] - fabsf(f1[i] - f2[i])) ? fabsf(f1[i] - f2[i]) : (nBoxSize[i] - fabsf(f1[i] - f2[i]));
  }

  return sqrtf(d[POS_X]*d[POS_X] + d[POS_Y]*d[POS_Y] + d[POS_Z]*d[POS_Z]);
}

float distInt(int f1[POS_MAX], int f2[POS_MAX]) {
  int d[POS_MAX];
  int i;
  for (i=0; i<POS_MAX; i++) {
    d[i] = abs(f1[i] - f2[i]);
    d[i] = d[i] > nBoxSize[i]/2 ? nBoxSize[i] - d[i] : d[i];
    //f[i] = (abs(f1[i] - f2[i]) < nBoxSize[i] - abs(f1[i] - f2[i])) ? abs(f1[i] - f2[i]) : (nBoxSize[i] - abs(f1[i] - f2[i]));
    //f[i] = f[i];
  }

  return sqrtf((float)(d[POS_X]*d[POS_X] + d[POS_Y]*d[POS_Y] + d[POS_Z]*d[POS_Z]));
}

float distBeadToVec(int beadID, int f1[POS_MAX]) {
    int d[POS_MAX];
    int i;
    for(i=0; i<POS_MAX;i++) {
        d[i] = abs(bead_info[beadID][i] - f1[i]);
        d[i] = d[i] > nBoxSize[i]/2 ? nBoxSize[i] - d[i] : d[i];
        //d[i] = (abs(bead_info[beadID][i] - f1[i]) < nBoxSize[i] - abs(bead_info[beadID][i] - f1[i])) ? abs(bead_info[beadID][i] - f1[i]) : (nBoxSize[i] - abs(bead_info[beadID][i] - f1[i]));
    }
    return sqrtf((float)(d[POS_X]*d[POS_X] + d[POS_Y]*d[POS_Y] + d[POS_Z]*d[POS_Z]));
}

float dist(int n1, int n2) {
  lInt d[POS_MAX];
  lInt i;

  for (i=0; i<POS_MAX; i++) {
    d[i] = abs(bead_info[n1][i] - bead_info[n2][i]);
    d[i] = d[i] > nBoxSize[i]/2 ? nBoxSize[i] - d[i] : d[i];
    //d[i] = (abs(bead_info[n1][i] - bead_info[n2][i]) < nBoxSize[i] - abs(bead_info[n1][i] - bead_info[n2][i])) ? abs(bead_info[n1][i] - bead_info[n2][i]) : (nBoxSize[i] - abs(bead_info[n1][i] - bead_info[n2][i]));
  }

  return sqrtf((float)(d[POS_X]*d[POS_X] + d[POS_Y]*d[POS_Y] + d[POS_Z]*d[POS_Z]));
}

int check_structure_topo(void){
  int i, j, k;//Looping variables
  int idx, idy;//Internal iterators for covalent bonds.
  int tmpR[POS_MAX];//Just a vector to store coordinates
  int bondPart;
  for(i=0; i < tot_beads; i++){
  idx = 0;
  bondPart = topo_info[i][idx];
  for(j=0; j<POS_MAX; j++){
    tmpR[j] = bead_info[i][j];
  }
  if(naTotLattice[LtIndV(tmpR)] == -1){
      printf("Lattice Position for bead %d is empty! Chain: %d\n", i, bead_info[i][BEAD_CHAINID]);
      return i+1;
  }
  if(i - naTotLattice[LtIndV(tmpR)] != 0){//This means there is a mismatch between where the bead is and where the lattice thinks the bead is
    printf("Bead position and lattice value not the same. Crashing\t\t");
    printf("B1:%d B2:%d\t C1:%d C2:%d\n", i, naTotLattice[LtIndV(tmpR)], bead_info[i][BEAD_CHAINID], bead_info[naTotLattice[LtIndV(tmpR)]][BEAD_CHAINID]);
    return i+1;
  }
  while(topo_info[i][idx] != -1 && idx < MAX_BONDS){
    bondPart = topo_info[i][idx];
    if(dist(i, bondPart) > 1.74 * linker_len[i][idx]){
      printf("Bad beads! %d\t(%d %d %d)\t\tTopo:(%d %d %d)\t\tLinkers:(%.5f\t%.5f\t%.5f)\n", i, bead_info[i][0], bead_info[i][1], bead_info[i][2], topo_info[i][0], topo_info[i][1], topo_info[i][2], (float)linker_len[i][0], (float)linker_len[i][1], (float)linker_len[i][2]);
      printf("\t\t\t\t\t-------------------->\t\t%f\tSHOULD BE\t%f\n", dist(i,bondPart), 1.74*(float)linker_len[i][idx]);
      printf("Bad beads! %d\t(%d %d %d)\t\tTopo:(%d %d %d)\t\tLinkers:(%.5f\t%.5f\t%.5f)\n\n", bondPart, bead_info[bondPart][0], bead_info[bondPart][1], bead_info[bondPart][2], topo_info[bondPart][0], topo_info[bondPart][1], topo_info[bondPart][2], (float)linker_len[bondPart][0], (float)linker_len[bondPart][1], (float)linker_len[bondPart][2]);
      return i+1;
    }
    idx++;
  }
  if(bead_info[i][BEAD_FACE] != -1){
    if(nBeadTypeIsSticker[bead_info[i][BEAD_TYPE]] != 1){
      printf("This bead -- %d -- should not have a bond. Crashing.\n", i );
      return i+1;
      }
      if(bead_info[i][BEAD_FACE] == i){
            printf("Self bonded.\n");
            return i+1;
      }
    if(i != bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE]){
      printf("Bad bond!\n\t%d %d %d %f\nCrashing.\n", i , bead_info[i][BEAD_FACE], bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE], fEnergy[bead_info[i][BEAD_TYPE]][bead_info[bead_info[i][BEAD_FACE]][BEAD_TYPE]][E_SC_SC]);
      return i+1;
    }
    if(dist(i, bead_info[i][BEAD_FACE]) > 1.74){
      printf("Bad bond! Distance is wrong\n\t%d %d %d\nCrashing.\n", i , bead_info[i][BEAD_FACE], bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE]);
      return i+1;
    }
    }
  }
    return 0;
}

void radial_distribution_tot(void){
  float x;  //For distance
  int i, j;
  float sph_norm = tot_beads * tot_beads * PI * 4.0 / 3.0 / nBoxSize[0] / nBoxSize[1] / nBoxSize[2]; //Volume of sphere constant, and number normalization

  //Initializing
  for (i=0;i<RDF_MAXBINS; i++){
      fRDF_TOT[i] = 0.0;
  }

  //Calculating where, and how many, pairs exist
  for (i=0; i < tot_beads; i++){
    for (j=i+1; j < tot_beads; j++){
      x = dist(i,j);
      if (x <= (float)nBoxSize[0]/2.0){// Obviously, the maximal radius can be N/2
        fRDF_TOT[(int)x * 4] += 2.0; //Adding a pair to that bin I am assuming for now that dr=1/4
      }
    }
  }

  //Normalization over a sphere
  for (i = 0; i < nBins_RDF; i++){
    x = 1.0 + (float)i*(3.0 + 3.0*(float)i);//Volume of shell between i and i+1
    x *= sph_norm;
    x /= 64.0; // Because (1/4)^(-3) is 64
    fRDF_TOT[i] /= x;
  }

}

void radial_distribution_split(void){//Only g(r) is normalized.

  float x;  //For distance
  int i, j, k;
  int resi, resj;
  int chaini, chainj;
  int myBin = 0;


  //Initializing

  for (j = 0; j < 7; j++){
    for (i = 0; i < RDF_MAXBINS; i++){
        ldRDF_ARR[j][i] = 0.;
    }
  }

  //Calculating where, and how many, pairs exist
  for (i=0; i < tot_beads; i++){
        resi = bead_info[i][BEAD_TYPE];
    for (j=i+1; j < tot_beads; j++){
      resj = bead_info[j][BEAD_TYPE];
      x = dist(i,j);// Calculate distance between beads
      //Calculate which bin:= floor(1/dr * x)
      myBin = (int)floor(4.*x);//I am assuming for now that dr=1/4
      //Note that the dist(i,j) ensures that no incorrect distances are calculated.
      ldRDF_ARR[0][myBin] += 2.0;  //Adding a pair to that bin
      if (resi == 0 && resj == 0){
          ldRDF_ARR[1][myBin] += 2.0;
        continue;
      }
      if ((resi == 0 && resj == 1) || (resj == 0 && resi == 1)){
          ldRDF_ARR[2][myBin] += 2.0;
        continue;
      }
      if ((resi == 0 && resj == 2) || (resj == 0 && resi == 2)){
          ldRDF_ARR[3][myBin] += 2.0;
        continue;
      }
      if ((resi == 1 && resj == 1)){
          ldRDF_ARR[4][myBin] += 2.0;
        continue;
      }
      if ((resi == 1 && resj == 2) || (resj == 1 && resi == 2)){
          ldRDF_ARR[5][myBin] += 2.0;
        continue;
      }
      if ((resi == 2 && resj == 2)){
          ldRDF_ARR[6][myBin] += 2.0;
      }
  }
  }



}

float vectorMag(const int MyVec[]){//Outputs the magnitude of the vector
  return sqrtf((float)(MyVec[0]*MyVec[0]+MyVec[1]*MyVec[1]+MyVec[2]*MyVec[2]));
}

void CalcGyrTensor(int ClusSize, int ClusIndex){
  //Calculate the components of the gyration tensor for a given cluster.
  //ClusSize is the size of the cluster -- obviously -- whereas ClusIndex tell us
  //where in naCluster the chain indecies are located. naCluster[ClusIndex][0-ClusSize] is all the chainID's I need
  //for the calculation
  //Remember that the Gyration Tensor is a 3x3 symmetric object so we only need 6 numbers.
    int i, k, j, j2;//Basic indecies for loops
    for(i=0;i<7;i++){ fGyrTensor[i] = 0.;}//Initializing
    int firstB, lastB;//Tracks the first and last bead of the given chain
    float tot_COM[POS_MAX] = {0.};//This is where we shall store the COM of the cluster.
    int NumRes = 0;//Tracks how many residues are in this cluster

  //The only thing one needs to be careful about is to take PBC into account; the rest is tedium.
  //Use good old differential geometry to map each coordinate to two new coordinates, and unpack at the end.
    float theta[POS_MAX] = {0.};
    float zeta[POS_MAX]  = {0.};//Extra coordinates for goodness
    float dumArg = 0.;//Just a dummy variable to be more efficient
    float dumArg2= 0.;//Another one
    // printf("Starting with COM\n");
    //Calculating the COM
    for(i=0; i<ClusSize; i++){
        firstB = chain_info[naCluster[ClusIndex][i]][CHAIN_START];
        lastB  = firstB + chain_info[naCluster[ClusIndex][i]][CHAIN_LENGTH];
       // printf("%d %d\n", firstB, lastB);
        //Just easier to track each chain like this
        for(k=firstB; k<lastB; k++){
            NumRes++;//Adding a residue to the total
            for(j=0; j<POS_MAX; j++){
                dumArg = 2.*PI*((float)bead_info[k][j]/(float)nBoxSize[j]);
            theta[j] += cosf(dumArg);
            zeta[j]  += sinf(dumArg);//Since I am taking the average just keep adding
            }
        }
       // printf("Next bead\n");
    }
    //printf("Done with COM\n");
    //Calculating average, and then COM
    for(j=0; j<POS_MAX; j++){
     theta[j] = theta[j]/(float)NumRes;
     zeta[j]  = zeta[j]/(float)NumRes;
     tot_COM[j] = atan2f(-theta[j],-zeta[j])+PI;
     tot_COM[j] = nBoxSize[j]*(tot_COM[j]/2./PI);
    }

  //Using the COM to calculate the Gyration Tensor
  // GyrTen_{ij} = 1/N sum_1^N (r_i-com_i)(r_j-com_j) so just take the sums and divide at the end
    for(i=0; i<ClusSize; i++){
        firstB = chain_info[naCluster[ClusIndex][i]][CHAIN_START];
        lastB  = firstB + chain_info[naCluster[ClusIndex][i]][CHAIN_LENGTH];
        //Just easier to track each chain like this
        for(k=firstB; k<lastB; k++){
            for(j=0; j<POS_MAX; j++){
                for(j2=j; j2<POS_MAX; j2++){
                dumArg = fabsf((float)bead_info[k][j] - tot_COM[j]) < (float)nBoxSize[j] - fabsf((float)bead_info[k][j] - tot_COM[j]) ? fabsf((float)bead_info[k][j] - tot_COM[j]) : (float)nBoxSize[j] - fabsf((float)bead_info[k][j] - tot_COM[j]);
                dumArg2 = fabsf((float)bead_info[k][j2] - tot_COM[j2]) < (float)nBoxSize[j2] - fabsf((float)bead_info[k][j2] - tot_COM[j2]) ? fabsf((float)bead_info[k][j2] - tot_COM[j2]) : (float)nBoxSize[j2] - fabsf((float)bead_info[k][j2] - tot_COM[j2]);
                    fGyrTensor[j + 3 * (j2 - j)] += dumArg * dumArg2;
                //0 = xx; 1 = yy; 2 = zz; 3 = xy; 6 = xz; 4 = yz; Need smarter indexing
                }
            }

        }
    }

    //Calculating the average;
    for(i=0;i<7;i++){ fGyrTensor[i] /= (float)NumRes; }//printf("%f\n",fGyrTensor[i]);}printf("\n");



   // exit(1);

}

void CalcTotGyrTensor(void){//Calculates the gyration tensor for the whole system
    int i, k, j, j2;//Basic indecies for loops
    for(i=0;i<7;i++){ fGyrTensor[i] = 0.;}//Initializing
    float tot_COM[POS_MAX] = {0.};//This is where we shall store the COM of the cluster.

    //The only thing one needs to be careful about is to take PBC into account; the rest is tedium.
  //Use good old differential geometry to map each coordinate to two new coordinates, and unpack at the end.
    float theta[POS_MAX] = {0.};
    float zeta[POS_MAX]  = {0.};//Extra coordinates for goodness
    float dumArg = 0.;//Just a dummy variable to be more efficient
    float dumArg2= 0.;//Another one

    for(i=0; i < tot_beads; i++){
        for(j=0;j<POS_MAX;j++){
            tot_COM[j]+=bead_info[i][j];
        }
    }

    //Calculating average, and then COM
    for(j=0; j<POS_MAX; j++){
     tot_COM[j] /= (float)tot_beads;
    }
  //printf("\n");
  //Using the COM to calculate the Gyration Tensor
  // GyrTen_{ij} = 1/N sum_1^N (r_i-com_i)(r_j-com_j) so just take the sums and divide at the end
    for(i=0; i < tot_beads; i++){
            for(j=0; j<POS_MAX; j++){
                dumArg = (float)bead_info[i][j] - tot_COM[j];
                for(j2=j; j2<POS_MAX; j2++){
                dumArg2 = (float)bead_info[i][j2] - tot_COM[j2];
                    fGyrTensor[j + 3 * (j2 - j)] += dumArg * dumArg2;
                //printf("%.2f\t", dumArg2);
                //0 = xx; 1 = yy; 2 = zz; 3 = xy; 6 = xz; 4 = yz; Need smarter indexing
                }
            }//printf("\n");
        }
        //exit(1);
    //Calculating the average;
    for(i=0;i<7;i++){ fGyrTensor[i] /= (float)tot_beads; }//printf("%f\n",fGyrTensor[i]);}printf("\n");
}

void CalcTotGyrRad(void){
  /*
  Only calculates the diagonals of the gyration tensor, and calculates the sum of the
  diagonals. Remember that Rg^2 = Tr(GyrTen) so we only need to calculate the diagonals, and then then sum.
  I shall borrow most of the code from above, and so read CalcGyrTensor for what's happening here.
  */
  int i, j;//Loop indecies
  float tot_COM[POS_MAX] = {0.};
  float dumArg = 0.;
  for(i=0;i<7;i++){//Initializing to 0
    fGyrTensor[i] = 0.;
  }

  for(i=0; i < tot_beads; i++){
      for(j=0;j<POS_MAX;j++){
          tot_COM[j]+=bead_info[i][j];
      }
  }
  for(j=0; j<POS_MAX; j++){//Dividing by number of beads
   tot_COM[j] /= (float)tot_beads;
  }

  for(i=0; i < tot_beads; i++){
          for(j=0; j<POS_MAX; j++){
              dumArg = (float)bead_info[i][j] - tot_COM[j];
              fGyrTensor[j] += dumArg * dumArg;
          }
      }

    fSysGyrRad += sqrtf((fGyrTensor[0] + fGyrTensor[1] + fGyrTensor[2]) / (float)tot_beads);//Adding to the total fSysGyrRad to be averaged at the end.
      nTotGyrRadCounter++;//Remembering that we have calculated the radius; this will be used to average the final value.

}

void avg_rdf_split(void){
  /*
  Calcutes the RDF and adds it all up so that it can be averaged out at the end of the run.
  */
  float x;  //For distance
  int i, j, k;
  int resi, resj;
  //int chaini, chainj;
  int myBin = 0;
  int array_pos;

  //Calculating where, and how many, pairs exist
  for (i=0; i < tot_beads; i++){
        resi = bead_info[i][BEAD_TYPE];
    for (j=i+1; j < tot_beads; j++){
      resj = bead_info[j][BEAD_TYPE];
      x = dist(i,j);
      //Note that dist(i,j) automatically ensures no distance is greater than (L/2)*sqrt(3)
      myBin = (int)floor(4.*x);//I am assuming for now that dr=1/4
      ldRDF_ARR[0][myBin] += 2.0;  //Adding a pair to that bin
      array_pos = resi == resj ? 1 + resi : 6 + resj - resi*(resi-9)/2;
      ldRDF_ARR[array_pos][myBin] += 2.0;
    }
  }
  nrdfCounter++;
  }

int check_linker_constraint(int beadID, int tmpR[]){
  //Check if the proposed new location for beadID is such that all the linkers are unbroken.
  int idx;//Iterator to loop over bond Partners
  int bondPartner;//It is what it is.

  idx = 0; bondPartner = topo_info[beadID][idx];//Initializing the two.
  while(idx < MAX_BONDS && topo_info[beadID][idx] != -1){//Keep going till we run out of partners
    bondPartner = topo_info[beadID][idx];
    if(distInt(bead_info[bondPartner], tmpR) > 1.74*(float)linker_len[beadID][idx]){
      return 0;//This means that we have broken one of the linkers.
    }
    idx++;
  }
  return 1;//This means that all linker constraints are satisfied.
}

int ShakeConstraint(int beadID, int tmpR[MAX_VALENCY][POS_MAX]){

  int curID = beadID;
  int idx, bPart;
  int j;
  int topIt = 0;
  int canI  = 1;

  while (curID != -1){
    for(j=0; j<POS_MAX; j++){
    bead_info[curID][j] = tmpR[topIt][j];//Moving
    }
    curID = topo_info[beadID][topIt++];
  }

  curID = beadID; topIt = 0;
  while(curID != -1 && canI == 1){
    idx = 0; bPart = topo_info[curID][idx];
    while(bPart != -1 && idx < MAX_BONDS){
      if(dist(curID, bPart) > 1.74*(float)linker_len[curID][idx]){
              canI = 0;
              break;
      }
      bPart = topo_info[curID][++idx];
    }
    curID = topo_info[beadID][topIt++];
  }
  curID = beadID; topIt = 0;
  while (curID != -1) {
      for (j=0; j<POS_MAX; j++) {
          bead_info[curID][j] = old_bead[curID][j];//Moving back
      }
      curID = topo_info[beadID][topIt++];
  }

  return canI;
}

void CenterMySystem(void){
  /*
  This function finds the geometric center of mass of the system, and then displaces the entire system such that the new COM is nBoxSize/2. Note that this function does't care about the periodic boundaries
  and just centers the system. This makes looking at droplets easier in VMD, and will be used in ascertaining structural quantities later.
  */

  int com_sys[POS_MAX];//Storing the COM of the system.
  int d_vct[POS_MAX];//The vector to move the system.
  int i, j ,k;//Looping iterators!

  //Initialize the vectors!
  for(j=0; j<POS_MAX; j++){
    com_sys[j] = 0;
    d_vct[j] = 0;
  }
  //Finding the center of mass.
  for(i=0; i < tot_beads; i++){
    for(j=0; j<POS_MAX; j++){
      com_sys[j]+=bead_info[i][j];
    }
  }
  //Dividing by the number of beads!
  for(j=0; j<POS_MAX; j++){
    com_sys[j]=floorf((float)com_sys[j]/(float)tot_beads);
  }
  //Finding the displacemement vector
  for(j=0; j<POS_MAX; j++){
    d_vct[j] = nBoxSize[j]/2 - com_sys[j];
  }
  int tmpR[POS_MAX];//Vectors to store coordinates before and after.
  //Firstly let's move all the beads, and empty the lattice out.
  for(i=0; i < tot_beads; i++){
    for(j=0; j<POS_MAX; j++){
      tmpR[j] = bead_info[i][j];
      bead_info[i][j] = (tmpR[j] + d_vct[j] + nBoxSize[j]) % nBoxSize[j];
    }
      naTotLattice[LtIndV(tmpR)] = -1;//Removing every bead from where it was.
  }
  //Now that all the beads are in the right place. Fill the lattice.
  for(i=0; i < tot_beads; i++){
    for(j=0; j<POS_MAX; j++){
      tmpR[j] = bead_info[i][j];
    }
      naTotLattice[LtIndV(tmpR)] = i;//Placing everything back.
  }

}
