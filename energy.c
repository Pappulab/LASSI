#include "global.h"
#include "energy.h"
#include "structure.h"

void Energy_Total_System() {
  int i;//Indecies
  // initialization
  for (i=0; i<MAX_E; i++) {
    faCurrEn[i] = 0.0;
  }
//printf("We have %d beads and %d chains", tot_beads, tot_chains);
  for (i=0; i < tot_beads; i++) {
    faCurrEn[E_OVLP] += Energy_Isotropic(i) / 2.;
    faCurrEn[E_SC_SC] += Energy_Anisotropic(i) / 2.;
    //printf("Done with bead %d\n", i);
  }

  for (i=1; i<MAX_E; i++) {
    faCurrEn[E_TOT] += faCurrEn[i];
  }
}

float Energy_Anisotropic(int beadID){//Calculates the SC-SC energy of the bead in question.
  float totEn = 0.0; //Storing total overlap energy
    if (bead_info[beadID][BEAD_FACE] != -1){
        totEn += fEnergy[bead_info[beadID][BEAD_TYPE]][bead_info[bead_info[beadID][BEAD_FACE]][BEAD_TYPE]][E_SC_SC];
        }
  return totEn;
}

float Energy_Isotropic(int beadID){//Calculate Contact and Overlap energy of bead i
    float totEn = 0.0; //Storing total overlap energy
    int i, j;//Indecies
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    int x, y, z; //Lattice indecies
    int BoxRad = 1;
    int secBi, resj;//Second bead index
    float xDis = 0.;//Distance between beads.
    int resi  = bead_info[beadID][BEAD_TYPE];
    totEn += nThermalization_Mode == 0 ? 0. : Energy_InitPotential(beadID);


    if (nBeadTypeCanOvlp[resi] == 0){
        return totEn;
    }//No need to do antying if there's no overlap cost.


//Going through possible neighbors For now going through (-1,1) neighbors so CONT is a small window
  for (j=0;j<POS_MAX;j++){
    tmpR[j]  = bead_info[beadID][j];
    }
  for (x = -BoxRad; x <= BoxRad; x++){
    for (y = -BoxRad; y <= BoxRad; y++){
      for (z = -BoxRad; z <= BoxRad; z++){
          tmpR2[0] = (tmpR[0] + x + nBoxSize[0]) % nBoxSize[0];
          tmpR2[1] = (tmpR[1] + y + nBoxSize[1]) % nBoxSize[1];
          tmpR2[2] = (tmpR[2] + z + nBoxSize[2]) % nBoxSize[2];
          secBi = naTotLattice[Lat_Ind_FromVec(tmpR2)];
          if (secBi != -1 && secBi != beadID){
            resj  = bead_info[secBi][BEAD_TYPE];
            xDis  = Dist_BeadToBead(beadID, secBi);
            //totEn += fEnergy[resi][resj][E_OVLP]/xDis/xDis;
            if (xDis <= 1.0) {
              totEn += fEnergy[resi][resj][E_OVLP] / 2.0;
            } else if (xDis <= 1.42) { // sqrt(2)
              totEn += fEnergy[resi][resj][E_OVLP] / 4.0;
            } else if (xDis <= 1.74) { // sqrt(3)
              totEn += fEnergy[resi][resj][E_OVLP] / 8.0;
            }
            /*else if (xDis <= fEnRad[resi][resj][E_CONT]*1.74){
              totEn += fEnergy[resi][resj][E_CONT];
            //}//This way contact and overlap do overlap
          }*/
          }
          //TODO: Add option for solvent interactions in the parfile
          /*if (secBi == -1 && fSolEnergy != 0.){
              totEn += fSolEnergy*(fCuTemp-fThetaTemp);
          }*/
        }
      }
    }


  return totEn;

}

float Energy_Of_Chain(int chainID){//Calculates the energy of the given chain
  float totEn = 0.0;
  int i;//Indecies

  for(i=chain_info[chainID][CHAIN_START];i<chain_info[chainID][CHAIN_START]+chain_info[chainID][CHAIN_LENGTH];i++){
    totEn += Energy_Anisotropic(i) + Energy_Isotropic(i);
  }

  return totEn;
}

float Energy_InitPotential(int beadID){

    int j;
    float totEn = 0.;
    int tmpR[POS_MAX];
    if (fCuTemp-fKT > 0.005) {
        switch (nThermalization_Mode) {
            case 1:
                    for (j = 0; j < POS_MAX; j++) {
                        tmpR[j] = bead_info[beadID][j];
                        tmpR[j] = tmpR[j] - nBoxSize[j] / 2;
                        totEn += (float) (tmpR[j] * tmpR[j]);
                    }
                    totEn = (fCuTemp - fKT) * totEn;

                break;

            case 2:
                    for (j = 0; j < POS_MAX; j++) {
                        tmpR[j] = bead_info[beadID][j];
                        tmpR[j] = tmpR[j] - nBoxSize[j] / 2;
                        totEn += (float) (tmpR[j] * tmpR[j]);
                    }
                    if (totEn >= 2500) {
                        totEn = fKT * ((float) tot_beads + totEn);
                    } else if (totEn <= 2000) {
                        totEn = fKT * ((float) tot_beads + 1. / (totEn + 0.02));
                    }
                break;
            default:
                totEn = 0.;
                break;
        }
    }
    else{
        nThermalization_Mode = 0;
    }

    return totEn;
}
