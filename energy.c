#include "global.h"
#include "energy.h"
#include "structure.h"

void Calc_Tot_Energy() {
  int i;//Indecies
  // initialization
  for (i=0; i<MAX_E; i++) {
    faCurrEn[i] = 0.0;
  }
//printf("We have %d beads and %d chains", tot_beads, tot_chains);
  for (i=0; i < tot_beads; i++) {
    faCurrEn[E_OVLP] += energy_cont_and_ovlp(i)/2.;
    faCurrEn[E_SC_SC] += energy_SC(i)/2.;
    //printf("Done with bead %d\n", i);
  }

  for (i=1; i<MAX_E; i++) {
    faCurrEn[E_TOT] += faCurrEn[i];
  }
}

float energy_SC(int i){//Calculates the SC-SC energy of the bead in question.
  float totEn = 0.0; //Storing total overlap energy
    if (bead_info[i][BEAD_FACE] != -1){
        totEn += fEnergy[bead_info[i][BEAD_TYPE]][bead_info[bead_info[i][BEAD_FACE]][BEAD_TYPE]][E_SC_SC];
        }
  return totEn;
}

float energy_cont_and_ovlp(int beadID){//Calculate Contact and Overlap energy of bead i
  float totEn = 0.0; //Storing total overlap energy
  int i, j;//Indecies
  int tmpR[POS_MAX], tmpR2[POS_MAX];//
  //TODO: Make an independent Indent function that does not clutter this function!
  if (Indent_Mode == 1){
      if (fCuTemp-fKT > 0.005) {
          for (j=0; j<POS_MAX; j++) {
              tmpR[j]  = bead_info[beadID][j];
              tmpR2[j] = tmpR[j] - nBoxSize[j] / 2;
              totEn   += (float) (tmpR2[j] * tmpR2[j]);
          }
          totEn = (fCuTemp - fKT) * totEn;
          //return totEn;
      }
  }
  if (Indent_Mode == 2){
      if (fCuTemp-fKT > 0.005) {
          for (j = 0; j < POS_MAX; j++) {
              tmpR[j] = bead_info[beadID][j];
              tmpR2[j] = tmpR[j] - nBoxSize[j] / 2;
              totEn += (float) (tmpR2[j] * tmpR2[j]);
          }
          if (totEn >= 2500) {
              //totEn = fKT*totEn;
              totEn = fKT * ((float) tot_beads + totEn);
              //return totEn;
          } else if (totEn <= 2000) {
              //totEn = fKT*totEn;
              totEn = fKT * ((float) tot_beads + 1. / (totEn + 0.02));
          }
      }
          else{
              totEn = 0.;
          }
  }

  int resi  = bead_info[beadID][BEAD_TYPE];
  if (TypeCanOvlp[resi] == 0){
    return totEn;
  }//No need to do antying if there's no overlap cost.
  int x, y, z; //Lattice indecies
  int BoxRad = 1;
  int secBi, resj;//Second bead index
  float xDis;//Distance between beads.

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
          secBi = naTotLattice[LtIndV(tmpR2)];
          if (secBi != -1 && secBi != beadID){
            resj  = bead_info[secBi][BEAD_TYPE];
            xDis  = dist(beadID, secBi);
            totEn += fEnergy[resi][resj][E_OVLP]/xDis/xDis;
            /*if (xDis <= 1.0) {
              totEn += fEnergy[resi][resj][E_OVLP] / 2.0;
            } else if (xDis <= 1.42) { // sqrt(2)
              totEn += fEnergy[resi][resj][E_OVLP] / 4.0;
            } else if (xDis <= 1.74) { // sqrt(3)
              totEn += fEnergy[resi][resj][E_OVLP] / 8.0;
            }*/
            /*else if (xDis <= fEnRad[resi][resj][E_CONT]*1.74){
              totEn += fEnergy[resi][resj][E_CONT];
            //}//This way contact and overlap do overlap
          }*/
          }
        }
      }
    }


  return totEn;

}

float energy_chain(int chainID){//Calculates the energy of the given chain
  float totEn = 0.0;
  int i;//Indecies

  for(i=chain_info[chainID][CHAIN_START];i<chain_info[chainID][CHAIN_START]+chain_info[chainID][CHAIN_LENGTH];i++){
    totEn += energy_SC(i) + energy_cont_and_ovlp(i);
  }

  return totEn;
}

