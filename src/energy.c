#include "global.h"
#include "energy.h"
#include "structure.h"

/// Energy_InitPotential calculates the biasing potential used in the thermalization/equilibration
/// The function calculates (T_current - T_final) and if < 0.001, sets nThermalization_Mode = 0.
/// \param beadID
/// \return The energy, given the nThermalization_Mode
float Energy_InitPotential(int beadID) {
    int j;
    float totEn = 0.;
    int tmpR[POS_MAX];
    if (fCuTemp - fKT > 0.005) {
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
                if (totEn >= 250.) {
                    totEn = (fCuTemp - fKT) * ( totEn);
                } else if (totEn <= 600.) {
                    totEn = (fCuTemp - fKT) * ((float) tot_beads + 1. / (totEn + 0.02));
                }
                else{
                    totEn=0.;
                }
                break;

            case 3:
                for (j = 0; j < POS_MAX; j++) {
                    tmpR[j] = bead_info[beadID][j];
                    tmpR[j] = tmpR[j] - nBoxSize[j] / 2;
                    totEn += (float) (tmpR[j] * tmpR[j]);
                }
                if (totEn >= 25*25) {
                    totEn = (fCuTemp) * totEn;
                }
                else {
                    totEn = 0.;
                }
                break;
            case 4:
                for (j = 0; j < POS_MAX; j++) {
                    tmpR[j] = bead_info[beadID][j];
                    tmpR[j] = tmpR[j] - nBoxSize[j] / 2;
                    totEn += (float) (tmpR[j] * tmpR[j]);
                }
                if (totEn <= 100.) {
                    totEn = (fCuTemp - fKT) * (totEn+0.2);
                }
                else{
                    totEn=0.;
                }
                break;
            case 5:
                for (j = 0; j < 1; j++) {
                    tmpR[j] = bead_info[beadID][j];
                    tmpR[j] = tmpR[j] - nBoxSize[j] / 2;
                    totEn += (float) (tmpR[j] * tmpR[j]);
                }
                totEn = (fCuTemp - fKT) * totEn;
                break;

            default:
                totEn = 0.;
                break;
        }
    } else {
        nThermalization_Mode = -1;
    }

    return totEn;
}

/// Energy_Anisotropic - returns the bond energy for beadID if it is bonded.
/// \param beadID
/// \return
float Energy_Anisotropic(int beadID) {//Calculates the SC-SC energy of the bead in question.
    float totEn = 0.0; //Storing total overlap energy
    if (bead_info[beadID][BEAD_FACE] != -1) {
        totEn += fEnergy[bead_info[beadID][BEAD_TYPE]][bead_info[bead_info[beadID][BEAD_FACE]][BEAD_TYPE]][E_SC_SC];
    }
    return totEn;
}

/// Energy_Anisotropic_Self - returns the bond energy for beadID only if it is bonded
/// to a bead in the same chain. Note that this does not account for double counting
/// \param beadID
/// \return
float Energy_Anisotropic_Self(int beadID) {//Calculates the SC-SC energy of the bead in question.
    float totEn = 0.0; //Storing total overlap energy
    int bP = bead_info[beadID][BEAD_FACE];
    if (bP != -1) {//This bead does have a bond
        if (bead_info[beadID][BEAD_CHAINID] == bead_info[bP][BEAD_CHAINID]) {
            totEn += fEnergy[bead_info[beadID][BEAD_TYPE]][bead_info[bP][BEAD_TYPE]][E_SC_SC];
        }
    }
    return totEn;
}

/// Energy_Anisotropic_For_Chain - returns the bond energy for beadID if it is bonded
/// If it is bonded to another in the same chain, divide the energy by 2.
/// This should only be used when calculating the energy of an entire chain!
/// \param beadID
/// \return
float Energy_Anisotropic_For_Chain(int beadID) {//Calculates the SC-SC energy of the bead in question.
    float totEn = 0.0; //Storing total overlap energy
    int bP = bead_info[beadID][BEAD_FACE];
    if (bP != -1) {//This bead does have a bond
        if (bead_info[beadID][BEAD_CHAINID] == bead_info[bP][BEAD_CHAINID]) {
            totEn += fEnergy[bead_info[beadID][BEAD_TYPE]][bead_info[bP][BEAD_TYPE]][E_SC_SC]/2.;
        }
        else{
            totEn += fEnergy[bead_info[beadID][BEAD_TYPE]][bead_info[bP][BEAD_TYPE]][E_SC_SC];
        }
    }
    return totEn;
}


/// Energy_Anisotropic_Contiguous_Range - returns the bond energy for beadID if it is bonded.
/// This function assumes that we get a beadID range [smallest_bead, largest_bead], and that
/// we will loop over every bead in this range. So if the bead partner is within this range,
/// we divide by 2 to account for double counting.
/// \param beadID
/// \return
float Energy_Anisotropic_Contiguous_Range(int beadID, int smallest_bead, int largest_bead){
    float totEn = 0.0; //Storing total overlap energy
    int bP = bead_info[beadID][BEAD_FACE];
    if (bP != -1) {//This bead does have a bond
        if (bead_info[beadID][BEAD_CHAINID] == bead_info[bP][BEAD_CHAINID]) {//Intra-molecular
            if (bP >= smallest_bead && bP <= largest_bead) {//Within subset
                totEn += fEnergy[bead_info[beadID][BEAD_TYPE]][bead_info[bP][BEAD_TYPE]][E_SC_SC] / 2.;
            }
            else{
                totEn += fEnergy[bead_info[beadID][BEAD_TYPE]][bead_info[bP][BEAD_TYPE]][E_SC_SC];
            }
        }
        else{
            totEn += fEnergy[bead_info[beadID][BEAD_TYPE]][bead_info[bP][BEAD_TYPE]][E_SC_SC];
        }
    }
    return totEn;
}

/// Energy_Anisotropic_With_List - returns the bond energy for beadID if it is bonded.
/// This function assumes that we get a list of beadID's of size list_size, but that the
/// last bead is at list_size-1. This also assumes that we will loop over the entire list
/// So if the bead partner is within this list, we divide by 2 to account for double counting.
/// \param beadID
/// \return
float Energy_Anisotropic_With_List(const int beadID, const int *bead_list, const int list_size){
    float totEn = 0.0; //Storing total overlap energy
    int bP = bead_info[beadID][BEAD_FACE];
    int check_bead = 0;
    int i;
    if (bP != -1) {//This bead does have a bond
        for (i=0; i<list_size; i++){//Checking if we should worry about double counting
            if (bP == bead_list[i]){
                check_bead = 1;
                break;
            }
        }
        if (check_bead == 1){
            totEn += fEnergy[bead_info[beadID][BEAD_TYPE]][bead_info[bP][BEAD_TYPE]][E_SC_SC] / 2.;
        }
        else{
            totEn += fEnergy[bead_info[beadID][BEAD_TYPE]][bead_info[bP][BEAD_TYPE]][E_SC_SC];
        }
    }
    return totEn;
}

/// Energy_Isotroptic_Old calculates the isotropic contribution to the energy by searching the 3^3-1 = 26 'neighbors'
/// The energy function  is like the BFM, where \f$\epislon$\f represents the overlap cost for total overlap, which is
/// forbidden explicitly in LASSI, so we have \f$\epislon/2$\f,\f$\epislon/4$\f and \f$\epislon/8$\f with increasing
/// distance.
/// \param beadID
/// \return
float Energy_Isotropic_Old(int beadID) {//Calculate Contact and Overlap energy of bead beadID
    float totEn = 0.0; //Storing total overlap energy
    int i, j;//Indecies
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    int x, y, z; //Lattice indecies
    int secBi, resj;//Second bead index
    float xDis = 0.;//Distance between beads.
    int resi = bead_info[beadID][BEAD_TYPE];
    totEn += nThermalization_Mode == -1 ? 0. : Energy_InitPotential(beadID);


    if (nBeadTypeCanOvlp[resi] == 0 && nBeadTypeCanCont[resi] == 0) {
        return totEn;
    }//No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont[resi] == 0 ? 1: 3;//No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++) {
        tmpR[j] = bead_info[beadID][j];
    }
    for (x = -BoxRad; x <= BoxRad; x++) {
        for (y = -BoxRad; y <= BoxRad; y++) {
            for (z = -BoxRad; z <= BoxRad; z++) {
                tmpR2[0] = (tmpR[0] + x + nBoxSize[0]) % nBoxSize[0];
                tmpR2[1] = (tmpR[1] + y + nBoxSize[1]) % nBoxSize[1];
                tmpR2[2] = (tmpR[2] + z + nBoxSize[2]) % nBoxSize[2];
                secBi = naTotLattice[Lat_Ind_FromVec(tmpR2)];
                if (secBi != -1 && secBi != beadID) {
                    resj = bead_info[secBi][BEAD_TYPE];
                    totEn += fEnergy[resi][resj][E_OVLP];
                    xDis = sqrtf((float)(x*x + y*y + z*z));
                    //xDis = Dist_BeadToBead(beadID, secBi);
                    if (xDis <= 1.0) {
                        totEn += fEnergy[resi][resj][E_OVLP] / 2.0;
                    } else if (xDis <= 1.42) { // sqrt(2)
                        totEn += fEnergy[resi][resj][E_OVLP] / 4.0;
                    } else if (xDis <= 1.74) { // sqrt(3)
                        totEn += fEnergy[resi][resj][E_OVLP] / 8.0;
                    }
                    else {//This way, contact is only outside ovlp
                        totEn += fEnergy[resi][resj][E_CONT] / xDis;
                  }
                }
            }
        }
    }

    return totEn;

}

/// Energy_Isotroptic calculates the isotropic contribution to the energy by searching the 3^3-1 = 26 'neighbors'
/// The energy function  is like the BFM, where \f$\epislon$\f represents the overlap cost for total overlap, which is
/// forbidden explicitly in LASSI, so we have \f$\epislon/2$\f,\f$\epislon/4$\f and \f$\epislon/8$\f with increasing
/// distance.
/// \param beadID
/// \return
float Energy_Isotropic(int beadID) {//Calculate Contact and Overlap energy of bead beadID
    float totEn = 0.0; //Storing total overlap energy
    int i, j;//Indecies
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    int x, y, z; //Lattice indecies
    int secBi, resj;//Second bead index
    float xDis = 0.;//Distance between beads.
    int resi = bead_info[beadID][BEAD_TYPE];
    totEn += nThermalization_Mode == -1 ? 0. : Energy_InitPotential(beadID);


    if (nBeadTypeCanOvlp[resi] == 0 && nBeadTypeCanCont[resi] == 0) {
        return totEn;
    }//No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont[resi] == 0 ? 1: 3;//No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++) {
        tmpR[j] = bead_info[beadID][j];
    }
    for (x = -BoxRad; x <= BoxRad; x++) {
        for (y = -BoxRad; y <= BoxRad; y++) {
            for (z = -BoxRad; z <= BoxRad; z++) {
                tmpR2[0] = (tmpR[0] + x + nBoxSize[0]) % nBoxSize[0];
                tmpR2[1] = (tmpR[1] + y + nBoxSize[1]) % nBoxSize[1];
                tmpR2[2] = (tmpR[2] + z + nBoxSize[2]) % nBoxSize[2];
                secBi = naTotLattice[Lat_Ind_FromVec(tmpR2)];
                if (secBi != -1 && secBi != beadID) {
                    resj = bead_info[secBi][BEAD_TYPE];
                    xDis = sqrtf((float)(x*x + y*y + z*z));
                    if (xDis <= 1.74) { // 1/r^3 potential
                        totEn += fEnergy[resi][resj][E_OVLP] / xDis / xDis / xDis;
                    }
                    // 1/r potential that goes till cube three
                    totEn += fEnergy[resi][resj][E_CONT] / xDis;
                }
            }
        }
    }

    return totEn;

}

/// Energy_Isotropic_Self - returns the bond energy for beadID only if it is interacting with a bead within the same
/// molecule. Note that it does not account for double counting
/// \param beadID
/// \return
float Energy_Isotropic_Self(int beadID) {//Calculate Contact and Overlap energy of beadID but only intra-molecular
    //interactions
    float totEn = 0.0; //Storing total overlap energy
    int i, j;//Indecies
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    int x, y, z; //Lattice indecies
    int secBi, resj;//Second bead index
    float xDis = 0.;//Distance between beads.
    int resi = bead_info[beadID][BEAD_TYPE];
    //totEn += nThermalization_Mode == -1 ? 0. : Energy_InitPotential(beadID);


    if (nBeadTypeCanOvlp[resi] == 0 && nBeadTypeCanCont[resi] == 0) {
        return totEn;
    }//No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont[resi] == 0 ? 1: 3;//No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++) {
        tmpR[j] = bead_info[beadID][j];
    }
    for (x = -BoxRad; x <= BoxRad; x++) {
        for (y = -BoxRad; y <= BoxRad; y++) {
            for (z = -BoxRad; z <= BoxRad; z++) {
                tmpR2[0] = (tmpR[0] + x + nBoxSize[0]) % nBoxSize[0];
                tmpR2[1] = (tmpR[1] + y + nBoxSize[1]) % nBoxSize[1];
                tmpR2[2] = (tmpR[2] + z + nBoxSize[2]) % nBoxSize[2];
                secBi = naTotLattice[Lat_Ind_FromVec(tmpR2)];
                if (secBi != -1 && secBi != beadID) {
                    if (bead_info[secBi][BEAD_CHAINID] == bead_info[beadID][BEAD_CHAINID]) {
                        resj = bead_info[secBi][BEAD_TYPE];
                        xDis = sqrtf((float) (x * x + y * y + z * z));
                        if (xDis <= 1.74) { // 1/r^3 potential
                            totEn += fEnergy[resi][resj][E_OVLP] ;/// xDis / xDis / xDis;
                        }
                        // 1/r potential that goes till cube three
                        totEn += fEnergy[resi][resj][E_CONT] / xDis;
                    }
                }
            }
        }
    }

    return totEn;

}


/// Energy_Isotropic_For_Chain - returns the isotropic interaction energy for beadID
/// If we have an intra-molecular interaction, divide the energy by 2.
/// This should only be used when calculating the energy of an entire chain!
/// \param beadID
/// \return
float Energy_Isotropic_For_Chain(int beadID) {//Calculate Contact and Overlap energy of beadID
    //Takes care of intra-molecular double counting
    float totEn = 0.0; //Storing total overlap energy
    int i, j;//Indecies
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    int x, y, z; //Lattice indecies
    int secBi, resj;//Second bead index
    float xDis = 0.;//Distance between beads.
    int resi = bead_info[beadID][BEAD_TYPE];
    totEn += nThermalization_Mode == -1 ? 0. : Energy_InitPotential(beadID);


    if (nBeadTypeCanOvlp[resi] == 0 && nBeadTypeCanCont[resi] == 0) {
        return totEn;
    }//No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont[resi] == 0 ? 1: 3;//No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++) {
        tmpR[j] = bead_info[beadID][j];
    }
    for (x = -BoxRad; x <= BoxRad; x++) {
        for (y = -BoxRad; y <= BoxRad; y++) {
            for (z = -BoxRad; z <= BoxRad; z++) {
                tmpR2[0] = (tmpR[0] + x + nBoxSize[0]) % nBoxSize[0];
                tmpR2[1] = (tmpR[1] + y + nBoxSize[1]) % nBoxSize[1];
                tmpR2[2] = (tmpR[2] + z + nBoxSize[2]) % nBoxSize[2];
                secBi = naTotLattice[Lat_Ind_FromVec(tmpR2)];
                if (secBi != -1 && secBi != beadID) {
                    resj = bead_info[secBi][BEAD_TYPE];
                    xDis = sqrtf((float) (x * x + y * y + z * z));
                    if (bead_info[secBi][BEAD_CHAINID] == bead_info[beadID][BEAD_CHAINID]) {//Intra-molecular
                        if (xDis <= 1.74) { // 1/r^3 potential
                            totEn += fEnergy[resi][resj][E_OVLP] /2.;// / xDis / xDis / xDis /2.;
                        }
                        // 1/r potential that goes till cube three
                        totEn += fEnergy[resi][resj][E_CONT] / xDis /2.;
                    }
                    else{//Inter-molecular
                        if (xDis <= 1.74) { // 1/r^3 potential
                            totEn += fEnergy[resi][resj][E_OVLP] ;// xDis / xDis / xDis;
                        }
                        // 1/r potential that goes till cube three
                        totEn += fEnergy[resi][resj][E_CONT] / xDis;
                    }
                }
            }
        }
    }

    return totEn;

}

/// Energy_Isotropic_Contiguous_Range - returns the isotropic interactions for beadID
/// This function assumes that we get a beadID range [smallest_bead, largest_bead], and that
/// we will loop over every bead in this range. So if the interactor is within this range,
/// we divide by 2 to account for double counting.
/// \param beadID
/// \return
float Energy_Isotropic_Contiguous_Range(int beadID, int smallest_bead, int largest_bead) {//Calculate Contact and Overlap energy of beadID
    //Takes care of intra-molecular double counting
    //If the bead is between smallest_bead and largest_bead, we divide by two.
    //This assumes that every bead between smallest_bead and largest_bead will be looped over
    float totEn = 0.0; //Storing total overlap energy
    int i, j;//Indecies
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    int x, y, z; //Lattice indecies
    int secBi, resj;//Second bead index
    float xDis = 0.;//Distance between beads.
    int resi = bead_info[beadID][BEAD_TYPE];
    totEn += nThermalization_Mode == -1 ? 0. : Energy_InitPotential(beadID);


    if (nBeadTypeCanOvlp[resi] == 0 && nBeadTypeCanCont[resi] == 0) {
        return totEn;
    }//No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont[resi] == 0 ? 1: 3;//No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++) {
        tmpR[j] = bead_info[beadID][j];
    }
    for (x = -BoxRad; x <= BoxRad; x++) {
        for (y = -BoxRad; y <= BoxRad; y++) {
            for (z = -BoxRad; z <= BoxRad; z++) {
                tmpR2[0] = (tmpR[0] + x + nBoxSize[0]) % nBoxSize[0];
                tmpR2[1] = (tmpR[1] + y + nBoxSize[1]) % nBoxSize[1];
                tmpR2[2] = (tmpR[2] + z + nBoxSize[2]) % nBoxSize[2];
                secBi = naTotLattice[Lat_Ind_FromVec(tmpR2)];
                if (secBi != -1 && secBi != beadID) {
                    resj = bead_info[secBi][BEAD_TYPE];
                    xDis = sqrtf((float) (x * x + y * y + z * z));
                    if (bead_info[secBi][BEAD_CHAINID] == bead_info[beadID][BEAD_CHAINID]) {//Intra-molecular
                        if (secBi >= smallest_bead && secBi <= largest_bead) {//Within subset
                            if (xDis <= 1.74) { // 1/r^3 potential
                                totEn += fEnergy[resi][resj][E_OVLP] /2.;// xDis / xDis / xDis / 2.;
                            }
                            // 1/r potential that goes till cube three
                            totEn += fEnergy[resi][resj][E_CONT] / xDis / 2.;
                        }
                        else{
                            if (xDis <= 1.74) { // 1/r^3 potential
                                totEn += fEnergy[resi][resj][E_OVLP] ;// xDis / xDis / xDis;
                            }
                            // 1/r potential that goes till cube three
                            totEn += fEnergy[resi][resj][E_CONT] / xDis;
                        }
                    }
                    else{//Inter-molecular
                        if (xDis <= 1.74) { // 1/r^3 potential
                            totEn += fEnergy[resi][resj][E_OVLP]; // / xDis / xDis / xDis;
                        }
                        // 1/r potential that goes till cube three
                        totEn += fEnergy[resi][resj][E_CONT] / xDis;
                    }
                }
            }
        }
    }

    return totEn;

}

/// Energy_Isotropic_With_List - returns the isotropic interactions for beadID.
/// This function assumes that we get a list of beadID's of size list_size, but that the
/// last bead is at list_size-1. This also assumes that we will loop over the entire list
/// So if the interactor is within this list, we divide by 2 to account for double counting.
/// \param beadID
/// \return
float Energy_Isotropic_With_List(const int beadID, const int *bead_list, const int list_size) {//Calculate Contact and Overlap energy of beadID
    //Takes care of intra-molecular double counting
    float totEn = 0.0; //Storing total overlap energy
    int i, j;//Indecies
    int bead_check = 0;
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    int x, y, z; //Lattice indecies
    int secBi, resj;//Second bead index
    float xDis = 0.;//Distance between beads.
    int resi = bead_info[beadID][BEAD_TYPE];
    totEn += nThermalization_Mode == -1 ? 0. : Energy_InitPotential(beadID);


    if (nBeadTypeCanOvlp[resi] == 0 && nBeadTypeCanCont[resi] == 0) {
        return totEn;
    }//No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont[resi] == 0 ? 1: 3;//No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++) {
        tmpR[j] = bead_info[beadID][j];
    }
    for (x = -BoxRad; x <= BoxRad; x++) {
        for (y = -BoxRad; y <= BoxRad; y++) {
            for (z = -BoxRad; z <= BoxRad; z++) {
                tmpR2[0] = (tmpR[0] + x + nBoxSize[0]) % nBoxSize[0];
                tmpR2[1] = (tmpR[1] + y + nBoxSize[1]) % nBoxSize[1];
                tmpR2[2] = (tmpR[2] + z + nBoxSize[2]) % nBoxSize[2];
                secBi = naTotLattice[Lat_Ind_FromVec(tmpR2)];
                if (secBi != -1 && secBi != beadID) {
                    resj = bead_info[secBi][BEAD_TYPE];
                    xDis = sqrtf((float) (x * x + y * y + z * z));

                    bead_check = 0;
                    for (i=0; i<list_size; i++){
                        if (secBi == bead_list[i]){
                            bead_check = 1;
                            break;
                        }
                    }
                    if (bead_check == 1) {//Intra-list
                        if (xDis <= 1.74) { // 1/r^3 potential
                            totEn += fEnergy[resi][resj][E_OVLP] /2.;// / xDis / xDis / xDis /2.;
                        }
                        // 1/r potential that goes till cube three
                        totEn += fEnergy[resi][resj][E_CONT] / xDis /2.;
                    }
                    else{//Inter-list
                        if (xDis <= 1.74) { // 1/r^3 potential
                            totEn += fEnergy[resi][resj][E_OVLP]; // / xDis / xDis / xDis;
                        }
                        // 1/r potential that goes till cube three
                        totEn += fEnergy[resi][resj][E_CONT] / xDis;
                    }
                }
            }
        }
    }

    return totEn;

}

//TODO: Change total energy calculation so that it properly decomposes the system energy into individual components.
/// Energy_Total_System - calculates the total energy of the system using the functions above.
/// Note the factor of 1/2 to account for double counting since all the energy contributions (for now)
/// are pair-wise. Furthermore, faCurrEn[] stores the two energies separately as well.
void Energy_Total_System(void) {
    int i;//Indecies
    // initialization
    for (i = 0; i < MAX_E; i++) {
        faCurrEn[i] = 0.0;
    }
//printf("We have %d beads and %d chains", tot_beads, tot_chains);
    for (i = 0; i < tot_beads; i++) {
        faCurrEn[E_OVLP] += Energy_Isotropic(i) / 2.;
        faCurrEn[E_SC_SC] += Energy_Anisotropic(i) / 2.;
        //printf("Done with bead %d\n", i);
    }

    for (i = 1; i < MAX_E; i++) {
        faCurrEn[E_TOT] += faCurrEn[i];
    }
}

/// Energy_Of_Chain - calculates the total energy of a molecule. Takes care of double-counting intra interactions
/// \param chainID - ID of the molecule to calculate the energy.
/// \return The total aniso + isotropic energy of this chain.
float Energy_Of_Chain(int chainID) {//Calculates the energy of the given chain.
    float totEn = 0.0;
    int i;//Looping index
    int fB = chain_info[chainID][CHAIN_START];
    int lB = chain_info[chainID][CHAIN_START] + chain_info[chainID][CHAIN_LENGTH];
    for (i = fB; i < lB; i++) {
        totEn += Energy_Anisotropic_For_Chain(i) + Energy_Isotropic_For_Chain(i);
    }

    return totEn;
}

/// Energy_Of_Chain_Self - only calculates the intra-molecular energy of this molecule
/// \param chainID - ID of the molecule to calculate the energy.
/// \return The total aniso + isotropic energy of this chain.
float Energy_Of_Chain_Self(int chainID) {//Calculates only intra-molecular interactions
    float totEn = 0.0;
    int i;//Looping index
    int fB = chain_info[chainID][CHAIN_START];
    int lB = chain_info[chainID][CHAIN_START] + chain_info[chainID][CHAIN_LENGTH];
    for (i = fB; i < lB; i++) {
        totEn += Energy_Anisotropic_Self(i) + Energy_Isotropic_Self(i);
    }

    return totEn/2.;
}

/// Energy_Of_Chain_OLD - calculates the total energy of a molecule. Double counts intra-molecular energies.
/// \param chainID - ID of the molecule to calculate the energy.
/// \return The total aniso + isotropic energy of this chain.
/// Note that this sub-routine is dumb and does not account for double counting within the same chain.
/// As such, this function should be used in both the old and new energy calculations so that the energy difference
/// is still correct.
float Energy_Of_Chain_OLD(int chainID) {//Calculates the energy of the given chain
    float totEn = 0.0;
    int i;//Looping index

    for (i = chain_info[chainID][CHAIN_START];
         i < chain_info[chainID][CHAIN_START] + chain_info[chainID][CHAIN_LENGTH]; i++) {
        totEn += Energy_Anisotropic(i) + Energy_Isotropic(i);
    }

    return totEn;
}
