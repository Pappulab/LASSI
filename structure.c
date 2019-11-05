#include "global.h"
#include "structure.h"

/// Lat_Ind_FromCoords - helper function to get the correct 1D index of this position
/// \param i
/// \param j
/// \param k
/// \return the 1D index location for (i,j,k) position
int Lat_Ind_FromCoords(int i, int j, int k) {//Lattice index from 3D to 1D array
    return i + nBoxSize[0] * (j + nBoxSize[1] * k);
}

/// Lat_Ind_FromVec - returns the 1D index given the array xArr
/// \param xArr
/// \return
int Lat_Ind_FromVec(int *xArr) {//Just the vector form of the above function. Easier to read sometimes.
    return xArr[POS_X] + nBoxSize[POS_X] * (xArr[POS_Y] + nBoxSize[POS_Y] * xArr[POS_Z]);
}

/// Lat_Ind_OfBead - returns the 1D index of this bead's location
/// \param beadID
/// \return
int Lat_Ind_OfBead(int beadID) {
    return bead_info[beadID][POS_X] +
           nBoxSize[POS_X] * (bead_info[beadID][POS_Y] + nBoxSize[POS_Y] * bead_info[beadID][POS_Z]);
}

//Note that the distance functions account for periodic boundaries

/// Dist_PointTotPoint_Float - euclidean distance between the vectors (arrays) f1 and f2 where f1 and f2 are floats.
/// \param f1
/// \param f2
/// \return
float Dist_PointTotPoint_Float(float *f1, float *f2) {
    float d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++) {
        d[i] = fabsf(f1[i] - f2[i]);
        d[i] = d[i] > (float) nBoxSize[i] / 2. ? (float) nBoxSize[i] - d[i] : d[i];
    }

    return sqrtf(d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]);
}

/// Dist_PointToPoint - euclidean distance between the vectors (arrays) f1 and f2 where f1 and f2 are integer points.
/// \param f1
/// \param f2
/// \return
float Dist_PointToPoint(int *f1, int *f2) {
    int d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++) {
        d[i] = abs(f1[i] - f2[i]);
        d[i] = d[i] > nBoxSize[i] / 2 ? nBoxSize[i] - d[i] : d[i];

    }

    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToPoint - euclidean distance between beadID and the vector f1.
/// \param beadID
/// \param f1
/// \return
float Dist_BeadToPoint(int beadID, int *f1) {
    int d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++) {
        d[i] = abs(bead_info[beadID][i] - f1[i]);
        d[i] = d[i] > nBoxSize[i] / 2 ? nBoxSize[i] - d[i] : d[i];
    }
    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToBead - euclidean distance between the two beads.
/// \param n1
/// \param n2
/// \return
float Dist_BeadToBead(int n1, int n2) {
    lInt d[POS_MAX];
    lInt i;

    for (i = 0; i < POS_MAX; i++) {
        d[i] = abs(bead_info[n1][i] - bead_info[n2][i]);
        d[i] = d[i] > nBoxSize[i] / 2 ? nBoxSize[i] - d[i] : d[i];
    }

    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Check_System_Structure -
/// \return 0 if everything is okay, beadID+1 if failed.
int Check_System_Structure(void) {
    int i, j;//Looping variables
    int idx;//Internal iterators for covalent bonds.
    int tmpR[POS_MAX];//Just a vector to store coordinates
    int bondPart;
    for (i = 0; i < tot_beads; i++) {
        idx = 0;
        bondPart = topo_info[i][idx];
        for (j = 0; j < POS_MAX; j++) {
            tmpR[j] = bead_info[i][j];
        }
        if (naTotLattice[Lat_Ind_FromVec(tmpR)] == -1) {
            printf("Lattice Position for bead %d is empty! Chain: %d\n", i, bead_info[i][BEAD_CHAINID]);
            return i + 1;
        }
        if (i - naTotLattice[Lat_Ind_FromVec(tmpR)] !=
            0) {//This means there is a mismatch between where the bead is and where the lattice thinks the bead is
            printf("Bead position and lattice value not the same. Crashing\t\t");
            printf("B1:%d B2:%d\t C1:%d C2:%d\n", i, naTotLattice[Lat_Ind_FromVec(tmpR)], bead_info[i][BEAD_CHAINID],
                   bead_info[naTotLattice[Lat_Ind_FromVec(
                           tmpR)]][BEAD_CHAINID]);
            return i + 1;
        }
        while (topo_info[i][idx] != -1 && idx < MAX_BONDS) {
            bondPart = topo_info[i][idx];
            if (Dist_BeadToBead(i, bondPart) > 1.74 * linker_len[i][idx]) {
                printf("Bad beads! %d\t(%d %d %d)\t\tTopo:(%d %d %d)\t\tLinkers:(%.5f\t%.5f\t%.5f)\n", i,
                       bead_info[i][0], bead_info[i][1], bead_info[i][2],
                       topo_info[i][0], topo_info[i][1], topo_info[i][2],
                       (float) linker_len[i][0], (float) linker_len[i][1], (float) linker_len[i][2]);
                printf("\t\t\t\t\t-------------------->\t\t%f\tSHOULD BE\t%f\n", Dist_BeadToBead(i, bondPart),
                       1.74 * (float) linker_len[i][idx]);
                printf("Bad beads! %d\t(%d %d %d)\t\tTopo:(%d %d %d)\t\tLinkers:(%.5f\t%.5f\t%.5f)\n\n", bondPart,
                       bead_info[bondPart][0], bead_info[bondPart][1], bead_info[bondPart][2],
                       topo_info[bondPart][0], topo_info[bondPart][1], topo_info[bondPart][2],
                       (float) linker_len[bondPart][0], (float) linker_len[bondPart][1],
                       (float) linker_len[bondPart][2]);
                return i + 1;
            }
            idx++;
        }
        if (bead_info[i][BEAD_FACE] != -1) {
            if (nBeadTypeIsSticker[bead_info[i][BEAD_TYPE]] != 1) {
                printf("This bead -- %d -- should not have a bond. Crashing.\n", i);
                return i + 1;
            }
            if (bead_info[i][BEAD_FACE] == i) {
                printf("Self bonded.\n");
                return i + 1;
            }
            if (i != bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE]) {
                printf("Bad bond!\n\t%d %d %d %f\nCrashing.\n", i, bead_info[i][BEAD_FACE],
                       bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE],
                       fEnergy[bead_info[i][BEAD_TYPE]][bead_info[bead_info[i][BEAD_FACE]][BEAD_TYPE]][E_SC_SC]);
                return i + 1;
            }
            if (Dist_BeadToBead(i, bead_info[i][BEAD_FACE]) > 1.74) {
                printf("Bad bond! Distance is wrong\n\t%d %d %d\n Distance is %f. Crashing.\n", i,
                       bead_info[i][BEAD_FACE], bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE],
                       Dist_BeadToBead(i, bead_info[i][BEAD_FACE]));
                return i + 1;
            }
        }
    }
    return 0;
}

/// Dist_VecMag - non periodic boundary euclidean magnitude of vector
/// \param f1
/// \return
float Dist_VecMag(const int *f1) {//Outputs the magnitude of the vector
    return sqrtf((float) (f1[0] * f1[0] + f1[1] * f1[1] + f1[2] * f1[2]));
}

/// GyrTensor_ClusterSpecific - calculates the total Gyration Tensor for a given cluster
/// \param ClusSize - the total size of the cluster.
/// \param ClusIndex - the index on naCluster where the cluster is stored.
/// THIS IS VERY OLD AND HASN'T BEEN LOOKED AT IN A WHILE
/// TODO: Update this for the new version
void GyrTensor_ClusterSpecific(int ClusSize, int ClusIndex) {
    //Calculate the components of the gyration tensor for a given cluster.
    //ClusSize is the size of the cluster -- obviously -- whereas ClusIndex tell us
    //where in naCluster the chain indecies are located. naCluster[ClusIndex][0-ClusSize] is all the chainID's I need
    //for the calculation
    //Remember that the Gyration Tensor is a 3x3 symmetric object so we only need 6 numbers.
    int i, k, j, j2;//Basic indecies for loops
    for (i = 0; i < 7; i++) { fGyrTensor[i] = 0.; }//Initializing
    int firstB, lastB;//Tracks the first and last bead of the given chain
    float tot_COM[POS_MAX] = {0.};//This is where we shall store the COM of the cluster.
    int NumRes = 0;//Tracks how many residues are in this cluster

    //The only thing one needs to be careful about is to take PBC into account; the rest is tedium.
    //Use good old differential geometry to map each coordinate to two new coordinates, and unpack at the end.
    float theta[POS_MAX] = {0.};
    float zeta[POS_MAX] = {0.};//Extra coordinates for goodness
    float dumArg = 0.;//Just a dummy variable to be more efficient
    float dumArg2 = 0.;//Another one
    // printf("Starting with COM\n");
    //Calculating the COM
    for (i = 0; i < ClusSize; i++) {
        firstB = chain_info[naCluster[ClusIndex][i]][CHAIN_START];
        lastB = firstB + chain_info[naCluster[ClusIndex][i]][CHAIN_LENGTH];
        // printf("%d %d\n", firstB, lastB);
        //Just easier to track each chain like this
        for (k = firstB; k < lastB; k++) {
            NumRes++;//Adding a residue to the total
            for (j = 0; j < POS_MAX; j++) {
                dumArg = 2. * PI * ((float) bead_info[k][j] / (float) nBoxSize[j]);
                theta[j] += cosf(dumArg);
                zeta[j] += sinf(dumArg);//Since I am taking the average just keep adding
            }
        }
        // printf("Next bead\n");
    }
    //printf("Done with COM\n");
    //Calculating average, and then COM
    for (j = 0; j < POS_MAX; j++) {
        theta[j] = theta[j] / (float) NumRes;
        zeta[j] = zeta[j] / (float) NumRes;
        tot_COM[j] = atan2f(-theta[j], -zeta[j]) + PI;
        tot_COM[j] = nBoxSize[j] * (tot_COM[j] / 2. / PI);
    }

    //Using the COM to calculate the Gyration Tensor
    // GyrTen_{ij} = 1/N sum_1^N (r_i-com_i)(r_j-com_j) so just take the sums and divide at the end
    for (i = 0; i < ClusSize; i++) {
        firstB = chain_info[naCluster[ClusIndex][i]][CHAIN_START];
        lastB = firstB + chain_info[naCluster[ClusIndex][i]][CHAIN_LENGTH];
        //Just easier to track each chain like this
        for (k = firstB; k < lastB; k++) {
            for (j = 0; j < POS_MAX; j++) {
                for (j2 = j; j2 < POS_MAX; j2++) {
                    dumArg = fabsf((float) bead_info[k][j] - tot_COM[j]) <
                             (float) nBoxSize[j] - fabsf((float) bead_info[k][j] - tot_COM[j]) ? fabsf(
                            (float) bead_info[k][j] - tot_COM[j]) : (float) nBoxSize[j] -
                                                                    fabsf((float) bead_info[k][j] - tot_COM[j]);
                    dumArg2 = fabsf((float) bead_info[k][j2] - tot_COM[j2]) <
                              (float) nBoxSize[j2] - fabsf((float) bead_info[k][j2] - tot_COM[j2]) ? fabsf(
                            (float) bead_info[k][j2] - tot_COM[j2]) : (float) nBoxSize[j2] -
                                                                      fabsf((float) bead_info[k][j2] - tot_COM[j2]);
                    fGyrTensor[j + 3 * (j2 - j)] += dumArg * dumArg2;
                    //0 = xx; 1 = yy; 2 = zz; 3 = xy; 6 = xz; 4 = yz; Need smarter indexing
                }
            }

        }
    }

    //Calculating the average;
    for (i = 0; i < 7; i++) { fGyrTensor[i] /= (float) NumRes; }//printf("%f\n",fGyrTensor[i]);}printf("\n");



    // exit(1);

}

/// GyrTensor_GyrRad - calculates the total Gyration Tensor of the system.
/// \param ClusSize - the total size of the cluster.
/// \param ClusIndex - the index on naCluster where the cluster is stored.
/// THIS IS VERY OLD AND HASN'T BEEN LOOKED AT IN A WHILE
/// TODO: Update this for the new version
void GyrTensor_GyrRad(void) {//Calculates the gyration tensor for the whole system
    int i, k, j, j2;//Basic indecies for loops
    for (i = 0; i < 7; i++) { fGyrTensor[i] = 0.; }//Initializing
    float tot_COM[POS_MAX] = {0.};//This is where we shall store the COM of the cluster.

    //The only thing one needs to be careful about is to take PBC into account; the rest is tedium.
    //Use good old differential geometry to map each coordinate to two new coordinates, and unpack at the end.
    float theta[POS_MAX] = {0.};
    float zeta[POS_MAX] = {0.};//Extra coordinates for goodness
    float dumArg = 0.;//Just a dummy variable to be more efficient
    float dumArg2 = 0.;//Another one

    for (i = 0; i < tot_beads; i++) {
        for (j = 0; j < POS_MAX; j++) {
            tot_COM[j] += bead_info[i][j];
        }
    }

    //Calculating average, and then COM
    for (j = 0; j < POS_MAX; j++) {
        tot_COM[j] /= (float) tot_beads;
    }
    //printf("\n");
    //Using the COM to calculate the Gyration Tensor
    // GyrTen_{ij} = 1/N sum_1^N (r_i-com_i)(r_j-com_j) so just take the sums and divide at the end
    for (i = 0; i < tot_beads; i++) {
        for (j = 0; j < POS_MAX; j++) {
            dumArg = (float) bead_info[i][j] - tot_COM[j];
            for (j2 = j; j2 < POS_MAX; j2++) {
                dumArg2 = (float) bead_info[i][j2] - tot_COM[j2];
                fGyrTensor[j + 3 * (j2 - j)] += dumArg * dumArg2;
                //printf("%.2f\t", dumArg2);
                //0 = xx; 1 = yy; 2 = zz; 3 = xy; 6 = xz; 4 = yz; Need smarter indexing
            }
        }//printf("\n");
    }
    //exit(1);
    //Calculating the average;
    for (i = 0; i < 7; i++) { fGyrTensor[i] /= (float) tot_beads; }//printf("%f\n",fGyrTensor[i]);}printf("\n");
}

/// GyrTensor_GyrRad_Avg - calculates the total radius of gyration of the system, while not being smart about the
/// periodic boundaries. This is used as a proxy to detect phase separation, but is a relic of the old formalism. The
/// RDF should be used in general. Although this can be used without the need for a non-interacting prior.
void GyrTensor_GyrRad_Avg(void) {
    /*
    Only calculates the diagonals of the gyration tensor, and calculates the sum of the
    diagonals. Remember that Rg^2 = Tr(GyrTen) so we only need to calculate the diagonals, and then then sum.
    I shall borrow most of the code from above, and so read GyrTensor_ClusterSpecific for what's happening here.
    */
    int i, j;//Loop indecies
    float tot_COM[POS_MAX] = {0.};
    float dumArg = 0.;
    for (i = 0; i < 7; i++) {//Initializing to 0
        fGyrTensor[i] = 0.;
    }

    for (i = 0; i < tot_beads; i++) {
        for (j = 0; j < POS_MAX; j++) {
            tot_COM[j] += bead_info[i][j];
        }
    }
    for (j = 0; j < POS_MAX; j++) {//Dividing by number of beads
        tot_COM[j] /= (float) tot_beads;
    }

    for (i = 0; i < tot_beads; i++) {
        for (j = 0; j < POS_MAX; j++) {
            dumArg = (float) bead_info[i][j] - tot_COM[j];
            fGyrTensor[j] += dumArg * dumArg;
        }
    }

    fSysGyrRad += sqrtf((fGyrTensor[0] + fGyrTensor[1] + fGyrTensor[2]) /
                        (float) tot_beads);//Adding to the total fSysGyrRad to be averaged at the end.
    nTotGyrRadCounter++;//Remembering that we have calculated the radius; this will be used to average the final value.

}

/// RDF_ComponentIndex - 1D index for the symmetric g_{ij} matrix given i and j.
/// \param i
/// \param j
/// \return The index of the array
/// The way it is set up, the indexing goes through the diagonal and then 0-1, 0-2, ... 0-N, 1-2, ... 1-N and so on
int RDF_ComponentIndex(const int i, const int j) {
    if (i > j) {
        return RDF_ComponentIndex(j, i);
    } else {
        return i == j ? 1 + i : nBeadTypes + j - (i * (3 + i - 2 * nBeadTypes)) / 2;
    }
}

/// RDFArr_Inde - 1D index for ld_TOTRDF_Arr which is used to globally store the different RDFs
/// \param run_cycle
/// \param rdf_comp
/// \param x_pos
/// \return 1D index for the totalRDFArray
int RDFArr_Index(const int run_cycle, const int rdf_comp, const int x_pos) {
    return x_pos + nRDF_TotBins * (rdf_comp + nRDF_TotComps * run_cycle);
}

/// RDF_ComponentWise_Avg - calculates the pair-distribution of the system where every bead acts as the center
/// of a radial histogram of pairs. Note that dr = 1/4 lattice units.
void RDF_ComponentWise_Avg(void) {
    /*
    Calculates the RDF and adds it all up so that it can be averaged out at the end of the run.
    */
    float x;  //For distance
    int i, j, k;
    int resi, resj;
    //int chaini, chainj;
    int myBin = 0;
    int array_pos;

    //Calculating where, and how many, pairs exist
    for (i = 0; i < tot_beads; i++) {
        resi = bead_info[i][BEAD_TYPE];
        for (j = i + 1; j < tot_beads; j++) {
            resj = bead_info[j][BEAD_TYPE];
            x = Dist_BeadToBead(i, j);
            //Note that Dist_BeadToBead(i,j) automatically ensures no distance is greater than (L/2)*sqrt(3)
            myBin = (int) floor(4. * x);//I am assuming for now that dr=1/4
            ldRDF_Arr[RDFArr_Index(0, 0, myBin)] += 2.0;//Adding a pair to that bin
            array_pos = RDF_ComponentIndex(resi, resj);
            ldRDF_Arr[RDFArr_Index(0, array_pos, myBin)] += 2.0;
        }
    }
    nRDFCounter++;
}

/// Check_LinkerConstraint - if I move beadID to tmpR, do I still satisfy the linker lengths for beadID?
/// \param beadID
/// \param tmpR
/// \return 1 means all is good, 0 means bad.
int Check_LinkerConstraint(int beadID, int *tmpR) {
    //Check if the proposed new location for beadID is such that all the linkers are unbroken.
    int idx;//Iterator to loop over bond Partners
    int bondPartner;//It is what it is.

    idx = 0;
    bondPartner = topo_info[beadID][idx];//Initializing the two.
    while (idx < MAX_BONDS && topo_info[beadID][idx] != -1) {//Keep going till we run out of partners
        bondPartner = topo_info[beadID][idx];
        if (Dist_PointToPoint(bead_info[bondPartner], tmpR) > 1.74 * (float) linker_len[beadID][idx]) {
            return 0;//This means that we have broken one of the linkers.
        }
        idx++;
    }
    return 1;//This means that all linker constraints are satisfied.
}

/// Check_MTLinkerConstraint - if I move beadID and all of it's covalent bonded beads to the locations stored in
/// tmpR[][], would all the linker constraints be satisfied?
/// \param beadID
/// \param tmpR
/// \return
int Check_MTLinkerConstraint(int beadID, int (*tmpR)[POS_MAX]) {

    int curID = beadID;
    int idx, bPart;
    int j;
    int topIt = 0;
    int canI = 1;

    while (curID != -1) {
        for (j = 0; j < POS_MAX; j++) {
            bead_info[curID][j] = tmpR[topIt][j];//Moving
        }
        curID = topo_info[beadID][topIt++];
    }

    curID = beadID;
    topIt = 0;
    while (curID != -1 && canI == 1) {
        idx = 0;
        bPart = topo_info[curID][idx];
        while (bPart != -1 && idx < MAX_BONDS) {
            if (Dist_BeadToBead(curID, bPart) > 1.74 * (float) linker_len[curID][idx]) {
                canI = 0;
                break;
            }
            bPart = topo_info[curID][++idx];
        }
        curID = topo_info[beadID][topIt++];
    }
    curID = beadID;
    topIt = 0;
    while (curID != -1) {
        for (j = 0; j < POS_MAX; j++) {
            bead_info[curID][j] = old_bead[curID][j];//Moving back
        }
        curID = topo_info[beadID][topIt++];
    }

    return canI;
}

