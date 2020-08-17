#include "global.h"
#include "initialize.h"
#include "structure.h"

/// Memory_Initialization_AtStart - allocates memory, and initializes the various global arrays.
void Memory_Initialization_AtStart(void) {
    int i, j;
    naTotLattice = malloc(nBoxSize[0] * nBoxSize[1] * nBoxSize[2] * sizeof(lInt));
    naClusHistList = malloc((1 + tot_chains) * sizeof(lLong));
    naChainCheckList = malloc((1 + tot_chains) * sizeof(lInt));
    fKT_Cycle = malloc((1 + nTot_CycleNum) * sizeof(float));
    ld_TOTCLUS_ARR = (lLDub **) malloc((nTot_CycleNum) * sizeof(lLDub));
    for (i = 0; i < nTot_CycleNum; i++) {
        ld_TOTCLUS_ARR[i] = (lLDub *) malloc((1 + tot_chains) * sizeof(lLDub));
        if (ld_TOTCLUS_ARR[i] == NULL) {
            printf("Malloc Failed for Cluster! Crashing. Probably ran out of memory. Reduce number of chains.\n");
            exit(1);
        }
    }
    ld_TOTGYRRAD_ARR = (lLDub **) malloc((nTot_CycleNum) * sizeof(lLDub));
    for (i = 0; i < nTot_CycleNum; i++) {
        ld_TOTGYRRAD_ARR[i] = (lLDub *) malloc((2) * sizeof(lLDub));
    }
    nRDF_TotComps = 2 + nBeadTypes + nBeadTypes * nBeadTypes;
    nRDF_TotComps /= 2;
    ld_TOTRDF_Arr = malloc((nTot_CycleNum * nRDF_TotComps * nRDF_TotBins) * sizeof(lLDub));
    ldRDF_Arr = malloc((nRDF_TotComps * nRDF_TotBins) * sizeof(lLDub));

    if (naTotLattice == NULL || naClusHistList == NULL ||
        naChainCheckList == NULL || ldRDF_Arr == NULL ||
        ld_TOTRDF_Arr == NULL) {
        printf("Malloc Failed! Crashing. Probably ran out of memory.\n");
        exit(1);
    } else {
        printf("Successfully allocated memory! Arrays initialized.\n");
    }
}

void Global_Array_Initialization_AtStart(void) {
    int i, j, k, xTemp, yTemp, zTemp;//Indecies
    int myCubeLen = 3;

    for (i = 0; i < nBoxSize[0] * nBoxSize[1] * nBoxSize[2]; i++) {//Initializing the lattice
        naTotLattice[i] = -1; //If -1, then there is no bead there.
    }
    i = 0;
    //Constructing the local array to search for bonding partners.
    for (xTemp = -myCubeLen; xTemp <= myCubeLen; xTemp++) {
        for (yTemp = -myCubeLen; yTemp <= myCubeLen; yTemp++) {
            for (zTemp = -myCubeLen; zTemp <= myCubeLen; zTemp++) {
                if ((xTemp != 0 || yTemp != 0 || zTemp != 0)
                    //&& (xTemp*xTemp + yTemp*yTemp + zTemp*zTemp > 3)
                    && (xTemp * xTemp + yTemp * yTemp + zTemp * zTemp <= 3)) {
                    LocalArr[i][0] = xTemp;
                    LocalArr[i][1] = yTemp;
                    LocalArr[i][2] = zTemp;
                    i++;
                }
            }
        }
    }
    //Initializing rotational orientational bias arrays.
    for (i = 0; i < MAX_VALENCY; i++) {
        for (j = 0; j < MAX_ROTSTATES; j++) {
            rot_trial[i][j] = -1;
        }
    }
    for (i = 0; i < MAX_ROTSTATES - 1; i++) {
        Rot_IndArr[i] = i;
    }
    //Initializing arrays required for cluster analyses.
    for (i = 0; i <= tot_chains; i++) {
        naChainCheckList[i] = 0;
        naClusHistList[i] = 0;
        for (j = 0; j <= tot_chains; j++) {
            naCluster[i][j] = -1;
        }
    }
    //Initializing arrays for pair-distribution calculations
    for (j = 0; j < nRDF_TotComps; j++) {
        for (i = 0; i < nRDF_TotBins; i++) {
            ldRDF_Arr[RDFArr_Index(0, j, i)] = 0.;
        }
    }
    for (k = 0; k < nTot_CycleNum; k++) {
        for (j = 0; j < nRDF_TotComps; j++) {
            for (i = 0; i < nRDF_TotBins; i++) {
                ld_TOTRDF_Arr[RDFArr_Index(k, j, i)] = 0.;
            }
        }
        ld_TOTGYRRAD_ARR[k][0] = 0.;
        ld_TOTGYRRAD_ARR[k][1] = 0.;
    }
    for (k = 0; k < nTot_CycleNum; k++) {
        for (i = 0; i <= tot_chains; i++) {
            ld_TOTCLUS_ARR[k][i] = 0.;
        }
    }
    //Setting counters
    fSysGyrRad = 0.;
    nTotGyrRadCounter = 0;
    nRDFCounter = 0;
    nTotClusCounter = 0;
    nLargestClusterRightNow = 0;

    //Checking which bead types interact rotationally and via overlap, separately.
    for (i = 0; i < MAX_AA; i++) {
        nBeadTypeIsSticker[i] = 0;//Assume beads don't rotationally interact.
        nBeadTypeCanOvlp[i] = 0;//Assume beads don't have an overlap cost
    }

    for (i = 0; i < MAX_AA; i++) {
        for (j = 0; j < MAX_AA; j++) {
            if (fEnergy[i][j][E_SC_SC] != 0.0) {//Seeing if this beadType rotationally interacts.
                nBeadTypeIsSticker[i] = 1;
            }
            if (fEnergy[i][j][E_OVLP] != 0.0) {//Seeing if this beadType rotationally interacts.
                nBeadTypeCanOvlp[i] = 1;
            }
            if (fEnergy[i][j][E_CONT] != 0.0) {//Seeing if this beadType rotationally interacts.
                nBeadTypeCanCont[i] = 1;
            }
        }
    }
    for (i = MV_NULL + 2; i < MAX_MV; i++) {
        fMCFreq[i] += fMCFreq[i - 1]; // Cumulative Frequencies
    }
    if (RotBias_Mode == 1) {
        fRot_Bias = expf(-fRot_Bias / fKT);
    }

    for (i = 0; i < nTot_CycleNum; i++) {
        fKT_Cycle[i] = fKT + (float) i * fdelta_temp;
    }
    ld_LogOfSmallestPossibleProb = logl((lLDub) 1. / (lLDub) RAND_MAX);
    printf("All setup has been completed!\n");
}

/// Reset_Global_Arrays - resets the various global counters and arrays for system analysis.
void Reset_Global_Arrays(void) {
    //Zero-ing out all the arrays used for data tracking!
    int i, j;
    for (i = 0; i <= tot_chains; i++) {
        naChainCheckList[i] = 0;
        naClusHistList[i] = 0;
        for (j = 0; j <= tot_chains; j++) {
            naCluster[i][j] = -1;
        }
    }
    //Initializing arrays for pair-distribution calculations
    for (j = 0; j < nRDF_TotComps; j++) {
        for (i = 0; i < nRDF_TotBins; i++) {
            ldRDF_Arr[RDFArr_Index(0, j, i)] = 0.;
        }
    }
    //Setting counters
    fSysGyrRad = 0.;
    nTotGyrRadCounter = 0;
    nRDFCounter = 0;
    nTotClusCounter = 0;
    nLargestClusterRightNow = 0;
}

/// Initial_Conditions_Simple - randomly place all the molecules.
void Initial_Conditions_Simple(void) {
    /*
    Generates initial conditions for a system with topology information. The idea is to see if the selected bead has been placed or not.
    If it has been placed, we move on to sequentially placing all of it's bonded beads within linker constraints.
    Therefore, each bead shall sprout forth its partners. Since I like goto statements, they have been used.
    If a bead has not been placed, it's first checked if one of its bonding partners might have been placed,
    and will thus act as an anchor. This means that for each chain we need to know which beads have already been placed,
    thus a temporary list is used.
    */

    int idx, idy;//Just internal counters for various things.
    int i, j, k;//Iterators for loops
    int bondPart;//Used to track the bond partner for a particular bead.
    int fB, lB;//Used to track the first and last bead of a particular chain for looping.
    int radUp, radLow;//Used as radii to generate coordinates within a linker-length sphere
    //Just remember tha lB-1 is the last bead in that chain, but in loops we go <lB.
    int TlCnt = 0;//Counter for trials around an anchor.
    int tmpR[POS_MAX], tmpR2[POS_MAX];//Arrays to store temporary coordinates;
    int temp_list[MAX_CHAINLEN];//Used to store already placed beads.
    int list_it = 0;//Iterator specifically for temp_list.

    for (i = 0; i < MAX_CHAINLEN; i++) {//initialize the list where -1 signifies emptiness.
        temp_list[i] = -1;
    }
    for (j = 0; j < POS_MAX; j++) {//Initialize the coordinate arrays to some numbers.
        tmpR[j] = rand() % nBoxSize[j];
        tmpR2[j] = tmpR[j];
    }

    for (k = 0; k < tot_chains; k++) {
        //This makes reading easier, honestly
        fB = chain_info[k][CHAIN_START];
        lB = fB + chain_info[k][CHAIN_LENGTH];

        //Reset the temp_list for each chain, and the iterator!
        for (i = 0; i < list_it; i++) {//initialize the list where -1 signifies emptiness.
            temp_list[i] = -1;
        }
        list_it = 0;

        for (i = fB; i < lB; i++) {//Going bead-bead in this chain.
            //Checking if the bead has already been placed.
            for (idy = 0; idy < list_it; idy++) {
                if (i == temp_list[idy]) {//This means this bead has already been placed.
                    goto SproutForth;
                }
            }
            //If we have reached here, this bead has not been placed before.
            idx = 0;
            bondPart = topo_info[i][idx];//Re-initialize these two.
            //Let's go through the bonded partners for this bead and see if one can be used as an anchor.
            while (topo_info[i][idx] != -1 && idx < MAX_BONDS) {//Again, -1 signifies that no bonded bead exists.
                bondPart = topo_info[i][idx];
                //Checking if thisBead has been placed already, and can be used as anchor.
                for (idy = 0; idy < list_it; idy++) {
                    if (bondPart == temp_list[idy]) {//Using bondPart as an anchor
                        goto FoundAnchor;//This is just so we don't keep going over the list needlessly and just pick the first one.
                    }
                }
                idx++;
                //If we got here, this means no bonded partner has already been placed, so no anchor.
                bondPart = -1;
            }

            FoundAnchor:
            if (bondPart != -1) {//We found an anchor.
                radUp = 2 * linker_len[i][idx] + 1;
                radLow = linker_len[i][idx];
                for (j = 0; j < POS_MAX; j++) {//Using thisBead as anchor
                    tmpR[j] = bead_info[bondPart][j];
                    tmpR2[j] = tmpR[j];
                }
                TlCnt = 0;//Reset this counter.
                while (naTotLattice[Lat_Ind_FromVec(tmpR2)] != -1 && TlCnt < 500000) {
                    for (j = 0; j < POS_MAX; j++) {
                        tmpR2[j] = (rand() % radUp) - radLow;
                        tmpR2[j] = (tmpR[j] + tmpR2[j] + nBoxSize[j]) % nBoxSize[j];
                    }
                    TlCnt++;
                }
                if (TlCnt == 500000) {
                    printf("Not enough space in the lattice. Crashing. Maybe try increasing max trials, "
                           "or make the box bigger!\t\n");
                    exit(1);
                }
                for (j = 0; j < POS_MAX; j++) {
                    bead_info[i][j] = tmpR2[j];
                }//Placing bead
                naTotLattice[Lat_Ind_FromVec(tmpR2)] = i;//Putting on lattice
                temp_list[list_it] = i;
                list_it++;//Remembering in hash list.
                //The bead has been placed around an appropriate anchor.
            } else {//There was no appropriate anchor. So we can put this bead wherever.
                for (j = 0; j < POS_MAX; j++) {
                    tmpR[j] = rand() % nBoxSize[j];
                }
                //This usually means this is the first bead in the chain.
                TlCnt = 0;
                while (naTotLattice[Lat_Ind_FromVec(tmpR)] != -1 &&
                       TlCnt < 500000) {//Keep looping till we find an empty lattice site.
                    TlCnt++;
                    for (j = 0; j < POS_MAX; j++) {//Generate a random point in the lattice.
                        tmpR[j] = rand() % nBoxSize[j];
                    }
                }
                if (TlCnt == 500000) {
                    printf("\n\nNot enough space in the lattice for first bead, bruh. Maybe try increasing max trials, "
                           "or make the box bigger!\n\n");
                    exit(1);
                }
                //Placing this bead wherever there is space, and using it as anchor.
                //Also updating the lattice, and temp_list
                for (j = 0; j < POS_MAX; j++) {
                    bead_info[i][j] = tmpR[j];
                }
                naTotLattice[Lat_Ind_FromVec(tmpR)] = i;//Placing on lattice!
                temp_list[list_it] = i;
                list_it++;//Remembering that this bead has been placed.
            }
            SproutForth:
            //Now going over the list of bondParts for i and sprouting them
            for (j = 0; j < POS_MAX; j++) {//Making the current bead an anchor.
                tmpR[j] = bead_info[i][j];
            }
            idx = 0;//Resets things!
            bondPart = topo_info[i][idx];//Resets things!
            while (topo_info[i][idx] != -1 && idx < MAX_BONDS) {//Again, -1 signifies that no bonded bead exists.
                bondPart = topo_info[i][idx];
                //If the bead has already been dealt with, move on to next potential bondPartner
                for (idy = 0; idy < list_it; idy++) {
                    if (bondPart == temp_list[idy]) {
                        goto SkipThisPartner;//Just so we don't keep looking further.
                    }
                }
                //Get the radii for placing this bead around bead i
                radUp = 2 * linker_len[i][idx] + 1;
                radLow = linker_len[i][idx];
                //Initializing new vector to be where bead i is.
                for (j = 0; j < POS_MAX; j++) {
                    tmpR2[j] = tmpR[j];
                }
                TlCnt = 0;//Reset this counter.
                while (naTotLattice[Lat_Ind_FromVec(tmpR2)] != -1 && TlCnt < 500000) {
                    for (j = 0; j < POS_MAX; j++) {
                        tmpR2[j] = (rand() % radUp) - radLow;
                        tmpR2[j] = (tmpR[j] + tmpR2[j] + nBoxSize[j]) % nBoxSize[j];
                    }
                    //printf("%d %d %d\n", tmpR2[0], tmpR2[1], tmpR2[2]);
                    TlCnt++;
                }
                if (TlCnt == 500000) {
                    printf("\n\nNot enough space in the lattice. Crashing. Maybe try increasing max trials, "
                           "or make the box bigger!\t %d\t%d %d\n\n", TlCnt, i, naTotLattice[Lat_Ind_FromVec(
                            tmpR2)]);
                    exit(1);
                }
                for (j = 0; j < POS_MAX; j++) {//Placing bead
                    bead_info[bondPart][j] = tmpR2[j];
                }

                naTotLattice[Lat_Ind_FromVec(tmpR2)] = bondPart;//Putting on lattice
                temp_list[list_it] = bondPart;
                list_it++;//Remembering in hash list.
                SkipThisPartner:
                idx++;
                bondPart = -1;
            }
            //Moving on to bead i+1
        }
        //Moving on to the next molecule!
    }
    //Initializing all beads with no physical bonds.
    for (i = 0; i < tot_beads; i++) {
        bead_info[i][BEAD_FACE] = -1;
    }
}

/// Initial_Conditions_FromFile - reads the restart file and sets the initial conditions for the system
void Initial_Conditions_FromFile(void) {
    printf("Reading file to generate initial configuration\n");
    char strLine[1000];
    int i, j;
    int tmp_beadinfo[BEADINFO_MAX];
    FILE *infile;
    infile = fopen(strRestartFile, "r");
    int tmpBeads = 0;
    while (fgets(strLine, sizeof(strLine), infile) != NULL) {
        tmpBeads++;
        if (strcmp(strLine, "ITEM: ATOMS id type mol x y z bP\n") == 0) {
            break;
        }
    }
    tmpBeads = 0;
    while (fgets(strLine, sizeof(strLine), infile) != NULL) {
        sscanf(strLine, "%d %d %d %d %d %d %d",
               &i,
               &tmp_beadinfo[BEAD_TYPE],
               &tmp_beadinfo[BEAD_CHAINID],
               &tmp_beadinfo[POS_X],
               &tmp_beadinfo[POS_Y],
               &tmp_beadinfo[POS_Z],
               &tmp_beadinfo[BEAD_FACE]);
        tmpBeads++;
        for (j = 0; j < BEADINFO_MAX; j++) {
            bead_info[i][j] = tmp_beadinfo[j];
        }
    }
    fclose(infile);
    if (tmpBeads != tot_beads) {
        printf("Given restart file does not have the same number of beads as the provided structure.");
        printf("Crashing!\n");
        exit(1);
    }
    for (i = 0; i < tot_beads; i++) {
        naTotLattice[Lat_Ind_OfBead(i)] = i;
    }
    printf("File read successfully and initial configuration set.\n");
}

/// Initial_Conditions_BreakBonds - glorified for loop to break every bond. Is used when a restart file with
/// thermalization steps.
void Initial_Conditions_BreakBonds(void) {
    int i;
    for (i = 0; i < tot_beads; i++) {
        bead_info[i][BEAD_FACE] = -1;
    }
    printf("Since system has thermalization cycle, deleting all physical bonds!\n");
}

/// Calculate_Rot_Bias - calculates the possible anisotropic Boltzmann weights at the beginning so that
/// we have a look-up table.
/// \param CurrentTemp
void Calculate_Rot_Bias(float CurrentTemp) {

    int i, j;
    for (i = 0; i < MAX_AA; i++) {
        for (j = 0; j < MAX_AA; j++) {
            dbias_bolt_fac[i][j] = expf(-fEnergy[i][j][E_SC_SC] / CurrentTemp);
        }
    }
    //TODO: Make sure that fRot_Bias can be used in the future to set a solvent 'anisotropy'
    fRot_Bias = expf(-f_globRotBias / CurrentTemp);

}

/// Temperature_Function - used to calculate what fCuTemp should (current temperature) based on how long the run has been
/// going.
/// \param mode - which of the four functions to use.
/// Note that the various modes are described below where \f$F(t)\f$ is the returned value, and t == nGen.
/// If mode = 0: \f$F(t) = fKT + \tanh(1.+ (nGen-nPre)/1250fPreKT\f$. A hyperbolic tangent to smoothly reduce temperature,
/// after nPreSteps
/// If mode = 1: \f$F(t) = fKT + 5\exp(-(nGen-nPre)/4nPre)\abs(sin((nGen-nPre)/nPre))\f$. An exponentially decaying
/// sinusoidal bouncing after nPreSteps
/// If mode = 2: \f$F(t) = fKT + 5\exp(-(nGen-10nPre)^2/10nPre^2)/fKT\f$. Gaussian like decay after nPreSteps.
/// If mode = 3: \f$F(t) = fKT + \exp(-4nGen/nPre)/fKT\f$. Exponential decay from the beginning.
/// \param nGen - basically 'time' or how many MC Steps have been done.
/// \return end_val - the current temperature.
float Temperature_Function(int mode, long nGen) {

    float x_val;
    float y_val;
    float end_val;

    switch (mode) {
        case 0:
            x_val = (float) (nPreSteps - nGen);
            x_val = x_val / 1250. / fPreKT;
            x_val = fPreKT * (tanhf(x_val) + 1.);
            end_val = fKT + x_val;

            break;

        case 1:
            x_val = (float) (nGen - nPreSteps) / (float) nPreSteps;
            y_val = -x_val;
            x_val = fabsf(sinf(x_val));
            y_val = expf(y_val / 4.);
            end_val = fKT + 5. * fKT * x_val * y_val;

            break;

        case 2:
            x_val = (float) (nGen - 10 * nPreSteps);
            x_val = x_val * x_val;
            y_val = (float) (nPreSteps * nPreSteps) * 10.;
            end_val = fKT + expf(-x_val / y_val) / fKT / 10.;

            break;

        case 3:
            x_val = -(float) (nGen);
            x_val = 4. * x_val;
            x_val = x_val / (float) (nPreSteps);
            end_val = fKT + expf(x_val);

            break;

        default:

            end_val = fKT;
            break;
    }

    return end_val;

}

