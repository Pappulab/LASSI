#include "initialize.h"
#include "global.h"
#include "structure.h"

/// Memory_Initialization_AtStart - allocates memory, and initializes the
/// various global arrays.
void Memory_Initialization_AtStart(void)
{
    char arr_name[100];

    strcpy(arr_name, "naTotLattice_glb");
    naTotLattice_glb = Create1DInt(naBoxSize_glb[0] * naBoxSize_glb[1] * naBoxSize_glb[2], arr_name);

    Memory_Allocate_NeighborLists();

    strcpy(arr_name, "laClusHistList_glb");
    laClusHistList_glb = Create1DLong((2 + tot_chains_glb), arr_name);

    strcpy(arr_name, "naChainCheckList_glb");
    naChainCheckList_glb = Create1DInt(2 + tot_chains_glb, arr_name);

    strcpy(arr_name, "fKTCycle");
    faKT_Cycle_glb = Create1DFloat(1 + nTotCycleNum_glb, arr_name);

    strcpy(arr_name, "naList_glb");
    naList_glb = Create1DInt(2 + tot_chains_glb, arr_name);

    strcpy(arr_name, "naClusterMatrix_glb");
    naClusterMatrix_glb = Create2DInt(tot_chains_glb + 2, tot_chains_glb + 2, arr_name);

    if (naReportFreqs_glb[REPORT_NETWORK] != 0)
        {
            strcpy(arr_name, "ldTotClusArr");
            ldaTOTCLUS_Arr_glb = Create2DLongdouble(2 + tot_chains_glb, nTotCycleNum_glb, arr_name);

            strcpy(arr_name, "ldMolClusArr");
            ldaMOLCLUS_Arr_glb = Create1DLongdouble(((tot_chains_glb + 1) * tot_chain_types_glb), arr_name);

            strcpy(arr_name, "ldTotMolClusArr");
            ldaTOTMOLCLUS_Arr_glb =
                Create1DLongdouble((tot_chains_glb * tot_chain_types_glb * nTotCycleNum_glb), arr_name);

            strcpy(arr_name, "ldTotGyrRad");
            ldaTOTRg_Arr_glb = Create2DLongdouble(2, nTotCycleNum_glb, arr_name);
        }
    nRDF_TotComps_glb = 2 + nBeadTypes_glb + nBeadTypes_glb * nBeadTypes_glb;
    nRDF_TotComps_glb /= 2;
    if (naReportFreqs_glb[REPORT_RDFTOT] != 0)
        {

            strcpy(arr_name, "ldTotRDFArr");
            ldaTOTRDF_Arr_glb = Create1DLongdouble((nTotCycleNum_glb * nRDF_TotComps_glb * nRDF_TotBins_glb), arr_name);

            strcpy(arr_name, "ldRDFArr");
            ldaRDF_Arr_glb = Create1DLongdouble((nRDF_TotComps_glb * nRDF_TotBins_glb), arr_name);
        }

    if (naReportFreqs_glb[REPORT_COMDEN] != 0)
        {
            nRadDen_TotComps_glb  = 2 * tot_chain_types_glb * (tot_chain_types_glb + 1);
            nRadDen_CompShift_glb = tot_chain_types_glb * (tot_chain_types_glb + 1);

            strcpy(arr_name, "ldRadDenArr");
            ldaRadDen_Arr_glb = Create1DLongdouble(nRadDen_TotComps_glb * nRDF_TotBins_glb, arr_name);

            strcpy(arr_name, "ldTotRadDenArr");
            ldaTOTRadDen_Arr_glb =
                Create1DLongdouble(nRadDen_TotComps_glb * nTotCycleNum_glb * nRDF_TotBins_glb, arr_name);
        }

    if (nTrajMode_glb == -1)
        {
            nTraj_FramesPerCycle_glb = nMCStepsPerCycle_glb / naReportFreqs_glb[REPORT_CONFIG] + 1;
            printf("\n***********************************************\n");
            printf("This feature is still experimental!\n");
            printf("\n***********************************************\n");

            strcpy(arr_name, "nTotTrajArr");
            naTOTTRAJ_Arr_glb = Create1DInt(nTraj_FramesPerCycle_glb * tot_beads_glb * BEADINFO_MAX, arr_name);
        }

    printf("Successfully allocated memory! Arrays initialized.\n");
}

void Memory_Allocate_NeighborLists(void)
{
    char arr_name[50];

    int num_of_points;

    num_of_points = 1 * 2 + 1;
    num_of_points = num_of_points * num_of_points * num_of_points;
    strcpy(arr_name, "naOldOvlpNeighs_glb");
    naOldOvlpNeighs_glb = Create1DInt(num_of_points, arr_name);
    strcpy(arr_name, "naNewOvlpNeighs_glb");
    naNewOvlpNeighs_glb = Create1DInt(num_of_points, arr_name);

    num_of_points = LARGEST_RADIUS * 2 + 1;
    num_of_points = num_of_points * num_of_points * num_of_points;
    strcpy(arr_name, "naOldContNeighs_glb");
    naOldContNeighs_glb = Create1DInt(num_of_points, arr_name);
    strcpy(arr_name, "naNewContNeighs_glb");
    naNewContNeighs_glb = Create1DInt(num_of_points, arr_name);
}

void Memory_VerifyMalloc(void)
{
    if (laClusHistList_glb == NULL)
        {
            printf("laClusHistList_glb malloc failed\n");
            exit(1);
        }
    if (naChainCheckList_glb == NULL)
        {
            printf("naChainCheckList_glb malloc failed\n");
            exit(1);
        }
    if (ldaRDF_Arr_glb == NULL && naReportFreqs_glb[REPORT_RDFTOT] != 0)
        {
            printf("ldaRDF_Arr_glb malloc failed\n");
            exit(1);
        }
    if (ldaTOTRDF_Arr_glb == NULL && naReportFreqs_glb[REPORT_RDFTOT] != 0)
        {
            printf("ldaTOTRDF_Arr_glb malloc failed\n");
            exit(1);
        }
    if (ldaRadDen_Arr_glb == NULL && naReportFreqs_glb[REPORT_COMDEN] != 0)
        {
            printf("ldaRadDen_Arr_glb malloc failed\n");
            exit(1);
        }
    if (ldaTOTRadDen_Arr_glb == NULL && naReportFreqs_glb[REPORT_COMDEN] != 0)
        {
            printf("ldaTOTRadDen_Arr_glb malloc failed\n");
            exit(1);
        }
    if (ldaMOLCLUS_Arr_glb == NULL && naReportFreqs_glb[REPORT_NETWORK] != 0)
        {
            printf("ldaMOLCLUS_Arr_glb malloc failed\n");
            exit(1);
        }
    if (ldaTOTMOLCLUS_Arr_glb == NULL && naReportFreqs_glb[REPORT_NETWORK] != 0)
        {
            printf("ldaTOTMOLCLUS_Arr_glb malloc failed\n");
            exit(1);
        }
    if ((naTOTTRAJ_Arr_glb == NULL) && (nTrajMode_glb == -1))
        {
            printf("naTOTTRAJ_Arr_glb malloc failed\n");
            exit(1);
        }
}

void Check_BeadTypewiseInteractions(void)
{
    int i, j;
    for (i = 0; i < MAX_AA; i++)
        {
            nBeadTypeIsSticker_glb[i] = 0; // Assume beads don't rotationally interact.
            nBeadTypeCanOvlp_glb[i]   = 0; // Assume beads don't have an overlap cost
            nBeadTypeCanCont_glb[i]   = 0; // Assume beads don't have contact costs.
            nBeadTypeCanFSol_glb[i]   = 0; // Assume beads have no stiffness
            nBeadTypeCanTInd_glb[i]   = 0; // Assume beads have no Temperature independent interactions
        }
    // Assume system has no interactions
    bSystemHasSCSC_glb = 0;
    bSystemHasOvlp_glb = 0;
    bSystemHasCont_glb = 0;
    bSystemHasFSol_glb = 0;
    bSystemHasTopo_glb = 0;

    // SC-SC Check
    for (i = 0; i < nBeadTypes_glb; i++)
        {
            for (j = 0; j < nBeadTypes_glb; j++)
                {
                    if (faEnergy_glb[i][j][E_SC_SC])
                        {
                            nBeadTypeIsSticker_glb[i] = 1;
                        }
                }
        }

    for (i = 0; i < nBeadTypes_glb; i++)
        {
            if (nBeadTypeIsSticker_glb[i])
                {
                    bSystemHasSCSC_glb = 1;
                }
        }

    // OVLP Check
    for (i = 0; i < nBeadTypes_glb; i++)
        {
            for (j = 0; j < nBeadTypes_glb; j++)
                {
                    if (faEnergy_glb[i][j][E_OVLP])
                        {
                            nBeadTypeCanOvlp_glb[i] = 1;
                        }
                }
        }
    for (i = 0; i < nBeadTypes_glb; i++)
        {
            if (nBeadTypeCanOvlp_glb[i])
                {
                    bSystemHasOvlp_glb = 1;
                }
        }

    // CONT Check
    for (i = 0; i < nBeadTypes_glb; i++)
        {
            for (j = 0; j < nBeadTypes_glb; j++)
                {
                    if (faEnergy_glb[i][j][E_CONT])
                        {
                            nBeadTypeCanCont_glb[i] = 1;
                        }
                }
        }
    for (i = 0; i < nBeadTypes_glb; i++)
        {
            if (nBeadTypeCanCont_glb[i])
                {
                    bSystemHasCont_glb = 1;
                }
        }

    // FSOL Check
    for (i = 0; i < nBeadTypes_glb; i++)
        {
            for (j = 0; j < nBeadTypes_glb; j++)
                {
                    if (faEnergy_glb[i][j][E_F_SOL])
                        {
                            nBeadTypeCanFSol_glb[i] = 1;
                        }
                }
        }
    for (i = 0; i < nBeadTypes_glb; i++)
        {
            if (nBeadTypeCanFSol_glb[i])
                {
                    bSystemHasFSol_glb = 1;
                }
        }

    // T_IND Check
    for (i = 0; i < nBeadTypes_glb; i++)
        {
            for (j = 0; j < nBeadTypes_glb; j++)
                {
                    if (faEnergy_glb[i][j][E_T_IND])
                        {
                            nBeadTypeCanTInd_glb[i] = 1;
                        }
                }
        }

    // STIFF Check
    for (i = 0; i < nBeadTypes_glb; i++)
        {
            for (j = 0; j < nBeadTypes_glb; j++)
                {
                    if (faEnergy_glb[i][j][E_STIFF])
                        {
                            bSystemHasTopo_glb = 1;
                        }
                }
        }
}

void Global_Array_Initialization_AtStart(void)
{
    int i, j, k, xTemp, yTemp, zTemp; // Indecies
    int myCubeLen = 3;

    for (i = 0; i < naBoxSize_glb[0] * naBoxSize_glb[1] * naBoxSize_glb[2]; i++)
        {                             // Initializing the lattice
            naTotLattice_glb[i] = -1; // If -1, then there is no bead there.
        }
    i = 0;
    // Constructing the local array to search for bonding partners.
    for (xTemp = -myCubeLen; xTemp <= myCubeLen; xTemp++)
        {
            for (yTemp = -myCubeLen; yTemp <= myCubeLen; yTemp++)
                {
                    for (zTemp = -myCubeLen; zTemp <= myCubeLen; zTemp++)
                        {
                            if ((xTemp != 0 || yTemp != 0 || zTemp != 0)
                                //&& (xTemp*xTemp + yTemp*yTemp + zTemp*zTemp > 3)
                                && (xTemp * xTemp + yTemp * yTemp + zTemp * zTemp <= 3))
                                {
                                    naLocalArr_glb[i][0] = xTemp;
                                    naLocalArr_glb[i][1] = yTemp;
                                    naLocalArr_glb[i][2] = zTemp;
                                    i++;
                                }
                        }
                }
        }

    // Initializing rotational orientational bias arrays.
    for (i = 0; i < MAX_VALENCY; i++)
        {
            for (j = 0; j < MAX_ROTSTATES; j++)
                {
                    naRotTrial_glb[i][j] = -1;
                }
        }
    for (i = 0; i < MAX_ROTSTATES - 1; i++)
        {
            naRot_IndArr_glb[i] = i;
        }
    // Initializing arrays required for cluster analyses.
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naChainCheckList_glb[i] = 0;
            laClusHistList_glb[i]   = 0;
            for (j = 0; j <= tot_chains_glb; j++)
                {
                    naClusterMatrix_glb[i][j] = -1;
                }
        }

    if (naReportFreqs_glb[REPORT_NETWORK] != 0)
        {
            for (i = 0; i < tot_chains_glb; i++)
                { // For MolWise Clus arrays
                    for (j = 0; j < tot_chain_types_glb; j++)
                        {
                            ldaMOLCLUS_Arr_glb[MolClusArr_Index(0, j, i)] = 0.;
                            for (k = 0; k < nTotCycleNum_glb; k++)
                                {
                                    ldaTOTMOLCLUS_Arr_glb[MolClusArr_Index(k, j, i)] = 0.;
                                }
                        }
                }
            for (k = 0; k < nTotCycleNum_glb; k++)
                { // For GyrRad
                    ldaTOTRg_Arr_glb[k][0] = 0.;
                    ldaTOTRg_Arr_glb[k][1] = 0.;
                }
            for (k = 0; k < nTotCycleNum_glb; k++)
                { // For ClusterDists
                    for (i = 0; i <= tot_chains_glb; i++)
                        {
                            ldaTOTCLUS_Arr_glb[k][i] = 0.;
                        }
                }
        }

    if (naReportFreqs_glb[REPORT_RDFTOT] != 0)
        {
            for (i = 0; i < nRDF_TotBins_glb; i++)
                { // For Radial distributions
                    for (j = 0; j < nRDF_TotComps_glb; j++)
                        { // For RDFs
                            ldaTOTRDF_Arr_glb[RDFArr_Index(0, j, i)] = 0.;
                            for (k = 0; k < nTotCycleNum_glb; k++)
                                {
                                    ldaTOTRDF_Arr_glb[RDFArr_Index(k, j, i)] = 0.;
                                }
                        }
                }
        }

    if (naReportFreqs_glb[REPORT_COMDEN] != 0)
        {
            for (j = 0; j < nRadDen_TotComps_glb; j++)
                { // For Density Dists wrt COM
                    for (i = 0; i < nRDF_TotBins_glb; i++)
                        { // For Radial distributions
                            ldaRadDen_Arr_glb[RadDenArr_Index(0, j, i)] = 0.;
                            for (k = 0; k < nTotCycleNum_glb; k++)
                                {
                                    ldaTOTRadDen_Arr_glb[RadDenArr_Index(k, j, i)] = 0.;
                                }
                        }
                }
        }

    // Setting counters
    Reset_Counters();

    // Checking which beads interaction how.

    Check_BeadTypewiseInteractions();

    for (i = MV_NULL + 2; i < MAX_MV; i++)
        {
            faMCFreq_glb[i] += faMCFreq_glb[i - 1]; // Cumulative Frequencies
        }
    for (i = 0; i < MAX_MV; i++)
        { // Zero out all the MCAccepts
            naMCAccepMat_glb[0][i] = 0;
            naMCAccepMat_glb[1][i] = 0;
        }

    for (i = 0; i < nTotCycleNum_glb; i++)
        {
            faKT_Cycle_glb[i] = fKT_glb + (float) i * fDeltaTemp_glb;
        }

    if (nTemp_inv_glb == 1)
        {
            for (i = 0; i < nTotCycleNum_glb; i++)
                {
                    faKT_Cycle_glb[i] = 1.f / faKT_Cycle_glb[i];
                }
        }

    fSquishRad_Sq_glb               = fSquishRad_glb * fSquishRad_glb;
    ldLogOfSmallestPossibleProb_glb = logl((lLDub) 1. / (lLDub) RAND_MAX);
    nLimitedClusterSize_glb         = (int) tot_chains_glb / 2;

    printf("All setup has been completed!\n");
}

void Reset_Counters(void)
{
    faSysGyrRad_glb             = 0.f;
    nTotGyrRadCounter_glb       = 0;
    nTotRDFCounter_glb          = 0;
    nTotRadDenCounter_glb       = 0;
    nTotClusCounter_glb         = 0;
    nLargestClusterRightNow_glb = 0;
    nTrajCurFrame_glb           = 0;
}

/// Reset_Global_Arrays - resets the various global counters and arrays for
/// system analysis.
void Reset_Global_Arrays(void)
{
    // Zero-ing out all the arrays used for data tracking!
    int i, j;
    for (i = 0; i <= tot_chains_glb; i++)
        {
            naChainCheckList_glb[i] = -1;
            laClusHistList_glb[i]   = 0;
            for (j = 0; j <= tot_chains_glb; j++)
                {
                    naClusterMatrix_glb[i][j] = -1;
                }
            if (naReportFreqs_glb[REPORT_NETWORK])
                {
                    for (j = 0; j < tot_chain_types_glb; j++)
                        {
                            ldaMOLCLUS_Arr_glb[MolClusArr_Index(0, j, i)] = 0.;
                        }
                }
        }
    // Initializing arrays for pair-distribution calculations
    if (naReportFreqs_glb[REPORT_RDFTOT])
        {
            for (j = 0; j < nRDF_TotComps_glb; j++)
                {
                    for (i = 0; i < nRDF_TotBins_glb; i++)
                        {
                            ldaRDF_Arr_glb[RDFArr_Index(0, j, i)] = 0.;
                        }
                }
        }
    if (naReportFreqs_glb[REPORT_COMDEN])
        {
            // Initalizing for density histograms wrt to the COM
            for (j = 0; j < nRadDen_TotComps_glb; j++)
                {
                    for (i = 0; i < nRDF_TotBins_glb; i++)
                        {
                            ldaRadDen_Arr_glb[RadDenArr_Index(0, j, i)] = 0.;
                        }
                }
        }
    // Setting counters
    Reset_Counters();
    // MCAccept arrays.
    for (i = 0; i < MAX_MV; i++)
        {
            naMCAccepMat_glb[0][i] = 0;
            naMCAccepMat_glb[1][i] = 0;
        }
}

/// Initial_Conditions_Simple - randomly place all the molecules.
void Initial_Conditions_Simple(void)
{
    /*
    Generates initial conditions for a system with topology information. The
    idea is to see if the selected bead has been placed or not. If it has been
    placed, we move on to sequentially placing all of it's bonded beads within
    linker constraints. Therefore, each bead shall sprout forth its partners.
    Since I like goto statements, they have been used. If a bead has not been
    placed, it's first checked if one of its bonding partners might have been
    placed, and will thus act as an anchor. This means that for each chain we
    need to know which beads have already been placed, thus a temporary list is
    used.
    */

    int idx, idy;      // Just internal counters for various things.
    int i, j, k;       // Iterators for loops
    int bondPart;      // Used to track the bond partner for a particular bead.
    int fB, lB;        // Used to track the first and last bead of a particular chain
                       // for looping.
    int radUp, radLow; // Used as radii to generate coordinates within a
                       // linker-length sphere
    // Just remember tha lB-1 is the last bead in that chain, but in loops we go
    // <lB.
    int TlCnt = 0;                     // Counter for trials around an anchor.
    int tmpR[POS_MAX], tmpR2[POS_MAX]; // Arrays to store temporary coordinates;
    int temp_list[MAX_CHAINLEN];       // Used to store already placed beads.
    int list_it = 0;                   // Iterator specifically for temp_list.

    for (i = 0; i < MAX_CHAINLEN; i++)
        { // initialize the list where -1 signifies emptiness.
            temp_list[i] = -1;
        }
    for (j = 0; j < POS_MAX; j++)
        { // Initialize the coordinate arrays to some numbers.
            tmpR[j]  = rand() % naBoxSize_glb[j];
            tmpR2[j] = tmpR[j];
        }

    for (k = 0; k < tot_chains_glb; k++)
        {
            // This makes reading easier, honestly
            fB = chain_info_glb[k][CHAIN_START];
            lB = fB + chain_info_glb[k][CHAIN_LENGTH];

            // Reset the temp_list for each chain, and the iterator!
            for (i = 0; i < list_it; i++)
                { // initialize the list where -1 signifies emptiness.
                    temp_list[i] = -1;
                }
            list_it = 0;

            for (i = fB; i < lB; i++)
                { // Going bead-bead in this chain.
                    // Checking if the bead has already been placed.
                    for (idy = 0; idy < list_it; idy++)
                        {
                            if (i == temp_list[idy])
                                { // This means this bead has already
                                  // been placed.
                                    goto SproutForth;
                                }
                        }
                    // If we have reached here, this bead has not been placed before.
                    idx      = 0;
                    bondPart = topo_info_glb[i][idx]; // Re-initialize these two.
                    // Let's go through the bonded partners for this bead and see if one
                    // can be used as an anchor.
                    while (topo_info_glb[i][idx] != -1 && idx < MAX_BONDS)
                        { // Again, -1 signifies that no bonded
                          // bead exists.
                            bondPart = topo_info_glb[i][idx];
                            // Checking if thisBead has been placed already, and can be used
                            // as anchor.
                            for (idy = 0; idy < list_it; idy++)
                                {
                                    if (bondPart == temp_list[idy])
                                        {                     // Using bondPart as an anchor
                                            goto FoundAnchor; // This is just so we don't keep going
                                                              // over the list needlessly and just
                                                              // pick the first one.
                                        }
                                }
                            idx++;
                            // If we got here, this means no bonded partner has already been
                            // placed, so no anchor.
                            bondPart = -1;
                        }

                FoundAnchor:
                    if (bondPart != -1)
                        { // We found an anchor.
                            radUp  = 2 * linker_len_glb[i][idx] + 1;
                            radLow = linker_len_glb[i][idx];
                            for (j = 0; j < POS_MAX; j++)
                                { // Using thisBead as anchor
                                    tmpR[j]  = bead_info_glb[bondPart][j];
                                    tmpR2[j] = tmpR[j];
                                }
                            TlCnt = 0; // Reset this counter.
                            while (naTotLattice_glb[Lat_Ind_FromVec(tmpR2)] != -1 && TlCnt < 5000000)
                                {
                                    for (j = 0; j < POS_MAX; j++)
                                        {
                                            tmpR2[j] = (rand() % radUp) - radLow;
                                            tmpR2[j] = (tmpR[j] + tmpR2[j] + naBoxSize_glb[j]) % naBoxSize_glb[j];
                                        }
                                    TlCnt++;
                                }
                            if (TlCnt == 5000000)
                                {
                                    printf("Not enough space in the lattice. Crashing. Maybe "
                                           "try increasing max trials, "
                                           "or make the box bigger!\t\n");
                                    exit(1);
                                }
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    bead_info_glb[i][j] = tmpR2[j];
                                }                                         // Placing bead
                            naTotLattice_glb[Lat_Ind_FromVec(tmpR2)] = i; // Putting on lattice
                            temp_list[list_it]                       = i;
                            list_it++; // Remembering in hash list.
                            // The bead has been placed around an appropriate anchor.
                        }
                    else
                        { // There was no appropriate anchor. So we can put this bead
                          // wherever.
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    tmpR[j] = rand() % naBoxSize_glb[j];
                                }
                            // This usually means this is the first bead in the chain.
                            TlCnt = 0;
                            while (naTotLattice_glb[Lat_Ind_FromVec(tmpR)] != -1 && TlCnt < 5000000)
                                { // Keep looping till we find
                                  // an empty lattice site.
                                    TlCnt++;
                                    for (j = 0; j < POS_MAX; j++)
                                        { // Generate a random point in the lattice.
                                            tmpR[j] = rand() % naBoxSize_glb[j];
                                        }
                                }
                            if (TlCnt == 5000000)
                                {
                                    printf("\n\nNot enough space in the lattice for first "
                                           "bead, bruh. Maybe try increasing max trials, "
                                           "or make the box bigger!\n\n");
                                    exit(1);
                                }
                            // Placing this bead wherever there is space, and using it as
                            // anchor. Also updating the lattice, and temp_list
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    bead_info_glb[i][j] = tmpR[j];
                                }
                            naTotLattice_glb[Lat_Ind_FromVec(tmpR)] = i; // Placing on lattice!
                            temp_list[list_it]                      = i;
                            list_it++; // Remembering that this bead has been placed.
                        }
                SproutForth:
                    // Now going over the list of bondParts for i and sprouting them
                    for (j = 0; j < POS_MAX; j++)
                        { // Making the current bead an anchor.
                            tmpR[j] = bead_info_glb[i][j];
                        }
                    idx      = 0;                     // Resets things!
                    bondPart = topo_info_glb[i][idx]; // Resets things!
                    while (topo_info_glb[i][idx] != -1 && idx < MAX_BONDS)
                        { // Again, -1 signifies that no bonded
                          // bead exists.
                            bondPart = topo_info_glb[i][idx];
                            // If the bead has already been dealt with, move on to next
                            // potential bondPartner
                            for (idy = 0; idy < list_it; idy++)
                                {
                                    if (bondPart == temp_list[idy])
                                        {
                                            goto SkipThisPartner; // Just so we don't keep looking
                                                                  // further.
                                        }
                                }
                            // Get the radii for placing this bead around bead i
                            radUp  = 2 * linker_len_glb[i][idx] + 1;
                            radLow = linker_len_glb[i][idx];
                            // Initializing new vector to be where bead i is.
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    tmpR2[j] = tmpR[j];
                                }
                            TlCnt = 0; // Reset this counter.
                            while (naTotLattice_glb[Lat_Ind_FromVec(tmpR2)] != -1 && TlCnt < 5000000)
                                {
                                    for (j = 0; j < POS_MAX; j++)
                                        {
                                            tmpR2[j] = (rand() % radUp) - radLow;
                                            tmpR2[j] = (tmpR[j] + tmpR2[j] + naBoxSize_glb[j]) % naBoxSize_glb[j];
                                        }
                                    // printf("%d %d %d\n", tmpR2[0], tmpR2[1], tmpR2[2]);
                                    TlCnt++;
                                }
                            if (TlCnt == 5000000)
                                {
                                    printf("\n\nNot enough space in the lattice. Crashing. "
                                           "Maybe try increasing max trials, "
                                           "or make the box bigger!\t %d\t%d %d\n\n",
                                           TlCnt, i, naTotLattice_glb[Lat_Ind_FromVec(tmpR2)]);
                                    exit(1);
                                }
                            for (j = 0; j < POS_MAX; j++)
                                { // Placing bead
                                    bead_info_glb[bondPart][j] = tmpR2[j];
                                }

                            naTotLattice_glb[Lat_Ind_FromVec(tmpR2)] = bondPart; // Putting on lattice
                            temp_list[list_it]                       = bondPart;
                            list_it++; // Remembering in hash list.
                        SkipThisPartner:
                            idx++;
                            bondPart = -1;
                        }
                    // Moving on to bead i+1
                }
            // Moving on to the next molecule!
        }
    // Initializing all beads with no physical bonds.
    for (i = 0; i < tot_beads_glb; i++)
        {
            bead_info_glb[i][BEAD_FACE] = -1;
        }
}

/// Initial_Conditions_FromFile - reads the restart file and sets the initial
/// conditions for the system
void Initial_Conditions_FromFile(void)
{
    printf("Reading file to generate initial configuration\n");
    char strLine[1000];
    int i, j;
    int tmp_beadinfo[BEADINFO_MAX];
    FILE* infile;
    infile       = fopen(strRestartFile_glb, "r");
    int tmpBeads = 0;
    while (fgets(strLine, sizeof(strLine), infile) != NULL)
        {
            tmpBeads++;
            if (strcmp(strLine, "ITEM: ATOMS id type mol x y z bP\n") == 0)
                {
                    break;
                }
        }
    tmpBeads = 0;
    while (fgets(strLine, sizeof(strLine), infile) != NULL)
        {
            sscanf(strLine, "%d %d %d %d %d %d %d", &i, &tmp_beadinfo[BEAD_TYPE], &tmp_beadinfo[BEAD_CHAINID],
                   &tmp_beadinfo[POS_X], &tmp_beadinfo[POS_Y], &tmp_beadinfo[POS_Z], &tmp_beadinfo[BEAD_FACE]);
            tmpBeads++;
            for (j = 0; j < BEADINFO_MAX; j++)
                {
                    bead_info_glb[i][j] = tmp_beadinfo[j];
                }
        }
    fclose(infile);
    if (tmpBeads != tot_beads_glb)
        {
            printf("Given restart file does not have the same number of beads as "
                   "the provided structure.");
            printf("Crashing!\n");
            exit(1);
        }
    for (i = 0; i < tot_beads_glb; i++)
        {
            naTotLattice_glb[Lat_Ind_OfBead(i)] = i;
        }
    printf("File read successfully and initial configuration set.\n");
}

/// Initial_Conditions_BreakBonds - glorified for loop to break every bond. Is
/// used when a restart file with thermalization steps.
void Initial_Conditions_BreakBonds(void)
{
    int i;
    for (i = 0; i < tot_beads_glb; i++)
        {
            bead_info_glb[i][BEAD_FACE] = -1;
        }
    printf("Since system has thermalization cycle, deleting all physical "
           "bonds!\n");
}

/// Calculate_Rot_Bias - calculates the possible anisotropic Boltzmann weights
/// at the beginning so that we have a look-up table. \param CurrentTemp
void Calculate_Rot_Bias(const float CurrentTemp)
{

    int i, j;
    for (i = 0; i < MAX_AA; i++)
        {
            for (j = 0; j < MAX_AA; j++)
                {
                    ldaBoltzFacNorm_glb[i][j] =
                        (lLDub) expl(-(lLDub) faEnergy_glb[i][j][E_SC_SC] / (lLDub) CurrentTemp);
                    // printf("%LE; ", ldaBoltzFacNorm_glb[i][j]);
                }
            // printf("\n");
        }
}

/// Temperature_Function - used to calculate what fCuTemp_glb should (current temperature) based on how long the
/// run has been going.
/// \param mode - which
/// of the four functions to use. Note that the various modes are described
/// below where $F(t)$ is the returned value, and t == nGen.
/// \param nGen - basically 'time' or how
/// many MC Steps have been done.
/// \return end_val - the current temperature.
float Temperature_Function(const int mode, const long nGen)
{

    float x_val;
    float end_val;

    switch (mode)
        {
            case 1:
                x_val   = -(float) (nGen);
                x_val   = 4.f * x_val;
                x_val   = x_val / (float) (nMCStepsForTherm_glb);
                end_val = fKT_glb + fPreKT_glb * expf(x_val);

                break;

            case 2:
                x_val   = -(float) (nGen);
                x_val   = fMC_TempRate_glb * x_val;
                x_val   = x_val / (float) (nMCStepsForTherm_glb);
                end_val = fKT_glb + expf(x_val);

                break;

            default:

                end_val = fKT_glb;
                break;
        }

    if (end_val - fKT_glb < 0.005)
        {
            puts("\n\n******************************\n");
            puts("Annealing Is Being Turned Off");
            puts("\n******************************\n\n");
            nAnnealing_Mode_glb = -1;

            if (nBiasPotential_CoupledToTemp_glb)
                {
                    BiasPotential_TurnOFF();
                }
        }

    return end_val;
}

void BiasPotential_TurnOFF(void)
{
    puts("******************************");
    puts("Bias Is Being Turned Off");
    puts("******************************");
    nBiasPotential_Mode_glb = -1;
}

/// Create1DInt create an array of ints of size of xDim.
/// \param xDim Length of array.
/// \param ArrName Dummy name of the array for debugging.
/// \return The pointer to the array.
int* Create1DInt(const size_t xDim, const char* ArrName)
{
    int* dumPtr;
#if DEBUG_BUILD
    printf("%s is being allocated ... ", ArrName);
#endif
    dumPtr = (int*) malloc(xDim * sizeof(int));

    if (dumPtr == NULL)
        {
            printf("Not enough memory!\n%s\nCrashing!\n", ArrName);
            exit(1);
        }
#if DEBUG_BUILD
    printf(" ... and allocated!\n");
#endif
    return dumPtr;
}

/// Create1DLong create an array of ints of size of xDim.
/// \param xDim Length of array.
/// \param ArrName Dummy name of the array for debugging.
/// \return The pointer to the array.
long* Create1DLong(const size_t xDim, const char* ArrName)
{
    long* dumPtr;
#if DEBUG_BUILD
    printf("%s is being allocated ... ", ArrName);
#endif
    dumPtr = (long*) malloc(xDim * sizeof(long));

    if (dumPtr == NULL)
        {
            printf("Not enough memory!\n%s\nCrashing!\n", ArrName);
            exit(1);
        }
#if DEBUG_BUILD
    printf(" ... and allocated!\n");
#endif
    return dumPtr;
}

/// Create1DLong create an array of ints of size of xDim.
/// \param xDim Length of array.
/// \param ArrName Dummy name of the array for debugging.
/// \return The pointer to the array.
float* Create1DFloat(const size_t xDim, const char* ArrName)
{
    float* dumPtr;
#if DEBUG_BUILD
    printf("%s is being allocated ... ", ArrName);
#endif
    dumPtr = (float*) malloc(xDim * sizeof(float));

    if (dumPtr == NULL)
        {
            printf("Not enough memory!\n%s\nCrashing!\n", ArrName);
            exit(1);
        }
#if DEBUG_BUILD
    printf(" ... and allocated!\n");
#endif
    return dumPtr;
}

/// Create1DLong create an array of ints of size of xDim.
/// \param xDim Length of array.
/// \param ArrName Dummy name of the array for debugging.
/// \return The pointer to the array.
long double* Create1DLongdouble(const size_t xDim, const char* ArrName)
{
    long double* dumPtr;
#if DEBUG_BUILD
    printf("%s is being allocated ... ", ArrName);
#endif
    dumPtr = (long double*) malloc(xDim * sizeof(long double));

    if (dumPtr == NULL)
        {
            printf("Not enough memory!\n%s\nCrashing!\n", ArrName);
            exit(1);
        }
#if DEBUG_BUILD
    printf(" ... and allocated!\n");
#endif
    return dumPtr;
}

/// Create2DInt creates a 2D array or arbitrary size of dimensions [yDim][xDim],
/// in C-style. We first allocate an array of pointers which have size yDim.
/// Call this yArr[] Using the first pointer, we then allocate a new array of
/// size xDim*yDim. Then for each of the pointers in yArr, we point to different
/// chunks of totAr; converting 2D indicies to 1D indicies. \param xDim: Size of
/// second index \param yDim: Size of first index \param ArrName: Dummy name for
/// the array, for debugging. \return Pointer to array-of-pointers, or a 2D
/// array.
int** Create2DInt(const size_t xDim, const size_t yDim, const char* ArrName)
{
#if DEBUG_BUILD
    printf("%s is being allocated ... ", ArrName);
#endif
    int** dumPtr;
    dumPtr = (int**) malloc(yDim * sizeof(int*));
    if (dumPtr == NULL)
        {
            printf("Not enough memory!\n%s\nCrashing!\n", ArrName);
            exit(1);
        }

    dumPtr[0] = (int*) malloc(xDim * yDim * sizeof(int));
    if (dumPtr[0] == NULL)
        {
            printf("Not enough memory!\n%s\nCrashing!\n", ArrName);
            exit(1);
        }

    size_t i;
    for (i = 1; i < yDim; i++)
        {
            dumPtr[i] = &dumPtr[0][i * xDim];
        }
#if DEBUG_BUILD
    printf(" ... and allocated!\n");
#endif
    return dumPtr;
}

/// Create2DInt creates a 2D array of arbitrary size of dimensions [yDim][xDim],
/// in C-style. The array is for long doubles We first allocate an array of
/// pointers which have size yDim. Call this yArr[] Using the first pointer, we
/// then allocate a new array of size xDim*yDim. Then for each of the pointers
/// in yArr, we point to different chunks of totAr; converting 2D indicies to 1D
/// indicies. \param xDim: Size of second index \param yDim: Size of first index
/// \param ArrName: Dummy name for the array, for debugging.
/// \return Pointer to array-of-pointers, or a 2D array.
long double** Create2DLongdouble(const size_t xDim, const size_t yDim, const char* ArrName)
{
#if DEBUG_BUILD
    printf("%s is being allocated ... ", ArrName);
#endif
    long double** dumPtr;
    dumPtr = (long double**) malloc(yDim * sizeof(long double*));
    if (dumPtr == NULL)
        {
            printf("Not enough memory!\n%s\nCrashing!\n", ArrName);
            exit(1);
        }

    dumPtr[0] = (long double*) malloc(xDim * yDim * sizeof(long double));
    if (dumPtr[0] == NULL)
        {
            printf("Not enough memory!\n%s\nCrashing!\n", ArrName);
            exit(1);
        }

    size_t i;
    for (i = 1; i < yDim; i++)
        {
            dumPtr[i] = &dumPtr[0][i * xDim];
        }
#if DEBUG_BUILD
    printf(" ... and allocated!\n");
#endif
    return dumPtr;
}
