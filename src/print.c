#include "print.h"
#include "initialize.h"
#include "cluster.h"
#include "energy.h"
#include "structure.h"

/// PrintToScreen_EnergyMatrix - fancy function to print symmetric matrices with labels
/// \param strTitle
/// \param nSeqEn
/// \param fArray
/// \param param
void PrintToScreen_EnergyMatrix(char* strTitle, int nSeqEn, float fArray[MAX_AA][MAX_AA][MAX_E], int param)
{
    const int nLen   = nSeqEn;
    const int outLen = (nLen + 3) * 5;
    int i, j;

    char sSectionHead[512];
    memset(sSectionHead, '-', outLen);
    sSectionHead[outLen] = '\0';

    printf("%s\n", sSectionHead);

    printf("| %-5s ", strTitle);
    for (i = 0; i < nLen; i++)
        {
            printf("%4d  ", i);
        }
    printf("%3s\n", "|");

    for (i = 0; i < nLen; i++)
        {
            printf("| %-5d ", i);
            for (j = 0; j < nLen; j++)
                {
                    printf("%5.2f ", fArray[i][j][param]);
                }
            printf("%3s\n", "|");
        }

    printf("%s\n", sSectionHead);
}

/// TrajArr_Index
/// \param beadID
/// \param nFrameNumber
/// \param beadProp
/// \return
long TrajArr_Index(const int beadID, const int nFrameNumber, const int beadProp)
{
    return beadProp + BEADINFO_MAX * (beadID + tot_beads_glb * nFrameNumber);
}

/// Write_ClusterDist - write the cluster histogram to a separate file. NOT USED
/// ANYMORE, but might find a use later Always appends to the file for this run.
/// Stopped using it because the IO load was slowing things down at our cluster
/// Would be a good way to gather proper statistics on clusters as the runs went
/// on \param filename \param nGen
void Write_ClusterDist(char* filename, long nGen)
{
    FILE* fp;
    int i;
    if (nGen == -1)
        {
            fp = fopen(filename, "w"); // overwrite
        }
    else
        {
            fp = fopen(filename, "a");
        }

    if (nGen == -1)
        { // title
            fprintf(fp, "#Step followed by histogram\n");
        }
    else
        {
            fprintf(fp, "#%ld\n", nGen);
            for (i = 0; i <= tot_chains_glb; i++)
                {
                    fprintf(fp, "%ld\t", laClusHistList_glb[i]);
                }
            fprintf(fp, "\n");
        }

    fclose(fp);
}

/// Write_GyrTen - write the total gyration tensor to a file. Remember that we
/// only have 6 things in a symmetric 3x3 tensor. Always appends to the file for
/// this run. Stopped using it because the IO load was slowing things down at
/// our cluster Would be a good way to gather proper statistics on the shapes
/// and such as the runs went on \param filename \param nGen
void Write_GyrTen(char* filename, long nGen)
{
    FILE* fp;
    int i;
    if (nGen == -1)
        {
            fp = fopen(filename, "w"); // overwrite
        }
    else
        {
            fp = fopen(filename, "a");
        }

    if (nGen == -1)
        { // title
            fprintf(fp, "#Step followed by Gyration Tensor Vals: xx yy zz xy yz xz\n");
        }
    else
        {
            fprintf(fp, "#%ld\n", nGen);
            for (i = 0; i < 7; i++)
                {
                    if (i == 5)
                        {
                            continue;
                        } // Not smart enough to make a better indexing function
                    fprintf(fp, "%f\t", faGyrTensor_glb[i]);
                }
            fprintf(fp, "\n");
        }

    fclose(fp);
}

/// FileIOUtil_WriteHeader_ForMCMove - write the header for the MC move acceptance file.
/// \param fileName
void FileIOUtil_WriteHeader_ForMCMove(const char* fileName)
{
    FILE* fp = fopen(fileName, "a");
    fprintf(fp, "#Steps, and Temp, are followed by (rej, acc) for each MC Move.\n");
    fprintf(fp, "#(nStep, T)             "
                "STROT                   "
                "LOCAL                   "
                "COLOCAL                 "
                "MTLOCAL                 "
                "SNAKE                   "
                "TRANS                   "
                "SM_CLST                 "
                "CLST                    "
                "PIVOT                   "
                "BRROT                   "
                "DBPVT                   "
                "PR_CLST                 "
                "\n");
    fclose(fp);
}

/// FileIO_WriteTo_MCMoveFile - writes the acceptance and rejection ratios of the various moves. Just keeps appending to
/// that file. Can be used to track the 'dynamics' during a simulation. ALSO: Zeros out all the acceptance arrays.
/// \param filename
/// \param nGen
/// \param fMCTemp
void FileIO_WriteTo_MCMoveFile(const char* filename, const long nGen, float const fMCTemp)
{
    FILE* fp = fopen(filename, "a+");
    int i;                                           // Iterator
    fprintf(fp, "%-10ld  %-10.2f  ", nGen, fMCTemp); // Step and Temp
    for (i = 1; i < MAX_MV; i++)
        {
            fprintf(fp, "%-10ld  %-10ld  ", naMCAccepMat_glb[0][i], naMCAccepMat_glb[1][i]);
        }
    fprintf(fp, "\n");

    for (i = 1; i < MAX_MV; i++)
        {
            naMCAccepMat_glb[0][i] = 0;
            naMCAccepMat_glb[1][i] = 0;
            // This way the print function will zero out the matrix every time
            // we print to a file!
        }

    fclose(fp);
}

/// Print_LogToScreen - prints current status of the run. Overall move
/// acceptance ratios, and energies. \param nGen \param run_it
void Print_LogToScreen(long nGen, const int run_it)
{
    ScreenIO_Print_Log_FullRun(nGen, run_it);
}

void FileIOUtil_WriteHeader_ForRDF(const char* strFileName)
{
    FILE* fp = fopen(strFileName, "w");

    fprintf(fp, "#Split RDFs. ALL-ALL; DIAGONALS (0-0, 1-1, ..) and then off-diagonals (0-1, 0-2, ..., 1-2, ...)\n");

    fclose(fp);
}

void FileIOUtil_WriteHeader_ForCOMDen(const char* strFileName)
{
    FILE* fp = fopen(strFileName, "w");

    fprintf(fp, "#Density Distribution from the COM outwards. Order is ChainType.\n");

    fclose(fp);
}

void FileIOUtil_WriteHeader_ForGyrRad(const char* strFileName)
{
    FILE* fp = fopen(strFileName, "w");

    fprintf(fp, "#Gyration Radius Of ENTIRE System & BoxSize.\n");

    fclose(fp);
}

void FileIOUtil_WriteHeader_ForCLUS(const char* strFileName)
{
    FILE* fp = fopen(strFileName, "w");

    fprintf(fp, "#Cluster size histograms. First is largest cluster, then n=1, n=2, and so on.\n");

    fclose(fp);
}

void FileIOUtil_WriteHeader_ForMolClus(const char* strFileName)
{
    FILE* fp = fopen(strFileName, "w");

    fprintf(fp, "#Cluster size histograms, split by molecule type. Each row is a ChainType, and each column is the "
                "cluster size starting from 1.\n");

    fclose(fp);
}

/// FileIOUtil_WriteHeader_ForEnergy - write the header for the energy file.
/// \param fileName
void FileIOUtil_WriteHeader_ForEnergy(const char* fileName)
{
    FILE* fp = fopen(fileName, "a+");
    fprintf(fp, "#"
                "STEP         "
                "TOT         "
                "OVLP        "
                "CONT        "
                "SC-SC       "
                "FSOL        "
                "T_IND       "
                "STIFF       "
                "BIAS        "
                "\n");
    fclose(fp);
}

/// FileIO_AppendEnergyTo_EnergyFile - write the decomposed energy of the system to a file.
/// Just appends to the file for this run.
/// \param fileNameStr
/// \param nGen
void FileIO_AppendEnergyTo_EnergyFile(const char* fileNameStr, const long nGen)
{
    FILE* fp = fopen(fileNameStr, "a");

    int i;
    fprintf(fp, "%-12ld", nGen);
    for (i = 0; i < (MAX_E); i++)
        {
            fprintf(fp, "  %-10.2E", faCurrEn_glb[i]);
        }
    fprintf(fp, "\n");

    fclose(fp);
}

/// FileIO_Trajectory_AppendFrame -
/// \param fileNameStr
/// \param run_it
/// \param nGen
void FileIO_Trajectory_AppendFrame(const char* fileNameStr, const int run_it, const long nGen)
{
    if (nTrajMode_glb == 1)
        {
            FileIOUtil_Traj_Bin_AppendFrame_ToFile(fileNameStr, 0);
        }
    else
        {
            if (run_it == -1) // Thermalization cycle
                {
                    FileIOUtil_Traj_Txt_AppendFrame_ToFile(fileNameStr, nGen);
                }
            else
                {
                    FileIOUtil_Traj_Txt_AppendFrame_ToFile(fileNameStr,
                                                           nGen + nMCStepsForTherm_glb + run_it * nMCStepsPerCycle_glb);
                }
        }
}

/// FileIOUtil_Traj_Txt_AppendFrame_ToFile -  write a LAMMPS style trajectory file for the system.
/// Appends to the file for this run.
/// \param filename
/// \param nGen
void FileIOUtil_Traj_Txt_AppendFrame_ToFile(const char* filename, const long nGen)
{
    // Writes the trajectory in LAMMPS format. To be viewed with VMD (hopefully). Read the LAMMPS dump documentation for
    // the actual format of the file
    FILE* fp = fopen(filename, "a"); // We will append to the file made for this run;

    int i; // Looping index
    fprintf(fp, "ITEM: TIMESTEP\n");
    fprintf(fp, "%ld\n", nGen); // Timestep

    fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fp, "%ld\n", (size_t) tot_beads_glb); // Total atom number

    fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n"); // BCs are always periodic for now
    fprintf(fp, "0 %d\n0 %d\n0 %d\n", naBoxSize_glb[0], naBoxSize_glb[1], naBoxSize_glb[2]); // Box dimensions

    fprintf(fp, "ITEM: ATOMS id type mol x y z bP\n"); // What we are printing

    for (i = 0; i < tot_beads_glb; i++)
        {
            fprintf(fp, "%d %d %d %d %d %d %d\n", i, bead_info_glb[i][BEAD_TYPE], bead_info_glb[i][BEAD_CHAINID],
                    bead_info_glb[i][POS_X], bead_info_glb[i][POS_Y], bead_info_glb[i][POS_Z],
                    bead_info_glb[i][BEAD_FACE]);
        }

    fclose(fp);
}

/// TrajUtil_SubselectDataFromBeadInfo - Assuming that sub_beads is a 2D array with enough space (nAtoms * 4),
/// we extract POS_X, POS_Y, POS_Z, and BEAD_FACE from each bead.
/// \param sub_beads
/// \param input_beads
void TrajUtil_SubselectDataFromBeadInfo(int** restrict sub_beads, int** restrict input_beads)
{
    int i;

    for (i = 0; i < tot_beads_glb; i++)
        {
            sub_beads[i][0] = input_beads[i][POS_X];
            sub_beads[i][1] = input_beads[i][POS_Y];
            sub_beads[i][2] = input_beads[i][POS_Z];
            sub_beads[i][3] = input_beads[i][BEAD_FACE];
        }
}

/// FileIOUtil_Traj_Bin_AppendFrame_ToFile - Appends the current frame in binary format. The format is flat,
/// and for every bead, in order of beadID, we write X, Y, Z, bP
/// \param filename
/// \param nGen
void FileIOUtil_Traj_Bin_AppendFrame_ToFile(const char* filename, const long nGen)
{
    const int nCrds = 4;

    int** subBeadInfo = Create2DInt(nCrds, tot_beads_glb, "TrajBinSubArr");

    TrajUtil_SubselectDataFromBeadInfo(subBeadInfo, bead_info_glb);

    FILE* fp = fopen(filename, "ab");

    const int nBeads = (int) tot_beads_glb;

    fwrite(&nBeads, sizeof(int), 1, fp);
    fwrite(subBeadInfo[0], sizeof(int), nCrds * tot_beads_glb, fp);
    fclose(fp);

    free(subBeadInfo[0]);
    free(subBeadInfo);
}

/// Save_Trajectory -  saves the current position onto the total array
/// \param nGen
/// \param curFrame
void Save_Trajectory(const long nGen, const long curFrame)
{

    int i, j;
    long ar_idx;
    for (i = 0; i < tot_beads_glb; i++)
        {
            for (j = 0; j < BEADINFO_MAX; j++)
                {
                    ar_idx                    = TrajArr_Index(i, curFrame, j);
                    naTOTTRAJ_Arr_glb[ar_idx] = bead_info_glb[i][j];
                }
        }
}

/// Write_Saved_Trajectory -  write a LAMMPS style trajectory file for the
/// system. This funcion writes all the stored frames for this temperature cycle
/// \param filename
/// \param nGen
void Write_Saved_Trajectory(char* filename, const int run_it)
{
    // Writes the trajectory in LAMMPS format. To be viewed with VMD
    // (Hopefully). Read the LAMMPS dump documentation for the actual formate of
    // the file
    FILE* fp;
    fp = fopen(filename, "a"); // We will append to the file made for this run

    int i, j, k; // Looping index
    int ar_idx;

    for (i = 0; i < nTrajCurFrame_glb; i++)
        {
            fprintf(fp, "ITEM: TIMESTEP\n");
            if (run_it >= 0)
                {
                    fprintf(fp, "%ld\n",
                            nMCStepsForTherm_glb + run_it * nMCStepsPerCycle_glb +
                                i * naReportFreqs_glb[REPORT_CONFIG]); // Timestep
                }
            else
                {
                    fprintf(fp, "%ld\n", i * naReportFreqs_glb[REPORT_CONFIG]); // Timestep
                }

            fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
            fprintf(fp, "%ld\n", (size_t) tot_beads_glb); // Total atom number

            fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n"); // BCs are always periodic for now
            fprintf(fp, "0 %d\n0 %d\n0 %d\n", naBoxSize_glb[0], naBoxSize_glb[1], naBoxSize_glb[2]); // Box dimensions

            fprintf(fp, "ITEM: ATOMS id type mol x y z bP\n"); // What we are printing

            for (j = 0; j < tot_beads_glb; j++)
                {
                    fprintf(fp, "%d", j);
                    ar_idx = TrajArr_Index(j, i, BEAD_TYPE);
                    fprintf(fp, " %d", naTOTTRAJ_Arr_glb[ar_idx]);
                    ar_idx = TrajArr_Index(j, i, BEAD_CHAINID);
                    fprintf(fp, " %d", naTOTTRAJ_Arr_glb[ar_idx]);
                    ar_idx = TrajArr_Index(j, i, POS_X);
                    fprintf(fp, " %d", naTOTTRAJ_Arr_glb[ar_idx]);
                    ar_idx = TrajArr_Index(j, i, POS_Y);
                    fprintf(fp, " %d", naTOTTRAJ_Arr_glb[ar_idx]);
                    ar_idx = TrajArr_Index(j, i, POS_Z);
                    fprintf(fp, " %d", naTOTTRAJ_Arr_glb[ar_idx]);
                    ar_idx = TrajArr_Index(j, i, BEAD_FACE);
                    fprintf(fp, " %d", naTOTTRAJ_Arr_glb[ar_idx]);
                    fprintf(fp, "\n");
                }
        }

    fclose(fp);
}

/// PrintToScreen_AllEnergyMatrices
void PrintToScreen_AllEnergyMatrices(void)
{

    const char lBrace[] = "<======      ";
    const char rBrace[] = "      ======>";

    printf("%s Energy Matrices %s\n", lBrace, rBrace);

    PrintToScreen_EnergyMatrix("OVLP ", nBeadTypes_glb, faEnergy_glb, E_OVLP);

    PrintToScreen_EnergyMatrix("CONT ", nBeadTypes_glb, faEnergy_glb, E_CONT);

    PrintToScreen_EnergyMatrix("SC_SC", nBeadTypes_glb, faEnergy_glb, E_SC_SC);

    PrintToScreen_EnergyMatrix("FSOL ", nBeadTypes_glb, faEnergy_glb, E_F_SOL);

    PrintToScreen_EnergyMatrix("STIFF", nBeadTypes_glb, faEnergy_glb, E_STIFF);
}

/// PrintToScreen_MCMoveFreqs
void PrintToScreen_MCMoveFreqs(void)
{
    const char* MoveName[MAX_MV];
    MoveName[MV_PIVOT]      = "Pivot           ";
    MoveName[MV_DBPVT]      = "Double Pivot    ";
    MoveName[MV_CLSTR]      = "Larger Cluster  ";
    MoveName[MV_SMCLSTR]    = "Smaller Cluster ";
    MoveName[MV_STROT]      = "Face Change     ";
    MoveName[MV_LOCAL]      = "Local           ";
    MoveName[MV_COLOCAL]    = "Co-local        ";
    MoveName[MV_MTLOCAL]    = "Shake           ";
    MoveName[MV_BRROT]      = "Rotate Branched ";
    MoveName[MV_SNAKE]      = "Slithering Snake";
    MoveName[MV_TRANS]      = "Translation     ";
    MoveName[MV_PR_SMCLSTR] = "Pr. Smal Cluster";
    MoveName[MV_PR_CLSTR]   = "Pr. Cluster     ";

    float freqMin = 1e10f;
    int i;
    for (i = MV_NULL + 1; i < MAX_MV; i++)
        {
            if (freqMin >= faMCFreq_glb[i] && faMCFreq_glb[i] != 0.)
                {
                    freqMin = faMCFreq_glb[i];
                }
        }

    for (i = 0; i < 30; i++)
        {
            printf("-");
        }
    printf("\n");
    printf("| MC Move Frequencies%9s\n", "|");
    for (i = MV_NULL + 1; i < MAX_MV; i++)
        {
            printf("| %-16s %5d %5s\n", MoveName[i], (int) ceilf(faMCFreq_glb[i] / freqMin), "|");
        }
    for (i = 0; i < 30; i++)
        {
            printf("-");
        }
}

void ScreenIOUtil_PrintTemperatures(void)
{
    printf("Thermalizing Temperature       = %.2f\n", fPreKT_glb);

    if (nTemp_inv_glb == 1)
        {
            puts("Temperature Inverted           = YES");
            printf("MC Temperatures: (First, Last) = (%.2f, %.2f)\n", fKT_glb,
                   1.f / (fKT_glb + (float) (nTotCycleNum_glb - 1) * fDeltaTemp_glb));
        }
    else
        {
            puts("Temperature Inverted           = NO");
            printf("MC Temperatures: (First, Last) = (%.2f, %.2f)\n", fKT_glb,
                   fKT_glb + (float) (nTotCycleNum_glb - 1) * fDeltaTemp_glb);
        }

    printf("Temperature Annealing Mode     = %d\n", nAnnealing_Mode_glb);
}

/// ScreenIO_Print_KeyFile - print the keyfile that was read in to the screen
void ScreenIO_Print_KeyFile(void)
{ // should be output-dependent (stdout, stderr, other files)

    int i;
    const char lBrace[] = "<======      ";
    const char rBrace[] = "      ======>";
    printf("%s System Settings %s\n", lBrace, rBrace);
    printf("Number of Bead Types = %d\n", nBeadTypes_glb);
    printf("Number of Beads      = %ld\n", (size_t) tot_beads_glb);
    printf("Number of Chains     = %ld\n", (size_t) tot_chains_glb);
    printf("Number of Components = %ld\n", (size_t) tot_chain_types_glb);
    printf("Lattice Dimensions   = %3d %3d %3d\n", naBoxSize_glb[0], naBoxSize_glb[1], naBoxSize_glb[2]);
    printf("Monomer Density      = %1.2e\n",
           (float) tot_beads_glb / (float) naBoxSize_glb[0] / (float) naBoxSize_glb[1] / (float) naBoxSize_glb[2]);
    printf("\n");

    PrintToScreen_AllEnergyMatrices();

    printf("\n");

    printf("%s Linker Info %s\n", lBrace, rBrace);
    printf("Linker length             = %.2f\n", fLinkerLength_glb);
    printf("Linker spring constant    = %.2f\n", fLinkerSprCon_glb);
    printf("Linker equilibrium length = %.2f\n", fLinkerEqLen_glb);
    printf("\n");

    printf("%s MC Setup %s\n", lBrace, rBrace);
    ScreenIOUtil_PrintTemperatures();

    printf("Indent Mode                    = %d\n", nBiasPotential_Mode_glb);
    //    printf("Rotational Bias Mode           = %d\n", RotBias_Mode_glb);
    printf("Number of MC Cycles            = %d\n", nTotCycleNum_glb);
    printf("Number of MC Steps/Cycle       = %e\n", (float) nMCStepsPerCycle_glb);
    printf("Number of Thermalizing Steps   = %e\n", (float) nMCStepsForTherm_glb);
    printf("RNG Seed                       = %d\n", nRNG_Seed_glb);
    printf("Clustering Mode                = %d\n", nClusteringMode_glb);

    PrintToScreen_MCMoveFreqs();
    printf("\n");
}

/// ScreenIO_Print_SystemEnergy
void ScreenIO_Print_SystemEnergy(void)
{

    int i;

    char sSectionHead[32];
    memset(sSectionHead, '-', 19);
    sSectionHead[19] = '\0';

    printf("%s\n", sSectionHead);
    puts("Energies          |");
    printf("Tot  : %-10.2e |\n", faCurrEn_glb[E_TOT]);
    printf("Ovlp : %-10.2e |\n", faCurrEn_glb[E_OVLP]);
    printf("Cont : %-10.2e |\n", faCurrEn_glb[E_CONT]);
    printf("Aniso: %-10.2e |\n", faCurrEn_glb[E_SC_SC]);
    printf("FSol : %-10.2e |\n", faCurrEn_glb[E_F_SOL]);
    printf("Stiff: %-10.2e |\n", faCurrEn_glb[E_STIFF]);
    printf("Bias : %-10.2e |\n", faCurrEn_glb[E_BIAS]);
    printf("%s\n", sSectionHead);
}

/// void ScreenIO_Print_AcceptanceRatios
void ScreenIO_Print_AcceptanceRatios(void)
{

    const char* MoveName[MAX_MV];
    MoveName[MV_PIVOT]      = "Pivot       ";
    MoveName[MV_DBPVT]      = "Double Pivot";
    MoveName[MV_CLSTR]      = "La Cluster  ";
    MoveName[MV_SMCLSTR]    = "Sm Cluster  ";
    MoveName[MV_STROT]      = "Face Change ";
    MoveName[MV_LOCAL]      = "Local       ";
    MoveName[MV_COLOCAL]    = "Co-local    ";
    MoveName[MV_MTLOCAL]    = "Multi Local ";
    MoveName[MV_BRROT]      = "Rot. Br.    ";
    MoveName[MV_SNAKE]      = "Sli. Snake  ";
    MoveName[MV_TRANS]      = "Translation ";
    MoveName[MV_PR_SMCLSTR] = "Pr. Sm. Cls.";
    MoveName[MV_PR_CLSTR]   = "Pr. Cls.    ";

    int i, j;
    float fAccRatio = 0.f;
    lLong nMoveSum  = 0L;

    char sSectionHead[32];
    memset(sSectionHead, '-', 22);
    sSectionHead[22] = '\0';

    printf("%s\n", sSectionHead);
    puts("Acceptance Ratios    |");
    for (i = 1; i < MAX_MV; i++)
        {
            nMoveSum = naMCAccepMat_glb[0][i] + naMCAccepMat_glb[1][i];
            if (nMoveSum)
                {
                    fAccRatio = 100.f * (float) naMCAccepMat_glb[1][i] / (float) nMoveSum;
                    printf("%-12s: %-7.2f|\n", MoveName[i], fAccRatio);
                }
            else
                {
                    printf("%-12s: %-7s|\n", MoveName[i], "NA");
                }
        }
    printf("%s\n", sSectionHead);
}

/// ScreenIO_Print_Log_Thermalization - print the log to the screen.
void ScreenIO_Print_Log_Thermalization(const long nGen)
{
    puts("****************************************");
    printf("Run Cycle: Thermalization\n");
    printf("Step     : %8.3e\n", (float) nGen);
    printf("MC Temp  : %8.3e\n", fCuTemp_glb);

    if (nTotClusCounter_glb > 0)
        {
            printf("Perc Phi : %8.3f\n",
                   ((float) nLargestClusterRightNow_glb) / ((float) nTotClusCounter_glb) / (float) tot_chains_glb);
        }

    ScreenIO_Print_SystemEnergy();
    ScreenIO_Print_AcceptanceRatios();
    puts("****************************************");
}

/// ScreenIO_Print_Log_FullRun - print the log to the screen.
void ScreenIO_Print_Log_FullRun(const long nGen, const int run_cycle)
{
    puts("****************************************");
    printf("Run Cycle: %d\n", run_cycle);
    printf("Step     : %8.3e\n", (float) nGen);
    printf("MC Temp  : %8.3e\n", fCuTemp_glb);

    if (nTotClusCounter_glb > 0)
        {
            printf("Perc Phi : %8.3f\n",
                   ((float) nLargestClusterRightNow_glb) / ((float) nTotClusCounter_glb) / (float) tot_chains_glb);
        }

    ScreenIO_Print_SystemEnergy();
    ScreenIO_Print_AcceptanceRatios();
    puts("****************************************");
}

/// ScreenIO_PrintStdErr_BeadPosition. Prints the beadID and it's position to STDERR.
/// \param beadID
void ScreenIO_PrintStdErr_BeadPosition(const int beadID)
{
    const int x = bead_info_glb[beadID][POS_X];
    const int y = bead_info_glb[beadID][POS_Y];
    const int z = bead_info_glb[beadID][POS_Z];
    fprintf(stderr, "bead: %d is a position: <%d, %d, %d>\n", beadID, x, y, z);
}

/// ScreenIO_Print_SanityCheckFailurePreamble. Prints the preamble to the error message for sanity check failures.
/// \param nGen Which MC Step this crash has occurred on.
/// \param run_cycle Which run_cycle this crash has occurred on.
void ScreenIO_Print_SanityCheckFailurePreamble(const long nGen, const int run_cycle)
{
    fputs("\n-------------------------------------------------------------------------------", stderr);
    fputs("\nERROR! Sanity check failed! Crashing!\n", stderr);
    fprintf(stderr, "Crash occurred at (run_cycle: %d; mc_step: %ld)\n", run_cycle, nGen);
    fputs("Snapshot of system saved to crash_snapshot.txt\n", stderr);
    fputs("CRASH TYPE: ", stderr);
}

/// ScreenIO_Print_SanityFail_BeadPosAndLattPos. Prints the crash information for when the bead position does not
/// correspond to the lattice representation. Prints the bead position, and what the lattice thinks is at that position.
/// \param badBead The bead for which to print the information.
void ScreenIO_Print_SanityFail_BeadPosAndLattPos(const int badBead)
{

    const int x = bead_info_glb[badBead][POS_X];
    const int y = bead_info_glb[badBead][POS_Y];
    const int z = bead_info_glb[badBead][POS_Z];
    fprintf(stderr, "bead: %d should be at: <%d, %d, %d>\n", badBead, x, y, z);

    const int latInd = Lat_Ind_OfBead(badBead);
    const int latVal = naTotLattice_glb[latInd];

    fprintf(stderr, "while\nlattice: <%d, %d, %d> (ind: %d) has: %d\n", x, y, z, latInd, latVal);
}

/// ScreenIO_Print_SanityFail_MolecularStructure. Prints the crash information for when the linkers for badBead are
/// broken. Prints out the chain badBead is in, all the bonded partners for badBead, where the beads are at, and the
/// distances between the beads.
/// \param badBead The bead for which to print the information.
void ScreenIO_Print_SanityFail_MolecularStructure(const int badBead)
{

    const int badChain = bead_info_glb[badBead][BEAD_CHAINID];

    fprintf(stderr, "chain: %d has broken structure\n", badChain);
    fprintf(stderr, "bead:  %d has broken linkers\n", badBead);
    ScreenIO_PrintStdErr_BeadPosition(badBead);
    fprintf(stderr, "\nbead:  %d has partners:\n", badBead);
    int idx      = 0;
    int bondPart = topo_info_glb[badBead][idx];
    while (bondPart != -1)
        {
            ScreenIO_PrintStdErr_BeadPosition(bondPart);
            idx++;
            bondPart = topo_info_glb[badBead][idx];
        }
    fputs("\n\nThe beads have the following distances and linker constraints:\n", stderr);

    idx             = 0;
    bondPart        = topo_info_glb[badBead][idx];
    float bdDist   = 0.;
    float linkCons = 0.;

    while (bondPart != -1)
        {

            bdDist   = Dist_BeadToBead(badBead, bondPart);
            linkCons = (float) linker_len_glb[badBead][idx] * LINKER_RSCALE;

            fprintf(stderr, "beads: (%d, %d) have dist: %5.3f with linker-constraint: %5.3f\n", badBead, bondPart,
                    bdDist, linkCons);
            idx++;
            bondPart = topo_info_glb[badBead][idx];
        }
}

/// ScreenIO_Print_SanityFail_SelfBond. Prints the crash information for when the a bead is self-bonded. Prints out the
/// id of the bead.
/// \param badBead The bead for which to print the information.
void ScreenIO_Print_SanityFail_SelfBond(const int badBead)
{
    fprintf(stderr, "bead: %d is bonded to itself!\n", badBead);
}

/// ScreenIO_Print_SanityFail_BeadBondSymmetry. Prints the crash information for when a bond is not symmetric. Prints
/// out the id of the bead that failed, which bead badBead is bonded to, and which bead the beadPartner is bonded to.
/// We also print the positions of the beads.
/// \param badBead The bead for which to print the information.
void ScreenIO_Print_SanityFail_BeadBondSymmetry(const int badBead)
{
    const int beadPart = bead_info_glb[badBead][BEAD_FACE];

    fprintf(stderr, "bead: %d is bonded to bead: %d\n", badBead, beadPart);

    const int partBond = bead_info_glb[beadPart][BEAD_FACE];

    fprintf(stderr, "while\nbead: %d is bonded to bead: %d\n", beadPart, partBond);

    ScreenIO_PrintStdErr_BeadPosition(badBead);
    ScreenIO_PrintStdErr_BeadPosition(beadPart);
}

/// ScreenIO_Print_SanityFail_BeadBondDistance. Prints the crash information for when the distance between bonded beads
/// is too large. Prints the relevant beads' ids, positions and the distance between them.
/// \param badBead The bead for which to print the information.
void ScreenIO_Print_SanityFail_BeadBondDistance(const int badBead)
{
    const int beadPart = bead_info_glb[badBead][BEAD_FACE];

    fprintf(stderr, "bead: %d is bonded to bead: %d\n", badBead, beadPart);

    const float bDist = Dist_BeadToBead(badBead, beadPart);
    fprintf(stderr, "while\ndistance: %5.3f\n", bDist);

    ScreenIO_PrintStdErr_BeadPosition(badBead);
    ScreenIO_PrintStdErr_BeadPosition(beadPart);
}

/// Write_RDF_ComponentWise - old implementation of printing the RDF, component
/// by component. Always appends to the file for this run. Stopped using it
/// because the IO load was slowing things down at our cluster Would be a good
/// way to gather proper statistics on clusters as the runs went on.
/// TODO: update this to work with the new indexing
/// \param filename
/// \param nGen
void Write_RDF_ComponentWise(char* filename, long nGen)
{
    FILE* fp;
    if (nGen == -1)
        {
            fp = fopen(filename, "w"); // overwrite
        }
    else
        {
            fp = fopen(filename, "a");
        }
    int i, k;
    if (nGen == -1)
        { // title
            fprintf(fp, "#Row-by-row split RDF. dr = 1/4\n");
        }
    else
        {
            for (k = 0; k < nRDF_TotComps_glb; k++)
                {
                    for (i = 0; i < nRDF_TotBins_glb; i++)
                        {
                            fprintf(fp, "%.5Lf\t", ldaRDF_Arr_glb[RDFArr_Index(0, k, i)]);
                        }
                    fprintf(fp, "\n");
                }
        }

    fclose(fp);
}

/// FileIO_WriteTo_TopFile - write a LAMMPS DATA file, which can be read in VMD to store
/// topological information for the system Only writes the topology once at the
/// beginning of the runs.
/// \param filename
void FileIO_WriteTo_TopFile(const char* filename)
{
    /*
    Writes a topology file for VMD. Since the trajectory is saved in the LAMMPS
    format, the topology file is also in the LAMMPS format. The format for this
    file is one which is used as input data for LAMMPS. This function copies the
    format given in read_data for LAMMPS:
    https://lammps.sandia.gov/doc/read_data.html
    */
    FILE* fp;
    fp = fopen(filename, "w"); // Just overwrite over if last file exists.
    int i, j, k;               // Loop iterators.
    int numBonds;              // Used to count total number of bonds!
    printf("Writing the topology file!\n");
    fprintf(fp, "LAMMPS Description\n");                   // The file must start with this.
    fprintf(fp, "\n");                                     // Empty line.
    fprintf(fp, "\t%ld\tatoms\n", (size_t) tot_beads_glb); // Listing total number of atoms
    numBonds = 0;
    for (i = 0; i < tot_beads_glb; i++)
        {
            for (j = 0; j < MAX_BONDS; j++)
                {
                    if (topo_info_glb[i][j] != -1)
                        { // There is a bond between i and topo_info_glb[i][j]
                            numBonds++;
                        }
                }
        }

    fprintf(fp, "\t%d\tbonds\n", numBonds);          // Listing number of bonds
    fprintf(fp, "\t0\tangles\n");                    // Systems don't have angle depenece yet
    fprintf(fp, "\t0\tdihedrals\n");                 // Systems don't have dihedrals  yet
    fprintf(fp, "\t0\timpropers\n");                 // Systems don't have imporopers yet
    fprintf(fp, "\n");                               // Empty line.
    fprintf(fp, "%d\tatom types\n", nBeadTypes_glb); // This many bead-types
    fprintf(fp, "1\tbond types\n");                  // System can have multiple bond-lengths
    fprintf(fp, "0\tangle types\n");                 // Systems don't have any angular forces yet
    fprintf(fp, "\n");                               // Empty line.
    fprintf(fp, " 0 %d xlo xhi\n", naBoxSize_glb[0]);
    fprintf(fp, " 0 %d ylo yhi\n", naBoxSize_glb[1]);
    fprintf(fp, " 0 %d zlo zhi\n", naBoxSize_glb[2]);
    fprintf(fp, "\n");       // Empty line.
    fprintf(fp, "Masses\n"); // These don't really mean anything
    fprintf(fp, "\n");       // Empty line.
    for (i = 0; i < nBeadTypes_glb; i++)
        {
            fprintf(fp, "%d\t1.0\n", i);
        }

    fprintf(fp, "\n");             // Empty line.
    fprintf(fp, "Atoms # bond\n"); // Signifying the beginning of atom coordinates.
    fprintf(fp, "\n");             // Empty line.
    for (i = 0; i < tot_beads_glb; i++)
        {
            fprintf(fp, "%d %d %d %d %d %d\n", i, bead_info_glb[i][BEAD_CHAINID], bead_info_glb[i][BEAD_TYPE],
                    bead_info_glb[i][POS_X], bead_info_glb[i][POS_Y], bead_info_glb[i][POS_Z]);
        }                   // Done with the coordinates
    fprintf(fp, "\n");      // Empty line.
    fprintf(fp, "Bonds\n"); // Signifying the beginning of atom coordinates.
    fprintf(fp, "\n");      // Empty line.
    k = 0;                  // This guy counts bondIDs
    for (i = 0; i < tot_beads_glb; i++)
        {
            for (j = 0; j < MAX_BONDS; j++)
                {
                    if (topo_info_glb[i][j] != -1)
                        { // There is a bond between i and topo_info_glb[i][j]
                            fprintf(fp, "%d %d %d %d\n", k, linker_len_glb[i][j] / (int) fLinkerLength_glb, i,
                                    topo_info_glb[i][j]);
                            k++;
                        }
                }
        }              // Done with the coordinates
    fprintf(fp, "\n"); // Empty line.

    fclose(fp);
}

/// Write_SysProp - writes the RDF, CLUS and GyrTen avg rrays to a file. OLD implementation that should not be used yet
/// Also, TODO: update this to work with the new indexing
/// \param filename
void Write_SysProp(char* filename)
{
    FILE* fp;
    int i, j;
    fp = fopen(filename, "w");
    fprintf(fp, "#Contains various averaged quatities.\n");
    // Gyration Radius
    fprintf(fp, "#Total Rg\n%f\t%f\n#Cluster Hist\n", faSysGyrRad_glb / (float) nTotGyrRadCounter_glb,
            (float) naBoxSize_glb[0] / 2.);
    // Cluster Histogram
    fprintf(fp, "%f\t", (float) nLargestClusterRightNow_glb / (float) nTotClusCounter_glb);
    for (i = 1; i <= tot_chains_glb; i++)
        {
            fprintf(fp, "%f\t", (float) laClusHistList_glb[i] / (float) nTotClusCounter_glb);
        }
    // Split RDFs
    fprintf(fp, "\n#Split RDFs. ALL-ALL; DIAGONALS and then from 0 onwards \n");
    for (j = 0; j < nRDF_TotComps_glb; j++)
        {
            for (i = 0; i < nRDF_TotBins_glb; i++)
                {
                    // fprintf(fp, "%LE\t", ldRDF_ARR[j][i] / (float)nTotRDFCounter_glb);
                    fprintf(fp, "%LE\t", ldaRDF_Arr_glb[RDFArr_Index(0, j, i)] / (float) nTotRDFCounter_glb);
                }
            fprintf(fp, "\n");
        }
    fprintf(fp, "\n#Done");
}

/// FileIO_WriteTo_RDFFile_TotFromGLB - writes the RDFTot file.
/// \param run_it Number of temperature cycles.
void FileIO_WriteTo_RDFFile_TotFromGLB(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_RDF.dat",
            strReportPrefix_glb); // Name Of the RDF Files
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Split RDFs. ALL-ALL; DIAGONALS and then from 0-0, 0-1, and onwards \n");
    for (i = 0; i < run_it; i++)
        {
            fprintf(fp, "#Run_Cycle = %d\n", i);
            for (j = 0; j < nRDF_TotComps_glb; j++)
                {
                    for (k = 0; k < nRDF_TotBins_glb; k++)
                        {
                            fprintf(fp, "%LE\t", ldaTOTRDF_Arr_glb[RDFArr_Index(i, j, k)]);
                        }
                    fprintf(fp, "\n");
                }
        }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_COMDenFile_TotFromGLB - writes the COMDenTot file.
/// \param run_it Number of temperature cycles.
void FileIO_WriteTo_COMDenFile_TotFromGLB(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_COMDen.dat",
            strReportPrefix_glb); // Name Of the COM Density Distribution
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Density distribution from the COM outwards. Order is ChainType.\n");
    for (i = 0; i < run_it; i++)
        {
            fprintf(fp, "#Run_Cycle = %d\n", i);
            for (j = 0; j < nRadDen_TotComps_glb; j++)
                {
                    for (k = 0; k < nRDF_TotBins_glb; k++)
                        {
                            fprintf(fp, "%LE\t", ldaTOTRadDen_Arr_glb[RadDenArr_Index(i, j, k)]);
                        }
                    fprintf(fp, "\n");
                }
        }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_ClusFile_TotFromGLB - writes the CLUS file.
/// \param run_it Number of temperature cycles.
void FileIO_WriteTo_ClusFile_TotFromGLB(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_CLUS.dat",
            strReportPrefix_glb); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Cluster Histograms: 1st column is largest cluster, and then clusters of size 1, 2, and so on\n");
    for (i = 0; i < run_it; i++)
        {
            fprintf(fp, "#Run_Cycle = %d\n", i);
            for (k = 0; k <= tot_chains_glb; k++)
                {
                    fprintf(fp, "%LE\t", ldaTOTCLUS_Arr_glb[i][k]);
                }
            fprintf(fp, "\n");
        }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_MolClus - writes the MolClus file
/// \param run_it Number of temperature cycles.
void FileIO_WriteTo_MolClusFile_TotFromGLB(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_MolClus.dat",
            strReportPrefix_glb); // Name Of the COM Density Distribution
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Cluster Histograms: 1st column is largest cluster, and then clusters of size 1, 2, and so on. Each "
                "row is a different moltype\n");
    for (i = 0; i < run_it; i++)
        {
            fprintf(fp, "#Run_Cycle = %d\n", i);
            for (j = 0; j < tot_chain_types_glb; j++)
                {
                    for (k = 0; k < tot_chains_glb; k++)
                        {
                            fprintf(fp, "%LE\t", ldaTOTMOLCLUS_Arr_glb[MolClusArr_Index(i, j, k)]);
                        }
                    fprintf(fp, "\n");
                }
        }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_GyrRad - writes the GR file.
/// \param run_it Number of temperature cycles.
void FileIO_WriteTo_GyrRadFile_TotFromGLB(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_GR.dat",
            strReportPrefix_glb); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Cluster Histograms: 1st column is largest cluster, and then "
                "clusters of size 1, 2, and so on\n");
    for (i = 0; i < run_it; i++)
        {
            fprintf(fp, "#Run_Cycle = %d\n", i);
            fprintf(fp, "%LE\t", ldaTOTRg_Arr_glb[i][0]);
            fprintf(fp, "%LE\n", ldaTOTRg_Arr_glb[i][1]);
        }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_Write_TotalSysProp_TotFromGLB - writes the RDF, CLUS and GyrTen arrays to their
/// respective files. These arrays store the data over the course of an ENTIRE
/// run, or over all the temp cycles.
/// \param run_it
void FileIO_Write_TotalSysProp_TotFromGLB(const int run_it)
{
    /* This function writes one large file with all the averaged values from
     * each run_cycle. run_it let's the function know how many cycles to write.
     * As opposed to Write_SysProp, each of the relevant averaged parameters
     * will be in their appropriately named files. The naming convention shall
     * be: $REPORTNAME_*.dat where * will be GR, CLUS and RDF for the different
     * measured quantities. Within each file, each run_cycle shall be written
     * one-by-one in the style of Write_SysProp
     */
    if (naReportFreqs_glb[REPORT_RDFTOT] != 0)
        {
            FileIO_WriteTo_RDFFile_TotFromGLB(run_it);
        }

    if (naReportFreqs_glb[REPORT_COMDEN] != 0)
        {
            FileIO_WriteTo_COMDenFile_TotFromGLB(run_it);
        }

    if (naReportFreqs_glb[REPORT_NETWORK] != 0)
        {
            FileIO_WriteTo_ClusFile_TotFromGLB(run_it);
            FileIO_WriteTo_MolClusFile_TotFromGLB(run_it);
            FileIO_WriteTo_GyrRadFile_TotFromGLB(run_it);
        }
}

/// FileIOUtil_CreateFile_Overwrite - Creates a new overwritten file with the given name.
/// \param fileName
void FileIOUtil_CreateFile_Overwrite(const char* fileName)
{
    FILE* fp = fopen(fileName, "w+");
    fclose(fp);
}

/// FileIOUtil_CreateFile_Overwrite - Creates a new overwritten file with the given name.
/// \param fileName
void FileIOUtil_CreateFile_Binary_Overwrite(const char* fileName)
{
    FILE* fp = fopen(fileName, "wb");
    fclose(fp);
}

/// FileIO_CreateRunningDataFiles - Creates the necessary files that are continuously written to over the course
/// of a simulation.
void FileIO_CreateRunningDataFiles(void)
{
    // Trajectory
    if (naReportFreqs_glb[REPORT_CONFIG])
        {
            sprintf(strFileTraj_glb, "%s_topo.lammpstrj", strReportPrefix_glb); // Name of the topology file
            FileIO_WriteTo_TopFile(strFileTraj_glb); // Write the topology file. Only need to write once

            if (nTrajMode_glb == 1)
                {
                    sprintf(strFileTraj_glb, "%s_EQ_trj.lassi", strReportPrefix_glb);
                    FileIOUtil_CreateFile_Binary_Overwrite(strFileTraj_glb); // Opens a new binary file.
                }
            else
                {
                    sprintf(strFileTraj_glb, "%s_EQ_trj.lammpstrj",
                            strReportPrefix_glb); // Naming convention for trajectory files.
                    FileIOUtil_CreateFile_Overwrite(
                        strFileTraj_glb); // This opens a new trajectory file; each run_it will have its own
                }
        }

    // Energy
    if (naReportFreqs_glb[REPORT_ENERGY])
        {
            sprintf(strFileEnergy_glb, "%s_energy.dat", strReportPrefix_glb);
            FileIOUtil_CreateFile_Overwrite(strFileEnergy_glb); // Open a new energy file; each run_it will have its own
            FileIOUtil_WriteHeader_ForEnergy(strFileEnergy_glb);
        }

    // MC Move Acceptance
    if (naReportFreqs_glb[REPORT_MCMOVE])
        {
            sprintf(strFileMCMove_glb, "%s_mcmove.dat", strReportPrefix_glb);
            FileIOUtil_CreateFile_Overwrite(strFileMCMove_glb); // Open a new MCInfo file; each run_it will have its own
            FileIOUtil_WriteHeader_ForMCMove(strFileMCMove_glb);
        }

    char strFileName[512];

    if (naReportFreqs_glb[REPORT_COMDEN])
        {
            sprintf(strFileName, "%s_COMDen.dat", strReportPrefix_glb);
            FileIOUtil_CreateFile_Overwrite(strFileName);
            FileIOUtil_WriteHeader_ForCOMDen(strFileName);
        }

    if (naReportFreqs_glb[REPORT_NETWORK])
        {
            sprintf(strFileName, "%s_CLUS.dat", strReportPrefix_glb);
            FileIOUtil_CreateFile_Overwrite(strFileName);
            FileIOUtil_WriteHeader_ForCLUS(strFileName);

            sprintf(strFileName, "%s_MolClus.dat", strReportPrefix_glb);
            FileIOUtil_CreateFile_Overwrite(strFileName);
            FileIOUtil_WriteHeader_ForMolClus(strFileName);

            sprintf(strFileName, "%s_GR.dat", strReportPrefix_glb);
            FileIOUtil_CreateFile_Overwrite(strFileName);
            FileIOUtil_WriteHeader_ForGyrRad(strFileName);
        }

    if (naReportFreqs_glb[REPORT_RDFTOT])
        {
            sprintf(strFileName, "%s_RDF.dat", strReportPrefix_glb);
            FileIOUtil_CreateFile_Overwrite(strFileName);
            FileIOUtil_WriteHeader_ForRDF(strFileName);
        }
}

/// ForPrinting_GetReportState. Calculates if the current step is a multiple of the reporting frequency
/// for this particular report.
/// Convenient wrapper around some modulo arithmetic.
/// \param nGen
/// \param thisReport
/// \return
char ForPrinting_GetReportState(const long nGen, const long thisReport)
{
    char dum_log = (char) ((nGen % thisReport) == 0);
    dum_log      = dum_log ? 1 : 0;
    return dum_log;
}

void FileIO_PrintCrashSnapshot(void)
{
    FILE* fp = fopen("crash_snapshot.txt", "w+");

    int i, x, y, z, bP;
    fputs("#CRASH SNAPSHOT\n#X, Y, Z, bondPartner\n", fp);
    for (i = 0; i < tot_beads_glb; i++)
        {
            x  = bead_info_glb[i][POS_X];
            y  = bead_info_glb[i][POS_Y];
            z  = bead_info_glb[i][POS_Z];
            bP = bead_info_glb[i][BEAD_FACE];
            fprintf(fp, "%d %d %d %d\n", x, y, z, bP);
        }

    fclose(fp);
}

/// DataPrinting_Thermalization - helper function for printing out data during the thermalization cycle.
/// No data analysis is performed during the thermalization procedure.
/// This function decides if, given the MCStep, it is time to print the following:
/// 1. The Log - to the screen.
/// 2. Trajectory - to the file (or saved if in total format)
/// 3. Energy - to the file.
/// 4. MCMove - to the file.
/// \param nGen
void DataPrinting_Thermalization(const long nGen)
{

    char cFlagForEnCal = 0;
    char cLogFlag      = 0;
    char cEnergyFlag   = 0;
    char cAccFlag      = 0;
    char cConfigFlag   = 0;

    if (naReportFreqs_glb[REPORT_LOG])
        {
            cLogFlag = ForPrinting_GetReportState(nGen, naReportFreqs_glb[REPORT_LOG]);
            if (cLogFlag)
                {
                    PerformRuntimeSanityChecks(nGen, -1);
                    Energy_Total_System();
                    cFlagForEnCal = 1;
                    ScreenIO_Print_Log_Thermalization(nGen);
                }
        }

    if (nGen)
        {
            if (naReportFreqs_glb[REPORT_CONFIG])
                {
                    cConfigFlag = ForPrinting_GetReportState(nGen, naReportFreqs_glb[REPORT_CONFIG]);
                    if (cConfigFlag)
                        {
                            FileIO_Trajectory_AppendFrame(strFileTraj_glb, -1, nGen);
                        }
                }

            if (naReportFreqs_glb[REPORT_ENERGY])
                {
                    // DO ENERGY FLAGS
                    cEnergyFlag = ForPrinting_GetReportState(nGen, naReportFreqs_glb[REPORT_ENERGY]);
                    if (cEnergyFlag)
                        {
                            if (! cFlagForEnCal)
                                {
                                    Energy_Total_System();
                                    cFlagForEnCal = 1;
                                }
                            FileIO_AppendEnergyTo_EnergyFile(strFileEnergy_glb, nGen);
                        }
                }

            if (naReportFreqs_glb[REPORT_MCMOVE])
                {
                    // DO MC_ACC PRINTING
                    cAccFlag = ForPrinting_GetReportState(nGen, naReportFreqs_glb[REPORT_MCMOVE]);
                    if (cAccFlag)
                        {
                            FileIO_WriteTo_MCMoveFile(strFileMCMove_glb, nGen, fCuTemp_glb);
                        }
                }
        }
}

/// DataPrinting_DuringRunCycles
/// \param nGen
/// \param run_it
void DataPrinting_DuringRunCycles(const long nGen, const int run_it)
{
    char cFlagForEnCal = 0;
    char cLogFlag      = 0;
    char cEnergyFlag   = 0;
    char cAccFlag      = 0;
    char cConfigFlag   = 0;

    if (naReportFreqs_glb[REPORT_LOG])
        {
            cLogFlag = ForPrinting_GetReportState(nGen, naReportFreqs_glb[REPORT_LOG]);
            if (cLogFlag)
                {
                    PerformRuntimeSanityChecks(nGen, run_it);
                    Energy_Total_System();
                    cFlagForEnCal = 1;
                    ScreenIO_Print_Log_FullRun(nGen, run_it);
                }
        }

    if (nGen)
        {
            if (naReportFreqs_glb[REPORT_CONFIG])
                {
                    cConfigFlag = ForPrinting_GetReportState(nGen, naReportFreqs_glb[REPORT_CONFIG]);
                    if (cConfigFlag)
                        {
                            FileIO_Trajectory_AppendFrame(strFileTraj_glb, run_it, nGen);
                        }
                }

            if (naReportFreqs_glb[REPORT_ENERGY])
                {
                    // DO ENERGY PRINTING
                    cEnergyFlag = ForPrinting_GetReportState(nGen, naReportFreqs_glb[REPORT_ENERGY]);
                    if (cEnergyFlag)
                        {
                            if (! cFlagForEnCal)
                                {
                                    Energy_Total_System();
                                    cFlagForEnCal = 1;
                                }
                            FileIO_AppendEnergyTo_EnergyFile(strFileEnergy_glb, nGen);
                        }
                }

            if (naReportFreqs_glb[REPORT_MCMOVE])
                {
                    // DO MC_ACC PRINTING
                    cAccFlag = ForPrinting_GetReportState(nGen, naReportFreqs_glb[REPORT_MCMOVE]);
                    if (cAccFlag)
                        {
                            FileIO_WriteTo_MCMoveFile(strFileMCMove_glb, nGen, fCuTemp_glb);
                        }
                }
        }
}

/// DataAnalysis_DuringRunCycles
/// \param nGen
/// \param run_it
void DataAnalysis_DuringRunCycles(const long nGen, const int run_it)
{

    //    char cRDF_flag  = 0;
    //    char cCOM_flag  = 0;
    //    char cCLUS_flag = 0;

    if (naReportFreqs_glb[REPORT_RDFTOT])
        { // SysProp is printed outside of this function in main.c, lol
            if (nGen % naReportFreqs_glb[REPORT_RDFTOT] == 0)
                {
                    // TODO: Update RDF Calculations
                    RDF_ComponentWise_Avg();
                }
        }
    if (naReportFreqs_glb[REPORT_COMDEN])
        { // SysProp is printed outside of this function in main.c, lol
            if (nGen % naReportFreqs_glb[REPORT_COMDEN] == 0)
                {
                    RadialDensityAnalysis_Perform_Analysis();
                }
        }
    if (naReportFreqs_glb[REPORT_NETWORK])
        { // SysProp is printed outside of this function in main.c, lol
            if (nGen % naReportFreqs_glb[REPORT_NETWORK] == 0)
                {
                    ClusAnalysis_Perform_Analysis(nClusteringMode_glb);
                    // TODO: Update Gyration calculations
                    GyrTensor_GyrRad_Avg();
                }
        }
}

/// FileIOUtil_PreCycle_Init
/// \param run_it
void FileIOUtil_PreCycle_Init(const int run_it)
{
    if (naReportFreqs_glb[REPORT_CONFIG])
        {
            if (nTrajMode_glb == 1)
                {
                    sprintf(strFileTraj_glb, "%s_trj.lassi", strReportPrefix_glb);
                    if (run_it == 0)
                        {
                            FileIOUtil_CreateFile_Binary_Overwrite(strFileTraj_glb); // Opens a new binary file.
                        }
                }
            else
                {
                    sprintf(strFileTraj_glb, "%s_trj.lammpstrj", strReportPrefix_glb);
                    if (run_it == 0)
                        {
                            FileIOUtil_CreateFile_Overwrite(strFileTraj_glb);
                        }
                }
        }

    if (naReportFreqs_glb[REPORT_ENERGY])
        {
            sprintf(strFileEnergy_glb, "%s_%d_energy.dat", strReportPrefix_glb, run_it);
            FileIOUtil_CreateFile_Overwrite(strFileEnergy_glb); // Open a new Energy file; each run_it will have its own
            FileIOUtil_WriteHeader_ForEnergy(strFileEnergy_glb);
        }
    if (naReportFreqs_glb[REPORT_MCMOVE])
        {
            sprintf(strFileMCMove_glb, "%s_%d_mcmove.dat", strReportPrefix_glb, run_it);
            FileIOUtil_CreateFile_Overwrite(strFileMCMove_glb); // Open a new MCInfo file; each run_it will have its own
            FileIOUtil_WriteHeader_ForMCMove(strFileMCMove_glb);
        }
}

/// FileIO_PostCycle_WriteSystemRestart
/// \param run_it
void FileIO_PostCycle_WriteSystemRestart(const int run_it)
{
    sprintf(strFileTraj_glb, "%s_%d_restart.lammpstrj", strReportPrefix_glb,
            run_it); // Naming convention for trajectory files.
    FileIOUtil_CreateFile_Overwrite(strFileTraj_glb);
    FileIOUtil_Traj_Txt_AppendFrame_ToFile(strFileTraj_glb, nMCStepsForTherm_glb + (run_it + 1) * nMCStepsPerCycle_glb);
}

/// FileIO_WriteRestart_ForThermalization
void FileIO_WriteRestart_ForThermalization(void)
{
    if (nMCStepsForTherm_glb)
        {
            sprintf(strFileTraj_glb, "%s_EQ_restart.lammpstrj",
                    strReportPrefix_glb); // Naming convention for trajectory files.
            FileIOUtil_CreateFile_Overwrite(strFileTraj_glb);
            FileIOUtil_Traj_Txt_AppendFrame_ToFile(strFileTraj_glb, nMCStepsForTherm_glb);
        }
}

/// CopyData_All - copies data from run_it specific data arrays to the overall
/// global data arrays that are printed later. Also note that the averaging, or
/// dividing by the frequency of acquisitions, occurs here.
/// \param run_it
void CopyData_All(const int run_it)
{
    if (naReportFreqs_glb[REPORT_RDFTOT])
        {
            CopyData_RDF(run_it);
        }
    if (naReportFreqs_glb[REPORT_COMDEN])
        {
            CopyData_COMDen(run_it);
        }
    if (naReportFreqs_glb[REPORT_NETWORK])
        {
            CopyData_Clus(run_it);
        }
}

/// CopyData_RDF: Copies ldaRDF_Arr_glb into ldaTOTRDF_Arr_glb
/// \param run_it: which run cycle we are on. Should be >= 0.
void CopyData_RDF(const int run_it)
{
    int i, j;
    const long double normFactor = 1.0 / (long double) nTotRDFCounter_glb;
    for (i = 0; i < nRDF_TotComps_glb; i++)
        {
            for (j = 0; j < nRDF_TotBins_glb; j++)
                {
                    ldaTOTRDF_Arr_glb[RDFArr_Index(run_it, i, j)] = ldaRDF_Arr_glb[RDFArr_Index(0, i, j)] * normFactor;
                }
        }
}

/// CopyData_RDF: Copies ldaRadDen_Arr_glb into ldaTOTRadDen_Arr_glb
/// \param run_it: which run cycle we are on. Should be >= 0.
void CopyData_COMDen(const int run_it)
{
    int i, j;
    const long double normFactor = 1.0 / (long double) nTotRadDenCounter_glb;
    for (i = 0; i < nRadDen_TotComps_glb; i++)
        {
            for (j = 0; j < nRDF_TotBins_glb; j++)
                {
                    ldaTOTRadDen_Arr_glb[RadDenArr_Index(run_it, i, j)] =
                        ldaRadDen_Arr_glb[RadDenArr_Index(0, i, j)] * normFactor;
                }
        }
}

/// CopyData_Clus - Copies ldMOLClUS_Arr into ldaTOTMOLCLUS_Arr_glb. Also copies laClusHistList_glb into
/// ldaTOTCLUS_Arr_glb. And stores the Gyration radius.
/// \param run_it: which run cycle we are on. Should be >= 0.
void CopyData_Clus(const int run_it)
{
    int i, j;
    const long double normFactor  = 1.0 / (long double) nTotClusCounter_glb;
    ldaTOTCLUS_Arr_glb[run_it][0] = (long double) nLargestClusterRightNow_glb * normFactor;
    for (i = 1; i <= tot_chains_glb; i++)
        {
            ldaTOTCLUS_Arr_glb[run_it][i] = (long double) laClusHistList_glb[i] * normFactor;
        }
    for (i = 0; i < tot_chains_glb; i++)
        {
            for (j = 0; j < tot_chain_types_glb; j++)
                {
                    ldaTOTMOLCLUS_Arr_glb[MolClusArr_Index(run_it, j, i)] =
                        ldaMOLCLUS_Arr_glb[MolClusArr_Index(0, j, i)] * normFactor;
                }
        }
    ldaTOTRg_Arr_glb[run_it][0] = (long double) faSysGyrRad_glb / (long double) nTotGyrRadCounter_glb;
    ldaTOTRg_Arr_glb[run_it][1] = (long double) naBoxSize_glb[0] / 2.0;
}

///
/// \param run_it
void FileIOUtil_AppendTo_GyrRadFile_ForRun(const int run_it)
{
    FILE* fp;
    char dumFile[512];
    sprintf(dumFile, "%s_GR.dat", strReportPrefix_glb); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "a");

    fprintf(fp, "#Run_Cycle = %d; Avg Over %d Samples\n", run_it, nTotClusCounter_glb);
    fprintf(fp, "%LE\t", (long double) faSysGyrRad_glb / (long double) nTotGyrRadCounter_glb);
    fprintf(fp, "%LE\n", (long double) naBoxSize_glb[0] / 2.0);

    fclose(fp);
}

///
/// \param run_it
void FileIOUtil_AppendTo_ClusFile_ForRun(const int run_it)
{
    FILE* fp;
    int i;
    char dumFile[512];
    sprintf(dumFile, "%s_CLUS.dat", strReportPrefix_glb); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "a");
    fprintf(fp, "#Run_Cycle = %d; Avg Over %d Samples\n", run_it, nTotClusCounter_glb);
    const double normFactor = 1.0 / (double) nTotClusCounter_glb;

    laClusHistList_glb[0] = nLargestClusterRightNow_glb;

    for (i = 0; i <= tot_chains_glb; i++)
        {
            fprintf(fp, "%E\t", (double) laClusHistList_glb[i] * normFactor);
        }
    fprintf(fp, "\n");
    fclose(fp);
}

///
/// \param run_it
void FileIOUtil_AppendTo_MolClusFile_ForRun(const int run_it)
{
    FILE* fp;
    int i, j;
    char dumFile[512];
    sprintf(dumFile, "%s_MolClus.dat", strReportPrefix_glb); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "a");
    fprintf(fp, "#Run_Cycle = %d; Avg Over %d Samples\n", run_it, nTotClusCounter_glb);

    const long double normFactor = 1.0 / (long double) nTotClusCounter_glb;

    for (j = 0; j < tot_chain_types_glb; j++)
        {
            for (i = 0; i < tot_chains_glb; i++)
                {
                    fprintf(fp, "%LE\t", ldaMOLCLUS_Arr_glb[MolClusArr_Index(0, j, i)] * normFactor);
                }
            fprintf(fp, "\n");
        }
    fclose(fp);
}

///
/// \param run_it
void FileIOUtil_AppendTo_RDFFile_ForRun(const int run_it)
{
    FILE* fp;
    int i, j;
    char dumFile[512];
    sprintf(dumFile, "%s_RDF.dat", strReportPrefix_glb); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "a");
    fprintf(fp, "#Run_Cycle = %d; Avg Over %d Samples\n", run_it, nTotRDFCounter_glb);

    const long double normFactor = 1.0 / (long double) nTotRDFCounter_glb;

    for (j = 0; j < nRDF_TotComps_glb; j++)
        {
            for (i = 0; i < nRDF_TotBins_glb; i++)
                {
                    fprintf(fp, "%LE\t", ldaRDF_Arr_glb[RDFArr_Index(0, j, i)] * normFactor);
                }
            fprintf(fp, "\n");
        }
    fclose(fp);
}

///
/// \param run_it
void FileIOUtil_AppendTo_COMDenFile_ForRun(const int run_it)
{
    FILE* fp;
    int i, j;
    char dumFile[512];
    sprintf(dumFile, "%s_COMDen.dat", strReportPrefix_glb);
    fp = fopen(dumFile, "a");
    fprintf(fp, "#Run_Cycle = %d; Avg Over %d Samples\n", run_it, nTotRadDenCounter_glb);

    const long double normFactor = 1.0 / (long double) nTotRadDenCounter_glb;

    for (j = 0; j < nRadDen_TotComps_glb; j++)
        {
            for (i = 0; i < nRDF_TotBins_glb; i++)
                {
                    fprintf(fp, "%LE\t", ldaRadDen_Arr_glb[RadDenArr_Index(0, j, i)] * normFactor);
                }
            fprintf(fp, "\n");
        }
    fclose(fp);
}

///
/// \param run_it
void FileIO_PostCycle_WriteCycleAvgData(const int run_it)
{
    if (naReportFreqs_glb[REPORT_RDFTOT])
        {
            FileIOUtil_AppendTo_RDFFile_ForRun(run_it);
        }
    if (naReportFreqs_glb[REPORT_COMDEN])
        {
            FileIOUtil_AppendTo_COMDenFile_ForRun(run_it);
        }
    if (naReportFreqs_glb[REPORT_NETWORK])
        {
            FileIOUtil_AppendTo_ClusFile_ForRun(run_it);
            FileIOUtil_AppendTo_MolClusFile_ForRun(run_it);
            FileIOUtil_AppendTo_GyrRadFile_ForRun(run_it);
        }
}
