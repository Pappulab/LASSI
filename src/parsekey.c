#include "parsekey.h"
#include "cluster.h"
#include "global.h"
#include "initialize.h"

int str2farr(char* strRaw, float fArray[MAX_AA]);

/// Parse_KeyFile - reads the parameter keyfile. I recommend using the Python
/// scripts to generate key files that are used in actual runs. Have to
/// painstakingly go through every different keyword that exists in the keyfile.
/// Everytime you want to add something to the keyfile, do it here and also the
/// python scripts so that everything is consistent.
/// \param filename
/// \return
int Parse_Keyfile(char* filename)
{
    FILE* infile;
    infile = fopen(filename, "r");
    int i;
    char strLine[1000];
    char strKeyword[1000];
    char strTemp[1000];
    char strEnergyFile[1000];
    char strStructFile[1000];
    int nStructFiletype = 0;
    char strTempArr[3][1000];
    strEnergyFile[0] = '\0';
    strStructFile[0] = '\0';

    int nLine = 0;

    int nErr = 0; // error code

    for (i = 0; i < MAX_MV; i++)
        {
            faMCFreq_glb[i] = 0.0f; // initialization; to be normalized
        }
    nBiasPotential_Mode_glb    = 0;
    nAnnealing_Mode_glb        = -1;
    while (fgets(strLine, sizeof(strLine), infile) != NULL)
        {
            nLine++;

            strKeyword[0] = '#';
            sscanf(strLine, "%s", strKeyword);

            if (strKeyword[0] == '#')
                { // if the first character is #,
                  // ignore the line
                }
            else
                {
                    for (i = 0; strLine[i] != '\0'; i++)
                        {
                            if (strLine[i] == '#')
                                {
                                    strLine[i] = '\0'; // ignore the content after #
                                    break;
                                }
                        }

                    if (strcmp(strKeyword, "BOX_SIZE") == 0)
                        {
                            sscanf(strLine, "%*s %s %s %s", strTempArr[0], strTempArr[1], strTempArr[2]);
                            if (strTempArr[1][0] < '0' || strTempArr[1][0] > '9')
                                {
                                    naBoxSize_glb[0] = atoi(strTempArr[0]);
                                    naBoxSize_glb[1] = naBoxSize_glb[0];
                                    naBoxSize_glb[2] = naBoxSize_glb[0];
                                }
                            else
                                {
                                    if (strTempArr[2][0] < '0' || strTempArr[2][0] > '9')
                                        {
                                            naBoxSize_glb[0] = atoi(strTempArr[0]);
                                            naBoxSize_glb[1] = atoi(strTempArr[1]);
                                            naBoxSize_glb[2] = 0; // error handling required
                                        }
                                    else
                                        {
                                            naBoxSize_glb[0] = atoi(strTempArr[0]);
                                            naBoxSize_glb[1] = atoi(strTempArr[1]);
                                            naBoxSize_glb[2] = atoi(strTempArr[2]);
                                        }
                                }
                            nRDF_TotBins_glb = naBoxSize_glb[0] * 4;
                        }
                    else if (strcmp(strKeyword, "MC_TEMP") == 0)
                        {
                            sscanf(strLine, "%*s %f", &fKT_glb);
                        }
                    else if (strcmp(strKeyword, "N_STEPS") == 0)
                        {
                            sscanf(strLine, "%*s %ld", &nMCStepsPerCycle_glb);
                        }
                    else if (strcmp(strKeyword, "PREEQ_STEPS") == 0)
                        {
                            sscanf(strLine, "%*s %ld", &nMCStepsForTherm_glb);
                        }
                    else if (strcmp(strKeyword, "PREEQ_TEMP") == 0)
                        {
                            sscanf(strLine, "%*s %f", &fPreKT_glb);
                        }
                    else if (strcmp(strKeyword, "MC_TEMP_RATE") == 0)
                        {
                            sscanf(strLine, "%*s %f", &fMC_TempRate_glb);
                        }
                    else if (strcmp(strKeyword, "MC_TEMP_MODE") == 0)
                        {
                            sscanf(strLine, "%*s %d", &nAnnealing_Mode_glb);
                        }
                    else if (strcmp(strKeyword, "MC_DELTA_TEMP") == 0)
                        {
                            sscanf(strLine, "%*s %f", &fDeltaTemp_glb);
                        }
                    else if (strcmp(strKeyword, "MC_INVERT_TEMP") == 0)
                        {
                            sscanf(strLine, "%*s %d", &nTemp_inv_glb);
                        }
                    else if (strcmp(strKeyword, "MC_CYCLE_NUM") == 0)
                        {
                            sscanf(strLine, "%*s %d", &nTotCycleNum_glb);
                        }
                    else if (strcmp(strKeyword, "MC_INDENT_MODE") == 0)
                        {
                            sscanf(strLine, "%*s %d", &nBiasPotential_Mode_glb);
                        }
                    else if (strcmp(strKeyword, "MV_TRANS_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_TRANS]);
                        }
                    else if (strcmp(strKeyword, "MV_CLSTR_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_CLSTR]);
                        }
                    else if (strcmp(strKeyword, "MV_SMCLSTR_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_SMCLSTR]);
                        }
                    else if (strcmp(strKeyword, "MV_STROT_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_STROT]);
                        }
                    else if (strcmp(strKeyword, "MV_LOCAL_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_LOCAL]);
                        }
                    else if (strcmp(strKeyword, "MV_SNAKE_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_SNAKE]);
                        }
                    else if (strcmp(strKeyword, "MV_DBPVT_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_DBPVT]);
                        }
                    else if (strcmp(strKeyword, "MV_COLOCAL_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_COLOCAL]);
                        }
                    else if (strcmp(strKeyword, "MV_MTLOCAL_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_MTLOCAL]);
                        }
                    else if (strcmp(strKeyword, "MV_PIVOT_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_PIVOT]);
                        }
                    else if (strcmp(strKeyword, "MV_BRROT_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_BRROT]);
                        }
                    else if (strcmp(strKeyword, "MV_PR_SMCLSTR") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_PR_SMCLSTR]);
                        }
                    else if (strcmp(strKeyword, "MV_PR_CLSTR") == 0)
                        {
                            sscanf(strLine, "%*s %f", &faMCFreq_glb[MV_PR_CLSTR]);
                        }
                    else if (strcmp(strKeyword, "RESTART_FILE") == 0)
                        {
                            sscanf(strLine, "%*s %s", strRestartFile_glb);
                        }
                    else if (strcmp(strKeyword, "STRUCT_FILETYPE") == 0)
                        {
                            sscanf(strLine, "%*s %d", &nStructFiletype);
                        }
                    else if (strcmp(strKeyword, "STRUCT_FILE") == 0)
                        {
                            sscanf(strLine, "%*s %s", strStructFile);
                        }
                    else if (strcmp(strKeyword, "ENERGY_FILE") == 0)
                        {
                            sscanf(strLine, "%*s %s", strEnergyFile);
                        }
                    else if (strcmp(strKeyword, "RANDOM_SEED") == 0)
                        {
                            sscanf(strLine, "%*s %d", &nRNG_Seed_glb);
                            nRNG_Seed_glb = nRNG_Seed_glb == 0 ? time(NULL) : nRNG_Seed_glb;
                        }
                    else if (strcmp(strKeyword, "REPORT_PREFIX") == 0)
                        {
                            sscanf(strLine, "%*s %s", strReportPrefix_glb);
                        }
                    else if (strcmp(strKeyword, "REPORT_LOG_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %ld", &naReportFreqs_glb[REPORT_LOG]);
                        }
                    else if (strcmp(strKeyword, "REPORT_ENERGY_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %ld", &naReportFreqs_glb[REPORT_ENERGY]);
                        }
                    else if (strcmp(strKeyword, "REPORT_CONFIG_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %ld", &naReportFreqs_glb[REPORT_CONFIG]);
                        }
                    else if (strcmp(strKeyword, "REPORT_MCMOVE_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %ld", &naReportFreqs_glb[REPORT_MCMOVE]);
                        }
                    else if (strcmp(strKeyword, "REPORT_NETWORK_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %ld", &naReportFreqs_glb[REPORT_NETWORK]);
                        }
                    else if (strcmp(strKeyword, "REPORT_RDFTOT_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %ld", &naReportFreqs_glb[REPORT_RDFTOT]);
                        }
                    else if (strcmp(strKeyword, "REPORT_COMDEN_FREQ") == 0)
                        {
                            sscanf(strLine, "%*s %ld", &naReportFreqs_glb[REPORT_COMDEN]);
                        }
                    else if (strcmp(strKeyword, "ANALYSIS_CLUSTER_MODE") == 0)
                        {
                            sscanf(strLine, "%*s %d", &nClusteringMode_glb);
                        }
                    else if (strcmp(strKeyword, "REPORT_CONFIG_MODE") == 0)
                        {
                            sscanf(strLine, "%*s %ld", &nTrajMode_glb);
                        }
                    else if (strcmp(strKeyword, "BIAS_INDENT_RADIUS") == 0)
                        {
                            sscanf(strLine, "%*s %f", &fSquishRad_glb);
                        }
                    else if (strcmp(strKeyword, "BIAS_INDENT_STIFF") == 0)
                        {
                            sscanf(strLine, "%*s %f", &fSquish_Stiff_glb);
                        }
                    else if (strcmp(strKeyword, "BIAS_KEEP_ON") == 0)
                        {
                            sscanf(strLine, "%*s %d", &nBiasPotential_KeepON_glb);
                        }
                    else if (strcmp(strKeyword, "BIAS_COUPLED_TO_TEMP") == 0)
                        {
                            sscanf(strLine, "%*s %d", &nBiasPotential_CoupledToTemp_glb);
                        }
                    else
                        {
                            fprintf(stderr, "ERROR: unable to parse line %d in %s.\n%s", nLine, filename, strLine);
                            exit(1);
                        }
                }
        }

    fclose(infile);

    nBiasPotential_KeepON_glb       = nBiasPotential_KeepON_glb == 1 ? 1 : 0;
    nBiasPotential_CoupledToTemp_glb = nBiasPotential_CoupledToTemp_glb == 1 ? 1 : 0;

    if (nBiasPotential_CoupledToTemp_glb && nBiasPotential_KeepON_glb)
        {
            fprintf(stderr, "KeyFile error <%s>: Cannot have BIAS_COUPLED_TO_TEMP=1 and BIAS_KEEP_ON=1\n", filename);
            exit(1);
        }

    float freq_tot = 0.0f;
    for (i = MV_NULL + 1; i < MAX_MV; i++)
        {
            freq_tot += faMCFreq_glb[i];
        }

    if (freq_tot == 0)
        {
            for (i = MV_NULL + 1; i < MAX_MV; i++)
                {
                    faMCFreq_glb[i] = 1.0f / (float) (MAX_MV - 1);
                }
        }
    else
        {
            for (i = MV_NULL + 1; i < MAX_MV; i++)
                {
                    faMCFreq_glb[i] /= freq_tot;
                }
        }

    if (strEnergyFile[0] != '\0')
        {
            nErr = Parse_EnergyFile(strEnergyFile);
        }

    if (strStructFile[0] != '\0')
        {

            Parse_StructureFile_CalcBeadsAndChains(strStructFile, &tot_beads_glb, &tot_chains_glb,
                                                   &tot_chain_types_glb);

            CreateBeadsAndChains(tot_beads_glb, tot_chains_glb);

            Parse_StructureFile(strStructFile);

            if (nStructFiletype == 0)
                {
                    bReadConf_glb = 0;
                }
            else if (nStructFiletype == 1)
                {
                    printf("Reading restart file provided to generate initial "
                           "configuration.\n");
                    bReadConf_glb = 1;
                }
            else
                {
                    fprintf(stderr, "ERROR: undefined value in STRUCT_FILETYPE of %s.\nCrashing!\n", filename);
                    exit(1);
                }
        }
    else
        {
            bReadConf_glb = -1;
        }
    if (bReadConf_glb == 1)
        {
            if (strRestartFile_glb[0] != '\0')
                {
                    printf("Restart file is %s\n", strRestartFile_glb);
                }
            else
                {
                    nErr = 5;
                    printf("No restart file provided.\n");
                }
        }

    return nErr;
}

void ForEnergyMatrix_FillWithZeros()
{
    int i, j, k;

    for (i = 0; i < MAX_AA; i++)
        {
            for (j = 0; j < MAX_AA; j++)
                {
                    for (k = 0; k < MAX_E; k++)
                        {
                            faEnergy_glb[i][j][k] = 0.f;
                        }
                }
        }
}

/// Parse_EnergyFile - reads the energy file
/// Have to painstakingly go through every different keyword that exists in the energy file. Reads all the matrices
/// \param strEnFile
/// \return
int Parse_EnergyFile(char* strEnFile)
{

    ForEnergyMatrix_FillWithZeros();

    int nRes = 0;

    FILE* infile;
    infile = fopen(strEnFile, "r");

    char strLine[256];
    int nFlag   = 0;
    char bOrder = 0;
    char strKey[256];
    int nRow;
    float fTemp[MAX_AA] = {0.f};
    int i, j;
    int nEntry = 0;

    while (fgets(strLine, sizeof(strLine), infile) != NULL)
        {
            if (strLine[0] == '#')
                { // All key-words start with '#'
                    sscanf(&strLine[1], "%s", strKey);
                    // STICKERS must be the first key-word
                    if (strcmp(strKey, "STICKERS") == 0)
                        {
                            nFlag = -1;
                        }
                    else if (bOrder == 0)
                        {
                            fprintf(stderr, "ERROR: #STICKERS is not the first entry in %s.\n", strEnFile);
                            nRes = 1;
                            break;
                        }
                    else
                        { // Checking for all other keywords
                            if (strcmp(strKey, "OVERLAP_POT") == 0)
                                {
                                    nFlag = 2 * (E_OVLP);
                                }
                            else if (strcmp(strKey, "CONTACT_POT") == 0)
                                {
                                    nFlag = 2 * (E_CONT);
                                }
                            else if (strcmp(strKey, "CONTACT_RAD") == 0)
                                {
                                    nFlag = 2 * (E_CONT) + 1;
                                }
                            else if (strcmp(strKey, "SC_SC_POT") == 0)
                                {
                                    nFlag = 2 * (E_SC_SC);
                                }
                            else if (strcmp(strKey, "FSOL_POT") == 0)
                                {
                                    nFlag = 2 * (E_F_SOL);
                                }
                            else if (strcmp(strKey, "T_IND_POT") == 0)
                                {
                                    nFlag = 2 * (E_T_IND);
                                }
                            else if (strcmp(strKey, "STIFF_POT") == 0)
                                {
                                    nFlag = 2 * (E_STIFF);
                                }
                            else if (strcmp(strKey, "LINKER_LENGTH") == 0)
                                {
                                    nFlag = -3;
                                }
                            else if (strcmp(strKey, "LINKER_SPRCON") == 0)
                                {
                                    nFlag = -4;
                                }
                            else if (strcmp(strKey, "LINKER_EQLEN") == 0)
                                {
                                    nFlag = -5;
                                }
                            else
                                {
                                    fprintf(stderr, "ERROR: irregular expression in %s: %s\n", strEnFile, strKey);
                                    nRes = 2;
                                    break;
                                }
                        }
                    // Should know which key-word we are dealing with.
                    // Next line should have numerical values. Either one number, or
                    // nBeadType numbers
                    nRow = 0;
                } // The line did not contain a keyword. Either empty or has numbers.
            else if (strcmp(strLine, "\r\n") && strcmp(strLine, "\n"))
                { // ignore empty lines
                    if (nFlag == -1)
                        { // sticker
                            sscanf(strLine, "%d", &nBeadTypes_glb);
                            if (nBeadTypes_glb > MAX_AA)
                                {
                                    fprintf(stderr, "ERROR: the number of AA types exceeds MAX_AA in %s.\n", strEnFile);
                                    nRes = 3;
                                    break;
                                }
                            bOrder = 1;
                        }
                    else if (nFlag == -3)
                        { // linker_length
                            sscanf(strLine, "%f", &fLinkerLength_glb);
                        }
                    else if (nFlag == -4)
                        { // linker_sprcon
                            sscanf(strLine, "%f", &fLinkerSprCon_glb);
                        }
                    else if (nFlag == -5)
                        { // linker_eqlen
                            sscanf(strLine, "%f", &fLinkerEqLen_glb);
                        }
                    else if (nFlag == 0)
                        {
                            fprintf(stderr, "ERROR: nFlag is not assigned in %s.\n", strEnFile);
                            nRes = 4;
                            break;
                        }
                    else
                        {                                      // Saving energy and radius matrices
                            nEntry = str2farr(strLine, fTemp); // Counting how many columns
                            if (nEntry != nBeadTypes_glb)
                                { // If columns not the same as sticker number
                                    if (nEntry == 1)
                                        { // If only 1-value, all-values in  that
                                          // matrix are this one
                                            for (i = 0; i < nBeadTypes_glb; i++)
                                                {
                                                    for (j = 0; j < nBeadTypes_glb; j++)
                                                        {
                                                            if (nFlag % 2 == 0)
                                                                { // energy
                                                                    faEnergy_glb[i][j][(int) (nFlag / 2)] = fTemp[0];
                                                                }
                                                            else
                                                                { // radius
                                                                    faEnRad_glb[i][j][(int) (nFlag / 2)] = fTemp[0];
                                                                }
                                                        }
                                                }
                                        }
                                    else
                                        {
                                            fprintf(stderr,
                                                    "ERROR: irregular expression in energy "
                                                    "matrices of %s.\n",
                                                    strEnFile);
                                            nRes = 4;
                                            break;
                                        }
                                }
                            else
                                { // Column number same as sticker-number; save the
                                  // matrix
                                    if (nFlag % 2 == 0)
                                        { // energy
                                            for (i = 0; i < nBeadTypes_glb; i++)
                                                {
                                                    faEnergy_glb[nRow][i][(int) (nFlag / 2)] = fTemp[i];
                                                }
                                        }
                                    else
                                        { // radius
                                            for (i = 0; i < nBeadTypes_glb; i++)
                                                {
                                                    faEnRad_glb[nRow][i][(int) (nFlag / 2)] = fTemp[i];
                                                }
                                        }
                                    nRow++;
                                }
                        }
                }
        }

    fclose(infile);

    return nRes;
}

/// str2farr - converts the string array to an array of floats. The nice thing
/// about C is the ease of string manipulation. used to generate matrices from
/// the energy file.
/// \param strRaw
/// \param fArray
/// \return Changes fArray to include be the array in
int str2farr(char strRaw[], float fArray[MAX_AA])
{
    int i         = 0;
    char* strTemp = strRaw;
    char* token;

    for (i = 0; i < MAX_AA; i++)
        {
            token = strtok_r(strTemp, " \t", &strTemp);
            if (token == NULL)
                {
                    break;
                }
            else
                {
                    fArray[i] = atof(token);
                }
        }

    return i;
}

/// Parse_StructureFile - reads the structure file filename.
/// The format of the structure file is in the function below. Should be easier
/// to use Python to generate the files, usually.
/// \param filename
void Parse_StructureFile(char* filename)
{
    /*
    This function reads in a structure file that also includes topology
    information. The format is: # The '#' is the commenting character. NEW{
    nCopiesOfMolecule
    #AtomID AtomType LinkerLengths BondedPartner
    }END
    Each molecule is within NEW{\0 ... \0}END which act as keywords. After NEW{,
    the next line contains the number of copies to be made of this molecule. It
    is also expected that whenever a bead is listed, ALL of its bonds follow in
    the next lines. E.G A molecule with 5 branches where each branch has 2 beads
    would be:
    ###
    NEW{
    nCopies
    0 0 1 1
    0 0 1 2
    0 0 1 3
    0 0 1 4
    0 0 1 5
    1 1 1 0
    1 1 1 6
    2 1 1 0
    2 1 1 7
    3 1 1 0
    3 1 1 8
    4 1 1 0
    4 1 1 9
    5 1 1 0
    5 1 1 10
    6 2 1 1
    7 2 1 2
    8 2 1 3
    9 2 1 4
    10 2 1 5
    }END
    ####
    As such, each molecule is independent and all molecule beads start from 0!
    For each listed bead, linker_len_glb[i][j] corresponds to the linker constraint
    between bead i, and bead (topo_info_glb[i][j]). After each molecule is read,
    nCopies copies are made, and THEN, the next line in the file is read.
    */
    FILE* inFile;
    inFile = fopen(filename, "r"); // Opening the file!

    char strLine[1000];    // Used to store each line from the input file.
    char strKeyword[1000]; // Used to convert strLine into specific keywords.
    int curID;             // Used to track the current beadID
    int curType, curLinker,
        curPartner;               // Used to track current
                                  // bead-type,linker-length(s),bond-partner.
    int nCopies;                  // Used to store how many copies to make.
    int nChainStart;              // Tracking which beadID each chain starts from.
    int i, j, k;                  // Just some iterators for looping
    int nCursor, nCursor2, nTemp; // Internal iterators to make sure beads and
                                  // their partners are not double counted.
    int nChainID,
        nChainType; // Internal iterator to count which chainID the chain is.
    int nFlag;      // Internal flag used to see which keyword is being read!
    int nBEADS;     // Internal counter to count unique beads in each chain!

    // Initialize the bead_info_glb and chain_info_glb
    for (i = 0; i < tot_beads_glb; i++)
        {
            for (j = 0; j < BEADINFO_MAX; j++)
                {
                    bead_info_glb[i][j] = -1;
                }
            for (j = 0; j < MAX_BONDS; j++)
                {
                    topo_info_glb[i][j]  = -1;
                    linker_len_glb[i][j] = -1;
                }
        }

    for (i = 0; i < tot_chains_glb; i++)
        {
            for (j = 0; j < CHAININFO_MAX; j++)
                {
                    chain_info_glb[i][j] = -1;
                }
        }

    // Initialization of the counters and iterators
    nCursor     = -1;
    nCursor2    = -1;
    nChainID    = -1;
    nChainType  = -1;
    nFlag       = -1;
    nTemp       = -1;
    nBEADS      = 0;
    nChainStart = 0;
    nCopies     = 0;
    while ((fgets(strLine, sizeof(strLine), inFile) != NULL) && (nFlag == -1))
        {
            // Keep reading the file until it ends or it's incorrectly formatted.

            strKeyword[0] = '#';               // This is the commenting character in all input files.
            sscanf(strLine, "%s", strKeyword); // Plain old read the line.

            if (strKeyword[0] != '#')
                { // If the line is not a comment, see what keyword it is
                    for (i = 0; strLine[i] != '\0'; i++)
                        { // Keep reading till end of that line
                            if (strLine[i] == '#')
                                {                      // This line now has a comment, so
                                                       // stop reading further
                                    strLine[i] = '\0'; // Just assume we have reached the end of the line
                                    break;
                                }
                        }
                    if (strcmp(strKeyword, "NEW{") == 0)
                        { // This signifies a new molecule type has been started
                            nFlag = 1;
                        }
                    if (strcmp(strKeyword, "}END") == 0)
                        { // This signifies a new molecule type has been started
                            nFlag = -1;
                        }
                }

            if (nFlag == 1)
                { // This signifies that a new molecule has started
                    // It's assumed that the next line contains the number of copies for
                    // this molecule.
                    nChainType++;
                    nChainTypeIsLinear_glb[nChainType] = 1; // Assume all chains are linear to begin with.
                    nChainStart += nBEADS;
                    nChainID++;
                    nBEADS = 0;
                    fgets(strLine, sizeof(strLine),
                          inFile); // Reading the next line, which has nMOLS
                    sscanf(strLine, "%d",
                           &nCopies); // Remembering how many copies top make.
                    while (fgets(strLine, sizeof(strLine), inFile) != NULL && nFlag == 1)
                        {
                            sscanf(strLine, "%s", strKeyword); // Just reading as a string
                            if (strcmp(strKeyword, "}END") == 0)
                                { // End of this molecule
                                    nFlag = -1;
                                    break;
                                }
                            sscanf(strLine, "%d %d %d %d", &curID, &curType, &curLinker, &curPartner);
                            curID += nChainStart; // Accounting for previously defined beads
                            if (bead_info_glb[curID][BEAD_TYPE] == -1)
                                { // This is to make sure that each bead is counted
                                  // once only even if it has many bonds.
                                    nBEADS++;
                                    bead_info_glb[curID][BEAD_TYPE]    = curType;
                                    bead_info_glb[curID][BEAD_CHAINID] = nChainID;
                                    nCursor = 0; // This is a counter for number of bonds, which
                                                 // should reset when you have a 'new' bead.
                                }

                            if (curPartner != -1)
                                { // This bead has a bonded partner
                                    if (nCursor > 1)
                                        { // This signifies that the chain is not
                                          // linear because a bead has more than two bonds
                                          // because indicies start at 0.
                                            nChainTypeIsLinear_glb[nChainType] = 0;
                                        }
                                    curPartner += nChainStart;                  // Accounts for all beads before,
                                                                                // like above.
                                    topo_info_glb[curID][nCursor] = curPartner; // Adding the ID of the partner
                                    linker_len_glb[curID][nCursor] =
                                        curLinker * (int) fLinkerLength_glb; // Adding the linker constraint.
                                    nCursor++;
                                }
                        }
                    chain_info_glb[nChainID][CHAIN_START]  = nChainStart;
                    chain_info_glb[nChainID][CHAIN_LENGTH] = nBEADS;
                    chain_info_glb[nChainID][CHAIN_TYPE]   = nChainType;
                    // We just fully store the first chain 'manually'. Now we just copy
                    // the chain nCopies times.
                    for (k = 1; k < nCopies; k++)
                        {                          // Now we just copy this molecule nMOL-1 times
                            nChainStart += nBEADS; // This accounts for chain lengths.
                            nChainID++;            // Going to the next chainID
                            chain_info_glb[nChainID][CHAIN_START]  = nChainStart;
                            chain_info_glb[nChainID][CHAIN_LENGTH] = nBEADS;
                            chain_info_glb[nChainID][CHAIN_TYPE]   = nChainType;

                            for (i = 0; i < nBEADS; i++)
                                {
                                    curID                              = i + nChainStart;
                                    nTemp                              = curID - nBEADS;
                                    bead_info_glb[curID][BEAD_TYPE]    = bead_info_glb[nTemp][BEAD_TYPE];
                                    bead_info_glb[curID][BEAD_CHAINID] = nChainID;
                                    for (j = 0; j < MAX_BONDS; j++)
                                        {
                                            linker_len_glb[curID][j] = linker_len_glb[nTemp][j];
                                            if (topo_info_glb[nTemp][j] != -1)
                                                {
                                                    topo_info_glb[curID][j] =
                                                        topo_info_glb[nTemp][j] + nBEADS; // Because we must account for
                                                                                          // chain lengths
                                                }
                                        }
                                }
                        }
                    nFlag = -1; // Resetting the flag.
                }

            // This is to make sure that the reading was done correctly.
            if (nFlag != -1)
                {
                    printf("Incorrectly formatted input structure file. I must crash "
                           ":(\n\n");
                    exit(1);
                }
        }
    fclose(inFile);
}

/// Parse_StructureFile_CalcBeadsAndChains - reads the structure-file filename,
/// and records the total number of beads, chains, bead-types, and chain-types.
/// \param filename: Full path of the file, or name of file if in the same directory.
/// \param n_bead_num: Stores how many total beads are in the structure file.
/// \param n_chain_num: Stores how many total chains are in the structure file.
/// \param n_chain_types: Stores how many different chain-types are in the file.
void Parse_StructureFile_CalcBeadsAndChains(char* filename, int* n_bead_num, int* n_chain_num,
                                            int* n_chain_types)
{
    int dum_beads       = 0;
    int dum_chains      = 0;
    int dum_chain_types = 0;
    int per_chain_num   = 0;
    int per_ch_bd_num   = 0;
    int errCode         = 0;
    int nFlag           = -1;
    char strLine[1000];
    char strKey[1000];

    int n_old_bd_id, n_new_bd_id;

    // So that this function always sets the values to 0.
    *n_chain_types = 0;
    *n_bead_num    = 0;
    *n_chain_num   = 0;

    FILE* inFile;
    inFile = fopen(filename, "r");

    while (fgets(strLine, sizeof(strLine), inFile) != NULL && errCode == 0)
        {
            sscanf(strLine, "%s", strKey);
            if (strKey[0] == '#')
                { // Ignore comments
                    continue;
                }

            if (strcmp(strKey, "NEW{") == 0)
                { // New molecule is starting
                    per_ch_bd_num = 0;
                    fgets(strLine, sizeof(strLine), inFile);
                    sscanf(strLine, "%s", strKey);
                    if (strKey[0] == '#')
                        {
                            errCode = 1;
                            break;
                        }

                    nFlag = sscanf(strLine, "%d", &per_chain_num);
                    if (nFlag != 1)
                        {
                            errCode = 1;
                            break;
                        }

                    dum_chain_types++;
                    dum_chains += per_chain_num;

                    n_new_bd_id   = 0;
                    n_old_bd_id   = -1;
                    per_ch_bd_num = 0;
                    while ((fgets(strLine, sizeof(strLine), inFile) != NULL) && (errCode == 0))
                        {
                            sscanf(strLine, "%s", strKey);
                            if (strKey[0] == '#')
                                {
                                    errCode = 1;
                                    break;
                                }
                            if (strcmp(strKey, "}END") == 0)
                                { // Molecule has ended
                                    dum_beads += per_ch_bd_num * per_chain_num;
                                    break;
                                }

                            nFlag = sscanf(strLine, "%d", &n_new_bd_id);
                            if (nFlag != 1)
                                {
                                    errCode = 1;
                                    break;
                                }

                            if (n_new_bd_id != n_old_bd_id)
                                {
                                    per_ch_bd_num++;
                                }
                            n_old_bd_id = n_new_bd_id;
                        }
                }
        }

    *n_bead_num    = dum_beads;
    *n_chain_num   = dum_chains;
    *n_chain_types = dum_chain_types;

    fclose(inFile);
}

/// CreateBeadsAndChains. Allocate memory for the given number of beads, chains,
/// and chain-types.
/// \param n_bead_num
/// \param n_chain_types
/// \param n_chain_num
void CreateBeadsAndChains(size_t n_bead_num, size_t n_chain_num)
{
    char strTemp[100];

    strcpy(strTemp, "Chain Info.");
    chain_info_glb = Create2DInt(CHAININFO_MAX, n_chain_num, strTemp);

    strcpy(strTemp, "Topo Info.");
    topo_info_glb = Create2DInt(MAX_BONDS, n_bead_num, strTemp);

    strcpy(strTemp, "Linker Len.");
    linker_len_glb = Create2DInt(MAX_BONDS, n_bead_num, strTemp);

    strcpy(strTemp, "Bead Info.");
    bead_info_glb = Create2DInt(BEADINFO_MAX, n_bead_num, strTemp);

    strcpy(strTemp, "Old Bead.");
    old_bead_glb = Create2DInt(BEADINFO_MAX, n_bead_num, strTemp);
}
