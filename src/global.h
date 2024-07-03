#ifndef _GLOBAL_H_ // include guard
#define _GLOBAL_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_AA         10
#define MAX_CHAINTYPES 10
#define MAX_CHAINLEN   500
#define MAX_BONDS      4
#define MAX_VALENCY    500

// energy parameters
#define E_TOT   0 // index zero must be assigned for total energy
#define E_OVLP  1
#define E_CONT  2
#define E_SC_SC 3
#define E_F_SOL 4
#define E_T_IND 5
#define E_STIFF 6
#define E_BIAS  7
#define MAX_E   8 // just for counting; must be the last number of this list

// MC move parameters
#define MV_NULL       0 // index zero indicates null move
#define MV_STROT      1
#define MV_LOCAL      2
#define MV_COLOCAL    3
#define MV_MTLOCAL    4
#define MV_SNAKE      5
#define MV_TRANS      6
#define MV_SMCLSTR    7
#define MV_CLSTR      8
#define MV_PIVOT      9
#define MV_BRROT      10
#define MV_DBPVT      11
#define MV_PR_SMCLSTR 12
#define MV_PR_CLSTR   13
#define MAX_MV        14 // just for counting; must be the last number of this list

// report parameters
#define REPORT_NULL    0 // index zero
#define REPORT_LOG     1
#define REPORT_ENERGY  2
#define REPORT_CONFIG  3
#define REPORT_MCMOVE  4
#define REPORT_NETWORK 5
#define REPORT_RDFTOT  6
#define REPORT_COMDEN  7
#define MAX_REPORT     8 // just for counting; must be the last number of this list

// auxiliary definitions
#define POS_X   0
#define POS_Y   1
#define POS_Z   2
#define POS_MAX 3

#define BEAD_CHAINID 3
#define BEAD_TYPE    4
#define BEAD_FACE    5
#define BEADINFO_MAX 6

#define CHAIN_TYPE    0
#define CHAIN_LENGTH  1
#define CHAIN_START   2
#define CHAININFO_MAX 3

#define MAX_ROTSTATES       27
#define CLUS_CONTACT_NEIGHS 27

typedef int lInt;
typedef long lLong;
typedef double lDub;
typedef long double lLDub;

// configurations and structural info

lInt** bead_info_glb;
lInt** old_bead_glb;
lInt** chain_info_glb;
lInt** topo_info_glb;
lInt** linker_len_glb;

int tot_beads_glb;
int tot_chains_glb;
int tot_chain_types_glb;
lInt nAnnealing_Mode_glb, nTemp_inv_glb;
lInt nBiasPotential_Mode_glb, RotBias_Mode_glb;
lInt nBiasPotential_KeepON_glb, nBiasPotential_CoupledToTemp_glb;

// system setup
lInt naBoxSize_glb[POS_MAX];
lInt naLocalArr_glb[MAX_ROTSTATES - 1][3]; // Used to quickly iterate over nearby points in an R-cube
lInt naRot_IndArr_glb[MAX_ROTSTATES - 1];
char bReadConf_glb;
char bSystemHasTopo_glb;
char bSystemHasCont_glb;
char bSystemHasOvlp_glb;
char bSystemHasFSol_glb;
char bSystemHasSCSC_glb;

// energy matrices for stickers
lInt nBeadTypes_glb;
float faEnergy_glb[MAX_AA][MAX_AA][MAX_E];
float faEnRad_glb[MAX_AA][MAX_AA][MAX_E];
#define LARGEST_RADIUS 3
lInt naRotTrial_glb[MAX_VALENCY][MAX_ROTSTATES]; // Used in orientational-bias MC moves
lLDub ldaBoltzFac_glb[MAX_ROTSTATES - 1];        // Used in orientational-bias
lLDub ldaBoltzNorm_glb[MAX_VALENCY];
lLDub ldaBoltzFacNorm_glb[MAX_AA][MAX_AA]; // For pre-calculating the factors.
lLDub ldLogOfSmallestPossibleProb_glb;     // Smallest probability possible logl(1/RAND_MAX)

float faCurrEn_glb[MAX_E]; // Vector for current energy

// Arrays to track certain topology and interaction information
lInt nBeadTypeIsSticker_glb[MAX_AA];         // Used to track if that beadType interacts via rotations.
lInt nChainTypeIsLinear_glb[MAX_CHAINTYPES]; // Used to track if this chainType is linear.
lInt nBeadTypeCanOvlp_glb[MAX_AA];           // Used to track if a certain beadType has an overlap cost.
lInt nBeadTypeCanCont_glb[MAX_AA];           // Used to track if a certain beadType has contact interactions
lInt nBeadTypeCanFSol_glb[MAX_AA];           // Used to track if a certain beadType has solvation energies
lInt nBeadTypeCanTInd_glb[MAX_AA]; // Used to track if a certain beadType has temperature independent solvation

float fLinkerLength_glb;
float fLinkerSprCon_glb;
float fLinkerEqLen_glb;

// MC setup
float fKT_glb, fPreKT_glb, fCuTemp_glb, fDeltaTemp_glb, fMC_TempRate_glb, fSquishRad_glb,
    fSquishRad_Sq_glb, fSquish_Stiff_glb;
float* faKT_Cycle_glb;
lLong nMCStepsPerCycle_glb, nMCStepsForTherm_glb;
float faMCFreq_glb[MAX_MV];
lInt nTotCycleNum_glb;

// random number generator nRNG_Seed_glb
lInt nRNG_Seed_glb;

// report-related
char strReportPrefix_glb[512];
char strFileEnergy_glb[512];
char strFileTraj_glb[512];
char strFileMCMove_glb[512];
char strFileSysProp_glb[512];
char strRestartFile_glb[512];
lLong naReportFreqs_glb[MAX_REPORT]; // Array to store report frequencies.
lLong nTrajMode_glb;
// Matrix to store acceptances and rejections 0: Rejected; 1: Accepted
// TODO: Have a more extensive way to record also where/when a particular move fails -- not just if it fails.
lLong naMCAccepMat_glb[2][MAX_MV];

// Cluster analysis
lInt nClusteringMode_glb;
lInt** naClusterMatrix_glb;
lInt* naList_glb;
lLong* laClusHistList_glb;
lInt* naChainCheckList_glb;
lInt nTotClusCounter_glb;
lLDub** ldaTOTCLUS_Arr_glb;
lLDub* ldaMOLCLUS_Arr_glb;
lLDub* ldaTOTMOLCLUS_Arr_glb;
lInt naTempR_glb[POS_MAX];
lInt nLargestClusterRightNow_glb;
lInt nLimitedClusterSize_glb;

// Neighbor search
lInt* naOldOvlpNeighs_glb;
lInt* naNewOvlpNeighs_glb;
lInt* naOldContNeighs_glb;
lInt* naNewContNeighs_glb;

// Radial Densities and PDFs
lLDub* ldaTOTRDF_Arr_glb;
lLDub* ldaTOTRadDen_Arr_glb;
lLDub* ldaRDF_Arr_glb;
lLDub* ldaRadDen_Arr_glb;
lInt nRDF_TotComps_glb;
lInt nTotRDFCounter_glb; // This counts how many times the RDF has been calculated for averaging at the end.
lInt nRDF_TotBins_glb;
lInt nRadDen_TotComps_glb;
lInt nRadDen_CompShift_glb;
lInt nTotRadDenCounter_glb; // This counts for the Radial Density histograms

// Gyration tensor
float faGyrTensor_glb[7]; // Gyration tensor
float faSysGyrRad_glb;    // Gyration radius of the system.
lLDub** ldaTOTRg_Arr_glb;
lInt nTotGyrRadCounter_glb; // Counter for total averaging

// Trajectory Saving
lInt* naTOTTRAJ_Arr_glb;
lLong nTraj_FramesPerCycle_glb;
lInt nTrajCurFrame_glb;

// Lattice To Remember Things
lInt* naTotLattice_glb;

#define LINKER_RSCALE 1.74f

#endif // _GLOBAL_H_
