#ifndef _GLOBAL_H_   // include guard
#define _GLOBAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define PI 3.14159

#define MAX_BEADS      53000
#define MAX_AA         6
#define MAX_CHAINTYPES 6
#define MAX_CHAINS     10000
#define MAX_CHAINLEN   100
#define MAX_BONDS      15
#define MAX_VALENCY    50

// energy parameters
#define E_TOT   0 // index zero must be assigned for total energy
#define E_OVLP  1
#define E_CONT  2
#define E_SC_SC 3
#define E_STIFF 4
#define MAX_E   5 // just for counting; must be the last number of this list

// MC move parameters
#define MV_NULL     0 // index zero indicates null move
#define MV_FACEC    1
#define MV_LOCAL    2
#define MV_COLOCAL  3
#define MV_SHAKE    4
#define MV_SNAKE    5
#define MV_TRANS    6
#define MV_SMCLSTR  7
#define MV_CLSTR    8
#define MV_PIVOT    9
#define MV_BRROT   10
#define MV_DBPVT   11
#define MAX_MV     12 // just for counting; must be the last number of this list

// report parameters
#define REPORT_NULL    0 // index zero
#define REPORT_LOG     1
#define REPORT_ENERGY  2
#define REPORT_CONFIG  3
#define REPORT_MCMOVE  4
#define REPORT_NETWORK 5
#define REPORT_RDFTOT  6
#define MAX_REPORT     7 // just for counting; must be the last number of this list

// auxiliary definitions
#define POS_X          0
#define POS_Y          1
#define POS_Z          2
#define POS_MAX        3

#define BEAD_CHAINID   3
#define BEAD_TYPE      4
#define BEAD_FACE      5
#define BEADINFO_MAX   6

#define CHAIN_TYPE     0
#define CHAIN_LENGTH   1
#define CHAIN_START    2
#define CHAININFO_MAX  3

#define RDF_MAXBINS     3000
#define TEMP_CYCLES_MAX 30
#define RDF_COMPS       22
#define MAX_ROTSTATES   27



// configurations
typedef int lInt;
typedef long lLong;
typedef double lDub;
typedef  long double lLDub;
lInt   bead_info[MAX_BEADS][BEADINFO_MAX];
lInt   old_bead[MAX_BEADS][BEADINFO_MAX]; //Redundant copy for MCSteps
lInt   linker_len[MAX_BEADS][MAX_BONDS];//Remember that this one is an INT, not float
lInt   topo_info[MAX_BEADS][MAX_BONDS];
lInt   chain_info[MAX_CHAINS][CHAININFO_MAX];
lInt   tot_beads;
lInt   tot_chains;
lInt   Temp_Mode;
lInt   Indent_Mode, RotBias_Mode;

// system setup
lInt   nBoxSize[POS_MAX];
lInt   LocalArr[MAX_ROTSTATES-1][3];//Used to quickly iterate over nearby points in a R-cube
lInt   Rot_IndArr[MAX_ROTSTATES-1];
char   bReadConf;

// energy matrices for stickers
lInt   nSeqEn;
float fEnergy[MAX_AA][MAX_AA][MAX_E];
float fEnRad[MAX_AA][MAX_AA][MAX_E];
lInt   rot_trial[MAX_VALENCY][MAX_ROTSTATES];//Used in orientational-bias MC moves
lDub  bolt_fac[MAX_ROTSTATES - 1];//Used in orientational-bias
lDub  bolt_norm[MAX_VALENCY];
lDub  dbias_bolt_fac[MAX_AA][MAX_AA];//For pre-calculating the factors.
float   faCurrEn[MAX_E]; //Vector for current energy

//Arrays to track certain topology and interaction information
lInt TypeCanRot[MAX_AA];//Used to track if that beadType interacts via rotations.
lInt TypeIsLinear[MAX_CHAINTYPES];//Used to track if this chainType is linear.
lInt TypeCanOvlp[MAX_AA];//Used to track if a certain beadType has an overlap cost.

float fLinkerLength;
float fLinkerSprCon;
float fLinkerEqLen;

// MC setup
float fKT, fPreKT, fCuTemp, fRot_Bias, f_globRotBias, fdelta_temp;
lLong  nSteps, nPreSteps;
float fMCFreq[MAX_MV];
lInt   nMCMaxTrials, nTot_CycleNum;

// random number generator seed
lInt   seed;

// report-related
char  strReportPrefix[100];
char fileEnergy[100];
char fileStruct[100];
char fileMCMove[100];
char fileSysProp[100];
lLong  nReport[MAX_REPORT];
//Matrix to store acceptances and rejections 0: Rejected; 1: Accepted
lLong   MCAccepMat[2][MAX_MV];


// Cluster analysis and MCMoves
lInt naCluster[MAX_CHAINS][MAX_CHAINS];
lInt naList[MAX_CHAINS];
lLong *naClusHistList;
lInt *naChainCheckList;
lInt *naChainCheckList2;
lInt nClusListCounter;
lLDub **ld_TOTCLUS_ARR;
lInt naTempR[POS_MAX];
lInt nLargestClusterRightNow;

//Radial Distribution Function
float fRDF_TOT[RDF_MAXBINS];
lLDub ldRDF_ARR[RDF_COMPS][RDF_MAXBINS];
lLDub ld_TOTRDF_ARR[TEMP_CYCLES_MAX][RDF_COMPS][RDF_MAXBINS];
lInt nrdfCounter;//This counts how many times the RDF has been calculated for averaging at the end.
lInt nBins_RDF;
float fGyrTensor[7];//Gyration tensor
float fSysGyrRad;//Gyration radius of the system.
lLDub ld_TOTGYRRAD_ARR[TEMP_CYCLES_MAX][2];
lInt nTotGyrRadCounter;//Counter for total averaging

//Lattice To Remember Things
lInt *naTotLattice;
//Random bookkeeping to reduce function overhead


#endif // _GLOBAL_H_
