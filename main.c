#include "global.h"
#include "parsekey.h"
#include "print.h"
#include "initialize.h"
#include "structure.h"
#include "energy.h"
#include "mcmove.h"
#include "cluster.h"

int main(int argc, char *argv[]) {
    char keyfile[100];
    char lBrace[] = "<======      ";
    char rBrace[] = "      ======>";
    if (argc == 3 && strcmp(argv[1], "-k") == 0) {
        strcpy(keyfile, argv[2]);
    } else {
        strcpy(keyfile, "param.key");
    }

    int errorcode;
    errorcode = Parse_Keyfile(keyfile);
    if (errorcode == 0) {
        printf("Key file %s was successfully parsed.\n\n", keyfile);
        Print_Key();
    } else {
        printf("ERROR: unable to parse key file %s.\n", keyfile);
        exit(1);
    }

    printf("%s Initializing %s\n", lBrace, rBrace);

    srand(RNG_Seed == 0 ? time(NULL) : RNG_Seed);


    clock_t tStart = clock();

    Memory_Initialization_AtStart();
    Global_Array_Initialization_AtStart();

    if (bReadConf == 0) {
        Initial_Conditions_Simple();
    } else if (bReadConf == 1) {
        Initial_Conditions_FromFile();
        if (nPreSteps > 0) {
            Initial_Conditions_BreakBonds();
        }
    } else {
        printf("Wrong initial conditions condition give. Crashing!\n");
        exit(1);
    }

    if (Check_System_Structure() == 0) {
        printf("Check structure sanity: OK\n");
    } else {
        printf("ERROR: wrong structure.\n");
        exit(1);
    }

    clock_t tEnd = clock();
    double elapsed_time = (tEnd - tStart) / (double) CLOCKS_PER_SEC;
    printf("Initialization done. %.2f sec elapsed.\n", elapsed_time);

    tStart = clock();

    long nGen;
    int nMCInfo; // MC accept/reject

    printf("%s Beginning MC Simulation %s\n", lBrace, rBrace);
    printf("******************************************\n\n");
    printf("_____________________\n");
    printf("Thermalizing system.\n");
    printf("---------------------\n\n");

    fCuTemp = fPreKT;
    Print_Data(-1, -1);//Initialization of files
    for (nGen = 0;
         nGen < nPreSteps; nGen++) {//Intentionally not performing any data acquisition in the thermalizing phase
        nMCInfo = MC_Step_Equil(fCuTemp);
        //printf("(%d,%d)\n", nMCInfo / 12, nMCInfo % 2);
        Print_Data(nGen, -1);
    }

    printf("____________________________\n");
    printf("System has been thermalized!\n");
    printf("----------------------------\n\n");
/*
The system has thermalized!
*/
    int run_cycle;

    for (run_cycle = 0; run_cycle < nTot_CycleNum; run_cycle++) {
        fKT = fKT_Cycle[run_cycle];//temp_cycle[run_cycle];
        Calculate_Rot_Bias(fKT_Cycle[run_cycle]);
        Print_Data(-1, run_cycle);
        for (nGen = 0; nGen <= nSteps; nGen++) {
            fCuTemp = Temperature_Function(Temp_Mode, nGen);
            nMCInfo = MC_Step(fCuTemp);
            //printf("(%d,%d)\n", nMCInfo / 12, nMCInfo % 2);
            Print_Data(nGen, run_cycle);
        }
        Temp_Mode = -1;
        Copy_Data(run_cycle);
        Reset_Global_Arrays();
    }

    Write_TotalSysProp(fileSysProp, run_cycle);

    tEnd = clock();
    elapsed_time = (tEnd - tStart) / (double) CLOCKS_PER_SEC;
    printf("____________________________________\n");
    printf("Simulation finished in %.2f minutes.\n", elapsed_time / 60.);
    printf("------------------------------------\n");
    printf("%s ENDING %s\n", lBrace, rBrace);
    printf("******************************************\n");

    return 0;
}
