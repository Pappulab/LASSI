#include "cluster.h"
#include "energy.h"
#include "global.h"
#include "initialize.h"
#include "mcmove.h"
#include "parsekey.h"
#include "print.h"
#include "structure.h"

int main(int argc, char* argv[])
{
    // Read in the system commands
    char keyfile[100];
    if (argc == 3 && strcmp(argv[1], "-k") == 0)
        {
            strcpy(keyfile, argv[2]);
        }
    else
        {
            strcpy(keyfile, "param.key");
        }
    // Text formatting helpers.
    char lBrace[] = "<======      ";
    char rBrace[] = "      ======>";
    // Read the param file!
    int errorcode;
    errorcode = Parse_Keyfile(keyfile);
    if (errorcode == 0)
        {
            printf("Key file %s was successfully parsed.\n\n", keyfile);
            srand(nRNG_Seed_glb);
            ScreenIO_Print_KeyFile();
        }
    else
        {
            printf("ERROR: unable to parse key file %s.\n", keyfile);
            exit(1);
        }

    printf("%s Initializing %s\n", lBrace, rBrace);

    clock_t tStart = clock();

    // Allocating memory and initializing all global arrays.
    Memory_Initialization_AtStart();
    Global_Array_Initialization_AtStart();

    // Reading in the structure file, and figuring out initial conditions.
    if (bReadConf_glb == 0)
        {
            Initial_Conditions_Simple();
        }
    else if (bReadConf_glb == 1)
        {
            Initial_Conditions_FromFile();
            if (nMCStepsForTherm_glb > 0)
                { // Remember that all thermalizing variants of
                  // MC moves cannot handle bonds.
                    Initial_Conditions_BreakBonds();
                }
        }
    else
        {
            printf("Wrong initial conditions condition given. Crashing!\n");
            exit(1);
        }

    // Performing a sanity check to see if all the beads and structures are
    // correct.


    puts("Performing initial sanity check.\n");
    PerformRuntimeSanityChecks(-1, -1);

    clock_t tEnd        = clock();
    double elapsed_time = (double) (tEnd - tStart) / (double) CLOCKS_PER_SEC;
    printf("Initialization done. %.2f sec elapsed.\n", elapsed_time);

    tStart = clock();

    long nGen;
    int nMCInfo; // MC accept/reject

    printf("%s Beginning MC Simulation %s\n", lBrace, rBrace);
    printf("******************************************\n\n");
    printf("_____________________\n");
    printf("Thermalizing system.\n");
    printf("---------------------\n\n");
    // Thermalizing the system.
    fCuTemp_glb = fPreKT_glb;
    // FILE Initialization
    FileIO_CreateRunningDataFiles();

    for (nGen = 0; nGen < nMCStepsForTherm_glb; nGen++)
        { // Intentionally not performing any data acquisition in the thermalization phase.
            nMCInfo = MC_Step_Therm(fCuTemp_glb);
            //        printf("(%d,%d)\n", nMCInfo / 12, nMCInfo % 2);
            DataPrinting_Thermalization(nGen);
        }
    /*
     * Post-thermalization
     * */
    FileIO_WriteRestart_ForThermalization();

    puts("____________________________");
    puts("System has been thermalized!");
    puts("----------------------------\n");

    if ((nAnnealing_Mode_glb == -1) && (nBiasPotential_CoupledToTemp_glb))
        {
            BiasPotential_TurnOFF();
        }

    /*
    The system has thermalized!
    */
    int run_cycle;
    // Going through the MC cycles.
    for (run_cycle = 0; run_cycle < nTotCycleNum_glb; run_cycle++)
        {
            /*
             * Pre run-cycle specific initialization.
             */
            fKT_glb = faKT_Cycle_glb[run_cycle];
            Calculate_Rot_Bias(fKT_glb);
            FileIOUtil_PreCycle_Init(run_cycle);
            for (nGen = 0; nGen < nMCStepsPerCycle_glb / 2; nGen++)
                {
                    fCuTemp_glb = nAnnealing_Mode_glb == -1 ? fKT_glb : Temperature_Function(nAnnealing_Mode_glb, nGen);
                    nMCInfo     = MC_Step(fCuTemp_glb);
                    //            printf("(%d,%d)\n", nMCInfo / 12, nMCInfo % 2);
                    DataPrinting_DuringRunCycles(nGen, run_cycle);
                }

            for (nGen = nMCStepsPerCycle_glb / 2; nGen <= nMCStepsPerCycle_glb; nGen++)
                {
                    fCuTemp_glb = nAnnealing_Mode_glb == -1 ? fKT_glb : Temperature_Function(nAnnealing_Mode_glb, nGen);
                    nMCInfo     = MC_Step(fCuTemp_glb);
                    //            printf("(%d,%d)\n", nMCInfo / 12, nMCInfo % 2);
                    DataPrinting_DuringRunCycles(nGen, run_cycle);
                    DataAnalysis_DuringRunCycles(nGen, run_cycle);
                }
            /*
             * Post run-cycle specific cleanup.
             */
            nAnnealing_Mode_glb        = -1;
            if ((! nBiasPotential_KeepON_glb) && (nBiasPotential_Mode_glb != -1))
                {
                    BiasPotential_TurnOFF();
                }

            FileIO_PostCycle_WriteSystemRestart(run_cycle);
            FileIO_PostCycle_WriteCycleAvgData(run_cycle);
            Reset_Global_Arrays();
        }

    tEnd         = clock();
    elapsed_time = (double) (tEnd - tStart) / (double) CLOCKS_PER_SEC;

    double elapsed_minutes = elapsed_time / 60.;
    double elapsed_hours   = elapsed_minutes / 60.;
    double elapsed_days    = elapsed_hours / 24.;

    printf("____________________________________\n");
    printf("Simulation finished! Timings:\n");
    printf("%.2lf mins.\n", elapsed_minutes);
    printf("%.2lf hours.\n", elapsed_hours);
    printf("%.2lf days.\n", elapsed_days);
    printf("------------------------------------\n");
    printf("%s ENDING %s\n", lBrace, rBrace);
    printf("******************************************\n");

    return 0;
}
