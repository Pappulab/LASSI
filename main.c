#include "global.h"
#include "parsekey.h"
#include "print.h"
#include "initialize.h"
#include "structure.h"
#include "energy.h"
#include "mcmove.h"
#include "cluster.h"

int main(int argc, char* argv[]) {
  char keyfile[100];
  char lBrace[] = "<======      ";
  char rBrace[] = "      ======>";
  if (argc == 3 && strcmp(argv[1], "-k") == 0) {
    strcpy(keyfile, argv[2]);
  } else {
    strcpy(keyfile, "param.key");
  }

  int errorcode;
  errorcode = parse_key(keyfile);
  if (errorcode == 0) {
    printf("Key file %s was successfully parsed.\n\n", keyfile);
    print_key();
  } else {
    printf("ERROR: unable to parse key file %s.\n", keyfile);
    exit(1);
  }

  printf("%s Initializing %s\n", lBrace, rBrace);

  srand(seed == 0 ? time(NULL) : seed);


  clock_t tStart = clock();

  MemoryAndBookkeepingInit();

  if (bReadConf == 0) {
    initialize_with_topo();
  } else if (bReadConf == -1) {
    printf("ERROR: no structure information.\n");
    exit(1);
  }

  if (check_structure_topo() == 0) {
    printf("Check structure sanity: OK\n");
  } else {
    printf("ERROR: wrong structure.\n");
    exit(1);
  }

  clock_t tEnd = clock();
  double elapsed_time = (tEnd-tStart)/(double)CLOCKS_PER_SEC;
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
  Print_Data(-1,-1);//Initialization of files
  for (nGen = 0; nGen < nPreSteps; nGen++) {//Intentionally not performing any data acquisition in the thermalizing phase
      nMCInfo = MC_Step_Equil(fCuTemp);
      Print_Data(nGen,-1);
  }

printf("____________________________\n");
printf("System has been thermalized!\n");
printf("----------------------------\n\n");
/*
The system has thermalized!
*/
  int run_cycle;
  float temp_cycle[10];
  for(run_cycle=0; run_cycle<10; run_cycle++){
      temp_cycle[run_cycle] = fKT*((float)run_cycle+0.1);
  }


  for(run_cycle = 0; run_cycle < 10; run_cycle++){
    fKT = temp_cycle[run_cycle];
    Calculate_Rot_Bias(temp_cycle[run_cycle]);
    Print_Data(-1, run_cycle);
  for(nGen = 0; nGen <= nSteps; nGen++){
    fCuTemp = Temperature_Function(Temp_Mode, nGen);
    nMCInfo = MC_Step(fCuTemp);
    Print_Data(nGen, run_cycle);
  }
  Temp_Mode = -1;
  write_SysProp(fileSysProp);
  Reset_Global_Arrays();
  }


  tEnd = clock();
  elapsed_time = (tEnd-tStart)/(double)CLOCKS_PER_SEC;
  printf("____________________________________\n");
  printf("Simulation finished in %.2f minutes.\n", elapsed_time/60.);
  printf("------------------------------------\n");
  printf("%s ENDING %s\n", lBrace, rBrace);
  printf("******************************************\n");

  return 0;
}
