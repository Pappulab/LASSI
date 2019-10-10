#include "global.h"
#include "print.h"
#include "cluster.h"
#include "structure.h"
#include "energy.h"

void print_matrix(char* strTitle,  int nSeqEn, float fArray[MAX_AA][MAX_AA][MAX_E], int param);

void write_network(char* filename, long nGen, int nT1, int nT2, float local_conc, float local_beadconc) {
  FILE *fp;

  if (nGen == -1) {
    fp = fopen(filename, "w"); // overwrite
  } else {
    fp = fopen(filename, "a");
  }

  if (nGen == -1) { // title
    fprintf(fp, "#step\tcont\tspec\tc_chain_loc\tc_chain_glo\tc_bead_local\tc_bead_global\n");
  } else {
    fprintf(fp, "%ld\t%.2f\t%.2f\t%.1e\t%.1e\t%.1e\t%.1e\n", nGen, (float)nT1/(float)tot_chains, (float)nT2 / (float)tot_chains, local_conc, (float)tot_chains / (float)nBoxSize[0] / (float)nBoxSize[1] / (float)nBoxSize[2], local_beadconc, (float)tot_beads / (float)nBoxSize[0] / (float)nBoxSize[1] / (float)nBoxSize[2]);
  }

  fclose(fp);
}

void write_cluster(char* filename, long nGen) {
  FILE *fp;
    int i;
  if (nGen == -1) {
    fp = fopen(filename, "w"); // overwrite
  } else {
    fp = fopen(filename, "a");
  }

  if (nGen == -1) { // title
    fprintf(fp, "#Step followed by histogram\n");
  } else {
    fprintf(fp, "#%ld\n",nGen);
    for(i = 0; i <= tot_chains; i++){
        fprintf(fp, "%ld\t", naClusHistList[i]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

void write_GyrTen(char* filename, long nGen) {
    FILE *fp;
    int i;
  if (nGen == -1) {
    fp = fopen(filename, "w"); // overwrite
  } else {
    fp = fopen(filename, "a");
  }

  if (nGen == -1) { // title
    fprintf(fp, "#Step followed by Gyration Tensor Vals: xx yy zz xy yz xz\n");
  } else {
    fprintf(fp, "#%ld\n",nGen);
    for(i = 0; i < 7; i++){
        if(i==5){continue;}//Not smart enough to make a better indexing function
        fprintf(fp, "%f\t", fGyrTensor[i]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

void write_mcmove(char* filename, long nGen, float fMCTemp) {
  FILE *fp;
  int i;//Iterator
  if (nGen == -1) {
    fp = fopen(filename, "w"); // overwrite
  } else {
    fp = fopen(filename, "a");
  }
  if(nGen == -1){
            fprintf(fp, "#Steps ,and Temp, are followed by (rej, acc) for each MC Move.\n");
            fprintf(fp, "#(nStep, T) TRANS CLSTER SM_CLSTER FACE LOCAL DBPVT PIVOT BR_ROT\n");
            for(i=0;i<MAX_MV;i++){
                MCAccepMat[0][i]=0;
                MCAccepMat[1][i]=0;
              //This way the print function will initialize the matrix at startup!
            }
      } else{
        printf("Acceptance Ratios:\n");
        fprintf(fp, "%ld\t%.2f\t", nGen, fMCTemp);//Step and Temp
      for(i=1;i<MAX_MV;i++){
        fprintf(fp, "%ld\t%ld\t", MCAccepMat[0][i], MCAccepMat[1][i]);
            printf("%.3f\t", 100.*(float)MCAccepMat[1][i]/((float)MCAccepMat[0][i]+1.+(float)
            MCAccepMat[1][i]));
          MCAccepMat[0][i]=0;
          MCAccepMat[1][i]=0;
        //This way the print function will zero out the matrix every time we print to a file!
      }
        fprintf(fp,"\n");
        printf("\n");
      }


  fclose(fp);
}

void print_log(long nGen, int run_it) {
  printf("Step       %.2e\n", (float)nGen);
  printf("Run Cycle: %d\n", run_it);
  //printf("MC Temp = %.3e;\tRot Bias Prob = %.3e /site;\n", fCuTemp, fRot_Bias);
  printf("MC Temp  = %.3e;\n", fCuTemp);
  printf("Total E  = %.3e;\tIso E = %.3e;\tAniso E = %.3e;\n", faCurrEn[E_TOT], faCurrEn[E_OVLP], faCurrEn[E_SC_SC]);
  printf("Percolation Parameter Is: %.3f\n", ((float)nLargestClusterRightNow) / ((float)nTotClusCounter + 0.0001) / (float)tot_chains);
  int i, j;
    printf("Acceptance Ratios:\n");
    for(i=1;i<MAX_MV;i++){
            printf("%.2f\t",  100. * (float) MCAccepMat[1][i] / ((float) MCAccepMat[0][i] + 0.00001 + (float)
                    MCAccepMat[1][i]));
    }
  printf("\n\n");
}

void write_energy(char* filename, long nGen) {
  FILE *fp;
  if (nGen == -1) {
    fp = fopen(filename, "w"); // overwrite
  } else {
    fp = fopen(filename, "a");
  }

  int i;
  if (nGen == -1) { // title
    fprintf(fp, "#step\ttotal\toverlap\tcontact\tSC-SC\tstiff\n");

  } else {
    fprintf(fp, "%ld", nGen);
    for (i=0; i<(MAX_E); i++) {
      fprintf(fp, "\t%.1f", faCurrEn[i]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

void write_trajectory(char* filename, long nGen){
//Writes the trajectory in LAMMPS format. To be viewed with VMD (Hopefully). Read the LAMMPS dump documentation for
//the actual formate of the file
  FILE *fp;
  if (nGen == -1){
    fp = fopen(filename, "w"); //This always overwrites a previos file
  }
  else{
    fp = fopen(filename, "a"); //We will append to the file made for this run

  int i;//Looping index
  fprintf(fp, "ITEM: TIMESTEP\n");
  fprintf(fp, "%ld\n", nGen);//Timestep

  fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
  fprintf(fp, "%d\n", tot_beads);//Total atom number

  fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");//BCs are always periodic for now
  fprintf(fp, "0 %d\n 0 %d\n 0 %d\n", nBoxSize[0], nBoxSize[1], nBoxSize[2]);//Box dimensions

  fprintf(fp, "ITEM: ATOMS id type mol x y z bP\n");//What we are printing

  for (i=0; i < tot_beads; i++){
    fprintf(fp, "%d %d %d %d %d %d %d\n", i, bead_info[i][BEAD_TYPE], bead_info[i][BEAD_CHAINID],
                              bead_info[i][POS_X], bead_info[i][POS_Y], bead_info[i][POS_Z],
                              bead_info[i][BEAD_FACE]);
                            }

  }

  fclose(fp);
}

void print_key(void) { // should be output-dependent (stdout, stderr, other files)

  int i;
  char lBrace[] = "<======      ";
  char rBrace[] = "      ======>";
  printf("%s System Settings %s\n", lBrace, rBrace);
  printf("Number of Components = %d\n", nSeqEn);
  printf("Number of Beads      = %d\n", tot_beads);
  printf("Number of Chains     = %d\n", tot_chains);
  printf("Box Size             = %d, %d, %d\n", nBoxSize[0], nBoxSize[1], nBoxSize[2]);
  printf("Monomer density      = %1.1e\n", (float)tot_beads / (float)nBoxSize[0] / (float)nBoxSize[1] / (float)nBoxSize[2]);
  printf("\n");

  printf("%s Energy Matrices %s\n", lBrace, rBrace);
  print_matrix("E_ovlp", nSeqEn, fEnergy, E_OVLP);
  print_matrix("E_cont", nSeqEn, fEnergy, E_CONT);
  print_matrix("R_cont", nSeqEn, fEnRad, E_CONT);
  print_matrix("SC_SC",  nSeqEn, fEnergy, E_SC_SC);
  printf("\n");

  printf("%s Linker Info %s\n", lBrace, rBrace);
  printf("Linker length             = %.2f\n", fLinkerLength);
  printf("Linker spring constant    = %.2f\n", fLinkerSprCon);
  printf("Linker equilibrium length = %.2f\n", fLinkerEqLen);
  printf("\n");

  printf("%s MC Setup %s\n", lBrace, rBrace);
  printf("MC Temperatures: (First, Last) = (%.2f, %.2f)\n", fKT, fKT+(float)(nTot_CycleNum-1)*fdelta_temp);
  printf("Temperature Mode               = %d\n", Temp_Mode);
  printf("Indent Mode                    = %d\n", Indent_Mode);
  printf("Rotational Bias Mode           = %d\n", RotBias_Mode);
  printf("Number of MC Cycles            = %e\n", (float)nTot_CycleNum);
  printf("Number of MC Steps/Cycle       = %e\n", (float)nSteps);
  printf("Thermalizing Temperature       = %.2f\n", fPreKT);
  printf("Number of Thermalizing Steps   = %e\n", (float)nPreSteps);
  printf("RNG Seed                       = %d\n", seed);
  char *MoveName[MAX_MV];
  MoveName[MV_PIVOT]   = "Pivot           ";
  MoveName[MV_DBPVT]   = "Double Pivot    ";
  MoveName[MV_CLSTR]   = "Larger Cluster  ";
  MoveName[MV_SMCLSTR] = "Smaller Cluster ";
  MoveName[MV_STROT]   = "Face Change     ";
  MoveName[MV_LOCAL]   = "Local           ";
  MoveName[MV_COLOCAL] = "Co-local        ";
  MoveName[MV_MTLOCAL]   = "Shake           ";
  MoveName[MV_BRROT]   = "Rotate Branched ";
  MoveName[MV_SNAKE]   = "Slithering Snake";
  MoveName[MV_TRANS]   = "Translation     ";
  float freqMin = 1e10;
  for (i=MV_NULL+1; i<MAX_MV; i++) {
      if(freqMin >= fMCFreq[i] && fMCFreq[i] != 0.){
        freqMin = fMCFreq[i];
      }
  }
  printf("MC Move Frequencies:\n");
  printf("--------------------------\n");
  for (i=MV_NULL+1; i<MAX_MV; i++) {
    printf("|%s  %d\n", MoveName[i], (int)ceilf(fMCFreq[i]/freqMin));
  }
  printf("--------------------------\n");
  printf("\n");
}

void print_matrix(char* strTitle, int nSeqEn, float fArray[MAX_AA][MAX_AA][MAX_E], int param) {
  int nLen;
  nLen = nSeqEn;
  int i, j;
  for (i=0; i<nLen+1; i++) {
    printf("--------");
  }
  printf("\n|%s", strTitle);
  for (i=0; i<nLen; i++) {
    printf("\t%d", i);
  }
  printf("      |\n");

  for (i=0; i<nLen; i++) {
    printf("|%d", i);
    for (j=0; j<nLen; j++) {
      printf("\t%.2f", fArray[i][j][param]);
    }
    printf("   |\n");
  }
  for (i=0; i<nLen+1; i++) {
    printf("--------");
  }
  printf("\n");
}

void write_rdftot(char* filename, long nGen){
  FILE *fp;
  if (nGen == -1) {
    fp = fopen(filename, "w"); // overwrite
  } else {
    fp = fopen(filename, "a");
  }
  int i;
  if (nGen == -1) { // title
    fprintf(fp, "#Row-by-row RDF. dr = 1\n");
  }
    else{
      for (i = 0; i < nBins_RDF;i++){
        fprintf(fp, "%.5f\t", fRDF_TOT[i]);
      }
      fprintf(fp, "\n");
    }

    fclose(fp);

}

void write_rdfsplit(char* filename, long nGen){
  FILE *fp;
  if (nGen == -1) {
    fp = fopen(filename, "w"); // overwrite
  } else {
    fp = fopen(filename, "a");
  }
  int i, k;
  if (nGen == -1) { // title
    fprintf(fp, "#Row-by-row split RDF. dr = 1/4\n");
  }
    else{
      for (k=0; k<1; k++){
        for (i = 0; i < nBins_RDF;i++){
          fprintf(fp, "%.5Lf\t", ldRDF_ARR[k][i]);
        }
        fprintf(fp, "\n");
      }
    }

    fclose(fp);

}

void write_topofile(char* filename){
      /*
      Writes a topology file for VMD. Since the trajectory is saved in the LAMMPS format,
      the topology file is also in the LAMMPS format. The format for this file is one which
      is used as input data for LAMMPS.
      This function copies the format given in read_data for LAMMPS: https://lammps.sandia.gov/doc/read_data.html
      */
  FILE *fp;
  fp = fopen(filename, "w");//Just overwrite over if last file exists.
  int i, j, k;//Loop iterators.
  int numBonds;//Used to count total number of bonds!
  printf("Writing the topology file!\n");
  fprintf(fp, "LAMMPS Description\n");//The file must start with this.
  fprintf(fp, "\n");//Empty line.
  fprintf(fp, "%d\tatoms\n", tot_beads);//Listing total number of atoms
  numBonds = 0;
  for(i=0; i < tot_beads; i++){
    for(j=0; j<MAX_BONDS; j++){
      if(topo_info[i][j] != -1){//There is a bond between i and topo_info[i][j]
        numBonds++;
      }
    }
  }

  fprintf(fp, "%d\tbonds\n", numBonds);//Listing number of bonds
  fprintf(fp, "0\tangles\n");//Systems don't have angle depenece yet
  fprintf(fp, "0\tdihedrals\n");//Systems don't have dihedrals  yet
  fprintf(fp, "0\timpropers\n");//Systems don't have imporopers yet
  fprintf(fp, "\n");//Empty line.
  fprintf(fp, "6\tatom types\n");//Just hard coded to have 6 sticker types for now
  fprintf(fp, "1\tbond types\n");//System can have multiple bond-lengths
  fprintf(fp, "0\tangle types\n");//Systems don't have any angular forces yet
  fprintf(fp, "\n");//Empty line.
  fprintf(fp, "0 %d xlo xhi\n", nBoxSize[0]);
  fprintf(fp, "0 %d ylo yhi\n", nBoxSize[1]);
  fprintf(fp, "0 %d zlo zhi\n", nBoxSize[2]);
  fprintf(fp, "\n");//Empty line.
  fprintf(fp, "Masses\n");//These don't really mean anything
  fprintf(fp, "\n");//Empty line.
  fprintf(fp, "0\t1.0\n1\t1.0\n2\t1.0\t3\t1.0\n4\t1.0\n5\t1.0\n");//Dummy masses for the 6 types.
  fprintf(fp, "\n");//Empty line.
  fprintf(fp, "Atoms\n");//Signifying the beginning of atom coordinates.
  fprintf(fp, "\n");//Empty line.
  for(i=0; i < tot_beads; i++){
    fprintf(fp, "%d %d %d 0.0 %d %d %d\n", i, bead_info[i][BEAD_CHAINID], bead_info[i][BEAD_TYPE], bead_info[i][POS_X], bead_info[i][POS_Y], bead_info[i][POS_Z]);
  }//Done with the coordinates
  fprintf(fp, "\n");//Empty line.
  fprintf(fp, "Bonds\n");//Signifying the beginning of atom coordinates.
  fprintf(fp, "\n");//Empty line.
  k=0;//This guy counts bondIDs
  for(i=0; i < tot_beads; i++){
    for(j=0; j<MAX_BONDS; j++){
      if(topo_info[i][j] != -1){//There is a bond between i and topo_info[i][j]
        fprintf(fp, "%d %d %d %d\n", k, linker_len[i][j] / (int)fLinkerLength, i, topo_info[i][j]);
        k++;
      }
    }
  }//Done with the coordinates
  fprintf(fp, "\n");//Empty line.

  fclose(fp);
}

void Write_SysProp(char* filename){

  FILE *fp;
  int i, j;
  fp = fopen(filename, "w");
  fprintf(fp, "#Contains various averaged quatities.\n");
  //Gyration Radius
  fprintf(fp, "#Total Rg\n%f\t%f\n#Cluster Hist\n", fSysGyrRad / (float)nTotGyrRadCounter, (float)nBoxSize[0] / 2.);
  //Cluster Histogram
  fprintf(fp, "%f\t", (float)nLargestClusterRightNow / (float)nTotClusCounter );
  for(i=1; i <= tot_chains; i++){
    fprintf(fp, "%f\t", (float)naClusHistList[i] / (float)nTotClusCounter);
  }
  //Split RDFs
  fprintf(fp, "\n#Split RDFs. ALL-ALL; DIAGONALS and then from 0 onwards \n");
  for (j=0; j<RDF_COMPS; j++){
    for (i=0; i<nBins_RDF; i++){
      fprintf(fp, "%LE\t", ldRDF_ARR[j][i] / (float)nRDFCounter);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n#Done");

  fclose(fp);

}

void Write_TotalSysProp(char* filename, int run_it){
    /* This function writes one large file with all the averaged values from each run_cycle. run_it let's the function
     * know how many cycles to write.
     * As opposed to Write_SysProp, each of the relevant averaged parameters will be in their appropriately named files.
     * The naming convention shall be: filename_*.dat where * will be GR, CLUS and RDF for the different measured
     * quantities. Within each file, each run_cycle shall be written one-by-one in the style of Write_SysProp
     */
    FILE *fp;
    int i,j,k;
    sprintf(filename, "%s_RDF.dat", strReportPrefix);//Name Of the RDF Files
    fp = fopen(filename, "w");
    fprintf(fp, "#Split RDFs. ALL-ALL; DIAGONALS and then from 0-0, 0-1, and onwards \n");
    for (i = 0; i < run_it; i++){
        fprintf(fp, "#Run_Cycle = %d\n", i);
        for (j = 0; j < RDF_COMPS ; j++){
            for (k = 0; k < nBins_RDF; k++){
                fprintf(fp, "%LE\t", ld_TOTRDF_ARR[i][j][k]);
            }
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "#Done");
    fclose(fp);

    sprintf(filename, "%s_CLUS.dat", strReportPrefix);//Name Of the ClusterHistogram Files
    fp = fopen(filename, "w");
    fprintf(fp, "#Cluster Histograms: 1st column is largest cluster, and then clusters of size 1, 2, and so on\n");
    for (i = 0; i < run_it; i++){
        fprintf(fp, "#Run_Cycle = %d\n", i);
            for (k = 0; k <= tot_chains; k++){
                fprintf(fp, "%LE\t", ld_TOTCLUS_ARR[i][k]);
            }
            fprintf(fp, "\n");
    }
    fprintf(fp, "#Done");
    fclose(fp);


    sprintf(filename, "%s_GR.dat", strReportPrefix);//Name Of the ClusterHistogram Files
    fp = fopen(filename, "w");
    fprintf(fp, "#Cluster Histograms: 1st column is largest cluster, and then clusters of size 1, 2, and so on\n");
    for (i = 0; i < run_it; i++){
        fprintf(fp, "#Run_Cycle = %d\n", i);
        fprintf(fp, "%LE\t", ld_TOTGYRRAD_ARR[i][0]);
        fprintf(fp, "%LE\n", ld_TOTGYRRAD_ARR[i][1]);
    }
    fprintf(fp, "#Done");
    fclose(fp);

}

void Print_Data(long nGen, int run_it){
    //This function handles all the data IO.
    int nFlagForEnCalc = 0; //Flag for total energy calculation
    //run_it == -1 corresponds to the thermalization sequence
    if (run_it ==-1){
        //Open the appropriate files before the thermalization sequence.
        if(nGen == -1) {
            if (nReport[REPORT_CONFIG] != 0) {
                sprintf(fileStruct, "%s_topo.lammpstrj", strReportPrefix);//Name of the topology file
                write_topofile(fileStruct);//Write the topology file. Only need to write once
                sprintf(fileStruct, "%s_trj.lammpstrj", strReportPrefix);//Naming convention for trajectory files.
                write_trajectory(fileStruct, -1);//This opens a new trajectory file; each run_it will have its own
            }
            if (nReport[REPORT_ENERGY] != 0) {
                sprintf(fileEnergy, "%s_energy.dat", strReportPrefix);
                write_energy(fileEnergy, -1);//Open a new energy file; each run_it will have its own
            }
            if (nReport[REPORT_MCMOVE] != 0) {
                sprintf(fileMCMove, "%s_mcmove.dat", strReportPrefix);
                write_mcmove(fileMCMove, -1, 0.0);//Open a new MCInfo file; each run_it will have its own
            }
            if (nReport[REPORT_NETWORK] != 0 || nReport[REPORT_RDFTOT] != 0) {
                sprintf(fileSysProp, "%s_SysProp.dat", strReportPrefix);//Naming convention or SysProp files; each run_it should have its own
            }
        }
        else{
            if (nReport[REPORT_LOG] != 0){
                if (nGen % nReport[REPORT_LOG] == 0){
                    if(Check_System_Structure() == 0){
                        printf("Structure is still a-okay!\n");
                    }
                    else{
                        printf("Molecular structure is inconsistent with initial structure.\nCRASHING\n\n");
                        exit(1);
                    }
                    Calc_Tot_Energy(); nFlagForEnCalc = 1;
                    print_log(nGen, run_it);
                }
            }
            if (nReport[REPORT_CONFIG] != 0){
                if (nGen % nReport[REPORT_CONFIG] == 0){
                    write_trajectory(fileStruct, nGen);
                }
            }
            if (nReport[REPORT_ENERGY] != 0){
                if (nGen % nReport[REPORT_ENERGY] == 0){
                    if(nFlagForEnCalc != 1){//Calculate the energy
                        Calc_Tot_Energy(); nFlagForEnCalc = 1;
                    }
                    write_energy(fileEnergy, nGen);
                }
            }
            if (nReport[REPORT_MCMOVE] != 0){
                if (nGen % nReport[REPORT_MCMOVE] == 0){
                    write_mcmove(fileMCMove, nGen, fCuTemp);
                }
            }
        }

    }
    if (run_it == 0 && nGen > 0){
        if (nReport[REPORT_LOG] != 0){
            if (nGen % nReport[REPORT_LOG] == 0){
                if(Check_System_Structure() == 0){
                    printf("Structure is still a-okay!\n");
                }
                else{
                    printf("Molecular structure is inconsistent with initial structure.\nCRASHING\n\n");
                    exit(1);
                }
                Calc_Tot_Energy(); nFlagForEnCalc = 1;
                print_log(nGen, run_it);
            }
        }
        if (nReport[REPORT_CONFIG] != 0){
            if (nGen % nReport[REPORT_CONFIG] == 0){
                write_trajectory(fileStruct, nGen+nPreSteps);
            }
        }
        if (nReport[REPORT_ENERGY] != 0){
            if (nGen % nReport[REPORT_ENERGY] == 0){
                if(nFlagForEnCalc != 1){//Calculate the energy
                    Calc_Tot_Energy(); nFlagForEnCalc = 1;
                }
                write_energy(fileEnergy, nGen+nPreSteps);
            }
        }
        if (nReport[REPORT_MCMOVE] != 0){
            if (nGen % nReport[REPORT_MCMOVE] == 0){
                write_mcmove(fileMCMove, nGen+nPreSteps, fCuTemp);
            }
        }
        if (nReport[REPORT_RDFTOT] != 0){//SysProp is printed outside of this function in main.c, lol
            if (nGen % nReport[REPORT_RDFTOT] == 0 && nGen > nSteps/2){
                avg_rdf_split();
            }
        }
        if (nReport[REPORT_NETWORK] != 0){//SysProp is printed outside of this function in main.c, lol
            if (nGen % nReport[REPORT_NETWORK] == 0 && nGen > nSteps/2){
                avg_clus_dist(naList);
                CalcTotGyrRad();
            }
        }
    }

    if (run_it > 0){
        if(nGen == -1) {
            if (nReport[REPORT_CONFIG] != 0) {
                sprintf(fileStruct, "%s_%d_trj.lammpstrj", strReportPrefix, run_it);//Naming convention for trajectory files.
                write_trajectory(fileStruct, -1);//This opens a new trajectory file; each run_it will have its own
                write_trajectory(fileStruct, 0);//End of previous run_it is initial conditions for this one!
                sprintf(fileStruct, "%s_trj.lammpstrj", strReportPrefix);//Returning back
            }
            if (nReport[REPORT_ENERGY] != 0) {
                sprintf(fileEnergy, "%s_%d_energy.dat", strReportPrefix, run_it);
                write_energy(fileEnergy, -1);//Open a new energy file; each run_it will have its own
                write_energy(fileEnergy, 0);//Last frame from previous run_it is first frame of this one
            }
            if (nReport[REPORT_MCMOVE] != 0) {
                sprintf(fileMCMove, "%s_%d_mcmove.dat", strReportPrefix, run_it);
                write_mcmove(fileMCMove, -1, 0.0);//Open a new MCInfo file; each run_it will have its own
            }
            if (nReport[REPORT_NETWORK] != 0 || nReport[REPORT_RDFTOT] != 0) {
                sprintf(fileSysProp, "%s_%d_SysProp.dat", strReportPrefix, run_it);//Naming convention or SysProp files; each run_it should have its own
            }
        }
        else{
            if (nReport[REPORT_LOG] != 0){
                if (nGen % nReport[REPORT_LOG] == 0){
                    if(Check_System_Structure() == 0){
                        printf("Structure is still a-okay!\n");
                    }
                    else{
                        printf("Molecular structure is inconsistent with initial structure.\nCRASHING\n\n");
                        exit(1);
                    }
                    Calc_Tot_Energy(); nFlagForEnCalc = 1;
                    print_log(nGen, run_it);
                }
            }
            if (nReport[REPORT_CONFIG] != 0){
                if (nGen % nReport[REPORT_CONFIG] == 0){
                    write_trajectory(fileStruct, nGen+(run_it*nSteps));
                }
            }
            if (nReport[REPORT_ENERGY] != 0){
                if (nGen % nReport[REPORT_ENERGY] == 0){
                    if(nFlagForEnCalc != 1){//Calculate the energy
                        Calc_Tot_Energy(); nFlagForEnCalc = 1;
                    }
                    write_energy(fileEnergy, nGen);
                }
            }
            if (nReport[REPORT_MCMOVE] != 0){
                if (nGen % nReport[REPORT_MCMOVE] == 0){
                    write_mcmove(fileMCMove, nGen, fCuTemp);
                }
            }
            if (nReport[REPORT_RDFTOT] != 0){//SysProp is printed outside of this function in main.c, lol
                if (nGen % nReport[REPORT_RDFTOT] == 0 && nGen > nSteps/2){
                    avg_rdf_split();
                }
            }
            if (nReport[REPORT_NETWORK] != 0){//SysProp is printed outside of this function in main.c, lol
                if (nGen % nReport[REPORT_NETWORK] == 0 && nGen > nSteps/2){
                    avg_clus_dist(naList);
                    CalcTotGyrRad();
                }
            }
        }
    }

}

void Copy_Data(int run_it){
    int i, j;
    for(i = 0; i < RDF_COMPS; i++){
        for (j = 0; j < RDF_MAXBINS; j++) {
            ld_TOTRDF_ARR[run_it][i][j] = ldRDF_ARR[i][j] / (long double)nRDFCounter;
        }
    }
    ld_TOTCLUS_ARR[run_it][0] = (long double)nLargestClusterRightNow / (long double)nTotClusCounter;
    for(i=1; i<=tot_chains; i++){
        ld_TOTCLUS_ARR[run_it][i] = (long double)naClusHistList[i] / (long double)nTotClusCounter;
    }
    ld_TOTGYRRAD_ARR[run_it][0] = (long double)fSysGyrRad / (long double)nTotGyrRadCounter;
    ld_TOTGYRRAD_ARR[run_it][1] = (long double)nBoxSize[0]/2.;
}