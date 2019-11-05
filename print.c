#include "global.h"
#include "print.h"
#include "cluster.h"
#include "structure.h"
#include "energy.h"

void Print_Matrix(char *strTitle, int nSeqEn, float fArray[MAX_AA][MAX_AA][MAX_E], int param);

/// Write_ClusterDist - write the cluster histogram to a separate file. NOT USED ANYMORE, but might find a use later
/// Always appends to the file for this run. Stopped using it because the IO load was slowing things down at our cluster
/// Would be a good way to gather proper statistics on clusters as the runs went on
/// \param filename
/// \param nGen
void Write_ClusterDist(char *filename, long nGen) {
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
        fprintf(fp, "#%ld\n", nGen);
        for (i = 0; i <= tot_chains; i++) {
            fprintf(fp, "%ld\t", naClusHistList[i]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/// Write_GyrTen - write the total gyration tensor to a file. Remember that we only have 6 things in a symmetric 3x3 tensor.
/// Always appends to the file for this run. Stopped using it because the IO load was slowing things down at our cluster
/// Would be a good way to gather proper statistics on the shapes and such as the runs went on
/// \param filename
/// \param nGen
void Write_GyrTen(char *filename, long nGen) {
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
        fprintf(fp, "#%ld\n", nGen);
        for (i = 0; i < 7; i++) {
            if (i == 5) { continue; }//Not smart enough to make a better indexing function
            fprintf(fp, "%f\t", fGyrTensor[i]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/// Write_MCMove - writes the acceptance and rejection ratios of the various moves.
/// Just keeps appending to that file. Can be used to track the 'dynamics' during a simulation.
/// \param filename
/// \param nGen
/// \param fMCTemp
void Write_MCMove(char *filename, long nGen, float fMCTemp) {
    FILE *fp;
    int i;//Iterator
    if (nGen == -1) {
        fp = fopen(filename, "w"); // overwrite
    } else {
        fp = fopen(filename, "a");
    }
    if (nGen == -1) {
        fprintf(fp, "#Steps ,and Temp, are followed by (rej, acc) for each MC Move.\n");
        fprintf(fp, "#(nStep, T) TRANS CLSTER SM_CLSTER FACE LOCAL DBPVT PIVOT BR_ROT\n");
        for (i = 0; i < MAX_MV; i++) {
            MCAccepMat[0][i] = 0;
            MCAccepMat[1][i] = 0;
            //This way the print function will initialize the matrix at startup!
        }
    } else {
        fprintf(fp, "%ld\t%.2f\t", nGen, fMCTemp);//Step and Temp
        for (i = 1; i < MAX_MV; i++) {
            fprintf(fp, "%ld\t%ld\t", MCAccepMat[0][i], MCAccepMat[1][i]);
            MCAccepMat[0][i] = 0;
            MCAccepMat[1][i] = 0;
            //This way the print function will zero out the matrix every time we print to a file!
        }
        fprintf(fp, "\n");
    }


    fclose(fp);
}

/// Print_LogToScreen - prints current status of the run. Overall move acceptance ratios, and energies.
/// \param nGen
/// \param run_it
void Print_LogToScreen(long nGen, int run_it) {
    printf("Step       %.2e\n", (float) nGen);
    printf("Run Cycle: %d\n", run_it);
    //printf("MC Temp = %.3e;\tRot Bias Prob = %.3e /site;\n", fCuTemp, fRot_Bias);
    printf("MC Temp  = %.3e;\n", fCuTemp);
    printf("Total E  = %.3e;\tIso E = %.3e;\tAniso E = %.3e;\n", faCurrEn[E_TOT], faCurrEn[E_OVLP], faCurrEn[E_SC_SC]);
    printf("Percolation Parameter Is: %.3f\n",
           ((float) nLargestClusterRightNow) / ((float) nTotClusCounter + 0.0001) / (float) tot_chains);
    int i, j;
    printf("Acceptance Ratios:\n");
    for (i = 1; i < MAX_MV; i++) {
        printf("%.2f\t", 100. * (float) MCAccepMat[1][i] / ((float) MCAccepMat[0][i] + 0.00001 + (float)
                MCAccepMat[1][i]));
    }
    printf("\n\n");
}

/// Write_Energy - write the decomposed energy of the system to a file.
/// Just appends to the file for this run.
/// \param filename
/// \param nGen
void Write_Energy(char *filename, long nGen) {
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
        for (i = 0; i < (MAX_E); i++) {
            fprintf(fp, "\t%.1f", faCurrEn[i]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/// Write_Trajectory -  write a LAMMPS style trajectory file for the system.
/// Appends to the file for this run.
/// \param filename
/// \param nGen
void Write_Trajectory(char *filename, long nGen) {
//Writes the trajectory in LAMMPS format. To be viewed with VMD (Hopefully). Read the LAMMPS dump documentation for
//the actual formate of the file
    FILE *fp;
    if (nGen == -1) {
        fp = fopen(filename, "w"); //This always overwrites a previos file
    } else {
        fp = fopen(filename, "a"); //We will append to the file made for this run

        int i;//Looping index
        fprintf(fp, "ITEM: TIMESTEP\n");
        fprintf(fp, "%ld\n", nGen);//Timestep

        fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
        fprintf(fp, "%d\n", tot_beads);//Total atom number

        fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");//BCs are always periodic for now
        fprintf(fp, "0 %d\n0 %d\n0 %d\n", nBoxSize[0], nBoxSize[1], nBoxSize[2]);//Box dimensions

        fprintf(fp, "ITEM: ATOMS id type mol x y z bP\n");//What we are printing

        for (i = 0; i < tot_beads; i++) {
            fprintf(fp, "%d %d %d %d %d %d %d\n", i, bead_info[i][BEAD_TYPE], bead_info[i][BEAD_CHAINID],
                    bead_info[i][POS_X], bead_info[i][POS_Y], bead_info[i][POS_Z],
                    bead_info[i][BEAD_FACE]);
        }

    }

    fclose(fp);
}

/// Print_Key - print the keyfile that was read in to the screen
void Print_Key(void) { // should be output-dependent (stdout, stderr, other files)

    int i;
    char lBrace[] = "<======      ";
    char rBrace[] = "      ======>";
    printf("%s System Settings %s\n", lBrace, rBrace);
    printf("Number of Components = %d\n", nBeadTypes);
    printf("Number of Beads      = %d\n", tot_beads);
    printf("Number of Chains     = %d\n", tot_chains);
    printf("Box Size             = %d, %d, %d\n", nBoxSize[0], nBoxSize[1], nBoxSize[2]);
    printf("Monomer density      = %1.1e\n",
           (float) tot_beads / (float) nBoxSize[0] / (float) nBoxSize[1] / (float) nBoxSize[2]);
    printf("\n");

    printf("%s Energy Matrices %s\n", lBrace, rBrace);
    Print_Matrix("E_ovlp", nBeadTypes, fEnergy, E_OVLP);
    Print_Matrix("E_cont", nBeadTypes, fEnergy, E_CONT);
    Print_Matrix("R_cont", nBeadTypes, fEnRad, E_CONT);
    Print_Matrix("SC_SC", nBeadTypes, fEnergy, E_SC_SC);
    printf("\n");

    printf("%s Linker Info %s\n", lBrace, rBrace);
    printf("Linker length             = %.2f\n", fLinkerLength);
    printf("Linker spring constant    = %.2f\n", fLinkerSprCon);
    printf("Linker equilibrium length = %.2f\n", fLinkerEqLen);
    printf("\n");

    printf("%s MC Setup %s\n", lBrace, rBrace);
    printf("MC Temperatures: (First, Last) = (%.2f, %.2f)\n", fKT, fKT + (float) (nTot_CycleNum - 1) * fdelta_temp);
    printf("Temperature Mode               = %d\n", Temp_Mode);
    printf("Indent Mode                    = %d\n", nThermalization_Mode);
    printf("Rotational Bias Mode           = %d\n", RotBias_Mode);
    printf("Number of MC Cycles            = %e\n", (float) nTot_CycleNum);
    printf("Number of MC Steps/Cycle       = %e\n", (float) nSteps);
    printf("Thermalizing Temperature       = %.2f\n", fPreKT);
    printf("Number of Thermalizing Steps   = %e\n", (float) nPreSteps);
    printf("RNG Seed                       = %d\n", RNG_Seed);
    char *MoveName[MAX_MV];
    MoveName[MV_PIVOT] = "Pivot           ";
    MoveName[MV_DBPVT] = "Double Pivot    ";
    MoveName[MV_CLSTR] = "Larger Cluster  ";
    MoveName[MV_SMCLSTR] = "Smaller Cluster ";
    MoveName[MV_STROT] = "Face Change     ";
    MoveName[MV_LOCAL] = "Local           ";
    MoveName[MV_COLOCAL] = "Co-local        ";
    MoveName[MV_MTLOCAL] = "Shake           ";
    MoveName[MV_BRROT] = "Rotate Branched ";
    MoveName[MV_SNAKE] = "Slithering Snake";
    MoveName[MV_TRANS] = "Translation     ";
    float freqMin = 1e10;
    for (i = MV_NULL + 1; i < MAX_MV; i++) {
        if (freqMin >= fMCFreq[i] && fMCFreq[i] != 0.) {
            freqMin = fMCFreq[i];
        }
    }
    printf("MC Move Frequencies:\n");
    printf("--------------------------\n");
    for (i = MV_NULL + 1; i < MAX_MV; i++) {
        printf("|%s  %d\n", MoveName[i], (int) ceilf(fMCFreq[i] / freqMin));
    }
    printf("--------------------------\n");
    printf("\n");
}

/// Print_Matrix - fancy function to print symmetric matrices with labels
/// \param strTitle
/// \param nSeqEn
/// \param fArray
/// \param param
void Print_Matrix(char *strTitle, int nSeqEn, float fArray[MAX_AA][MAX_AA][MAX_E], int param) {
    int nLen;
    nLen = nSeqEn;
    int i, j;
    for (i = 0; i < nLen + 1; i++) {
        printf("--------");
    }
    printf("\n|%s", strTitle);
    for (i = 0; i < nLen; i++) {
        printf("\t%d", i);
    }
    printf("      |\n");

    for (i = 0; i < nLen; i++) {
        printf("|%d", i);
        for (j = 0; j < nLen; j++) {
            printf("\t%.2f", fArray[i][j][param]);
        }
        printf("   |\n");
    }
    for (i = 0; i < nLen + 1; i++) {
        printf("--------");
    }
    printf("\n");
}

/// Write_RDF_ComponentWise - old implementation of printing the RDF, component by component.
/// Always appends to the file for this run. Stopped using it because the IO load was slowing things down at our cluster
/// Would be a good way to gather proper statistics on clusters as the runs went on. TODO: update this to work with
/// the new indexing
/// \param filename
/// \param nGen
void Write_RDF_ComponentWise(char *filename, long nGen) {
    FILE *fp;
    if (nGen == -1) {
        fp = fopen(filename, "w"); // overwrite
    } else {
        fp = fopen(filename, "a");
    }
    int i, k;
    if (nGen == -1) { // title
        fprintf(fp, "#Row-by-row split RDF. dr = 1/4\n");
    } else {
        for (k = 0; k < nRDF_TotComps; k++) {
            for (i = 0; i < nRDF_TotBins; i++) {
                fprintf(fp, "%.5Lf\t", ldRDF_Arr[RDFArr_Index(0, k, i)]);
            }
            fprintf(fp, "\n");
        }
    }

    fclose(fp);

}

/// Write_TopFile - write a LAMMPS DATA file, which can be read in VMD to store topological information for the system
/// Only writes the topology once at the beginning of the runs.
/// \param filename
void Write_TopFile(char *filename) {
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
    for (i = 0; i < tot_beads; i++) {
        for (j = 0; j < MAX_BONDS; j++) {
            if (topo_info[i][j] != -1) {//There is a bond between i and topo_info[i][j]
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
    fprintf(fp, "1\t1.0\n2\t1.0\n3\t1.0\n4\t1.0\n5\t1.0\n6\t1.0\n");//Dummy masses for the 6 types.
    fprintf(fp, "\n");//Empty line.
    fprintf(fp, "Atoms\n");//Signifying the beginning of atom coordinates.
    fprintf(fp, "\n");//Empty line.
    for (i = 0; i < tot_beads; i++) {
        fprintf(fp, "%d %d %d 0.0 %d %d %d\n", i, bead_info[i][BEAD_CHAINID], bead_info[i][BEAD_TYPE],
                bead_info[i][POS_X], bead_info[i][POS_Y], bead_info[i][POS_Z]);
    }//Done with the coordinates
    fprintf(fp, "\n");//Empty line.
    fprintf(fp, "Bonds\n");//Signifying the beginning of atom coordinates.
    fprintf(fp, "\n");//Empty line.
    k = 0;//This guy counts bondIDs
    for (i = 0; i < tot_beads; i++) {
        for (j = 0; j < MAX_BONDS; j++) {
            if (topo_info[i][j] != -1) {//There is a bond between i and topo_info[i][j]
                fprintf(fp, "%d %d %d %d\n", k, linker_len[i][j] / (int) fLinkerLength, i, topo_info[i][j]);
                k++;
            }
        }
    }//Done with the coordinates
    fprintf(fp, "\n");//Empty line.

    fclose(fp);
}

/// Write_SysProp - writes the RDF, CLUS and GyrTen avg rrays to a file. OLD implementation that should not be used yet
/// Also, TODO: update this to work with the new indexing
/// \param filename
void Write_SysProp(char *filename) {

    FILE *fp;
    int i, j;
    fp = fopen(filename, "w");
    fprintf(fp, "#Contains various averaged quatities.\n");
    //Gyration Radius
    fprintf(fp, "#Total Rg\n%f\t%f\n#Cluster Hist\n", fSysGyrRad / (float) nTotGyrRadCounter, (float) nBoxSize[0] / 2.);
    //Cluster Histogram
    fprintf(fp, "%f\t", (float) nLargestClusterRightNow / (float) nTotClusCounter);
    for (i = 1; i <= tot_chains; i++) {
        fprintf(fp, "%f\t", (float) naClusHistList[i] / (float) nTotClusCounter);
    }
    //Split RDFs
    fprintf(fp, "\n#Split RDFs. ALL-ALL; DIAGONALS and then from 0 onwards \n");
    for (j = 0; j < nRDF_TotComps; j++) {
        for (i = 0; i < nRDF_TotBins; i++) {
            //fprintf(fp, "%LE\t", ldRDF_ARR[j][i] / (float)nRDFCounter);
            fprintf(fp, "%LE\t", ldRDF_Arr[RDFArr_Index(0, j, i)] / (float) nRDFCounter);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n#Done");

    fclose(fp);

}

/// Write_TotalSysProp - writes the RDF, CLUS and GyrTen arrays to their respective files.
/// These arrays store the data over the course of an ENTIRE run, or over all the temp cycles.
/// \param filename
/// \param run_it
void Write_TotalSysProp(char *filename, int run_it) {
    /* This function writes one large file with all the averaged values from each run_cycle. run_it let's the function
     * know how many cycles to write.
     * As opposed to Write_SysProp, each of the relevant averaged parameters will be in their appropriately named files.
     * The naming convention shall be: filename_*.dat where * will be GR, CLUS and RDF for the different measured
     * quantities. Within each file, each run_cycle shall be written one-by-one in the style of Write_SysProp
     */
    FILE *fp;
    int i, j, k;
    sprintf(filename, "%s_RDF.dat", strReportPrefix);//Name Of the RDF Files
    fp = fopen(filename, "w");
    fprintf(fp, "#Split RDFs. ALL-ALL; DIAGONALS and then from 0-0, 0-1, and onwards \n");
    for (i = 0; i < run_it; i++) {
        fprintf(fp, "#Run_Cycle = %d\n", i);
        for (j = 0; j < nRDF_TotComps; j++) {
            for (k = 0; k < nRDF_TotBins; k++) {
                fprintf(fp, "%LE\t", ld_TOTRDF_Arr[RDFArr_Index(i, j, k)]);
                //fprintf(fp, "%LE\t", ld_TOTRDF_ARR[i][j][k]);
            }
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "#Done");
    fclose(fp);

    sprintf(filename, "%s_CLUS.dat", strReportPrefix);//Name Of the ClusterHistogram Files
    fp = fopen(filename, "w");
    fprintf(fp, "#Cluster Histograms: 1st column is largest cluster, and then clusters of size 1, 2, and so on\n");
    for (i = 0; i < run_it; i++) {
        fprintf(fp, "#Run_Cycle = %d\n", i);
        for (k = 0; k <= tot_chains; k++) {
            fprintf(fp, "%LE\t", ld_TOTCLUS_ARR[i][k]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "#Done");
    fclose(fp);


    sprintf(filename, "%s_GR.dat", strReportPrefix);//Name Of the ClusterHistogram Files
    fp = fopen(filename, "w");
    fprintf(fp, "#Cluster Histograms: 1st column is largest cluster, and then clusters of size 1, 2, and so on\n");
    for (i = 0; i < run_it; i++) {
        fprintf(fp, "#Run_Cycle = %d\n", i);
        fprintf(fp, "%LE\t", ld_TOTGYRRAD_ARR[i][0]);
        fprintf(fp, "%LE\n", ld_TOTGYRRAD_ARR[i][1]);
    }
    fprintf(fp, "#Done");
    fclose(fp);

}

/// Print_Data - helper function that decides given nGen and run_it which things to print.
/// Usually don't print much during the thermalization but could change that here
/// \param nGen
/// \param run_it
void Print_Data(long nGen, int run_it) {
    //This function handles all the data IO.
    int nFlagForEnCalc = 0; //Flag for total energy calculation
    //run_it == -1 corresponds to the thermalization sequence
    if (run_it == -1) {
        //Open the appropriate files before the thermalization sequence.
        if (nGen == -1) {
            if (nReport[REPORT_CONFIG] != 0) {
                sprintf(fileStruct, "%s_topo.lammpstrj", strReportPrefix);//Name of the topology file
                Write_TopFile(fileStruct);//Write the topology file. Only need to write once
                sprintf(fileStruct, "%s_trj.lammpstrj", strReportPrefix);//Naming convention for trajectory files.
                Write_Trajectory(fileStruct, -1);//This opens a new trajectory file; each run_it will have its own
            }
            if (nReport[REPORT_ENERGY] != 0) {
                sprintf(fileEnergy, "%s_energy.dat", strReportPrefix);
                Write_Energy(fileEnergy, -1);//Open a new energy file; each run_it will have its own
            }
            if (nReport[REPORT_MCMOVE] != 0) {
                sprintf(fileMCMove, "%s_mcmove.dat", strReportPrefix);
                Write_MCMove(fileMCMove, -1, 0.0);//Open a new MCInfo file; each run_it will have its own
            }
            if (nReport[REPORT_NETWORK] != 0 || nReport[REPORT_RDFTOT] != 0) {
                sprintf(fileSysProp, "%s_SysProp.dat",
                        strReportPrefix);//Naming convention or SysProp files; each run_it should have its own
            }
        } else {
            if (nReport[REPORT_LOG] != 0) {
                if (nGen % nReport[REPORT_LOG] == 0) {
                    if (Check_System_Structure() == 0) {
                        printf("Structure is still a-okay!\n");
                    } else {
                        printf("Molecular structure is inconsistent with initial structure.\nCRASHING\n\n");
                        exit(1);
                    }
                    Energy_Total_System();
                    nFlagForEnCalc = 1;
                    Print_LogToScreen(nGen, run_it);
                }
            }
            if (nReport[REPORT_CONFIG] != 0) {
                if (nGen % nReport[REPORT_CONFIG] == 0) {
                    Write_Trajectory(fileStruct, nGen);
                }
            }
            if (nReport[REPORT_ENERGY] != 0) {
                if (nGen % nReport[REPORT_ENERGY] == 0) {
                    if (nFlagForEnCalc != 1) {//Calculate the energy
                        Energy_Total_System();
                        nFlagForEnCalc = 1;
                    }
                    Write_Energy(fileEnergy, nGen);
                }
            }
            if (nReport[REPORT_MCMOVE] != 0) {
                if (nGen % nReport[REPORT_MCMOVE] == 0) {
                    Write_MCMove(fileMCMove, nGen, fCuTemp);
                }
            }
        }

    }
    if (run_it == 0 && nGen > 0) {
        if (nReport[REPORT_LOG] != 0) {
            if (nGen % nReport[REPORT_LOG] == 0) {
                if (Check_System_Structure() == 0) {
                    printf("Structure is still a-okay!\n");
                } else {
                    printf("Molecular structure is inconsistent with initial structure.\nCRASHING\n\n");
                    exit(1);
                }
                Energy_Total_System();
                nFlagForEnCalc = 1;
                Print_LogToScreen(nGen, run_it);
            }
        }
        if (nReport[REPORT_CONFIG] != 0) {
            if (nGen % nReport[REPORT_CONFIG] == 0) {
                Write_Trajectory(fileStruct, nGen + nPreSteps);
            }
        }
        if (nReport[REPORT_ENERGY] != 0) {
            if (nGen % nReport[REPORT_ENERGY] == 0) {
                if (nFlagForEnCalc != 1) {//Calculate the energy
                    Energy_Total_System();
                    nFlagForEnCalc = 1;
                }
                Write_Energy(fileEnergy, nGen + nPreSteps);
            }
        }
        if (nReport[REPORT_MCMOVE] != 0) {
            if (nGen % nReport[REPORT_MCMOVE] == 0) {
                Write_MCMove(fileMCMove, nGen + nPreSteps, fCuTemp);
            }
        }
        if (nReport[REPORT_RDFTOT] != 0) {//SysProp is printed outside of this function in main.c, lol
            if (nGen % nReport[REPORT_RDFTOT] == 0 && nGen > nSteps / 2) {
                RDF_ComponentWise_Avg();
            }
        }
        if (nReport[REPORT_NETWORK] != 0) {//SysProp is printed outside of this function in main.c, lol
            if (nGen % nReport[REPORT_NETWORK] == 0 && nGen > nSteps / 2) {
                Clus_Distribution_Avg();
                GyrTensor_GyrRad_Avg();
            }
        }
    }

    if (run_it > 0) {
        if (nGen == -1) {
            if (nReport[REPORT_CONFIG] != 0) {
                sprintf(fileStruct, "%s_%d_trj.lammpstrj", strReportPrefix,
                        run_it);//Naming convention for trajectory files.
                Write_Trajectory(fileStruct, -1);//This opens a new trajectory file; each run_it will have its own
                Write_Trajectory(fileStruct, 0);//End of previous run_it is initial conditions for this one!
                sprintf(fileStruct, "%s_trj.lammpstrj", strReportPrefix);//Returning back
            }
            if (nReport[REPORT_ENERGY] != 0) {
                sprintf(fileEnergy, "%s_%d_energy.dat", strReportPrefix, run_it);
                Write_Energy(fileEnergy, -1);//Open a new energy file; each run_it will have its own
                Write_Energy(fileEnergy, 0);//Last frame from previous run_it is first frame of this one
            }
            if (nReport[REPORT_MCMOVE] != 0) {
                sprintf(fileMCMove, "%s_%d_mcmove.dat", strReportPrefix, run_it);
                Write_MCMove(fileMCMove, -1, 0.0);//Open a new MCInfo file; each run_it will have its own
            }
            if (nReport[REPORT_NETWORK] != 0 || nReport[REPORT_RDFTOT] != 0) {
                sprintf(fileSysProp, "%s_%d_SysProp.dat", strReportPrefix,
                        run_it);//Naming convention or SysProp files; each run_it should have its own
            }
        } else {
            if (nReport[REPORT_LOG] != 0) {
                if (nGen % nReport[REPORT_LOG] == 0) {
                    if (Check_System_Structure() == 0) {
                        printf("Structure is still a-okay!\n");
                    } else {
                        printf("Molecular structure is inconsistent with initial structure.\nCRASHING\n\n");
                        exit(1);
                    }
                    Energy_Total_System();
                    nFlagForEnCalc = 1;
                    Print_LogToScreen(nGen, run_it);
                }
            }
            if (nReport[REPORT_CONFIG] != 0) {
                if (nGen % nReport[REPORT_CONFIG] == 0) {
                    Write_Trajectory(fileStruct, nGen + (run_it * nSteps));
                }
            }
            if (nReport[REPORT_ENERGY] != 0) {
                if (nGen % nReport[REPORT_ENERGY] == 0) {
                    if (nFlagForEnCalc != 1) {//Calculate the energy
                        Energy_Total_System();
                        nFlagForEnCalc = 1;
                    }
                    Write_Energy(fileEnergy, nGen);
                }
            }
            if (nReport[REPORT_MCMOVE] != 0) {
                if (nGen % nReport[REPORT_MCMOVE] == 0) {
                    Write_MCMove(fileMCMove, nGen, fCuTemp);
                }
            }
            if (nReport[REPORT_RDFTOT] != 0) {//SysProp is printed outside of this function in main.c, lol
                if (nGen % nReport[REPORT_RDFTOT] == 0 && nGen > nSteps / 2) {
                    RDF_ComponentWise_Avg();
                }
            }
            if (nReport[REPORT_NETWORK] != 0) {//SysProp is printed outside of this function in main.c, lol
                if (nGen % nReport[REPORT_NETWORK] == 0 && nGen > nSteps / 2) {
                    Clus_Distribution_Avg();
                    GyrTensor_GyrRad_Avg();
                }
            }
        }
    }

}

/// Copy_Data - copies data from run_it specific data arrays to the overall global data arrays that are printed later.
/// Also note that the averaging, or dividing by the frequency of acquisitions, occurs here.
/// \param run_it
void Copy_Data(int run_it) {
    int i, j;
    for (i = 0; i < nRDF_TotComps; i++) {
        for (j = 0; j < nRDF_TotBins; j++) {
            ld_TOTRDF_Arr[RDFArr_Index(run_it, i, j)] =
                    ldRDF_Arr[RDFArr_Index(0, i, j)] / (long double) nRDFCounter;

        }
    }
    ld_TOTCLUS_ARR[run_it][0] = (long double) nLargestClusterRightNow / (long double) nTotClusCounter;
    for (i = 1; i <= tot_chains; i++) {
        ld_TOTCLUS_ARR[run_it][i] = (long double) naClusHistList[i] / (long double) nTotClusCounter;
    }
    ld_TOTGYRRAD_ARR[run_it][0] = (long double) fSysGyrRad / (long double) nTotGyrRadCounter;
    ld_TOTGYRRAD_ARR[run_it][1] = (long double) nBoxSize[0] / 2.;
}