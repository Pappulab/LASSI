from __future__ import division
import numpy as np
import scipy as sp
import os
import shutil
import subprocess as sproc
import time
import pickle
import matplotlib.pyplot as plt
__author__ = 'Furqan Dar'

def MKDirCatch(this_dir):
    try:
        os.mkdir(this_dir)
    except OSError as myErr:
        if myErr[0] != 17:
            print("There is something wrong!")
            raise

def Generate_Dir_Tree(systems_list, linker_lengths, box_size, runs_per_condition, path_of_tree):
    for struc in systems_list:
        struc_dir = path_of_tree + struc
        MKDirCatch(struc_dir)
        for linLen in linker_lengths:
            lin_dir = struc_dir + linLen
            MKDirCatch(lin_dir)
            for boxSize in box_size:
                box_dir = lin_dir + boxSize + '/'
                MKDirCatch(box_dir)
                for run_con in ['NoInt/', 'WInt/']:
                    con_dir = box_dir + run_con
                    MKDirCatch(con_dir)
                    for run_it in range(runs_per_condition):
                        run_dir = con_dir + str(run_it + 1)
                        MKDirCatch(run_dir)

def Read_ParamFile(file_name):
    dum_dict = {};
    dum_list = [];
    with open(file_name) as pFile:
        totFile = pFile.readlines()
        for a_line in totFile:
            # print a_line[0]
            if a_line[0] == '#' or a_line[0] == '\n':
                continue
            s_line = a_line[:-1].split(" ")
            s_line = [a_key for a_key in s_line if a_key != '']
            s_line = s_line[:2]
            dum_list.append(s_line[0])
            dum_dict[s_line[0]] = s_line[1]
    return [dum_list, dum_dict]

def Write_Param_InTree(systems_list, linker_lengths, box_size, runs_per_condition, path_of_tree, param_Arr):
    dKeys = param_Arr[0]
    dDict = param_Arr[1]
    for struc in systems_list:
        struc_dir = path_of_tree + struc
        for linLen in linker_lengths:
            lin_dir = struc_dir + linLen
            for boxSize in box_size:
                dDict['BOX_SIZE'] = boxSize
                box_dir = lin_dir + boxSize + '/'
                for run_con in ['NoInt/', 'WInt/']:
                    con_dir = box_dir + run_con
                    for run_it in range(runs_per_condition):
                        run_dir = con_dir + str(run_it + 1)
                        dDict['RANDOM_SEED'] = str(run_it + 1)
                        with open(run_dir + '/param.key', 'w+') as pFile:
                            for a_key in dKeys:
                                N_spcs = 25 - len(a_key)
                                pFile.write(a_key + ' ' * N_spcs + dDict[a_key] + '\n')

def WriteEnergy_InTree(systems_list, linker_lengths, box_size, runs_per_condition, path_of_tree, enInt_filename,
                       enNoInt_filename):
    for struc in systems_list:
        struc_dir = path_of_tree + struc
        for linLen in linker_lengths:
            lin_dir = struc_dir + linLen
            for boxSize in box_size:
                box_dir = lin_dir + boxSize + '/'
                for run_con, cur_file in zip(['NoInt/', 'WInt/'], [enNoInt_filename, enInt_filename]):
                    con_dir = box_dir + run_con
                    orig_file = cur_file
                    for run_it in range(runs_per_condition):
                        run_dir = con_dir + str(run_it + 1)
                        this_file = run_dir + '/energy.prm'
                        shutil.copyfile(orig_file, this_file)

def WriteStruc_InTree(systems_list, linker_lengths, box_size, runs_per_condition, path_of_tree, struc_filename):
    for struc in systems_list:
        struc_dir = path_of_tree + struc
        for linLen in linker_lengths:
            lin_dir = struc_dir + linLen
            for boxSize in box_size:
                box_dir = lin_dir + boxSize + '/'
                for run_con in ['NoInt/', 'WInt/']:
                    con_dir = box_dir + run_con
                    for run_it in range(runs_per_condition):
                        run_dir = con_dir + str(run_it + 1)
                        this_file = run_dir + '/structure.prm'
                        shutil.copyfile(struc_filename, this_file)

def SubmitJobs_ToQueue(systems_list, linker_lengths, box_size, runs_per_condition, path_of_tree, qsub_command,
                       path_to_LASSI):
    for s_ID, struc in enumerate(systems_list):
        struc_dir = path_of_tree + struc
        for l_ID, linLen in enumerate(linker_lengths):
            lin_dir = struc_dir + linLen
            for b_ID, boxSize in enumerate(box_size):
                box_dir = lin_dir + boxSize + '/'
                for c_ID, run_con in enumerate(['NoInt/', 'WInt/']):
                    con_dir = box_dir + run_con
                    for run_it in range(runs_per_condition):
                        run_dir = con_dir + str(run_it + 1)
                        run_num = ['I', s_ID, b_ID, run_con[0], run_it]
                        run_num = [str(a_val) for a_val in run_num]
                        run_num = "_".join(run_num)
                        run_command = qsub_command + ' ' + run_num + ' ' + path_to_LASSI
                        os.chdir(run_dir)
                        # print run_command.split(" ")
                        try:
                            ret_code = sproc.check_output(run_command, shell=True, stderr=sproc.STDOUT)
                            if ret_code < 0:
                                print("The submission failed! Signal: ", -retcode)
                            else:
                                print(ret_code)
                        except OSError as myErr:
                            print("Couldn't submit because ", myErr)
                            raise
                        time.sleep(2)
    # Return back to starting directory
    os.chdir(path_of_tree)
    os.chdir('../')

def sGen_Implicit_Linear(Mol_Num, Bead_Num, Bead_Type, Linker_Len):
    if not type(Mol_Num) == type(Bead_Num) == type(Bead_Type) \
           == type(Linker_Len) == int:
        print("Every paramater needs to be an integer!")
        return
    try:
        struc_len = (Bead_Num - 1) * 2
        dumArr = np.zeros((struc_len, 4), dtype=int)

        stickList = np.array([int(np.ceil(i / 2)) for i in range(struc_len)])
        bondList = np.array([int(np.ceil((aBead + 1) / 2)) if aBead % 2 == 0 else int(np.ceil((aBead - 1) / 2)) \
                             for aBead in range(struc_len)])

        dumArr.T[0] = stickList
        dumArr.T[1] = Bead_Type
        dumArr.T[2] = Linker_Len
        dumArr.T[3] = bondList
    except TypeError:
        print("Did you format the molecule correctly?")

    return [Mol_Num, dumArr]

def sGen_Monomer(Mol_Num, Bead_Type):
    if not type(Mol_Num) == type(Bead_Type) == int:
        print("Every paramater needs to be an integer!")
        return
    try:
        dumArr = np.zeros((1, 4), dtype=int)

        dumArr.T[0] = 0
        dumArr.T[1] = Bead_Type
        dumArr.T[2] = 1
        dumArr.T[3] = -1
    except TypeError:
        print("Did you format the molecule correctly?")

    return [Mol_Num, dumArr]

def sGen_Dimer(Mol_Num, Bead_Type1, Bead_Type2, Linker_Len):
    if not type(Mol_Num) == type(Bead_Type2) == type(Bead_Type1) \
           == type(Linker_Len) == int:
        print("Every paramater needs to be an integer!")
        return
    try:
        struc_len = (2 - 1) * 2
        dumArr = np.zeros((struc_len, 4), dtype=int)

        stickList = np.array([int(np.ceil(i / 2)) for i in range(struc_len)])
        bondList = np.array([int(np.ceil((aBead + 1) / 2)) if aBead % 2 == 0 else int(np.ceil((aBead - 1) / 2)) \
                             for aBead in range(struc_len)])

        dumArr.T[0] = [0, 1]
        dumArr.T[1] = [Bead_Type1, Bead_Type2]
        dumArr.T[2] = Linker_Len
        dumArr.T[3] = [1, 0]
    except TypeError:
        print("Did you format the molecule correctly?")

    return [Mol_Num, dumArr]

def sGen_Implicit_SymmBranched(Mol_Num, Branch_Num, Beads_PerBranch, Bead_Type, Hub_Type, Linker_Len, Hub_Linker):
    if not type(Mol_Num) == type(Branch_Num) == type(Bead_Type) \
           == type(Linker_Len) == type(Beads_PerBranch) == type(Hub_Type) == int:
        print("Every paramater needs to be an integer!")
        return
    try:
        if Branch_Num == 2:
            dumArr = sGen_Implicit_Linear(-1, Beads_PerBranch * 2 + 1, Bead_Type, Linker_Len)[1]
            dumArr[Beads_PerBranch * 2 - 1][1] = Hub_Type
            dumArr[Beads_PerBranch * 2][1] = Hub_Type
        else:
            struc_len = Branch_Num * (1 + Beads_PerBranch * 2 - 1)
            dumArr = np.zeros((struc_len, 4), dtype=int)

            stick_list = [[0 for i in range(Branch_Num)]]
            stick_list.append([int(np.ceil(i / 2)) for i in range(1, Beads_PerBranch * 2)])
            for aBranch in range(Branch_Num - 1):
                stick_list.append([stick_list[-1][-1] + int(np.ceil(i / 2)) for i in range(1, Beads_PerBranch * 2)])
            stick_list = [aBead for aList in stick_list for aBead in aList]
            stick_list = np.array(stick_list)

            bond_list = [[i*Beads_PerBranch + 1 for i in range(Branch_Num)]]
            bond_list.append([int(np.ceil((aBead + 1) / 2)) if aBead % 2 == 0 else int(np.ceil((aBead - 1) / 2)) \
                              for aBead in range(1, Beads_PerBranch * 2)])
            for aBranch in range(Branch_Num - 1):
                dum_list = [1 + bond_list[-1][-1] + int(np.ceil((aBead + 1) / 2)) if aBead % 2 == 0 \
                                else 1 + bond_list[-1][-1] + int(np.ceil((aBead - 1) / 2)) \
                            for aBead in range(1, Beads_PerBranch * 2)]
                dum_list[0] = 0
                bond_list.append(dum_list)
            bond_list = [aBead for aList in bond_list for aBead in aList]
            bond_list = np.array(bond_list)

            dumArr.T[0] = stick_list
            dumArr.T[1] = Bead_Type;
            dumArr.T[1][:Branch_Num] = Hub_Type
            dumArr.T[2] = Linker_Len
            dumArr.T[2][:Branch_Num] = Hub_Linker
            dumArr.T[2][Branch_Num::Beads_PerBranch*2-1] = Hub_Linker
            dumArr.T[3] = bond_list
    except TypeError:
        print("Did you format the molecule correctly?")
    return [Mol_Num, dumArr]

def sGen_Explicit_Linear(Mol_Num, Bead_Num, Bead_Type, Linker_Len, Linker_Type):
    if not type(Mol_Num) == type(Bead_Num) == type(Bead_Type) == type(Linker_Len) == type(Linker_Type) == int:
        print("Every paramater needs to be an integer!")
        return
    try:
        if Linker_Len < 2:
            print("Explicit linker molecules must have at least a linker length"
                  "of 2 lattice units! Doing nothing.")
            return
        else:
            tot_beads = Bead_Num + (Bead_Num - 1) * (Linker_Len - 1)
            struc_len = (tot_beads - 1) * 2
            dumArr = np.zeros((struc_len, 4), dtype=int)

            stickList = np.array([int(np.ceil(i / 2)) for i in range(struc_len)])
            bondList = np.array([int(np.ceil((aBead + 1) / 2)) if aBead % 2 == 0 else int(np.ceil((aBead - 1) / 2)) \
                                 for aBead in range(struc_len)])
            dumArr.T[0] = stickList
            dumArr.T[1] = Linker_Type

            dumArr.T[1][0] = Bead_Type;
            dumArr.T[1][-1] = Bead_Type;
            dumArr.T[1][1 + 2 * (Linker_Len - 1)::1 + 2 * (Linker_Len - 1)] = Bead_Type
            dumArr.T[1][1 + 2 * (Linker_Len - 1) + 1::1 + 2 * (Linker_Len - 1)] = Bead_Type

            dumArr.T[2] = 1
            dumArr.T[3] = bondList
    except TypeError:
        print("Did you format the molecule correctly?")

    return [Mol_Num, dumArr]

def sGen_WriteStructuresToFile(StrucList, file_name):
    with open(file_name, "w+") as strucFile:
        for a_Mol in StrucList:
            mol_num = a_Mol[0];
            mol_struc = a_Mol[1];
            beads_per_mol = len(np.unique(mol_struc.T[0]))
            strucFile.write("#New Molecule Type:= {:} beads per molecule\n".format(beads_per_mol))
            strucFile.write("NEW{\n")
            strucFile.write(str(mol_num) + "\n")
            for a_line in mol_struc:
                this_line = [str(a_num) + "\t" for a_num in a_line]
                this_line.append("\n")
                strucFile.write("".join(this_line))
            strucFile.write("}END\n")
    #print("Wrote structures to file: {:}".format(file_name))

def Index_RDF(compA, compB, totComp):
    assert type(compA) == type(compB) == type(totComp) == int
    if compA > compB:
        return Index_RDF(compB, compA, totComp)
    elif compA == compB == -1:
        return 0
    else:
        return 1 + compA if compA == compB else totComp + compB - int((compA * (3 + compA - 2 * totComp)) / 2)

def Index_COMDen(compA, compB, totMols):
    assert type(compA) == type(compB) == type(totMols) == int
    if compA < 0:
        return compB
    else:
        return totMols + compB + totMols * compA

def MeanAndError(theArray, theAxis):
    return np.mean(theArray, axis=theAxis), np.std(theArray, axis=theAxis)

def Get_Cubic_Lattice_Normalization(boxSize):
    dum_list = []
    with open('/home/fdar/Work/LASSI/PythonCodes/MyModules/Lat_RadialNorm_dr4.txt', 'r') as nFile:
        for i, line in enumerate(nFile):
            if i == boxSize-10+1:
                dum_list = np.fromstring(line,sep=' ')[1:]
    return dum_list

def RDFVol_F1(xVal, boxSize):
    return np.arctan(np.sqrt(-2. + (4. * (xVal ** 2.) / (boxSize ** 2.))))

def RDFVol_F2(xVal, boxSize):
    return (8. * xVal * np.arctan((2. * xVal * (-3. + (4. * (xVal ** 2.)) / (boxSize ** 2.))) /
                                  (boxSize * np.sqrt(-2. + (4. * (xVal ** 2.)) / (boxSize ** 2.)) *
                                   (1. + (4. * (xVal ** 2.)) / (boxSize ** 2.))))) / boxSize

def RDFVolumeElement(xAr, boxSize):
    '''
    Correctly calculates the spherical volume element in a cubic lattice with periodic boundary conditions
    '''
    new_x_Ar = np.ones(len(xAr));
    for xId, xVal in enumerate(xAr):
        if xVal <= boxSize / 2.:
            new_x_Ar[xId] = 4. * np.pi * (xVal ** 2.);
        if xVal > boxSize / 2. and xVal <= np.sqrt(2.) * boxSize / 2.:
            new_x_Ar[xId] = 2. * np.pi * (3. * boxSize - 4. * xVal) * xVal
        if xVal > np.sqrt(2.) * boxSize / 2. and xVal <= np.sqrt(3.) * boxSize / 2.:
            new_x_Ar[xId] = 2. * boxSize * (
                    -12. * RDFVol_F1(xVal, boxSize) + RDFVol_F2(xVal, boxSize) + 3. * np.pi) * xVal
    return new_x_Ar

def AddPlot_WTieLines(plObj, comp1, comp2, plColor, erKW=None, plKW=None):
    c1, c1_er = comp1;
    c2, c2_er = comp2;

    if erKW == None:
        erKW = {'capthick':2, 'capsize':2.5, 'elinewidth':2, 'alpha':0.7, 'fmt':'o', 'markersize':10}
    if plKW == None:
        plKW = {'lw':2, 'markersize':10, 'alpha':0.7}


    plObj.errorbar(c1[1:], c2[1:], xerr=c1_er, yerr=c2_er, color=plColor, **erKW)

    plObj.plot([c1[1], c1[0]], [c2[1], c2[0]], '-', color=plColor, marker='o', **plKW)
    plObj.plot([c1[0]], [c2[0]], 'o', color=plColor, **plKW)
    plObj.plot([c1[2], c1[0]], [c2[2], c2[0]], '-', color=plColor, marker='s', **plKW)

def AddPlot_WTieLinesDashed(plObj, comp1, comp2, plColor, erKW=None, plKW=None):
    c1, c1_er = comp1;
    c2, c2_er = comp2;

    if erKW == None:
        erKW = {'capthick':2, 'capsize':2.5, 'elinewidth':2, 'alpha':0.7, 'fmt':'o', 'markersize':10}
    if plKW == None:
        plKW = {'lw':2, 'markersize':10, 'alpha':0.7}

    plObj.errorbar(c1[1:], c2[1:], xerr=c1_er, yerr=c2_er, color=plColor, **erKW)

    plObj.plot([c1[1], c1[0]], [c2[1], c2[0]], '--', color=plColor, **plKW)
    plObj.plot([c1[0]], [c2[0]], 'o', color=plColor, **plKW)
    plObj.plot([c1[2], c1[0]], [c2[2], c2[0]], '--', color=plColor, **plKW)

def AddPlot_WTieLines_Explicit(plObj, myXAr, myXEr, myYAr, myYEr, con1_0, con2_0, plColor, erKW, plKW):
    plObj.errorbar(myXAr, myYAr, xerr=myXEr, yerr=myYEr, color=plColor, **erKW)

    plObj.plot([myXAr[0], con1_0], [myYAr[0], con2_0], '-', color=plColor, **plKW)
    plObj.plot([con1_0], [con2_0], 'o', color=plColor, **plKW)
    plObj.plot([myXAr[1], con1_0], [myYAr[1], con2_0], '-', color=plColor, **plKW)

def AddPlot_WOTieLines(plObj, comp1, comp2, plColor, erKW=None):
    c1, c1_er = comp1;
    c2, c2_er = comp2;

    if erKW == None:
        erKW = {'capthick':2, 'capsize':2.5, 'elinewidth':2, 'alpha':0.7, 'fmt':'o', 'markersize':10}


    plObj.errorbar(c1[1:], c2[1:], xerr=c1_er, yerr=c2_er, color=plColor, **erKW)

def AddPlot_WOTieLines_EXP(plObj, myXAr, myXEr, myYAr, myYEr, con1_0, con2_0, plColor, erKW, plKW):
    plObj.errorbar(myXAr, myYAr, xerr=myXEr, yerr=myYEr, color=plColor, **erKW)

    #plObj.plot([myXAr[0], con1_0], [myYAr[0], con2_0], '-', color=plColor, **plKW)
    #plObj.plot([con1_0], [con2_0], 'o', color=plColor, **plKW)
    #plObj.plot([myXAr[1], con1_0], [myYAr[1], con2_0], '-', color=plColor, **plKW)

def ForPlot_Get_OneComp_TieLineConcs(tot_tie_data, temp_id, box_id, comp1, comp2):

    dat_comp1  = tot_tie_data[box_id][comp1][temp_id]
    dat_comp2  = tot_tie_data[box_id][comp2][temp_id]

    cons_comp1 = np.array([dat_comp1[0],
                           dat_comp1[1], dat_comp1[3]])
    err_comp1  = np.array([dat_comp1[2], dat_comp1[4]])
    cons_comp2 = np.array([dat_comp2[0],
                           dat_comp2[1], dat_comp2[3]])
    err_comp2  = np.array([dat_comp2[2], dat_comp2[4]])

    return [[cons_comp1, err_comp1], [cons_comp2, err_comp2]]

def ForPlot_Get_TwoComp_TieLineConcs(tot_tie_data, temp_id, box_id, num_id, comp1, comp2):

    dat_comp1  = tot_tie_data[temp_id][box_id][num_id][comp1]
    dat_comp2  = tot_tie_data[temp_id][box_id][num_id][comp2]

    cons_comp1 = np.array([dat_comp1[0],
                           dat_comp1[1], dat_comp1[3]])
    err_comp1  = np.array([dat_comp1[2], dat_comp1[4]])
    cons_comp2 = np.array([dat_comp2[0],
                           dat_comp2[1], dat_comp2[3]])
    err_comp2  = np.array([dat_comp2[2], dat_comp2[4]])

    return [[cons_comp1, err_comp1], [cons_comp2, err_comp2]]

def ForPlot_Get_Ortho_TieLineConcs(tot_tie_data, temp_id, box_id, comp1, comp2):

    dat_comp1  = tot_tie_data[temp_id][comp1].T[box_id]
    dat_comp2  = tot_tie_data[temp_id][comp2].T[box_id]

    cons_comp1 = np.array([dat_comp1[0],
                           dat_comp1[1], dat_comp1[3]])
    err_comp1  = np.array([dat_comp1[2], dat_comp1[4]])
    cons_comp2 = np.array([dat_comp2[0],
                           dat_comp2[1], dat_comp2[3]])
    err_comp2  = np.array([dat_comp2[2], dat_comp2[4]])

    return [[cons_comp1, err_comp1], [cons_comp2, err_comp2]]

def ForPlot_GenLogConcs(comp1, comp2):

    comp1_log = np.log10(comp1[0]); comp1_log_err = comp1[1]/comp1[0][1:]/np.log(10.)
    comp2_log = np.log10(comp2[0]); comp2_log_err = comp2[1]/comp2[0][1:]/np.log(10.)

    return [[comp1_log, comp1_log_err], [comp2_log, comp2_log_err]]

def ForPlot_GenRatioConcs(comp1, comp2):

    tot_con = comp1[0]+comp2[0]
    tot_con_err = np.sqrt(comp1[1]**2. + comp2[1]**2.)
    tot_con_log = np.log10(tot_con)
    tot_con_log_err = tot_con_err/tot_con[1:]/np.log(10.)

    con_ratio = comp1[0]/comp2[0]
    con_ratio_err = con_ratio[1:] * np.sqrt((comp1[1]/comp1[0][1:])**2. + (comp2[1]/comp2[0][1:])**2.)
    con_ratio_log = np.log10(con_ratio)
    con_ratio_log_err = con_ratio_err/con_ratio[1:]/np.log(10.)

    return [[tot_con_log, tot_con_log_err], [con_ratio_log, con_ratio_log_err]]

def AddPlot_SplitConcs(plObj, comp1, comp2, plColor, erKW=None):
    c1, c1_er = comp1;
    c2, c2_er = comp2;

    if erKW == None:
        erKW = {'capthick':2, 'capsize':2.5, 'elinewidth':2, 'alpha':0.7, 'markersize':10}

    plObj.errorbar(c1[1], c2[1], xerr=c1_er[0], yerr=c2_er[0], fmt = 'o', fillstyle='right',  color=plColor, **erKW)
    plObj.errorbar(c1[2], c2[2], xerr=c1_er[1], yerr=c2_er[1], fmt = '^', fillstyle='left',   color=plColor, **erKW)

def AddPlot_Ortho_TotLin_SplitConcs(plObj, this_dir, which_comp=0, bID_List=None, plColor='k', erKW=None):
    if bID_List == None:
        bID_List = range(7,22,2)
    for bID in bID_List:
        this_dat = this_dir[bID]
        plDat = this_dat[which_comp]
        comp1 = plDat[0]
        comp2 = plDat[1]
        AddPlot_SplitConcs(plObj, comp1, comp2, plColor, erKW)

def AddPlot_TotLin_SplitConcs(plObj, this_dir, which_comp=0, bID_List=None, sID_List=None, plColor='k', erKW=None):
    if bID_List == None:
        bID_List = [6,7,8,9,10]
    if sID_List == None:
        sID_List = range(7,18)
    for bID in bID_List:
        for sID in sID_List:
            this_dat = this_dir[bID][sID]
            plDat = this_dat[which_comp]
            comp1 = plDat[0]
            comp2 = plDat[1]
            AddPlot_SplitConcs(plObj, comp1, comp2, plColor, erKW)

def AddPlot_WTieLinesSplit(plObj, comp1, comp2, plColor, erKW=None, plKW=None):
    c1, c1_er = comp1;
    c2, c2_er = comp2;

    if erKW == None:
        erKW = {'capthick':2.5, 'capsize':2.5, 'elinewidth':2.5, 'alpha':0.7, 'markersize':10}
    if plKW == None:
        plKW = {'lw':3, 'markersize':10, 'alpha':0.7}


    plObj.errorbar(c1[1], c2[1], xerr=c1_er[0], yerr=c2_er[0], fmt='o', color=plColor, fillstyle='right', **erKW)
    plObj.errorbar(c1[2], c2[1], xerr=c1_er[1], yerr=c2_er[1], fmt='^', color=plColor, **erKW)

    plObj.plot([c1[1], c1[0]], [c2[1], c2[0]], '-', color=plColor, **plKW)
    plObj.plot([c1[0]], [c2[0]], 'o', color=plColor, **plKW)
    plObj.plot([c1[2], c1[0]], [c2[2], c2[0]], '-', color=plColor, **plKW)

def ForPlot_OneComp_GenData(this_box, eps_list, cDiffCut=0.2):
    relDiff = this_box[-1]
    upTO    = np.argwhere(relDiff < cDiffCut)
    if len(upTO) == 0:
        upTO = len(relDiff)
    else:
        upTO = upTO[0][0]
    cLo = this_box[3][:upTO]; cLoE = this_box[4][:upTO]
    cHi = this_box[1][:upTO]; cHiE = this_box[2][:upTO]
    return [cLo, cLoE, cHi, cHiE, eps_list[:upTO]]

def ForPlot_OneComp_GenLogData(this_box, eps_list, cDiffCut=0.2):
    relDiff = this_box[-1]
    upTO    = np.argwhere(relDiff < cDiffCut)
    if len(upTO) == 0:
        upTO = len(relDiff)
    else:
        upTO = upTO[0][0]
    cLo = this_box[3][:upTO]; cLoE = this_box[4][:upTO]/cLo/np.log(10)
    cHi = this_box[1][:upTO]; cHiE = this_box[2][:upTO]/cHi/np.log(10)
    return [np.log10(cLo), cLoE, np.log10(cHi), cHiE, eps_list[:upTO]]

def AddPlot_OneComp_Binodal(plObj, this_box, erKW=None):
    cLo = this_box[0]; cLoE = this_box[1];
    cHi = this_box[2]; cHiE = this_box[3];
    eps_list = this_box[-1]
    if erKW == None:
        erKW = {'fmt':'o-', 'capthick':2.5, 'capsize':2.5, 'elinewidth':2.5, 'alpha':0.7, 'markersize':10, 'lw':2}

    plObj.errorbar(cLo, eps_list, xerr=cLoE, fillstyle='left', color='b', **erKW)
    plObj.errorbar(cHi, eps_list, xerr=cHiE, fillstyle='right', color='r', **erKW)

class Generate_Linear_Molecule:
    def __init__(self):
        self.Blocks = []
        self.Linkers = []
        self.InterBlockLinkers = []
        self.MolStructure = []

    def Add_Block(self, bead_list, linker_list):
        assert len(bead_list) == (len(linker_list) + 1), \
            "number of bead and linkers does not generate a correct molecule!"
        self.Blocks.append(bead_list)
        self.Linkers.append(linker_list)

    def Add_InterBlockLinker(self, linker_len):
        self.InterBlockLinkers.append(linker_len)

    def Repeat_Block(self, bead_list, linker_list, inter_linker, repeat_num=5):
        assert len(bead_list) == (len(linker_list) + 1), \
            "number of bead and linkers does not generate a correct molecule!"
        for a_block in range(repeat_num):
            self.Blocks.append(bead_list)
            self.Linkers.append(linker_list)
            if a_block < repeat_num - 1:
                self.Add_InterBlockLinker(inter_linker)

    def Form_Structure(self):
        block_number = len(self.Blocks)
        inter_number = len(self.InterBlockLinkers)
        assert inter_number == (block_number - 1), \
            "number of blocks and inter-block-linkers does not generate a correct structure!"

        tot_linker_list = []
        for block_ID, a_block in enumerate(self.Linkers):
            if block_ID > 0:
                tot_linker_list.append(self.InterBlockLinkers[block_ID - 1])
            for lin_ID, a_lin in enumerate(a_block):
                tot_linker_list.append(a_lin)

        tot_bead_list = []
        for block_ID, a_block in enumerate(self.Blocks):
            for bead_ID, a_bead in enumerate(a_block):
                tot_bead_list.append(a_bead)
        num_linkers = len(tot_linker_list)
        num_beads = len(tot_bead_list)
        assert num_beads == (num_linkers + 1), 'the linker and bead numbers are not consistent!'
        tot_bead_id = 0;
        tot_struc = []
        for bd_ID, a_bead in enumerate(tot_bead_list):
            if bd_ID == 0:
                if num_beads == block_number == 1:
                    tot_struc.append([tot_bead_id, a_bead, -1, -1])
                else:
                    tot_struc.append([tot_bead_id, a_bead, tot_linker_list[bd_ID], tot_bead_id + 1])
                    tot_bead_id += 1
            elif bd_ID == num_beads - 1:
                tot_struc.append([tot_bead_id, a_bead, tot_linker_list[bd_ID - 1], tot_bead_id - 1])
                tot_bead_id += 1
            else:
                tot_struc.append([tot_bead_id, a_bead, tot_linker_list[bd_ID - 1], tot_bead_id - 1])
                tot_struc.append([tot_bead_id, a_bead, tot_linker_list[bd_ID], tot_bead_id + 1])
                tot_bead_id += 1
        self.MolStructure = np.array(tot_struc)

class Generate_Energy_File:
    def __init__(self, tot_stick_num):
        self.St_Num = tot_stick_num
        self.Ovlp_En = np.zeros((tot_stick_num, tot_stick_num))
        self.Cont_En = np.zeros((tot_stick_num, tot_stick_num))
        self.Sti_En = np.zeros((tot_stick_num, tot_stick_num))
        self.Tot_File = []

    def Add_OvlpEnergy_Between(self, energy, a_pair):
        comp1 = a_pair[0];
        comp2 = a_pair[1];
        self.Ovlp_En[comp1, comp2] = energy
        self.Ovlp_En[comp2, comp1] = energy

    def Add_ContEnergy_Between(self, energy, a_pair):
        comp1 = a_pair[0];
        comp2 = a_pair[1];
        self.Cont_En[comp1, comp2] = energy
        self.Cont_En[comp2, comp1] = energy

    def Add_StiEnergy_Between(self, energy, a_pair):
        comp1 = a_pair[0];
        comp2 = a_pair[1];
        self.Sti_En[comp1, comp2] = energy
        self.Sti_En[comp2, comp1] = energy

    def Add_Energy_Between(self, energy, a_pair, en_mode="S"):
        if en_mode == "O":
            self.Add_OvlpEnergy_Between(energy, a_pair)
        if en_mode == "C":
            self.Add_ContEnergy_Between(energy, a_pair)
        if en_mode == "S":
            self.Add_StiEnergy_Between(energy, a_pair)

    def Add_Energy_For(self, energy, comp1, en_mode="S"):
        for comp2 in range(self.St_Num):
            if comp1 != comp2:
                self.Add_Energy_Between(energy, comp1, comp2, en_mode)

    def Write_Matrix(self, this_mat):
        unique_ints = np.unique(this_mat)
        num_of_ints = len(unique_ints)
        write_mat = [];
        if num_of_ints == 1:
            write_mat.append(str(this_mat[0, 0]))
        else:
            for aRow in this_mat:
                this_line = []
                for aNum in aRow:
                    this_line.append(str(aNum))
                this_line = " ".join(this_line)
                write_mat.append(this_line)
        return write_mat

    def Form_Energy_File(self):
        tot_file = []
        tot_file.append("#STICKERS")
        tot_file.append(str(self.St_Num))
        tot_file.append("")

        tot_file.append("#OVERLAP_POT")
        for aline in self.Write_Matrix(self.Ovlp_En):
            tot_file.append(aline)
        tot_file.append("")

        tot_file.append("#CONTACT_POT")
        for aline in self.Write_Matrix(self.Cont_En):
            tot_file.append(aline)
        tot_file.append("")

        tot_file.append("#CONTACT_RAD")
        tot_file.append(str(0.0))
        tot_file.append("")

        tot_file.append("#SC_SC_POT")
        for aline in self.Write_Matrix(self.Sti_En):
            tot_file.append(aline)
        tot_file.append("")

        tot_file.append("#LINKER_LENGTH")
        tot_file.append(str(1.0))
        tot_file.append("")

        tot_file.append("#LINKER_SPRCON")
        tot_file.append(str(0.0))
        tot_file.append("")

        tot_file.append("#LINKER_EQLEN")
        tot_file.append(str(1.0))
        tot_file.append("")
        self.Tot_File = tot_file

    def Print_File(self):
        for a_line in self.Tot_File:
            print a_line

    def Write_File_To(self, path_to_file):
        with open(path_to_file, "w") as myFile:
            for aline in self.Tot_File:
                myFile.write(aline + "\n")


class Sim_Setup:
    def __init__(self, struc_Ar, lin_Ar, param_file):
        self.GlobalParamFile = param_file
        self.CurrentDir = os.getcwd() + '/'
        self.SysInfo = {}
        structures = struc_Ar
        if len(lin_Ar) != len(struc_Ar):
            print("Setting all linker lengths to 1.0 lattice units")
            self.linkers = ['1.0' for i in struc_Ar]
        else:
            linkers = lin_Ar
        #print("Initialized a LASSI setup with the following structures:")
        for aStr, aLin in zip(structures, linkers):
            self.SysInfo[aStr] = {}
            self.SysInfo[aStr]['Linker Length'] = aLin
            #print(" " * 2 + "{:}:= Linker length {:} lattice units.".format(aStr, aLin))
            self.SysInfo[aStr]['Structure'] = []
        self.AddParamFiles_ForAll(param_file)
        self.SetRunName_ForAll()
        self.SetNumberOfRuns_ForAll(2)

    def Set_SimulationPath(self, end_path_name):
        self.SimulationPath = end_path_name
        print("Simulations shall be done in dir: {:}".format(self.SimulationPath))

    def Set_QSUB_Command(self, qsub_command):
        self.QSubCommand = qsub_command
        self.QSubIter = 0;

    def Set_QSUB_Queues(self, WIntQueue, NoIntQueue):
        self.QSUB_WIntQ = WIntQueue
        self.QSUB_NoIntQ = NoIntQueue

    def Set_PathToLASSI(self, path_to_lassi):
        self.PathToLASSI = path_to_lassi

    def Reset_QSUB_Iter(self):
        self.QSubIter = 0;

    def Read_StrucFileFor(self, sysName):
        struc_file = self.SysInfo[sysName]['Structure File']
        dum_struc_tot = []
        with open(struc_file) as sFile:
            tot_file = sFile.readlines()
            tot_file = iter(tot_file)
            for a_line in tot_file:
                if a_line[0] == '#' or a_line[0] == '\n':
                    continue
                if a_line[:4] == 'NEW{':
                    dum_struc_in = []
                    a_line = tot_file.next()
                    num_mol = int(a_line[:-1])
                    a_line = tot_file.next()
                    while (a_line[:4] != "}END"):
                        a_line = a_line[:-1].split()
                        this_line = [int(aVal) for aVal in a_line]
                        a_line = tot_file.next()
                        dum_struc_in.append(this_line)
                    dum_struc_tot.append([num_mol, np.array(dum_struc_in)])
        return dum_struc_tot

    def AddNewSystem(self, sysName, lin_len):
        if sysName in self.SysInfo.keys():
            print("This structure name already exists! Doing nothing.")
        else:
            self.SysInfo[sysName] = {}
            self.SysInfo[sysName]['Linker Length'] = lin_len
            self.SysInfo[sysName]['Structure'] = []
            print("Added {:}:= Linker length {:} lattice units.".format(sysName, lin_len))
            self.AddParamFileTo(sysName, self.GlobalParamFile)
            self.SetRunNameFor(sysName, sysName)
            self.SetNumberOfRunsFor(sysName, 2)

    def AddStrucFileTo(self, sysName, file_name):
        try:
            self.SysInfo[sysName]['Structure File'] = file_name
            print("{:} has structure file: {:}".format(sysName, file_name))
            self.SysInfo[sysName]['Key File'][1]['STRUCT_FILE'] = self.CurrentDir + file_name
            self.SysInfo[sysName]['Structure'] = self.Read_StrucFileFor(sysName)
            self.CalcNonZeroRDFComps(sysName)
            self.SysInfo[sysName]['Tot Molecules'] = self.CalcTotalMolecules(sysName)
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def AddEnergyFileTo(self, sysName, fileName_NoInt, fileName_WInt):
        try:
            with open(self.CurrentDir + fileName_WInt) as WIntFile:
                pass
            with open(self.CurrentDir + fileName_NoInt) as NoIntFile:
                pass
        except IOError:
            print("Do these files exist?")
            return
        try:
            self.SysInfo[sysName]['Int Energy File'] = self.CurrentDir + fileName_WInt
            self.SysInfo[sysName]['NoInt Energy File'] = self.CurrentDir + fileName_NoInt
            #print("{:} has energy files: Int:{:}  NoInt:{:}".format(sysName, self.SysInfo[sysName]['Int Energy File'],
            #                                                        self.SysInfo[sysName]['NoInt Energy File']))
            with open(self.CurrentDir + fileName_WInt) as WIntFile:
                WIntFile.readline()  # First line is a comment
                tot_stickers = int(WIntFile.readline())
            self.SysInfo[sysName]['Tot Bead Types'] = tot_stickers
        except KeyError:
            print("Failed! Did you type the correct system name?")
            return

    def AddEnergyFiles_ForAll(self, fileName_NoInt, fileName_WInt):
        for aSys in self.SysInfo.keys():
            self.AddEnergyFileTo(aSys, fileName_NoInt, fileName_WInt)

    def AddParamFileTo(self, sysName, fileName):
        self.SysInfo[sysName]['Key File'] = Read_ParamFile(fileName)

    def AddParamFiles_ForAll(self, fileName):
        for aSys in self.SysInfo.keys():
            self.AddParamFileTo(aSys, fileName)

    def PrintParamsFor(self, sysName):
        try:
            dumKeys = self.SysInfo[sysName]['Key File']
            dumDict = dumKeys[1]
            dumKeys = dumKeys[0]
            for aKey in dumKeys:
                dum_spaces = 25 - len(aKey)
                print('{:}{:}{:}'.format(aKey, ' ' * dum_spaces, dumDict[aKey]))
        except KeyError as myErr:
            if myErr[0] == 'Key File':
                print("Did you import the file?")
            else:
                print("Failed! Did you type the correct system name?")

    def ResetStruc(self, sysName):
        self.SysInfo[sysName]['Structure'] = []

    def AddStruc_ImpLinear(self, sysName, mol_num, bead_num, bead_type, lin_len):
        try:
            self.SysInfo[sysName]['Structure'].append(sGen_Implicit_Linear(mol_num, bead_num, bead_type, lin_len))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_Monomer(self, sysName, mol_num, bead_type):
        try:
            self.SysInfo[sysName]['Structure'].append(sGen_Monomer(mol_num, bead_type))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_Dimer(self, sysName, mol_num, bead_type1, bead_type2, lin_len):
        try:
            self.SysInfo[sysName]['Structure'].append(sGen_Dimer(mol_num, bead_type1, bead_type2, lin_len))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_ExpLinear(self, sysName, mol_num, bead_num, bead_type, lin_len, lin_type):
        try:
            self.SysInfo[sysName]['Structure'].append(
                sGen_Explicit_Linear(mol_num, bead_num, bead_type, lin_len, lin_type))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_ImpBranSymm(self, sysName, mol_num, branch_num, bead_per_branch, bead_type, hub_type, lin_len, hub_len):
        try:
            self.SysInfo[sysName]['Structure'].append(
                sGen_Implicit_SymmBranched(mol_num, branch_num, bead_per_branch, bead_type, hub_type, lin_len, hub_len))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def Write_StructureFileFor(self, sysName, file_name):
        try:
            sGen_WriteStructuresToFile(self.SysInfo[sysName]['Structure'], file_name)
            self.AddStrucFileTo(sysName, file_name)
            self.CalcNonZeroRDFComps(sysName)
            self.SysInfo[sysName]['Tot Molecules'] = self.CalcTotalMolecules(sysName)
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def Write_StructureFiles_ForAll(self):
        for aSys in self.SysInfo.keys():
            sGen_WriteStructuresToFile(self.SysInfo[aSys]['Structure'], aSys + '_struc.prm')
            self.AddStrucFileTo(sysName, aSys + '_struc.prm')

    def SetTemperaturesFor(self, sysName, init_temp, final_temp, temp_steps, therm_temp):
        assert (type(temp_steps) == int) and (type(init_temp) == type(final_temp) == float)
        assert temp_steps > 0. and init_temp > 0. and final_temp > 0.
        delta_temp = (final_temp - init_temp) / (temp_steps - 1)
        try:
            key_dum = self.SysInfo[sysName]['Key File'][1]
            key_dum['MC_TEMP'] = init_temp
            key_dum['MC_CYCLE_NUM'] = temp_steps
            key_dum['MC_DELTA_TEMP'] = delta_temp
            key_dum['PREEQ_TEMP'] = therm_temp
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetTemperatures_Inverted_For(self, sysName, init_temp, delta_temp, temp_steps, therm_temp=1000.0):
        assert (type(temp_steps) == int) and (type(init_temp) == type(delta_temp) == float)
        assert temp_steps > 0 and init_temp > 0. and delta_temp*(temp_steps-1) + init_temp > 0.
        try:
            key_dum = self.SysInfo[sysName]['Key File'][1]
            key_dum['MC_TEMP'] = init_temp
            key_dum['MC_CYCLE_NUM'] = temp_steps
            key_dum['MC_DELTA_TEMP'] = delta_temp
            key_dum['MC_INVERT_TEMP'] = 1
            key_dum['PREEQ_TEMP'] = therm_temp
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetTemperatures_ForAll(self, init_temp, final_temp, temp_steps, therm_temp):
        for aSys in self.SysInfo.keys():
            self.SetTemperaturesFor(aSys, init_temp, final_temp, temp_steps, therm_temp)

    def SetMCStepsFor(self, sysName, therm_steps, run_steps):
        try:
            key_dum = self.SysInfo[sysName]['Key File'][1]
            key_dum['N_STEPS'] = int(run_steps)
            key_dum['PREEQ_STEPS'] = int(therm_steps)
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetMCSteps_ForAll(self, therm_steps, run_steps):
        for aSys in self.SysInfo.keys():
            self.SetMCStepsFor(aSys, therm_steps, run_steps)

    def CalcBoxSizeArr(self, low_con, high_con, tot_beads, tot_boxes):
        if not low_con < high_con < 0.0:
            print("Low conc should be lower than high conc!"
                  "Both numbers need to be negative. (10^low_con, 10^high_con).")
            return
        dum_li = np.linspace(low_con, high_con, tot_boxes)
        dum_ar = 10. ** dum_li
        dum_ar = tot_beads / dum_ar
        dum_ar = np.array(dum_ar ** (1. / 3.), dtype=int)
        dum_ar_s = np.sort(dum_ar)
        return dum_ar_s

    def CalcTotBeadsFor(self, sysName):
        try:
            struc_ar = self.SysInfo[sysName]['Structure']
            tot_beads = []
            for a_mol in struc_ar:
                num_mol = a_mol[0]
                struc_list = a_mol[1].T[0]
                num_beads_per = len(np.unique(struc_list))
                tot_beads.append(num_beads_per * num_mol)
            return np.sum(np.array(tot_beads))
        except KeyError:
            print("Failed! Did you type the correct system name?")
            return

    def CalcNonZeroRDFComps(self, sysName):
        try:
            struc_list = self.SysInfo[sysName]['Structure']
        except KeyError:
            print("Did you type the system name correctly?")
            return
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            return

        bead_type_list = []
        for a_mol in struc_list:
            bead_types = a_mol[1].T[1]
            bead_type_list.append(np.unique(bead_types))
        bead_type_list = [aType for a_mol in bead_type_list for aType in a_mol]
        bead_types = np.array(bead_type_list)
        bead_types = bead_types.flatten()
        bead_types = np.unique(bead_types)
        bead_type_list = [[-1, -1]]  # For the total
        for i in range(len(bead_types)):
            for j in range(i, len(bead_types)):
                bead_type_list.append([bead_types[i], bead_types[j]])
        self.SysInfo[sysName]['Comp List'] = np.array(bead_type_list, dtype=int)

    def CalcTotalMolecules(self, sysName):
        try:
            struc_ar = self.SysInfo[sysName]['Structure']
            tot_mols = []
            for a_mol in struc_ar:
                num_mol = a_mol[0]
                tot_mols.append(num_mol)
            return np.sum(np.array(tot_mols))
        except KeyError:
            print("Failed! Did you type the correct system name?")
            return

    def SetBoxSizesFor(self, sysName, low_con, high_con, tot_boxes):
        try:
            tot_beads = self.CalcTotBeadsFor(sysName)
            self.SysInfo[sysName]['Boxes'] = self.CalcBoxSizeArr(low_con, high_con, tot_beads, tot_boxes)
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetBoxSizes_ForAll(self, low_con, high_con, tot_boxes):
        for aSys in self.SysInfo.keys():
            self.SetBoxSizesFor(aSys, low_con, high_con, tot_boxes)

    def SetBoxSizes_To(self, sysName):
        try:
            tot_ar = self.SysInfo[sysName]['Boxes']
        except KeyError:
            print("Failed! Did you type the correct system name? "
                  "Or have you setup the boxes for this system?")
            return
        for aSys in self.SysInfo.keys():
            if aSys != sysName:
                self.SysInfo[aSys]['Boxes'] = tot_ar

    def SetRunNameFor(self, sysName, runName):
        try:
            self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX'] = runName
        except KeyError:
            print("Failed! Did you type the correct system name?"
                  " Does this system exist?")

    def SetRunName_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.SetRunNameFor(aSys, aSys)

    def SetNumberOfRunsFor(self, sysName, run_num):
        try:
            self.SysInfo[sysName]['Runs'] = run_num
        except KeyError:
            print("Failed! Did you type the correct system name?"
                  " Does this system exist?")

    def SetNumberOfRuns_ForAll(self, run_num):
        for aSys in self.SysInfo.keys():
            self.SetNumberOfRunsFor(aSys, run_num)

    def MakeDirs_NoIntFor(self, sysName):
        try:
            dum_dir = self.SimulationPath
            MKDirCatch(dum_dir)
        except AttributeError:
            print("You have not set the simulation path yet!")
            return
        try:
            self.SysInfo[sysName]
            dum_dir += sysName + '/'
            MKDirCatch(dum_dir)
        except KeyError:
            print("Did you write the correct system, "
                  "and does it exist?")
            return
        dum_dir += str(self.SysInfo[sysName]['Linker Length']) + '/'
        MKDirCatch(dum_dir)
        try:
            for a_box in self.SysInfo[sysName]['Boxes']:
                this_dir = dum_dir + str(a_box) + '/'
                MKDirCatch(this_dir)
                this_dir += 'NoInt/'
                MKDirCatch(this_dir)
                for a_run in range(self.SysInfo[sysName]['Runs']):
                    final_dir = this_dir + str(a_run + 1) + '/'
                    MKDirCatch(final_dir)
        except KeyError:
            print("Have you setup the boxes yet?")

    def MakeDirs_WIntFor(self, sysName):
        try:
            dum_dir = self.SimulationPath
            MKDirCatch(dum_dir)
        except AttributeError:
            print("You have not set the simulation path yet!")
            return
        try:
            self.SysInfo[sysName]
            dum_dir += sysName + '/'
            MKDirCatch(dum_dir)
        except KeyError:
            print("Did you write the correct system, "
                  "and does it exist?")
            return
        dum_dir += str(self.SysInfo[sysName]['Linker Length']) + '/'
        MKDirCatch(dum_dir)
        try:
            for a_box in self.SysInfo[sysName]['Boxes']:
                this_dir = dum_dir + str(a_box) + '/'
                MKDirCatch(this_dir)
                this_dir += 'WInt/'
                MKDirCatch(this_dir)
                for a_run in range(self.SysInfo[sysName]['Runs']):
                    final_dir = this_dir + str(a_run + 1) + '/'
                    MKDirCatch(final_dir)
        except KeyError:
            print("Have you setup the boxes yet?")

    def MakeDirs_For(self, sysName):
        self.MakeDirs_WIntFor(sysName)
        self.MakeDirs_NoIntFor(sysName)

    def MakeDirs_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.MakeDirs_For(aSys)

    def MakeDirs_ForAll_NoInt(self):
        for aSys in self.SysInfo.keys():
            self.MakeDirs_NoIntFor(aSys)

    def MakeDirs_ForAll_WInt(self):
        for aSys in self.SysInfo.keys():
            self.MakeDirs_WIntFor(aSys)

    def Write_ParamsWIntFor(self, sysName):
        try:
            dum_dir = self.SimulationPath
        except AttributeError:
            print("You have not set the simulation path yet!")
            return
        try:
            self.SysInfo[sysName]
            dum_dir += sysName + '/'
            param_copy = self.SysInfo[sysName]['Key File'][:]
        except KeyError:
            print("Did you write the correct system, "
                  "and does it exist?")
            return
        dum_dir += str(self.SysInfo[sysName]['Linker Length']) + '/'
        try:
            for a_box in self.SysInfo[sysName]['Boxes']:
                param_copy[1]['BOX_SIZE'] = a_box
                this_dir = dum_dir + str(a_box) + '/'
                this_dir += 'WInt/'
                param_copy[1]['ENERGY_FILE'] = self.SysInfo[sysName]['Int Energy File']
                for a_run in range(self.SysInfo[sysName]['Runs']):
                    final_dir = this_dir + str(a_run + 1) + '/'
                    param_copy[1]['RANDOM_SEED'] = 0
                    with open(final_dir + 'param.key', 'w+') as pFile:
                        for a_key in param_copy[0]:
                            N_spcs = 25 - len(a_key)
                            # print(type(param_copy[1][a_key]))
                            pFile.write(a_key + ' ' * N_spcs + str(param_copy[1][a_key]) + '\n')
        except KeyError:
            print("Have you setup the boxes yet?")

    def Write_ParamsNoIntFor(self, sysName):
        try:
            dum_dir = self.SimulationPath
        except AttributeError:
            print("You have not set the simulation path yet!")
            return
        try:
            self.SysInfo[sysName]
            dum_dir += sysName + '/'
            param_copy = self.SysInfo[sysName]['Key File'][:]
        except KeyError:
            print("Did you write the correct system, "
                  "and does it exist?")
            return
        dum_dir += str(self.SysInfo[sysName]['Linker Length']) + '/'
        try:
            for a_box in self.SysInfo[sysName]['Boxes']:
                param_copy[1]['BOX_SIZE'] = a_box
                this_dir = dum_dir + str(a_box) + '/'
                this_dir += 'NoInt/'
                param_copy[1]['ENERGY_FILE'] = self.SysInfo[sysName]['NoInt Energy File']
                param_copy[1]['MV_STROT_FREQ'] = '0.0';
                param_copy[1]['MV_COLOCAL_FREQ'] = '0.0';
                param_copy[1]['MV_SMCLSTR_FREQ'] = '0.0';
                param_copy[1]['MV_CLSTR_FREQ'] = '0.0';
                param_copy[1]['REPORT_NETWORK_FREQ'] = '0';
                for a_run in range(self.SysInfo[sysName]['Runs']):
                    final_dir = this_dir + str(a_run + 1) + '/'
                    param_copy[1]['RANDOM_SEED'] = 0
                    with open(final_dir + 'param.key', 'w+') as pFile:
                        for a_key in param_copy[0]:
                            N_spcs = 25 - len(a_key)
                            pFile.write(a_key + ' ' * N_spcs + str(param_copy[1][a_key]) + '\n')
        except KeyError:
            print("Have you setup the boxes yet?")

    def Write_ParamsFor(self, sysName):
        self.Write_ParamsWIntFor(sysName)
        self.Write_ParamsNoIntFor(sysName)

    def Write_ParamsWInt_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.Write_ParamsWIntFor(aSys)

    def Write_ParamsNoInt_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.Write_ParamsNoIntFor(aSys)

    def Write_Params_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.Write_ParamsWIntFor(aSys)
            self.Write_ParamsNoIntFor(aSys)

    def SubmitWIntJobs_ToQueueFor(self, sysName, WIntQueue):
        try:
            dum_dir = self.SimulationPath
        except AttributeError:
            print("You have not set the simulation path yet!")
            return
        try:
            self.SysInfo[sysName]
            dum_dir += sysName + '/'
        except KeyError:
            print("Did you write the correct system, "
                  "and does it exist?")
            return
        dum_dir += str(self.SysInfo[sysName]['Linker Length']) + '/'
        try:
            for b_ID, a_box in enumerate(self.SysInfo[sysName]['Boxes']):
                this_dir = dum_dir + str(a_box) + '/'
                this_dir += 'WInt/'
                for a_run in range(self.SysInfo[sysName]['Runs']):
                    final_dir = this_dir + str(a_run + 1) + '/'
                    run_dir = final_dir
                    run_num = ['W' + sysName[:2], b_ID, self.QSubIter]
                    run_num = [str(a_val) for a_val in run_num]
                    run_num = "_".join(run_num)
                    run_command = self.QSubCommand + ' ' + WIntQueue + ' ' + run_num + ' /.' + self.PathToLASSI
                    try:
                        os.chdir(run_dir)
                    except IOError:
                        print("Did you make all the directories?")
                        return
                    try:
                        ret_code = sproc.check_output(run_command, shell=True, stderr=sproc.STDOUT)
                        if ret_code < 0:
                            print("The submission failed! Signal: ", -retcode)
                        else:
                            print(ret_code)
                    except OSError as myErr:
                        print("Couldn't submit because ", myErr)
                        raise
                    time.sleep(2)
                    self.QSubIter += 1
        except KeyError:
            print("Have you setup the boxes yet?")
        os.chdir(self.CurrentDir)

    def SubmitNoIntJobs_ToQueueFor(self, sysName, NoIntQueue):
        try:
            dum_dir = self.SimulationPath
        except AttributeError:
            print("You have not set the simulation path yet!")
            return
        try:
            self.SysInfo[sysName]
            dum_dir += sysName + '/'
        except KeyError:
            print("Did you write the correct system, "
                  "and does it exist?")
            return
        dum_dir += str(self.SysInfo[sysName]['Linker Length']) + '/'
        try:
            for b_ID, a_box in enumerate(self.SysInfo[sysName]['Boxes']):
                this_dir = dum_dir + str(a_box) + '/'
                this_dir += 'NoInt/'
                for a_run in range(self.SysInfo[sysName]['Runs']):
                    final_dir = this_dir + str(a_run + 1) + '/'
                    run_dir = final_dir
                    run_num = ['N' + sysName[:2], b_ID, self.QSubIter]
                    run_num = [str(a_val) for a_val in run_num]
                    run_num = "_".join(run_num)
                    run_command = self.QSubCommand + ' ' + NoIntQueue + ' ' + run_num + ' /.' + self.PathToLASSI
                    try:
                        os.chdir(run_dir)
                    except IOError:
                        print("Did you make all the directories?")
                        return
                    try:
                        ret_code = sproc.check_output(run_command, shell=True, stderr=sproc.STDOUT)
                        if ret_code < 0:
                            print("The submission failed! Signal: ", -retcode)
                        else:
                            print(ret_code)
                    except OSError as myErr:
                        print("Couldn't submit because ", myErr)
                        raise
                    time.sleep(2)
                    self.QSubIter += 1
        except KeyError:
            print("Have you setup the boxes yet?")
        os.chdir(self.CurrentDir)

    def SubmitWIntJobs_ForAll(self, WIntQueue):
        for aSys in self.SysInfo.keys():
            self.SubmitWIntJobs_ToQueueFor(aSys, WIntQueue)

    def SubmitNoIntJobs_ForAll(self, NoIntQueue):
        for aSys in self.SysInfo.keys():
            self.SubmitNoIntJobs_ToQueueFor(aSys, NoIntQueue)

    def SubmitJobs_ForAll(self):
        try:
            self.QSubCommand
            self.QSUB_WIntQ
            self.QSUB_NoIntQ
        except AttributeError:
            print("You need to set up the qsub command, and which queues to use!")
            return
        self.Reset_QSUB_Iter()
        for aSys in self.SysInfo.keys():
            self.SubmitNoIntJobs_ToQueueFor(aSys, self.QSUB_NoIntQ)
            self.SubmitWIntJobs_ToQueueFor(aSys, self.QSUB_WIntQ)

    def SetAndRun_For(self, sysName):
        self.MakeDirs_For(sysName)
        self.Write_ParamsFor(sysName)
        self.Reset_QSUB_Iter()
        self.SubmitWIntJobs_ToQueueFor(sysName, self.QSUB_WIntQ)
        self.SubmitWIntJobs_ToQueueFor(sysName, self.QSUB_NoIntQ)

class Analysis(Sim_Setup):
    def __init__(self, SimSetupInstance):
        try:
            self.SysInfo = SimSetupInstance.SysInfo
            self.GlobalParamFile = SimSetupInstance.GlobalParamFile
            self.SimulationPath = SimSetupInstance.SimulationPath
            self.CurrentDir = SimSetupInstance.CurrentDir
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise

    def Read_PDFFileFor(self, sysName, file_name, box_size, tot_temps):
        try:
            dum_dat = np.loadtxt(file_name)
            comp_list = self.SysInfo[sysName]['Comp List']
            tot_comps = self.SysInfo[sysName]['Tot Bead Types']
            tot_comps_pos = 1 + Index_RDF(tot_comps - 2, tot_comps - 1, tot_comps)
        except IOError:
            print("Does the RDF file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        actual_comps = len(comp_list)
        ret_dat = np.zeros((actual_comps, tot_temps, box_size * 4))
        for a_temp in range(tot_temps):
            for compID, compPair in enumerate(comp_list):
                compA = int(compPair[0]);
                compB = int(compPair[1])
                ret_dat[compID][a_temp] = dum_dat[tot_comps_pos * a_temp + Index_RDF(compA, compB, tot_comps)]
        return ret_dat

    def Read_CLUSFileFor(self, file_name):
        try:
            dum_dat = np.loadtxt(file_name)
        except IOError:
            print("Does the RDF file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Read_COMDenFileFor(self, file_name):
        try:
            dum_dat = np.loadtxt(file_name)
        except IOError:
            print("Does the RDF file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Read_MolClusFileFor(self, file_name):
        try:
            dum_dat = np.loadtxt(file_name)
        except IOError:
            print("Does the RDF file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Read_GYRADFileFor(self, file_name):
        try:
            dum_dat = np.loadtxt(file_name)
        except IOError:
            print("Does the RDF file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Save_NoIntPDF_For(self, sysName):
        try:
            dum_dat = self.Collect_NoInt_PDFs_For(sysName)
            file_name = sysName + '_N_PDF.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_NoInt_PDF'] = file_name
            pass
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntPDF_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_PDFs_For(sysName)
            file_name = sysName + '_W_PDF.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_WInt_PDF'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntClus_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_CLUS_For(sysName)
            file_name = sysName + '_CLUS.b'
            np.save(file_name, dum_dat)
            self.SysInfo[sysName]['Raw_CLUS'] = file_name + '.npy'
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntCOMDens_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_COMDens_For(sysName)
            file_name = sysName + '_W_COMDen.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_WInt_COMDen'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_NoIntCOMDens_For(self, sysName):
        try:
            dum_dat = self.Collect_NoInt_COMDens_For(sysName)
            file_name = sysName + '_N_COMDen.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_NoInt_COMDen'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntMolClus_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_MolClus_For(sysName)
            file_name = sysName + '_W_MolClus.b'
            np.save(file_name, dum_dat)
            self.SysInfo[sysName]['Raw_WInt_MolClus'] = file_name + '.npy'
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_GyrRad_For(self, sysName):
        try:
            dum_dat = self.Collect_GYRAD_For(sysName)
            file_name = sysName + '_GyrRad.b'
            np.save(file_name, dum_dat)
            self.SysInfo[sysName]['Raw_GyrRad'] = file_name + '.npy'
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Set_PDFFileNames_For(self, sysName):
        self.SysInfo[sysName]['Raw_NoInt_PDF'] = sysName + '_N_PDF.b'
        self.SysInfo[sysName]['Raw_WInt_PDF'] = sysName + '_W_PDF.b'

    def Collect_NoInt_PDFs_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir + str(a_box) + '/NoInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir + str(a_run + 1) + '/' + run_name + '_RDF.dat'
                tot_raw_dat_temp.append(self.Read_PDFFileFor(sysName, file_name, a_box, tot_temps))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat

    def Collect_WInt_PDFs_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir + str(a_box) + '/WInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir + str(a_run + 1) + '/' + run_name + '_RDF.dat'
                tot_raw_dat_temp.append(self.Read_PDFFileFor(sysName, file_name, a_box, tot_temps))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat

    def Collect_WInt_CLUS_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir + str(a_box) + '/WInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir + str(a_run + 1) + '/' + run_name + '_CLUS.dat'
                tot_raw_dat_temp.append(self.Read_CLUSFileFor(file_name))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat

    def Collect_WInt_COMDens_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir + str(a_box) + '/WInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir + str(a_run + 1) + '/' + run_name + '_COMDen.dat'
                tot_raw_dat_temp.append(self.Read_COMDenFileFor(file_name))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat

    def Collect_NoInt_COMDens_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir + str(a_box) + '/NoInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir + str(a_run + 1) + '/' + run_name + '_COMDen.dat'
                tot_raw_dat_temp.append(self.Read_COMDenFileFor(file_name))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat

    def Collect_WInt_MolClus_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir + str(a_box) + '/WInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir + str(a_run + 1) + '/' + run_name + '_MolClus.dat'
                tot_raw_dat_temp.append(self.Read_MolClusFileFor(file_name))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat

    def Collect_GYRAD_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir + str(a_box) + '/WInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir + str(a_run + 1) + '/' + run_name + '_GR.dat'
                tot_raw_dat_temp.append(self.Read_GYRADFileFor(file_name))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat

    def Set_OrdParamFileNames_For(self, sysName):
        self.SysInfo[sysName]['Perc File'] = sysName + '_phic.c.npy'
        self.SysInfo[sysName]['RhoBar File'] = sysName + '_rhobar.c.npy'

    def Gen_pRDFs_For(self, sysName, my_WInt=None, my_NoInt=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            comp_list = this_sys['Comp List']
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
            self.Set_PDFFileNames_For(sysName)
        except KeyError:
            print("Is this the correct system?")
            return
        if my_WInt != None:
            this_sys['Raw_WInt_PDF'] = my_WInt
        if my_NoInt != None:
            this_sys['Raw_NoInt_PDF'] = my_NoInt
        nFile = this_sys['Raw_NoInt_PDF']
        iFile = this_sys['Raw_WInt_PDF']
        with open(iFile) as dumFile:
            i_pdfs = pickle.load(dumFile)
        with open(nFile) as dumFile:
            n_pdfs = pickle.load(dumFile)

        tot_rdf_arr = []
        for boxID, a_box in enumerate(box_ar):
            x_ar = np.arange(0., a_box, 0.25)
            n_ar = np.array(n_pdfs[boxID])
            i_ar = np.array(i_pdfs[boxID])
            n_run_ave, n_run_err = MeanAndError(n_ar, 2);
            n_temp_ave, n_temp_err = MeanAndError(n_run_ave, 0)
            i_run_ave, i_run_err = MeanAndError(i_ar, 0);
            zero_points = np.argwhere(n_temp_ave[0] == 0.).T[0]
            x_ar = np.delete(x_ar, zero_points)
            rdf_box_ar = []
            for aComp, compPair in enumerate(comp_list):
                n_dum = n_temp_ave[aComp]
                n_dum = np.delete(n_dum, zero_points)
                rdf_temp_ar = []
                for a_temp in range(tot_temps):
                    i_dum = i_run_ave[aComp][a_temp]
                    i_dum = np.delete(i_dum, zero_points)
                    p_rdf_ar = i_dum / n_dum
                    rdf_temp_ar.append(p_rdf_ar)
                rdf_box_ar.append(rdf_temp_ar)
            tot_rdf_arr.append([x_ar, rdf_box_ar])
        rFile = sysName + '_CalcRDF.b'
        this_sys['RDF File'] = rFile
        with open(rFile, "wb+") as rdf_file:
            pickle.dump(tot_rdf_arr, rdf_file)

    def Gen_RhoBar_For(self, sysName, my_rdfFile=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            comp_list = this_sys['Comp List']
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
            this_sys['RDF File'] = sysName + '_CalcRDF.b'
        except KeyError:
            print("Is this the correct system?")
            return
        if my_rdfFile != None:
            this_sys['RDF File'] = my_rdfFile
        rdf_file = this_sys['RDF File']
        with open(rdf_file) as rFile:
            tot_rdfs = pickle.load(rFile)

        rho_bar_mat = []
        for box_size, a_box in zip(box_ar, tot_rdfs):
            x_ar = a_box[0]
            pRDFs = a_box[1]
            rho_bar_comp = []
            for aComp in pRDFs:
                rho_bar_temp = []
                for a_temp in aComp:
                    g_of_r = a_temp;
                    my_func = np.abs(g_of_r - 1.)
                    volume_norm = RDFVolumeElement(x_ar, box_size)
                    my_func_norm = my_func * volume_norm
                    rho_bar = np.trapz(my_func_norm, x_ar) / (box_size ** 3.)
                    # rho_bar = np.trapz(g_of_r*my_func_norm,x_ar)/np.trapz(g_of_r*volume_norm,x_ar)
                    rho_bar_temp.append(rho_bar)
                rho_bar_comp.append(rho_bar_temp)
            rho_bar_mat.append(rho_bar_comp)
        rho_bar_mat = np.array(rho_bar_mat)
        rhoFile = sysName + '_rhobar.c'
        np.save(rhoFile, rho_bar_mat)
        rhoFile += '.npy'
        this_sys['RhoBar File'] = rhoFile

    def Gen_PhiC_For(self, sysName, my_clusFile=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            mol_num = this_sys['Tot Molecules']
            this_sys['Raw_CLUS'] = sysName + '_CLUS.b.npy'
        except KeyError:
            print("Is this the correct system?")
            return
        if my_clusFile != None:
            this_sys['Raw_CLUS'] = my_clusFile
        clus_file = this_sys['Raw_CLUS']
        tot_clus = np.load(clus_file)
        phi_c_mat = []
        for a_box in tot_clus:
            phi_run_avg = np.mean(a_box, axis=0)
            phi_run_err = np.std(a_box, axis=0)
            phi_c_temp = []
            for a_temp_ave, a_temp_err in zip(phi_run_avg, phi_run_err):
                phi_c_val = a_temp_ave[0]
                phi_c_err = a_temp_err[0]
                phi_c_temp.append([phi_c_val, phi_c_err])
            phi_c_mat.append(phi_c_temp)
        phi_c_mat = np.array(phi_c_mat)
        phiFile = sysName + '_phic.c'
        np.save(phiFile, phi_c_mat)
        phiFile += '.npy'
        this_sys['Perc File'] = phiFile

    def Gen_GyRad_For(self, sysName, my_gyradFile=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            mol_num = this_sys['Tot Molecules']
            this_sys['Raw_GyrRad'] = sysName + '_GyrRad.b.npy'
        except KeyError:
            print("Is this the correct system?")
            return
        if my_gyradFile != None:
            this_sys['Raw_GyrRad'] = my_gyradFile
        gyr_file = this_sys['Raw_GyrRad']
        tot_gyr = np.load(gyr_file)
        gyr_mat = []
        for a_box in tot_gyr:
            gyr_run_avg = np.mean(a_box, axis=0)
            gyr_run_err = np.std(a_box, axis=0)
            gyr_temp = []
            for a_temp_ave, a_temp_err in zip(gyr_run_avg, gyr_run_err):
                gyr_val = a_temp_ave[1]/a_temp_ave[0]
                gyr_err = a_temp_ave[1]*a_temp_err[0]/(gyr_val**2.)
                gyr_temp.append([gyr_val, gyr_err])
            gyr_mat.append(gyr_temp)
        gyr_mat = np.array(gyr_mat)
        gyrFile = sysName + '_gyr.c'
        np.save(gyrFile, gyr_mat)
        gyrFile += '.npy'
        this_sys['GYR File'] = gyrFile

    def Gen_CorrectDen_For(self, sysName, my_WInt=None, my_NoInt=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            tot_mol_types = len(this_sys['Structure'])
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
            comden_comps = tot_mol_types*(tot_mol_types+1)
        except KeyError:
            print("Is this the correct system?")
            return
        if my_WInt != None:
            this_sys['Raw_WInt_COMDen'] = my_WInt
        else:
            this_sys['Raw_WInt_COMDen'] = sysName + '_W_COMDen.b'
        if my_NoInt != None:
            this_sys['Raw_NoInt_COMDen'] = my_NoInt
        else:
            this_sys['Raw_NoInt_COMDen'] = sysName + '_N_COMDen.b'

        nFile = this_sys['Raw_NoInt_COMDen']
        iFile = this_sys['Raw_WInt_COMDen']

        with open(iFile) as dumFile:
            i_pdfs = pickle.load(dumFile)
        with open(nFile) as dumFile:
            n_pdfs = pickle.load(dumFile)

        per_box_den = []
        for boxID, a_box in enumerate(box_ar):
            x_ar = np.arange(0., a_box, 0.25)
            n_ar = np.array(n_pdfs[boxID])
            i_ar = np.array(i_pdfs[boxID])
            i_run_ave, i_run_err = MeanAndError(i_ar, 0);
            n_run_ave, n_run_err = MeanAndError(n_ar, 0);

            per_comp_den = []
            for compID in range(comden_comps):
                comp_IAr = i_run_ave[compID::comden_comps]
                comp_NAr = n_run_ave[compID::comden_comps] #Can average over all temperature for NAr
                n_temp_ave, n_temp_err = MeanAndError(comp_NAr, 0)
                zero_points = np.argwhere(n_temp_ave == 0).T[0]
                comp_xAr = np.delete(x_ar, zero_points)
                dum_nAr = np.delete(n_temp_ave, zero_points)
                per_temp_den = []
                for tempID, aTemp in enumerate(comp_IAr):
                    dum_iAr = np.delete(aTemp, zero_points)
                    per_temp_den.append(dum_iAr)
                per_comp_den.append([comp_xAr, dum_nAr, per_temp_den])
            per_box_den.append(per_comp_den)


        rFile = sysName + '_CorrectDen.b'
        this_sys['CorrectDen File'] = rFile
        with open(rFile, "wb+") as den_file:
            pickle.dump(per_box_den, den_file)

    def Gen_TieLineData_For(self, sysName, my_corrDenFile=None, idx_upto=13, idx_upfrom=-13, idx_uptill=-3):
        this_sys = self.SysInfo[sysName]
        if my_corrDenFile != None:
            this_sys['CorrectDen File'] = my_corrDenFile
        else:
            this_sys['CorrectDen File'] = sysName + '_CorrectDen.b'
        dFile = this_sys['CorrectDen File'];
        with open(dFile) as dumFile:
            tot_den_data = pickle.load(dumFile)


        box_ar        = this_sys['Boxes']
        tot_struc     = this_sys['Structure']
        tot_mol_types = len(tot_struc)
        tot_poss_den  = tot_mol_types*(tot_mol_types+1)
        tot_temps     = this_sys['Key File'][1]['MC_CYCLE_NUM']
        tot_num_ar = []
        for additional_mol in tot_struc[:]:
            this_num = additional_mol[0]
            this_num = this_num
            tot_num_ar.append(this_num)
        tot_num_ar = np.array(tot_num_ar).T

        #The tie-line data are organized as follows:
        # Temperature, BoxSize, Molecule Number, Component, and then
        # [highCon, highCon_std, lowCon, lowCon_std, errSum/deltaCon]

        tie_line_data = np.zeros((len(box_ar), tot_poss_den, tot_temps, 6))
        for bID, aBox in enumerate(box_ar):
            inverse_vol = 1./(aBox**3.)
            for compID in range(tot_poss_den):
                this_comp_den = tot_num_ar[compID % tot_mol_types]*inverse_vol
                #this_comp_xAr = tot_den_data[bID][compID][0][1:]
                this_comp_nAr = tot_den_data[bID][compID][1][1:]
                for tempID in range(tot_temps):
                    this_comp_iAr = tot_den_data[bID][compID][2][tempID][1:]
                    this_comp_yAr = this_comp_den*this_comp_iAr/this_comp_nAr
                    highCon, highConEr = MeanAndError(this_comp_yAr[:idx_upto], 0)
                    lowCon, lowConEr   = MeanAndError(this_comp_yAr[idx_upfrom:idx_uptill], 0)
                    deltaCon           = np.abs(highCon-lowCon)
                    totErr             = np.sqrt(highConEr**2. + lowConEr**2.)

                    tie_line_data[bID][compID][tempID][0] = this_comp_den
                    tie_line_data[bID][compID][tempID][1] = highCon
                    tie_line_data[bID][compID][tempID][2] = highConEr
                    tie_line_data[bID][compID][tempID][3] = lowCon
                    tie_line_data[bID][compID][tempID][4] = lowConEr
                    tie_line_data[bID][compID][tempID][5] = totErr/deltaCon
        tieFile = sysName + '_TieLines.c'
        np.save(tieFile, tie_line_data)
        tieFile += '.npy'
        this_sys['TieLines File'] = tieFile

    def Gen_Binodal_For(self, sysName, my_TieFile=None):
        this_sys = self.SysInfo[sysName]
        if my_TieFile != None:
            this_sys['TieLines File'] = my_TieFile
        else:
            this_sys['TieLines File'] = sysName + '_TieLines.c.npy'
        dFile = this_sys['TieLines File'];
        tot_tie_data = np.load(dFile)

        box_ar        = this_sys['Boxes']
        tot_boxes     = len(box_ar)
        tot_struc     = this_sys['Structure']
        tot_mol_types = len(tot_struc)
        tot_temps     = this_sys['Key File'][1]['MC_CYCLE_NUM']
        #The binodal data are organized as follows:
        # Temperature, BoxSize, Molecule Number, Component, and then
        # [highCon, highCon_std, lowCon, lowCon_std, errSum/deltaCon]
        # This function adds up all the components in the system to have one master binodal

        binodal_data = np.zeros((tot_boxes, 6, tot_temps))
        for bxID in range(tot_boxes):
            this_box = tot_tie_data[bxID][:tot_mol_types].T
            dum_box = np.ones((6, tot_temps))
            for lID in [0, 1, 3]:
                dum_box[lID] = np.sum(this_box[lID], axis=1)
            for lID in [2, 4]:
                dum_box[lID] = np.sqrt(np.sum(this_box[lID] ** 2., axis=1))
            dum_box[5] = np.abs(np.log10(dum_box[1]) - np.log10(dum_box[3] + 1e-12))
            binodal_data[bxID] = dum_box

        binFile = sysName + '_Binodal.c'
        np.save(binFile, binodal_data)
        binFile += '.npy'
        this_sys['Binodal File'] = binFile

    def Set_PRDFsFileName_For(self, sysName):
        self.SysInfo[sysName]['RDF File'] = sysName + '_CalcRDF.b'

    def Gen_RhoBarInterp_For(self, sysName, of_comp=0, temp_scale=2., inter_points=50):
        try:
            this_sys = self.SysInfo[sysName]
            key_dum = this_sys['Key File'][1]
        except KeyError:
            print("Is this the correct system?")

        init_temp = key_dum['MC_TEMP']
        temp_steps = key_dum['MC_CYCLE_NUM']
        delta_temp = key_dum['MC_DELTA_TEMP']
        mol_num = this_sys['Tot Molecules']
        box_ar = [[a_box for a_box in this_sys['Boxes']]] * temp_steps
        temp_ar = [[init_temp + i * delta_temp for i in range(temp_steps)]] * len(box_ar[0])
        box_ar = np.array(box_ar)
        temp_ar = np.array(temp_ar) / temp_scale
        conc_ar = np.log10(mol_num / (box_ar ** 3.))
        conc_ar = conc_ar.flatten()
        temp_ar = temp_ar.T.flatten()
        conc_grid = np.linspace(conc_ar.min(), conc_ar.max(), inter_points)
        temp_grid = np.linspace(temp_ar.min(), temp_ar.max(), inter_points)
        conc_grid, temp_grid = np.meshgrid(conc_grid, temp_grid)
        this_sys['RhoBar File'] = sysName + '_rhobar.c.npy'
        rho_barFile = this_sys['RhoBar File']
        rho_bar = np.load(rho_barFile)
        this_rho = rho_bar[:, of_comp, :].T
        this_rho = this_rho.flatten()

        rho_func = sp.interpolate.Rbf(conc_ar, temp_ar, this_rho, function='linear')
        rho_grid = rho_func(conc_grid, temp_grid)
        return [conc_grid, temp_grid, rho_grid]

    def Gen_PhiCInterp_For(self, sysName, temp_scale=2., inter_points=50):
        try:
            this_sys = self.SysInfo[sysName]
            key_dum = this_sys['Key File'][1]
        except KeyError:
            print("Is this the correct system?")

        init_temp = key_dum['MC_TEMP']
        temp_steps = key_dum['MC_CYCLE_NUM']
        delta_temp = key_dum['MC_DELTA_TEMP']
        mol_num = this_sys['Tot Molecules']
        box_ar = [[a_box for a_box in this_sys['Boxes']]] * temp_steps
        temp_ar = [[init_temp + i * delta_temp for i in range(temp_steps)]] * len(box_ar[0])
        box_ar = np.array(box_ar)
        temp_ar = np.array(temp_ar) / temp_scale
        conc_ar = np.log10(mol_num / (box_ar ** 3.))
        conc_ar = conc_ar.flatten()
        temp_ar = temp_ar.T.flatten()
        conc_grid = np.linspace(conc_ar.min(), conc_ar.max(), inter_points)
        temp_grid = np.linspace(temp_ar.min(), temp_ar.max(), inter_points)
        conc_grid, temp_grid = np.meshgrid(conc_grid, temp_grid)
        this_sys['Perc File'] = sysName + '_phic.c.npy'
        phi_cFile = this_sys['Perc File']
        phi_c = np.load(phi_cFile)
        phi_c = phi_c[:, :, 0].T
        phi_c = phi_c.flatten() / mol_num

        phi_func = sp.interpolate.Rbf(conc_ar, temp_ar, phi_c, function='linear')
        phi_grid = phi_func(conc_grid, temp_grid)
        return [conc_grid, temp_grid, phi_grid]

    def Gen_GyrRadInterp_For(self, sysName, temp_scale=2., inter_points=50):
        try:
            this_sys = self.SysInfo[sysName]
            key_dum = this_sys['Key File'][1]
        except KeyError:
            print("Is this the correct system?")

        init_temp = key_dum['MC_TEMP']
        temp_steps = key_dum['MC_CYCLE_NUM']
        delta_temp = key_dum['MC_DELTA_TEMP']
        mol_num = this_sys['Tot Molecules']
        box_ar = [[a_box for a_box in this_sys['Boxes']]] * temp_steps
        temp_ar = [[init_temp + i * delta_temp for i in range(temp_steps)]] * len(box_ar[0])
        box_ar = np.array(box_ar)
        temp_ar = np.array(temp_ar) / temp_scale
        conc_ar = np.log10(mol_num / (box_ar ** 3.))
        conc_ar = conc_ar.flatten()
        temp_ar = temp_ar.T.flatten()
        conc_grid = np.linspace(conc_ar.min(), conc_ar.max(), inter_points)
        temp_grid = np.linspace(temp_ar.min(), temp_ar.max(), inter_points)
        conc_grid, temp_grid = np.meshgrid(conc_grid, temp_grid)

        gyrRad_File = this_sys['GYR File']
        gyrRad = np.load(gyrRad_File)
        gyrRad = gyrRad[:, :, 0].T
        gyrRad = gyrRad.flatten()

        gyrRad_func = sp.interpolate.Rbf(conc_ar, temp_ar, gyrRad, function='linear')
        gyrRad_grid = gyrRad_func(conc_grid, temp_grid)
        return [conc_grid, temp_grid, gyrRad_grid]

    def Gen_OrderParamInterp_For(self, sysName, of_comp=0, temp_scale=2., inter_points=50):
        try:
            this_sys = self.SysInfo[sysName]
            key_dum = this_sys['Key File'][1]
        except KeyError:
            print("Is this the correct system?")

        init_temp = key_dum['MC_TEMP']
        temp_steps = key_dum['MC_CYCLE_NUM']
        delta_temp = key_dum['MC_DELTA_TEMP']
        mol_num = this_sys['Tot Molecules']
        box_ar = [[a_box for a_box in this_sys['Boxes']]] * temp_steps
        temp_ar = [[init_temp + i * delta_temp for i in range(temp_steps)]] * len(box_ar[0])
        box_ar = np.array(box_ar)
        temp_ar = np.array(temp_ar) / temp_scale
        conc_ar = np.log10(mol_num / (box_ar ** 3.))
        conc_ar = conc_ar.flatten()
        temp_ar = temp_ar.T.flatten()
        conc_grid = np.linspace(conc_ar.min(), conc_ar.max(), inter_points)
        temp_grid = np.linspace(temp_ar.min(), temp_ar.max(), inter_points)
        conc_grid, temp_grid = np.meshgrid(conc_grid, temp_grid)

        rho_barFile = this_sys['RhoBar File']
        rho_bar = np.load(rho_barFile)
        this_rho = rho_bar[:, of_comp, :].T
        this_rho = this_rho.flatten()
        phi_cFile = this_sys['Perc File']
        phi_c = np.load(phi_cFile)
        phi_c = phi_c[:, :, 0].T
        phi_c = phi_c.flatten() / mol_num

        rho_func = sp.interpolate.Rbf(conc_ar, temp_ar, this_rho, function='linear')
        rho_grid = rho_func(conc_grid, temp_grid)
        phi_func = sp.interpolate.Rbf(conc_ar, temp_ar, phi_c, function='linear')
        phi_grid = phi_func(conc_grid, temp_grid)
        return [conc_grid, temp_grid, rho_grid, phi_grid]

    def Store_OrderParamInterp_For(self, sysName, of_comp=0, temp_scale=2., inter_points=50):
        try:
            this_sys = self.SysInfo[sysName]
        except KeyError:
            print("Is this the correct system?")
        this_sys['Order Params'] = self.Gen_OrderParamInterp_For(sysName, of_comp, temp_scale, inter_points)


class TwoComp_Setup:
    def __init__(self, param_file):
        self.GlobalParamFile = param_file
        self.CurrentDir = os.getcwd() + '/'
        self.Set_SimulationPath(self.CurrentDir)
        self.SysInfo = {}

    def Set_SimulationPath(self, end_path_name):
        self.SimulationPath = end_path_name
        print("Simulations shall be done in dir: {:}".format(self.SimulationPath))

    def Set_QSUB_Command(self, qsub_command):
        self.QSubCommand = qsub_command
        self.Reset_QSUB_Iter()

    def Set_QSUB_Queues(self, WIntQueue, NoIntQueue):
        self.QSUB_WIntQ = WIntQueue
        self.QSUB_NoIntQ = NoIntQueue

    def Set_PathToLASSI(self, path_to_lassi):
        self.PathToLASSI = path_to_lassi

    def Reset_QSUB_Iter(self):
        self.QSubIter = 0;

    def Read_StrucFileFor(self, sysName):
        struc_file = self.SysInfo[sysName]['Structure File']
        dum_struc_tot = []
        with open(struc_file) as sFile:
            tot_file = sFile.readlines()
            tot_file = iter(tot_file)
            for a_line in tot_file:
                if a_line[0] == '#' or a_line[0] == '\n':
                    continue
                if a_line[:4] == 'NEW{':
                    dum_struc_in = []
                    a_line = tot_file.next()
                    num_mol = int(a_line[:-1])
                    a_line = tot_file.next()
                    while (a_line[:4] != "}END"):
                        a_line = a_line[:-1].split()
                        this_line = [int(aVal) for aVal in a_line]
                        a_line = tot_file.next()
                        dum_struc_in.append(this_line)
                    dum_struc_tot.append([num_mol, np.array(dum_struc_in)])
        return dum_struc_tot

    def AddNewSystem(self, sysName, lin_len, mol1_min, mol2_min, totmol_max):
        if sysName in self.SysInfo.keys():
            print("This structure name already exists! Doing nothing.")
        else:
            self.SysInfo[sysName] = {}
            self.SysInfo[sysName]['Linker Length'] = lin_len
            self.SysInfo[sysName]['Structure'] = []
            self.SysInfo[sysName]['Mol Mins']  = np.array([mol1_min, mol2_min])
            self.SysInfo[sysName]['Mol Max']   = totmol_max
            self.AddParamFileTo(sysName, self.GlobalParamFile)
            self.SetRunNameFor(sysName, sysName)
            self.SetNumberOfRunsFor(sysName, 2)

    def AddStrucFileTo(self, sysName, file_name):
        try:
            self.SysInfo[sysName]['Structure File'] = file_name
            print("{:} has structure file: {:}".format(sysName, file_name))
            self.SysInfo[sysName]['Key File'][1]['STRUCT_FILE'] = self.CurrentDir + file_name
            self.SysInfo[sysName]['Structure'] = self.Read_StrucFileFor(sysName)
            if len(self.SysInfo[sysName]['Structure'] < 2):
                raise KeyError
            self.CalcNonZeroRDFComps(sysName)
            # self.SysInfo[sysName]['Tot Molecules'] = self.CalcTotalMolecules(sysName)
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def AddEnergyFileTo(self, sysName, fileName_NoInt, fileName_WInt):
        try:
            with open(self.CurrentDir + fileName_WInt) as WIntFile:
                pass
            with open(self.CurrentDir + fileName_NoInt) as NoIntFile:
                pass
        except IOError:
            print("Do these files exist?")
            return
        try:
            self.SysInfo[sysName]['Int Energy File'] = self.CurrentDir + fileName_WInt
            self.SysInfo[sysName]['NoInt Energy File'] = self.CurrentDir + fileName_NoInt
            # print("{:} has energy files: Int:{:}  NoInt:{:}".format(sysName, self.SysInfo[sysName]['Int Energy File'],
            #                                                         self.SysInfo[sysName]['NoInt Energy File']))
            with open(self.CurrentDir + fileName_WInt) as WIntFile:
                WIntFile.readline()  # First line is a comment
                tot_stickers = int(WIntFile.readline())
            self.SysInfo[sysName]['Tot Bead Types'] = tot_stickers
        except KeyError:
            print("Failed! Did you type the correct system name?")
            return

    def AddEnergyFiles_ForAll(self, fileName_NoInt, fileName_WInt):
        for aSys in self.SysInfo.keys():
            self.AddEnergyFileTo(aSys, fileName_NoInt, fileName_WInt)

    def AddParamFileTo(self, sysName, fileName):
        self.SysInfo[sysName]['Key File'] = Read_ParamFile(fileName)

    def AddParamFiles_ForAll(self, fileName):
        for aSys in self.SysInfo.keys():
            self.AddParamFileTo(aSys, fileName)

    def PrintParamsFor(self, sysName):
        try:
            dumKeys = self.SysInfo[sysName]['Key File']
            dumDict = dumKeys[1]
            dumKeys = dumKeys[0]
            for aKey in dumKeys:
                dum_spaces = 25 - len(aKey)
                print('{:}{:}{:}'.format(aKey, ' ' * dum_spaces, dumDict[aKey]))
        except KeyError as myErr:
            if myErr[0] == 'Key File':
                print("Did you import the file?")
            else:
                print("Failed! Did you type the correct system name?")

    def ResetStruc(self, sysName):
        self.SysInfo[sysName]['Structure'] = []

    def AddStruc_ImpLinear(self, sysName, mol_num, bead_num, bead_type, lin_len):
        try:
            self.SysInfo[sysName]['Structure'].append(sGen_Implicit_Linear(mol_num, bead_num, bead_type, lin_len))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_Monomer(self, sysName, mol_num, bead_type):
        try:
            self.SysInfo[sysName]['Structure'].append(sGen_Monomer(mol_num, bead_type))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_Dimer(self, sysName, mol_num, bead_type1, bead_type2, lin_len):
        try:
            self.SysInfo[sysName]['Structure'].append(sGen_Dimer(mol_num, bead_type1, bead_type2, lin_len))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_ExpLinear(self, sysName, mol_num, bead_num, bead_type, lin_len, lin_type):
        try:
            self.SysInfo[sysName]['Structure'].append(
                sGen_Explicit_Linear(mol_num, bead_num, bead_type, lin_len, lin_type))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_ImpBranSymm(self, sysName, mol_num, branch_num, bead_per_branch, bead_type, hub_type, lin_len, hub_len):
        try:
            self.SysInfo[sysName]['Structure'].append(
                sGen_Implicit_SymmBranched(mol_num, branch_num, bead_per_branch, bead_type, hub_type, lin_len, hub_len))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def Write_TwoCompStructureFileFor(self, sysName, file_prefix, mol1_num, mol2_num):
        dummyStruc = self.SysInfo[sysName]['Structure'][:]
        dummyStruc[0][0] = mol1_num
        dummyStruc[1][0] = mol2_num
        file_name = file_prefix + str(mol1_num) + '_' + str(mol2_num) + '.prm'
        sGen_WriteStructuresToFile(dummyStruc, file_name)

    def Write_StructureFilesFor(self, sysName):
        this_sys_num_list = self.SysInfo[sysName]['Molecule Numbers']
        MKDirCatch(sysName + '_Structures/')
        file_prefix = sysName + '_Structures/'
        for mol1, mol2 in this_sys_num_list:
            self.Write_TwoCompStructureFileFor(sysName, file_prefix, mol1, mol2)
        self.CalcNonZeroRDFComps(sysName)
        self.CalcCOMDenComps(sysName)
        self.Set_StructureFolderFor(sysName)

    def Set_StructureFolderFor(self, sysName):
        file_prefix = sysName + '_Structures/'
        self.SysInfo[sysName]['Structure Folder'] = self.CurrentDir + file_prefix;

    def Write_StructureFiles_ForAll(self):
        for aSys in self.SysInfo.keys():
            Write_StructureFilesFor(aSys, aSys)

    def SetTemperaturesFor(self, sysName, init_temp, final_temp, temp_steps, therm_temp):
        assert (type(temp_steps) == int) and (type(init_temp) == type(final_temp) == float)
        assert temp_steps > 0. and init_temp > 0. and final_temp > 0.
        if temp_steps > 1:
            delta_temp = (final_temp - init_temp) / (temp_steps - 1)
        else:
            delta_temp = 0.0
        try:
            key_dum = self.SysInfo[sysName]['Key File'][1]
            key_dum['MC_TEMP'] = init_temp
            key_dum['MC_CYCLE_NUM'] = temp_steps
            key_dum['MC_DELTA_TEMP'] = delta_temp
            key_dum['PREEQ_TEMP'] = therm_temp
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetTemperatures_ForAll(self, init_temp, final_temp, temp_steps, therm_temp):
        for aSys in self.SysInfo.keys():
            self.SetTemperaturesFor(aSys, init_temp, final_temp, temp_steps, therm_temp)

    def SetMCStepsFor(self, sysName, therm_steps, run_steps):
        try:
            key_dum = self.SysInfo[sysName]['Key File'][1]
            key_dum['N_STEPS'] = int(run_steps)
            key_dum['PREEQ_STEPS'] = int(therm_steps)
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetMCSteps_ForAll(self, therm_steps, run_steps):
        for aSys in self.SysInfo.keys():
            self.SetMCStepsFor(aSys, therm_steps, run_steps)

    def CalcMolNumArr(self, sysName, tot_bin_num):
        assert tot_bin_num % 2 == 1  # We should always have an odd number input
        actual_bins = (tot_bin_num) / 2 + 1
        num_min = self.SysInfo[sysName]['Mol Mins'].max()
        num_max = self.SysInfo[sysName]['Mol Max']
        n1Mols = np.ceil(10.0 ** np.linspace(np.log10(num_min), np.log10(num_max / 2.), actual_bins))
        n1Mols = np.array(n1Mols, dtype=int)
        mol_ar = np.sort(np.append(n1Mols, num_max - n1Mols[-2::-1]))
        self.SysInfo[sysName]['Molecule Numbers'] = np.array([mol_ar, num_max - mol_ar]).T

    def CalcMolNum_Linear_Arr(self, sysName, tot_bin_num):
        assert tot_bin_num % 2 == 1  # We should always have an odd number input
        num_min = self.SysInfo[sysName]['Mol Mins'].max()
        num_max = self.SysInfo[sysName]['Mol Max']
        n1Mols = np.linspace(num_min, num_max-num_min, tot_bin_num)
        n1Mols = np.array(n1Mols, dtype=int)
        n2Mols = num_max - n1Mols
        self.SysInfo[sysName]['Molecule Numbers'] = np.array([n1Mols, n2Mols]).T

    def CalcHybridMolNum_Arr(self, sysName, tot_bin_num):
        assert tot_bin_num % 2 == 1  # We should always have an odd number input
        #First find linear
        num_min = self.SysInfo[sysName]['Mol Mins'].max()
        num_max = self.SysInfo[sysName]['Mol Max']
        linMols  = np.linspace(num_min, num_max - num_min, tot_bin_num+1)
        linMols  = np.array(linMols, dtype=int)
        #Now log-space
        actual_bins = (tot_bin_num) / 2 + 1
        logMols = np.ceil(10.0 ** np.linspace(np.log10(num_min), np.log10(num_max / 2.), actual_bins))
        logMols = np.array(logMols, dtype=int)
        logMols = np.sort(np.append(logMols, num_max - logMols[-2::-1]))
        #Combine the two
        totmAr = np.unique(np.append(linMols, logMols))
        totmAr = np.array([totmAr, totmAr[::-1]]).T
        self.SysInfo[sysName]['Molecule Numbers'] = totmAr[::2]

    def CalcBoxSizeArr(self, low_con, high_con, tot_beads, tot_boxes):
        if not low_con < high_con < 0.0:
            print("Low conc should be lower than high conc!"
                  "Both numbers need to be negative. (10^low_con, 10^high_con).")
            return
        dum_li = np.linspace(low_con, high_con, tot_boxes)
        dum_ar = 10. ** dum_li
        dum_ar = tot_beads / dum_ar
        dum_ar = np.array(dum_ar ** (1. / 3.), dtype=int)
        dum_ar_s = np.sort(dum_ar)
        return dum_ar_s

    def CalcBoxSize_Linear_Arr(self, low_con, high_con, tot_beads, tot_boxes):
        if not 0.0 < low_con < high_con < 1.0:
            print("Low conc should be lower than high conc!"
                  "Both numbers need to be in the interval (0,1).")
            return
        dum_ar = np.linspace(low_con, high_con, tot_boxes)
        dum_ar = tot_beads / dum_ar
        dum_ar = np.array(dum_ar ** (1. / 3.), dtype=int)
        dum_ar_s = np.sort(dum_ar)
        return dum_ar_s

    def CalcMaxBeadsFor(self, sysName):
        # Figure out the maximum possible number of beads given the structure and molecule numbers.
        try:
            struc_ar = self.SysInfo[sysName]['Structure']
            mol_nums = [];
            bead_nums = [];
            for a_mol in struc_ar:
                num_mol = a_mol[0]
                struc_list = a_mol[1].T[0]
                num_beads_per = len(np.unique(struc_list))
                bead_nums.append(num_beads_per)
                mol_nums.append(num_mol)
            mol_nums = np.array(mol_nums)
            bead_nums = np.array(bead_nums)
            largest_mol = np.argmax(bead_nums)
            smaller_mol = (largest_mol + 1) % 2
            mol_num_arr = self.SysInfo[sysName]['Molecule Numbers']
            largest_num = np.argmax(mol_num_arr.T[largest_mol])
            smaller_num = mol_num_arr[largest_num][smaller_mol]
            largest_num = mol_num_arr[largest_num][largest_mol]
            mol_nums[largest_mol] = largest_num;
            mol_nums[smaller_mol] = smaller_num;

            return np.dot(mol_nums, bead_nums)
        except KeyError:
            print("Failed! Did you type the correct system name?")
            raise

    def CalcNonZeroRDFComps(self, sysName):
        try:
            struc_list = self.SysInfo[sysName]['Structure']
        except KeyError:
            print("Did you type the system name correctly?")
            raise

        bead_type_list = []
        for a_mol in struc_list:
            bead_types = a_mol[1].T[1]
            bead_type_list.append(np.unique(bead_types))
        bead_type_list = [aType for a_mol in bead_type_list for aType in a_mol]
        bead_types = np.array(bead_type_list)
        bead_types = bead_types.flatten()
        bead_types = np.unique(bead_types)
        bead_type_list = [[-1, -1]]  # For the total
        for i in range(len(bead_types)):
            for j in range(i, len(bead_types)):
                bead_type_list.append([bead_types[i], bead_types[j]])
        self.SysInfo[sysName]['Comp List'] = np.array(bead_type_list, dtype=int)

    def CalcCOMDenComps(self, sysName):
        try:
            struc_list = self.SysInfo[sysName]['Structure']
        except KeyError:
            print("Did you type the system name correctly?")
            raise

        num_mols = len(struc_list) #Total number of molecules
        dum_comp_list = [];
        for i in range(num_mols):
            dum_comp_list.append([-1, i])
        for i in range(num_mols):
            for j in range(num_mols):
                dum_comp_list.append([i, j])
        self.SysInfo[sysName]['COMDen List'] = np.array(dum_comp_list, dtype=int)

    def SetBoxSizesFor(self, sysName, low_con, high_con, tot_boxes):
        try:
            tot_beads = self.CalcMaxBeadsFor(sysName)
            self.SysInfo[sysName]['Boxes'] = self.CalcBoxSizeArr(low_con, high_con, tot_beads, tot_boxes)
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetLinearBoxSizesFor(self, sysName, low_con, high_con, tot_boxes):
        try:
            tot_beads = self.CalcMaxBeadsFor(sysName)
            self.SysInfo[sysName]['Boxes'] = self.CalcBoxSize_Linear_Arr(low_con, high_con, tot_beads, tot_boxes)
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetHybridBoxSizesFor(self, sysName, low_con, high_con, tot_boxes):
        try:
            tot_beads = self.CalcMaxBeadsFor(sysName)
            boxAr1    = self.CalcBoxSize_Linear_Arr(10.**low_con, 10.**high_con, tot_beads, tot_boxes)
            boxAr2    = self.CalcBoxSizeArr(low_con, high_con, tot_beads, tot_boxes)
            tot_boxAr = np.unique(np.append(boxAr1, boxAr2).flatten())
            self.SysInfo[sysName]['Boxes'] = tot_boxAr[::2]
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetBoxSizes_ForAll(self, low_con, high_con, tot_boxes):
        for aSys in self.SysInfo.keys():
            self.SetBoxSizesFor(aSys, low_con, high_con, tot_boxes)

    def SetLinearBoxSizes_ForAll(self, low_con, high_con, tot_boxes):
        for aSys in self.SysInfo.keys():
            self.SetLinearBoxSizesFor(aSys, low_con, high_con, tot_boxes)

    def SetBoxSizes_To(self, sysName):
        try:
            tot_ar = self.SysInfo[sysName]['Boxes'][:]
        except KeyError:
            print("Failed! Did you type the correct system name? "
                  "Or have you setup the boxes for this system?")
            return
        for aSys in self.SysInfo.keys():
            if aSys != sysName:
                self.SysInfo[aSys]['Boxes'] = tot_ar[:]

    def ShowSamplingGridFor(self, sysName, view_lims):
        plt.figure(figsize=[6,6])
        mol_num = self.SysInfo[sysName]['Molecule Numbers']
        box_arr = self.SysInfo[sysName]['Boxes']
        run_num = self.SysInfo[sysName]['Runs']
        N1_arr  = mol_num.T[0]
        N2_arr  = mol_num.T[1]
        print(r'A total of {:} independent conditions! With {:} runs per condition, there are a total of {:} independent simulations (WInt & NoInt).'.format(
            len(N1_arr) * len(box_arr), run_num, len(N1_arr) * len(box_arr) * run_num * 2))

        for aBox in box_arr:
            plt.plot(np.log10(N1_arr/(aBox**3.)), np.log10(N2_arr/(aBox**3.)), 'ro', alpha=0.9)
        plt.plot(np.linspace(-10,0,10),np.linspace(-10,0,10), 'k--',lw=3,alpha=0.5)
        plt.xlim(view_lims, 0)
        plt.ylim(view_lims, 0)
        plt.grid(lw=2,alpha=0.7)
        plt.show()

    def ShowLinearSamplingGridFor(self, sysName, view_lims):
        plt.figure(figsize=[6,6])
        mol_num = self.SysInfo[sysName]['Molecule Numbers']
        box_arr = self.SysInfo[sysName]['Boxes']
        run_num = self.SysInfo[sysName]['Runs']
        N1_arr  = mol_num.T[0]
        N2_arr  = mol_num.T[1]
        print(r'A total of {:} independent conditions! With {:} runs per condition, there are a total of {:} independent simulations (WInt & NoInt).'.format(
            len(N1_arr) * len(box_arr), run_num, len(N1_arr) * len(box_arr) * run_num * 2))

        for aBox in box_arr:
            plt.plot((N1_arr/(aBox**3.)), (N2_arr/(aBox**3.)), 'ro', alpha=0.9)
        plt.plot(np.linspace(0,1,10),np.linspace(0,1,10), 'k--',lw=3,alpha=0.5)
        plt.xlim(0.0, view_lims)
        plt.ylim(0.0, view_lims)
        plt.grid(lw=2,alpha=0.7)
        plt.show()

    def ShowHybridSamplingGridFor(self, sysName, lower_lims, upper_lims):
        myFig = plt.figure(figsize=[12, 6])
        myAx  = myFig.subplots(1,2)
        mol_num = self.SysInfo[sysName]['Molecule Numbers']
        box_arr = self.SysInfo[sysName]['Boxes']
        run_num = self.SysInfo[sysName]['Runs']
        N1_arr = mol_num.T[0]
        N2_arr = mol_num.T[1]
        print(
            r'A total of {:} independent conditions! With {:} runs per condition, there are a total of {:} independent simulations (WInt & NoInt).'.format(
                len(N1_arr) * len(box_arr), run_num, len(N1_arr) * len(box_arr) * run_num * 2))

        for aBox in box_arr:
            myAx[0].loglog((N1_arr / (aBox ** 3.)), (N2_arr / (aBox ** 3.)), 'ro', alpha=0.9)
            myAx[1].plot((N1_arr / (aBox ** 3.)), (N2_arr / (aBox ** 3.)), 'ro', alpha=0.9)
        [anAx.plot(np.linspace(lower_lims*0.9, upper_lims*1.1, 10), np.linspace(lower_lims*0.9, upper_lims*1.1, 10), 'k--', lw=3, alpha=0.5) for anAx in myAx]
        [anAx.set_xlim(lower_lims, upper_lims)]
        [anAx.set_ylim(lower_lims, upper_lims)]
        [anAx.grid(lw=2, alpha=0.7) for anAx in myAx]
        plt.tight_layout()
        plt.show()

    def SetRunNameFor(self, sysName, runName):
        try:
            self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX'] = runName
        except KeyError:
            print("Failed! Did you type the correct system name?"
                  " Does this system exist?")

    def SetRunName_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.SetRunNameFor(aSys, aSys)

    def SetNumberOfRunsFor(self, sysName, run_num):
        try:
            self.SysInfo[sysName]['Runs'] = run_num
        except KeyError:
            print("Failed! Did you type the correct system name?"
                  " Does this system exist?")

    def SetNumberOfRuns_ForAll(self, run_num):
        for aSys in self.SysInfo.keys():
            self.SetNumberOfRunsFor(aSys, run_num)

    def MakeDirs_NoIntFor(self, sysName):

        dum_dir = self.SimulationPath
        dum_dir += 'PhaseDiags/'
        MKDirCatch(dum_dir)
        dum_dir += sysName + '/'
        MKDirCatch(dum_dir)

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']

        for b_ID, a_box in enumerate(box_arr):
            this_dir = dum_dir + str(a_box) + '/'
            MKDirCatch(this_dir)
            for a_con in num_arr:
                con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/'
                MKDirCatch(con_dir)
                con_dir += 'NoInt/'
                MKDirCatch(con_dir)
                for a_run in range(num_of_runs):
                    final_dir = con_dir + str(a_run + 1) + '/'
                    MKDirCatch(final_dir)

    def MakeDirs_WIntFor(self, sysName):

        dum_dir = self.SimulationPath
        dum_dir += 'PhaseDiags/'
        MKDirCatch(dum_dir)
        dum_dir += sysName + '/'
        MKDirCatch(dum_dir)

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']

        for b_ID, a_box in enumerate(box_arr):
            this_dir = dum_dir + str(a_box) + '/'
            MKDirCatch(this_dir)
            for a_con in num_arr:
                con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/'
                MKDirCatch(con_dir)
                con_dir += 'WInt/'
                MKDirCatch(con_dir)
                for a_run in range(num_of_runs):
                    final_dir = con_dir + str(a_run + 1) + '/'
                    MKDirCatch(final_dir)

    def MakeDirs_For(self, sysName):
        self.MakeDirs_WIntFor(sysName)
        self.MakeDirs_NoIntFor(sysName)

    def MakeDirs_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.MakeDirs_For(aSys)

    def MakeDirs_ForAll_NoInt(self):
        for aSys in self.SysInfo.keys():
            self.MakeDirs_NoIntFor(aSys)

    def MakeDirs_ForAll_WInt(self):
        for aSys in self.SysInfo.keys():
            self.MakeDirs_WIntFor(aSys)

    def Write_ParamsWIntFor(self, sysName):
        dum_dir = self.SimulationPath+'PhaseDiags/'+sysName+'/'

        param_copy = self.SysInfo[sysName]['Key File'][:]

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']
        param_copy[1]['ENERGY_FILE'] = self.SysInfo[sysName]['Int Energy File']
        struc_folder = self.SysInfo[sysName]['Structure Folder']

        for b_ID, a_box in enumerate(box_arr):
            param_copy[1]['BOX_SIZE'] = a_box
            this_dir = dum_dir + str(a_box) + '/'
            for a_con in num_arr:
                con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
                param_copy[1]['STRUCT_FILE'] = struc_folder + str(a_con[0]) + '_' + str(a_con[1]) + '.prm'
                for a_run in range(num_of_runs):
                    param_copy[1]['RANDOM_SEED'] = 0
                    final_dir = con_dir + str(a_run + 1) + '/'
                    with open(final_dir + 'param.key', 'w+') as pFile:
                        for a_key in param_copy[0]:
                            N_spcs = 25 - len(a_key)
                            pFile.write(a_key + ' ' * N_spcs + str(param_copy[1][a_key]) + '\n')

    def Write_ParamsNoIntFor(self, sysName):
        dum_dir = self.SimulationPath+'PhaseDiags/'+sysName+'/'

        param_copy = self.SysInfo[sysName]['Key File'][:]

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']
        param_copy[1]['ENERGY_FILE'] = self.SysInfo[sysName]['NoInt Energy File']
        param_copy[1]['MV_STROT_FREQ'] = '0.0';
        param_copy[1]['MV_COLOCAL_FREQ'] = '0.0';
        param_copy[1]['MV_SMCLSTR_FREQ'] = '0.0';
        param_copy[1]['MV_CLSTR_FREQ'] = '0.0';
        param_copy[1]['REPORT_NETWORK_FREQ'] = '0';
        struc_folder = self.SysInfo[sysName]['Structure Folder']

        for b_ID, a_box in enumerate(box_arr):
            param_copy[1]['BOX_SIZE'] = a_box
            this_dir = dum_dir + str(a_box) + '/'
            for a_con in num_arr:
                con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/NoInt/'
                param_copy[1]['STRUCT_FILE'] = struc_folder + str(a_con[0]) + '_' + str(a_con[1]) + '.prm'
                for a_run in range(num_of_runs):
                    param_copy[1]['RANDOM_SEED'] = 0
                    final_dir = con_dir + str(a_run + 1) + '/'
                    with open(final_dir + 'param.key', 'w+') as pFile:
                        for a_key in param_copy[0]:
                            N_spcs = 25 - len(a_key)
                            pFile.write(a_key + ' ' * N_spcs + str(param_copy[1][a_key]) + '\n')

    def Write_ParamsFor(self, sysName):
        self.Write_ParamsWIntFor(sysName)
        self.Write_ParamsNoIntFor(sysName)

    def Write_ParamsWInt_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.Write_ParamsWIntFor(aSys)

    def Write_ParamsNoInt_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.Write_ParamsNoIntFor(aSys)

    def Write_Params_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.Write_ParamsWIntFor(aSys)
            self.Write_ParamsNoIntFor(aSys)

    def SubmitWIntJobs_ToQueueFor(self, sysName, WIntQueue):
        dum_dir = self.SimulationPath+'PhaseDiags/'+sysName+'/'

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']

        for b_ID, a_box in enumerate(box_arr[::]):
            this_dir = dum_dir + str(a_box) + '/'
            for a_con in num_arr:
                con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
                for a_run in range(num_of_runs):
                    final_dir = con_dir + str(a_run + 1) + '/'
                    run_dir = final_dir
                    run_num = ['W' + sysName[:2], b_ID, self.QSubIter]
                    run_num = [str(a_val) for a_val in run_num]
                    run_num = "_".join(run_num)
                    run_command = self.QSubCommand + ' ' + WIntQueue + ' ' + run_num + ' /.' + self.PathToLASSI
                    try:
                        os.chdir(run_dir)
                    except IOError:
                        print("Did you make all the directories?")
                        raise
                    try:
                        ret_code = sproc.check_output(run_command, shell=True, stderr=sproc.STDOUT)
                        if ret_code < 0:
                            print("The submission failed! Signal: ", -retcode)
                        else:
                            print(ret_code)
                    except OSError as myErr:
                        print("Couldn't submit because ", myErr)
                        raise
                    time.sleep(2)
                    self.QSubIter += 1
        os.chdir(self.CurrentDir)

    def SubmitNoIntJobs_ToQueueFor(self, sysName, NoIntQueue):
        dum_dir = self.SimulationPath+'PhaseDiags/'+sysName+'/'

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']

        for b_ID, a_box in enumerate(box_arr[::]):
            this_dir = dum_dir + str(a_box) + '/'
            for a_con in num_arr:
                con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/NoInt/'
                for a_run in range(num_of_runs):
                    final_dir = con_dir + str(a_run + 1) + '/'
                    run_dir = final_dir
                    run_num = ['N' + sysName[:2], b_ID, self.QSubIter]
                    run_num = [str(a_val) for a_val in run_num]
                    run_num = "_".join(run_num)
                    run_command = self.QSubCommand + ' ' + NoIntQueue + ' ' + run_num + ' /.' + self.PathToLASSI
                    try:
                        os.chdir(run_dir)
                    except IOError:
                        print("Did you make all the directories?")
                        raise
                    try:
                        ret_code = sproc.check_output(run_command, shell=True, stderr=sproc.STDOUT)
                        if ret_code < 0:
                            print("The submission failed! Signal: ", -retcode)
                        else:
                            print(ret_code)
                    except OSError as myErr:
                        print("Couldn't submit because ", myErr)
                        raise
                    time.sleep(2)
                    self.QSubIter += 1

        os.chdir(self.CurrentDir)

    def SubmitWIntJobs_ForAll(self, WIntQueue):
        for aSys in self.SysInfo.keys():
            self.SubmitWIntJobs_ToQueueFor(aSys, WIntQueue)

    def SubmitNoIntJobs_ForAll(self, NoIntQueue):
        for aSys in self.SysInfo.keys():
            self.SubmitWIntJobs_ToQueueFor(aSys, NoIntQueue)

    def SubmitJobs_ForAll(self):
        try:
            self.QSubCommand
            self.QSUB_WIntQ
            self.QSUB_NoIntQ
        except AttributeError:
            print("You need to set up the qsub command, and which queues to use!")
            return
        self.Reset_QSUB_Iter()
        for aSys in self.SysInfo.keys():
            self.SubmitNoIntJobs_ToQueueFor(aSys, self.QSUB_NoIntQ)
            self.SubmitWIntJobs_ToQueueFor(aSys, self.QSUB_WIntQ)

    def SetAndRun_For(self, sysName):
        self.MakeDirs_For(sysName)
        self.Write_ParamsFor(sysName)
        self.Reset_QSUB_Iter()
        self.SubmitWIntJobs_ToQueueFor(sysName, self.QSUB_WIntQ)
        self.SubmitWIntJobs_ToQueueFor(sysName, self.QSUB_NoIntQ)

class TwoComp_Analysis(TwoComp_Setup):
    def __init__(self, SimSetupInstance):
        try:
            self.SysInfo = SimSetupInstance.SysInfo
            self.GlobalParamFile = SimSetupInstance.GlobalParamFile
            self.SimulationPath = SimSetupInstance.SimulationPath
            self.CurrentDir = SimSetupInstance.CurrentDir
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise

    def Read_PDFFileFor(self, sysName, file_name, box_size, tot_temps):
        try:
            dum_dat = np.loadtxt(file_name)
            comp_list = self.SysInfo[sysName]['Comp List']
            tot_comps = self.SysInfo[sysName]['Tot Bead Types']
            tot_comps_pos = 1 + Index_RDF(tot_comps - 2, tot_comps - 1, tot_comps)
        except IOError:
            print("Does the PDF file exist? Did these runs finish and output data?")
            print("Failed for "+file_name)
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        actual_comps = len(comp_list)
        ret_dat = np.zeros((actual_comps, tot_temps, box_size * 4))
        for a_temp in range(tot_temps):
            for compID, compPair in enumerate(comp_list):
                compA = int(compPair[0]);
                compB = int(compPair[1])
                ret_dat[compID][a_temp] = dum_dat[tot_comps_pos * a_temp + Index_RDF(compA, compB, tot_comps)]
        return ret_dat

    def Read_CLUSFileFor(self, file_name):
        try:
            dum_dat = np.loadtxt(file_name)
        except IOError:
            print("Does the CLUS file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Read_COMDenFileFor(self, sysName, file_name, box_size, tot_temps):
        try:
            dum_dat   = np.loadtxt(file_name)
            comp_list = self.SysInfo[sysName]['COMDen List']
            num_mols = len(np.unique(comp_list.T[0]))-1
            tot_comps = num_mols*(num_mols+1)
        except IOError:
            print("Does the PDF file exist? Did these runs finish and output data?")
            print("Failed for " + file_name)
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        ret_dat = np.zeros((tot_comps, tot_temps, box_size * 4))
        for a_temp in range(tot_temps):
            for compID, compPair in enumerate(comp_list):
                compA = int(compPair[0]);
                compB = int(compPair[1])
                ret_dat[compID][a_temp] = dum_dat[tot_comps * a_temp + Index_COMDen(compA, compB, num_mols)]
        return ret_dat

    def Read_MolClusFileFor(self, file_name):
        try:
            dum_dat = np.loadtxt(file_name)
        except IOError:
            print("Does the MolClus file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Read_GYRADFileFor(self, file_name):
        try:
            dum_dat = np.loadtxt(file_name)
        except IOError:
            print("Does the GYRRAD file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Save_NoIntPDF_For(self, sysName):
        try:
            dum_dat = self.Collect_NoInt_PDFs_For(sysName)
            file_name = sysName + '_N_PDF.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_NoInt_PDF'] = file_name
            pass
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntPDF_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_PDFs_For(sysName)
            file_name = sysName + '_W_PDF.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_WInt_PDF'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntCOMDen_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_COMDens_For(sysName)
            file_name = sysName + '_W_COMDen.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_WInt_COMDen'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_NoIntCOMDen_For(self, sysName):
        try:
            dum_dat = self.Collect_NoInt_COMDens_For(sysName)
            file_name = sysName + '_N_COMDen.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_NoInt_COMDen'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntMolClus_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_MolClus_For(sysName)
            file_name = sysName + '_W_MolClus.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_WInt_MolClus'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntClus_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_CLUS_For(sysName)
            file_name = sysName + '_CLUS.b'
            np.save(file_name, dum_dat)
            self.SysInfo[sysName]['Raw_CLUS'] = file_name + '.npy'
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Set_Auto_PDFFileNames_For(self, sysName):
        self.SysInfo[sysName]['Raw_NoInt_PDF'] = sysName + '_N_PDF.b'
        self.SysInfo[sysName]['Raw_WInt_PDF'] = sysName + '_W_PDF.b'

    def Collect_NoInt_PDFs_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'PhaseDiags/'+sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        tot_raw_dat = []
        for b_ID, a_box in enumerate(box_arr):
            box_dir = over_dir + str(a_box)+'/'
            tot_raw_per_box = []
            for a_con in num_arr:
                con_dir = box_dir + str(a_con[0])+'_'+str(a_con[1])+'/NoInt/'
                tot_raw_dat_temp = []
                for a_run in range(run_num):
                    file_name = con_dir + str(a_run + 1) + '/' + run_name + '_RDF.dat'
                    tot_raw_dat_temp.append(self.Read_PDFFileFor(sysName, file_name, a_box, tot_temps))
                tot_raw_per_box.append(tot_raw_dat_temp)
            tot_raw_dat.append(tot_raw_per_box)
        return tot_raw_dat

    def Collect_WInt_PDFs_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'PhaseDiags/' + sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        tot_raw_dat = []
        for b_ID, a_box in enumerate(box_arr):
            box_dir = over_dir + str(a_box) + '/'
            tot_raw_per_box = []
            for a_con in num_arr:
                con_dir = box_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
                tot_raw_dat_temp = []
                for a_run in range(run_num):
                    file_name = con_dir + str(a_run + 1) + '/' + run_name + '_RDF.dat'
                    tot_raw_dat_temp.append(self.Read_PDFFileFor(sysName, file_name, a_box, tot_temps))
                tot_raw_per_box.append(tot_raw_dat_temp)
            tot_raw_dat.append(tot_raw_per_box)
        return tot_raw_dat

    def Collect_WInt_COMDens_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'PhaseDiags/' + sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        tot_raw_dat = []
        for b_ID, a_box in enumerate(box_arr):
            box_dir = over_dir + str(a_box) + '/'
            tot_raw_per_box = []
            for a_con in num_arr:
                con_dir = box_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
                tot_raw_dat_temp = []
                for a_run in range(run_num):
                    file_name = con_dir + str(a_run + 1) + '/' + run_name + '_COMDen.dat'
                    tot_raw_dat_temp.append(self.Read_COMDenFileFor(sysName, file_name, a_box, tot_temps))
                tot_raw_per_box.append(tot_raw_dat_temp)
            tot_raw_dat.append(tot_raw_per_box)
        return tot_raw_dat

    def Collect_NoInt_COMDens_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'PhaseDiags/' + sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        tot_raw_dat = []
        for b_ID, a_box in enumerate(box_arr):
            box_dir = over_dir + str(a_box) + '/'
            tot_raw_per_box = []
            for a_con in num_arr:
                con_dir = box_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/NoInt/'
                tot_raw_dat_temp = []
                for a_run in range(run_num):
                    file_name = con_dir + str(a_run + 1) + '/' + run_name + '_COMDen.dat'
                    tot_raw_dat_temp.append(self.Read_COMDenFileFor(sysName, file_name, a_box, tot_temps))
                tot_raw_per_box.append(tot_raw_dat_temp)
            tot_raw_dat.append(tot_raw_per_box)
        return tot_raw_dat

    def Collect_WInt_MolClus_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'PhaseDiags/' + sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        tot_raw_dat = []
        for b_ID, a_box in enumerate(box_arr):
            box_dir = over_dir + str(a_box) + '/'
            tot_raw_per_box = []
            for a_con in num_arr:
                con_dir = box_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
                tot_raw_dat_temp = []
                for a_run in range(run_num):
                    file_name = con_dir + str(a_run + 1) + '/' + run_name + '_MolClus.dat'
                    tot_raw_dat_temp.append(self.Read_MolClusFileFor(file_name))
                tot_raw_per_box.append(tot_raw_dat_temp)
            tot_raw_dat.append(tot_raw_per_box)
        return tot_raw_dat

    def Collect_WInt_CLUS_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'PhaseDiags/' + sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        tot_raw_dat = []
        for b_ID, a_box in enumerate(box_arr):
            box_dir = over_dir + str(a_box) + '/'
            tot_raw_per_box = []
            for a_con in num_arr:
                con_dir = box_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
                tot_raw_dat_temp = []
                for a_run in range(run_num):
                    file_name = con_dir + str(a_run + 1) + '/' + run_name + '_CLUS.dat'
                    tot_raw_dat_temp.append(self.Read_CLUSFileFor(file_name))
                tot_raw_per_box.append(tot_raw_dat_temp)
            tot_raw_dat.append(tot_raw_per_box)
        return tot_raw_dat

    def Gen_pRDFs_For(self, sysName, my_WInt=None, my_NoInt=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            self.Set_Auto_PDFFileNames_For(sysName)
        except KeyError:
            print("Is this the correct system?")
            return
        if my_WInt != None:
            this_sys['Raw_WInt_PDF'] = my_WInt
        if my_NoInt != None:
            this_sys['Raw_NoInt_PDF'] = my_NoInt
        nFile = this_sys['Raw_NoInt_PDF']
        iFile = this_sys['Raw_WInt_PDF']
        with open(iFile) as dumFile:
            i_pdfs = pickle.load(dumFile)
        with open(nFile) as dumFile:
            n_pdfs = pickle.load(dumFile)

        tot_rdf_arr = []
        for boxID, a_box in enumerate(box_ar):
            x_ar = np.arange(0, a_box, 0.25)
            n_ar = np.array(n_pdfs[boxID]);
            i_ar = np.array(i_pdfs[boxID]);
            n_run_ave, n_run_err = MeanAndError(n_ar, 3);
            n_temp_ave, n_temp_err = MeanAndError(n_run_ave, 1);
            n_temp_err = np.sqrt(n_temp_err**2. + np.sum(n_run_err**2.,axis=1))
            i_run_ave, i_run_err = MeanAndError(i_ar, 1);
            n_sup_ave_for_zeros = np.mean(np.mean(n_temp_ave, axis = 1), axis=0);
            zero_points = np.argwhere(n_sup_ave_for_zeros == 0.).T[0]
            x_ar = np.delete(x_ar, zero_points)
            rdf_per_num = []
            for numID, a_num in enumerate(zip(n_temp_ave, i_run_ave)):
                rdf_per_comp = []
                for compID, aComp in enumerate(zip(a_num[0], a_num[1])):
                    n_dum = np.delete(aComp[0], zero_points)
                    rdf_per_temp = []
                    for tempID, aTemp in enumerate(aComp[1]):
                        i_dum = np.delete(aTemp, zero_points)
                        p_rdf_val = i_dum / n_dum
                        rdf_per_temp.append(p_rdf_val)
                    rdf_per_comp.append(rdf_per_temp)
                rdf_per_num.append(rdf_per_comp)
            tot_rdf_arr.append([x_ar, rdf_per_num])
        rFile = sysName + '_CalcRDF.b'
        this_sys['RDF File'] = rFile
        with open(rFile, "wb+") as rdf_file:
            pickle.dump(tot_rdf_arr, rdf_file)

    def Gen_RhoBar_For(self, sysName, my_rdfFile=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
            comp_list = this_sys['Comp List']
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
            this_sys['RDF File'] = sysName + '_CalcRDF.b'
        except KeyError:
            print("Is this the correct system?")
            return
        if my_rdfFile != None:
            this_sys['RDF File'] = my_rdfFile
        rdf_file = this_sys['RDF File']
        with open(rdf_file) as rFile:
            tot_rdfs = pickle.load(rFile)

        rho_bar_mat = np.zeros((len(box_ar), len(num_arr),  tot_temps, len(comp_list)))
        for boxID, a_box in enumerate(tot_rdfs):
            x_ar = a_box[0]
            pRDFs = a_box[1]
            box_size = box_ar[boxID]
            volume_norm = RDFVolumeElement(x_ar, box_size)
            for numID, aNum in enumerate(pRDFs):
                for compID, aComp in enumerate(aNum):
                    for tempID, aTemp in enumerate(aComp):
                        g_of_r = aTemp;
                        my_func = np.abs(g_of_r - 1.)
                        my_func_norm = my_func * volume_norm
                        rho_bar = np.trapz(my_func_norm, x_ar) / (box_size ** 3.)
                        # rho_bar = np.trapz(g_of_r*my_func_norm,x_ar)/np.trapz(g_of_r*volume_norm,x_ar)
                        rho_bar_mat[boxID][numID][tempID][compID] = rho_bar;
        rhoFile = sysName + '_rhobar.c'
        np.save(rhoFile, rho_bar_mat)
        rhoFile += '.npy'
        this_sys['RhoBar File'] = rhoFile

    def Gen_PhiC_For(self, sysName, my_clusFile=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            num_ar = self.SysInfo[sysName]['Molecule Numbers'];
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
            this_sys['Raw_CLUS'] = sysName + '_CLUS.b.npy'
        except KeyError:
            print("Is this the correct system?")
            return
        if my_clusFile != None:
            this_sys['Raw_CLUS'] = my_clusFile
        clus_file = this_sys['Raw_CLUS']
        tot_clus = np.load(clus_file)

        clus_run_avg = np.mean(tot_clus, axis=2)
        clus_run_err = np.std(tot_clus, axis=2)
        num_mols = float(num_ar[0,0]+num_ar[0,1])
        CMols    = np.linspace(1,num_mols,num_mols)
        phi_c_mat = np.zeros((len(box_ar), len(num_ar),  tot_temps, 5, 2))
        for boxID, aBox in enumerate(box_ar):
            for numID, aNum in enumerate(num_ar):
                for tempID in range(tot_temps):
                    clus_raw_m = clus_run_avg[boxID,numID,tempID]
                    clus_raw_e = clus_run_err[boxID,numID,tempID]
                    this_phi = phi_c_mat[boxID,numID,tempID]

                    this_phi[0,0] = clus_raw_m[0]/num_mols
                    this_phi[0,1] = clus_raw_e[0]/num_mols

                    CDist_Sum = np.sum(clus_raw_m[1:])
                    CDist_Norm = clus_raw_m[1:]/CDist_Sum
                    CErr_Norm = clus_raw_e[1:]/CDist_Sum
                    sigma_mean = np.dot(CMols,CDist_Norm)
                    sigma_m_err = np.sqrt(np.sum((CDist_Norm * CErr_Norm)**2.))
                    sigma_var = np.dot(CMols*CMols, CDist_Norm) - sigma_mean**2.
                    sigma_v_err1 = np.sqrt(np.sum((CMols*CMols*CErr_Norm)**2.))
                    sigma_v_err2 = 2.*sigma_m_err*sigma_mean
                    sigma_var_err = np.sqrt(sigma_v_err1**2. + sigma_v_err2**2.)
                    this_phi[1,0] = sigma_mean
                    this_phi[1,1] = sigma_m_err
                    this_phi[2,0] = sigma_var
                    this_phi[2,1] = sigma_var_err

                    CDist_Sum = np.sum(clus_raw_m[1:] * CMols)
                    CDist_Norm = clus_raw_m[1:] * CMols / CDist_Sum
                    CErr_Norm  = clus_raw_e[1:] * CMols / CDist_Sum
                    sigma_mean = np.dot(CMols, CDist_Norm)
                    sigma_m_err = np.sqrt(np.sum((CDist_Norm * CErr_Norm) ** 2.))
                    sigma_var = np.dot(CMols * CMols, CDist_Norm) - sigma_mean ** 2.
                    sigma_v_err1 = np.sqrt(np.sum((CMols * CMols * CErr_Norm) ** 2.))
                    sigma_v_err2 = 2. * sigma_m_err * sigma_mean
                    sigma_var_err = np.sqrt(sigma_v_err1 ** 2. + sigma_v_err2 ** 2.)
                    this_phi[3, 0] = sigma_mean
                    this_phi[3, 1] = sigma_m_err
                    this_phi[4, 0] = sigma_var
                    this_phi[4, 1] = sigma_var_err

        phiFile = sysName + '_phic.c'
        np.save(phiFile, phi_c_mat)
        phiFile += '.npy'
        this_sys['Perc File'] = phiFile

    def Gen_CorrectDen_For(self, sysName, my_WInt=None, my_NoInt=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            tot_mol_types = len(this_sys['Structure'])
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
            comden_comps = tot_mol_types*(tot_mol_types+1)
        except KeyError:
            print("Is this the correct system?")
            return
        if my_WInt != None:
            this_sys['Raw_WInt_COMDen'] = my_WInt
        else:
            this_sys['Raw_WInt_COMDen'] = sysName + '_W_COMDen.b'
        if my_NoInt != None:
            this_sys['Raw_NoInt_COMDen'] = my_NoInt
        else:
            this_sys['Raw_NoInt_COMDen'] = sysName + '_N_COMDen.b'

        nFile = this_sys['Raw_NoInt_COMDen']
        iFile = this_sys['Raw_WInt_COMDen']

        with open(iFile) as dumFile:
            i_pdfs = pickle.load(dumFile)
        with open(nFile) as dumFile:
            n_pdfs = pickle.load(dumFile)

        per_box_den = []
        for boxID, a_box in enumerate(box_ar):
            x_ar = np.arange(0., a_box, 0.25)
            n_ar = np.array(n_pdfs[boxID])
            i_ar = np.array(i_pdfs[boxID])
            i_run_ave, i_run_err = MeanAndError(i_ar, 1);
            n_run_ave, n_run_err = MeanAndError(n_ar, 1);

            per_num_den  = [];
            for numID, a_num in enumerate(zip(n_run_ave, i_run_ave)):
                num_NAr = a_num[0];
                num_WAr = a_num[1];
                per_comp_den = []
                for compID in range(comden_comps):
                    comp_NAr = num_NAr[compID::comden_comps]
                    comp_IAr = num_WAr[compID::comden_comps]
                    n_temp_ave, n_temp_err = MeanAndError(comp_NAr, 0)
                    zero_points = np.argwhere(n_temp_ave == 0).T[0]
                    comp_xAr = np.delete(x_ar, zero_points)
                    dum_nAr = np.delete(n_temp_ave, zero_points)
                    per_temp_den = []
                    for tempID, aTemp in enumerate(comp_IAr):
                        dum_iAr = np.delete(aTemp, zero_points)
                        per_temp_den.append(dum_iAr)
                    per_comp_den.append([comp_xAr, dum_nAr, per_temp_den])
                per_num_den.append(per_comp_den)
            per_box_den.append(per_num_den)

        rFile = sysName + '_CorrectDen.b'
        this_sys['CorrectDen File'] = rFile
        with open(rFile, "wb+") as den_file:
            pickle.dump(per_box_den, den_file)

    def Gen_TieLineData_For(self, sysName, my_corrDenFile=None, idx_upto=13, idx_upfrom=-13, idx_uptill=-3):
        this_sys = self.SysInfo[sysName]
        if my_corrDenFile != None:
            this_sys['CorrectDen File'] = my_corrDenFile
        else:
            this_sys['CorrectDen File'] = sysName + '_CorrectDen.b'
        dFile = this_sys['CorrectDen File'];
        with open(dFile) as dumFile:
            tot_den_data = pickle.load(dumFile)


        box_ar        = this_sys['Boxes']
        tot_struc     = this_sys['Structure']
        tot_mol_types = len(tot_struc)
        tot_poss_den  = tot_mol_types*(tot_mol_types+1)
        tot_temps     = this_sys['Key File'][1]['MC_CYCLE_NUM']
        temp_num_ar   = this_sys['Molecule Numbers'].T #This only has molecule numbers from the two sampled molecules
        #but ignores the rest, which we should manually add.
        tot_num_ar = []
        tot_num_ar.append(temp_num_ar[0])
        tot_num_ar.append(temp_num_ar[1])
        for additional_mol in tot_struc[2:]:
            this_num = additional_mol[0]
            this_num = this_num*np.ones(len(temp_num_ar[0]))
            tot_num_ar.append(this_num)
        tot_num_ar = np.array(tot_num_ar).T

        #The tie-line data are organized as follows:
        # Temperature, BoxSize, Molecule Number, Component, and then
        # [highCon, highCon_std, lowCon, lowCon_std, errSum/deltaCon]

        tie_line_data = np.zeros((tot_temps, len(box_ar), len(tot_num_ar), tot_poss_den, 6))
        for bID, aBox in enumerate(box_ar):
            inverse_vol = 1./(aBox**3.)
            for numID, aNum in enumerate(tot_num_ar):
                for compID in range(tot_poss_den):
                    this_comp_den = aNum[compID % tot_mol_types]*inverse_vol
                    #this_comp_xAr = tot_den_data[bID][numID][compID][0][1:]
                    this_comp_nAr = tot_den_data[bID][numID][compID][1][1:]
                    for tempID in range(tot_temps):
                        this_comp_iAr = tot_den_data[bID][numID][compID][2][tempID][1:]
                        this_comp_yAr = this_comp_den*this_comp_iAr/this_comp_nAr
                        highCon, highConEr = MeanAndError(this_comp_yAr[:idx_upto], 0)
                        lowCon, lowConEr   = MeanAndError(this_comp_yAr[idx_upfrom:idx_uptill], 0)
                        deltaCon           = np.abs(highCon-lowCon)
                        totErr             = np.sqrt(highConEr**2. + lowConEr**2.)

                        tie_line_data[tempID][bID][numID][compID][0] = this_comp_den
                        tie_line_data[tempID][bID][numID][compID][1] = highCon
                        tie_line_data[tempID][bID][numID][compID][2] = highConEr
                        tie_line_data[tempID][bID][numID][compID][3] = lowCon
                        tie_line_data[tempID][bID][numID][compID][4] = lowConEr
                        tie_line_data[tempID][bID][numID][compID][5] = totErr/deltaCon
        tieFile = sysName + '_TieLines.c'
        np.save(tieFile, tie_line_data)
        tieFile += '.npy'
        this_sys['TieLines File'] = tieFile

    def Gen_InterpolationGrids_For(self, sysName, num_points=75):
        this_sys   = self.SysInfo[sysName]
        mol_num_ar = this_sys['Molecule Numbers'][:]
        box_ar     = this_sys['Boxes'][:]
        cGrid1s    = [np.log10((mol_num_ar.T[0] + mol_num_ar.T[1]) / ((aBox) ** 3.)) for aBox in box_ar]
        cGrid1s    = np.array(cGrid1s)
        cGrid2s    = [np.log10(np.array(mol_num_ar.T[0], dtype=float) / mol_num_ar.T[1]) for aBox in box_ar]
        cGrid2s    = np.array(cGrid2s)
        conc_grids = np.linspace(cGrid1s.min(), cGrid1s.max(), num_points)
        rat_grids  = np.linspace(cGrid2s.min(), cGrid2s.max(), num_points)
        conc_grids, rat_grids = np.meshgrid(conc_grids, rat_grids)
        cGrid1  = [np.log10(mol_num_ar.T[0] / ((aBox) ** 3.)) for aBox in box_ar]
        cGrid1  = np.array(cGrid1)
        cGrid2  = [np.log10(mol_num_ar.T[1] / ((aBox) ** 3.)) for aBox in box_ar]
        cGrid2  = np.array(cGrid2)
        c1_grid = np.linspace(cGrid1.min(), cGrid1.max(), num_points)
        c2_grid = c1_grid[::-1]
        c1_grid, c2_grid = np.meshgrid(c1_grid, c2_grid)
        xGrid = conc_grids + rat_grids - np.log10(1. + (10. ** rat_grids));
        yGrid = conc_grids - np.log10(1. + (10. ** rat_grids))

        return [conc_grids, rat_grids, xGrid, yGrid]

    def Gen_PlottingDataFor(self, sysName, comp1, comp2, my_TieLineFile=None):
        this_sys = self.SysInfo[sysName]
        if my_TieLineFile != None:
            this_sys['TieLines File'] = my_TieLineFile
        else:
            this_sys['TieLines File'] = sysName + '_TieLines.c.npy'

        tie_file = this_sys['TieLines File']

        tot_tie_line_data = np.load(tie_file)

        tot_plot_data = []
        for tempID, aTemp in enumerate(tot_tie_line_data):
            per_box_data = []
            for boxID, aBox in enumerate(aTemp):
                per_num_data = []
                for numID, aNum in enumerate(aBox):
                    cons_lin = ForPlot_Get_TwoComp_TieLineConcs(tot_tie_line_data, tempID, boxID, numID, comp1, comp2)
                    cons_log = ForPlot_GenLogConcs(cons_lin[0], cons_lin[1])
                    cons_rat = ForPlot_GenRatioConcs(cons_lin[0], cons_lin[1])
                    per_num_data.append([cons_lin, cons_log, cons_rat])
                per_box_data.append(per_num_data)
            tot_plot_data.append(per_box_data)
        #tot_plot_data = np.array(tot_plot_data)

        this_sys['Plotting Data'] = tot_plot_data


class TwoComp_Ortho_Setup:
    def __init__(self, param_file):
        self.GlobalParamFile = param_file
        self.CurrentDir = os.getcwd() + '/'
        self.Set_SimulationPath(self.CurrentDir)
        self.SysInfo = {}

    def Set_SimulationPath(self, end_path_name):
        self.SimulationPath = end_path_name
        print("Simulations shall be done in dir: {:}".format(self.SimulationPath))

    def Set_QSUB_Command(self, qsub_command):
        self.QSubCommand = qsub_command
        self.Reset_QSUB_Iter()

    def Set_QSUB_Queues(self, WIntQueue, NoIntQueue):
        self.QSUB_WIntQ = WIntQueue
        self.QSUB_NoIntQ = NoIntQueue

    def Set_PathToLASSI(self, path_to_lassi):
        self.PathToLASSI = path_to_lassi

    def Reset_QSUB_Iter(self):
        self.QSubIter = 0;

    def Read_StrucFileFor(self, sysName):
        struc_file = self.SysInfo[sysName]['Structure File']
        dum_struc_tot = []
        with open(struc_file) as sFile:
            tot_file = sFile.readlines()
            tot_file = iter(tot_file)
            for a_line in tot_file:
                if a_line[0] == '#' or a_line[0] == '\n':
                    continue
                if a_line[:4] == 'NEW{':
                    dum_struc_in = []
                    a_line = tot_file.next()
                    num_mol = int(a_line[:-1])
                    a_line = tot_file.next()
                    while (a_line[:4] != "}END"):
                        a_line = a_line[:-1].split()
                        this_line = [int(aVal) for aVal in a_line]
                        a_line = tot_file.next()
                        dum_struc_in.append(this_line)
                    dum_struc_tot.append([num_mol, np.array(dum_struc_in)])
        return dum_struc_tot

    def AddNewSystem(self, sysName, lin_len, mol1_min, mol2_min, totmol_max):
        if sysName in self.SysInfo.keys():
            print("This structure name already exists! Doing nothing.")
        else:
            self.SysInfo[sysName] = {}
            self.SysInfo[sysName]['Linker Length'] = lin_len
            self.SysInfo[sysName]['Structure'] = []
            self.SysInfo[sysName]['Mol Mins']  = np.array([mol1_min, mol2_min])
            self.SysInfo[sysName]['Mol Max']   = totmol_max
            self.AddParamFileTo(sysName, self.GlobalParamFile)
            self.SetRunNameFor(sysName, sysName)
            self.SetNumberOfRunsFor(sysName, 2)

    def AddStrucFileTo(self, sysName, file_name):
        try:
            self.SysInfo[sysName]['Structure File'] = file_name
            print("{:} has structure file: {:}".format(sysName, file_name))
            self.SysInfo[sysName]['Key File'][1]['STRUCT_FILE'] = self.CurrentDir + file_name
            self.SysInfo[sysName]['Structure'] = self.Read_StrucFileFor(sysName)
            if len(self.SysInfo[sysName]['Structure'] < 2):
                raise KeyError
            self.CalcNonZeroRDFComps(sysName)
            # self.SysInfo[sysName]['Tot Molecules'] = self.CalcTotalMolecules(sysName)
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def AddEnergyFileTo(self, sysName, fileName_NoInt, fileName_WInt):
        try:
            with open(self.CurrentDir + fileName_WInt) as WIntFile:
                pass
            with open(self.CurrentDir + fileName_NoInt) as NoIntFile:
                pass
        except IOError:
            print("Do these files exist?")
            return
        try:
            self.SysInfo[sysName]['Int Energy File'] = self.CurrentDir + fileName_WInt
            self.SysInfo[sysName]['NoInt Energy File'] = self.CurrentDir + fileName_NoInt
            # print("{:} has energy files: Int:{:}  NoInt:{:}".format(sysName, self.SysInfo[sysName]['Int Energy File'],
            #                                                         self.SysInfo[sysName]['NoInt Energy File']))
            with open(self.CurrentDir + fileName_WInt) as WIntFile:
                WIntFile.readline()  # First line is a comment
                tot_stickers = int(WIntFile.readline())
            self.SysInfo[sysName]['Tot Bead Types'] = tot_stickers
        except KeyError:
            print("Failed! Did you type the correct system name?")
            return

    def AddEnergyFiles_ForAll(self, fileName_NoInt, fileName_WInt):
        for aSys in self.SysInfo.keys():
            self.AddEnergyFileTo(aSys, fileName_NoInt, fileName_WInt)

    def AddParamFileTo(self, sysName, fileName):
        self.SysInfo[sysName]['Key File'] = Read_ParamFile(fileName)

    def AddParamFiles_ForAll(self, fileName):
        for aSys in self.SysInfo.keys():
            self.AddParamFileTo(aSys, fileName)

    def PrintParamsFor(self, sysName):
        try:
            dumKeys = self.SysInfo[sysName]['Key File']
            dumDict = dumKeys[1]
            dumKeys = dumKeys[0]
            for aKey in dumKeys:
                dum_spaces = 25 - len(aKey)
                print('{:}{:}{:}'.format(aKey, ' ' * dum_spaces, dumDict[aKey]))
        except KeyError as myErr:
            if myErr[0] == 'Key File':
                print("Did you import the file?")
            else:
                print("Failed! Did you type the correct system name?")

    def ResetStruc(self, sysName):
        self.SysInfo[sysName]['Structure'] = []

    def AddStruc_ImpLinear(self, sysName, mol_num, bead_num, bead_type, lin_len):
        try:
            self.SysInfo[sysName]['Structure'].append(sGen_Implicit_Linear(mol_num, bead_num, bead_type, lin_len))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_Monomer(self, sysName, mol_num, bead_type):
        try:
            self.SysInfo[sysName]['Structure'].append(sGen_Monomer(mol_num, bead_type))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_Dimer(self, sysName, mol_num, bead_type1, bead_type2, lin_len):
        try:
            self.SysInfo[sysName]['Structure'].append(sGen_Dimer(mol_num, bead_type1, bead_type2, lin_len))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_ExpLinear(self, sysName, mol_num, bead_num, bead_type, lin_len, lin_type):
        try:
            self.SysInfo[sysName]['Structure'].append(
                sGen_Explicit_Linear(mol_num, bead_num, bead_type, lin_len, lin_type))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_ImpBranSymm(self, sysName, mol_num, branch_num, bead_per_branch, bead_type, hub_type, lin_len, hub_len):
        try:
            self.SysInfo[sysName]['Structure'].append(
                sGen_Implicit_SymmBranched(mol_num, branch_num, bead_per_branch, bead_type, hub_type, lin_len, hub_len))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def Write_TwoCompStructureFileFor(self, sysName, file_prefix, mol1_num, mol2_num):
        dummyStruc = self.SysInfo[sysName]['Structure'][:]
        dummyStruc[0][0] = mol1_num
        dummyStruc[1][0] = mol2_num
        file_name = file_prefix + str(mol1_num) + '_' + str(mol2_num) + '.prm'
        sGen_WriteStructuresToFile(dummyStruc, file_name)

    def Write_StructureFilesFor(self, sysName):
        this_sys_num_list = self.SysInfo[sysName]['Molecule Numbers']
        MKDirCatch(sysName + '_Structures/')
        file_prefix = sysName + '_Structures/'
        for mol1, mol2 in this_sys_num_list:
            self.Write_TwoCompStructureFileFor(sysName, file_prefix, mol1, mol2)
        for mol2, mol1 in this_sys_num_list:
            self.Write_TwoCompStructureFileFor(sysName, file_prefix, mol1, mol2)

        self.CalcNonZeroRDFComps(sysName)
        self.CalcCOMDenComps(sysName)
        self.Set_StructureFolderFor(sysName)

    def Set_StructureFolderFor(self, sysName):
        file_prefix = sysName + '_Structures/'
        self.SysInfo[sysName]['Structure Folder'] = self.CurrentDir + file_prefix;

    def Write_StructureFiles_ForAll(self):
        for aSys in self.SysInfo.keys():
            Write_StructureFilesFor(aSys, aSys)

    def SetTemperaturesFor(self, sysName, init_temp, final_temp, temp_steps, therm_temp):
        assert (type(temp_steps) == int) and (type(init_temp) == type(final_temp) == float)
        assert temp_steps > 0. and init_temp > 0. and final_temp > 0.
        if temp_steps > 1:
            delta_temp = (final_temp - init_temp) / (temp_steps - 1)
        else:
            delta_temp = 0.0
        try:
            key_dum = self.SysInfo[sysName]['Key File'][1]
            key_dum['MC_TEMP'] = init_temp
            key_dum['MC_CYCLE_NUM'] = temp_steps
            key_dum['MC_DELTA_TEMP'] = delta_temp
            key_dum['PREEQ_TEMP'] = therm_temp
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetTemperatures_ForAll(self, init_temp, final_temp, temp_steps, therm_temp):
        for aSys in self.SysInfo.keys():
            self.SetTemperaturesFor(aSys, init_temp, final_temp, temp_steps, therm_temp)

    def SetMCStepsFor(self, sysName, therm_steps, run_steps):
        try:
            key_dum = self.SysInfo[sysName]['Key File'][1]
            key_dum['N_STEPS'] = int(run_steps)
            key_dum['PREEQ_STEPS'] = int(therm_steps)
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetMCSteps_ForAll(self, therm_steps, run_steps):
        for aSys in self.SysInfo.keys():
            self.SetMCStepsFor(aSys, therm_steps, run_steps)

    def CalcMolNumArr(self, sysName, tot_bin_num):
        assert tot_bin_num % 2 == 1  # We should always have an odd number input
        actual_bins = (tot_bin_num) / 2 + 1
        num_min = self.SysInfo[sysName]['Mol Mins'].max()
        num_max = self.SysInfo[sysName]['Mol Max']
        n1Mols = np.ceil(10.0 ** np.linspace(np.log10(num_min), np.log10(num_max / 2.), actual_bins))
        n1Mols = np.array(n1Mols, dtype=int)
        mol_ar = np.sort(np.append(n1Mols, num_max - n1Mols[-2::-1]))
        self.SysInfo[sysName]['Molecule Numbers'] = np.array([mol_ar, num_max - mol_ar]).T

    def CalcMolNum_Linear_Arr(self, sysName, tot_bin_num):
        assert tot_bin_num % 2 == 1  # We should always have an odd number input
        num_min = self.SysInfo[sysName]['Mol Mins'].max()
        num_max = self.SysInfo[sysName]['Mol Max']
        n1Mols = np.linspace(num_min, num_max-num_min, tot_bin_num)
        n1Mols = np.array(n1Mols, dtype=int)
        n2Mols = num_max - n1Mols
        self.SysInfo[sysName]['Molecule Numbers'] = np.array([n1Mols, n2Mols]).T

    def CalcHybridMolNum_Arr(self, sysName, tot_bin_num):
        assert tot_bin_num % 2 == 1  # We should always have an odd number input
        #First find linear
        num_min = self.SysInfo[sysName]['Mol Mins'].max()
        num_max = self.SysInfo[sysName]['Mol Max']
        linMols  = np.linspace(num_min, num_max - num_min, tot_bin_num+1)
        linMols  = np.array(linMols, dtype=int)
        #Now log-space
        actual_bins = (tot_bin_num) / 2 + 1
        logMols = np.ceil(10.0 ** np.linspace(np.log10(num_min), np.log10(num_max/2), actual_bins))

        logMols = np.array(logMols, dtype=int)
        logMols = np.sort(np.append(logMols, num_max - logMols[-2::-1]))
        #Combine the two
        totmAr = np.unique(np.append(linMols, logMols))
        totmAr = np.array([totmAr, totmAr[::-1]]).T
        self.SysInfo[sysName]['Molecule Numbers'] = totmAr[::2]

    def SetBoxSizesFor(self, sysName, target_con, low_con, high_con, mol_num_boxes):
        self.CalcHybridMolNum_Arr(sysName, mol_num_boxes)
        if not 0.0 < low_con < high_con < 1.0:
            print("Low conc should be lower than high conc!"
                  "Both numbers need to be in the interval (0,1).")
            return
        raw_molecule_ar = self.SysInfo[sysName]['Molecule Numbers'][:]

        tot_new_box = []
        tot_conc_ar = []
        for aPair in raw_molecule_ar:
            this_box_raw = aPair[0]/target_con
            this_box_raw = np.ceil(this_box_raw)**(1./3.)
            this_box_raw = int(this_box_raw)
            tot_new_box.append(this_box_raw)
            tot_conc_ar.append((aPair[0]+aPair[1])/(this_box_raw**3.))
        tot_conc_ar = np.array(tot_conc_ar)
        tot_new_box = np.array(tot_new_box)
        boxes_uni, ids_uni = np.unique(tot_new_box, return_index=True)
        mols_uni = raw_molecule_ar[ids_uni]
        conc_uni = tot_conc_ar[ids_uni]
        idx_start = np.argwhere(conc_uni >= high_con).T[0]
        idx_end   = np.argwhere(conc_uni <= low_con).T[0]
        if len(idx_start) > 0:
            idx_start = idx_start[-1]
        else:
            idx_start = 0;
        if len(idx_end) > 0:
            idx_end = idx_end[0]
        else:
            idx_end = len(conc_uni)

        new_box_ar = boxes_uni[idx_start:idx_end]
        new_mol_ar = mols_uni[idx_start:idx_end]
        self.SysInfo[sysName]['Boxes'] = new_box_ar[:]
        self.SysInfo[sysName]['Molecule Numbers'] = new_mol_ar[:]

    def CalcMaxBeadsFor(self, sysName):
        # Figure out the maximum possible number of beads given the structure and molecule numbers.
        try:
            struc_ar = self.SysInfo[sysName]['Structure']
            mol_nums = [];
            bead_nums = [];
            for a_mol in struc_ar:
                num_mol = a_mol[0]
                struc_list = a_mol[1].T[0]
                num_beads_per = len(np.unique(struc_list))
                bead_nums.append(num_beads_per)
                mol_nums.append(num_mol)
            mol_nums = np.array(mol_nums)
            bead_nums = np.array(bead_nums)
            largest_mol = np.argmax(bead_nums)
            smaller_mol = (largest_mol + 1) % 2
            mol_num_arr = self.SysInfo[sysName]['Molecule Numbers']
            largest_num = np.argmax(mol_num_arr.T[largest_mol])
            smaller_num = mol_num_arr[largest_num][smaller_mol]
            largest_num = mol_num_arr[largest_num][largest_mol]
            mol_nums[largest_mol] = largest_num;
            mol_nums[smaller_mol] = smaller_num;

            return np.dot(mol_nums, bead_nums)
        except KeyError:
            print("Failed! Did you type the correct system name?")
            raise

    def CalcNonZeroRDFComps(self, sysName):
        try:
            struc_list = self.SysInfo[sysName]['Structure']
        except KeyError:
            print("Did you type the system name correctly?")
            raise

        bead_type_list = []
        for a_mol in struc_list:
            bead_types = a_mol[1].T[1]
            bead_type_list.append(np.unique(bead_types))
        bead_type_list = [aType for a_mol in bead_type_list for aType in a_mol]
        bead_types = np.array(bead_type_list)
        bead_types = bead_types.flatten()
        bead_types = np.unique(bead_types)
        bead_type_list = [[-1, -1]]  # For the total
        for i in range(len(bead_types)):
            for j in range(i, len(bead_types)):
                bead_type_list.append([bead_types[i], bead_types[j]])
        self.SysInfo[sysName]['Comp List'] = np.array(bead_type_list, dtype=int)

    def CalcCOMDenComps(self, sysName):
        try:
            struc_list = self.SysInfo[sysName]['Structure']
        except KeyError:
            print("Did you type the system name correctly?")
            raise

        num_mols = len(struc_list) #Total number of molecules
        dum_comp_list = [];
        for i in range(num_mols):
            dum_comp_list.append([-1, i])
        for i in range(num_mols):
            for j in range(num_mols):
                dum_comp_list.append([i, j])
        self.SysInfo[sysName]['COMDen List'] = np.array(dum_comp_list, dtype=int)

    def SetBoxSizes_ForAll(self, target_con, low_con, high_con, mol_num_boxes):
        for aSys in self.SysInfo.keys():
            self.SetBoxSizesFor(aSys, target_con, low_con, high_con, mol_num_boxes)

    def SetBoxSizesTo(self, thisSys, SysToCopy):
        box_ar = self.SysInfo[SysToCopy]['Boxes'][:]+0
        mol_ar = self.SysInfo[SysToCopy]['Molecule Numbers'][:] + 0
        self.SysInfo[thisSys]['Boxes'] = box_ar[:]
        self.SysInfo[thisSys]['Molecule Numbers'] = mol_ar[:]

    def ShowHybridSamplingGridFor(self, sysName, lower_lims, upper_lims):
        myFig = plt.figure(figsize=[12, 6])
        myAx  = myFig.subplots(1,2)
        mol_num = self.SysInfo[sysName]['Molecule Numbers']
        box_arr = self.SysInfo[sysName]['Boxes']
        run_num = self.SysInfo[sysName]['Runs']
        N1_arr = mol_num.T[0]
        N2_arr = mol_num.T[1]
        print(
            r'A total of {:} independent conditions! With {:} runs per condition, there are a total of {:} independent simulations (WInt & NoInt).'.format(
                len(box_arr)*2-1, run_num, (len(box_arr)*2-1) * run_num * 2))

        for aBox, n1, n2 in zip(box_arr, N1_arr, N2_arr):
            myAx[0].plot(np.log10(n1 / (aBox ** 3.)), np.log10(n2 / (aBox ** 3.)), 'ro', alpha=0.9)
            myAx[1].plot((n1 / (aBox ** 3.)), (n2 / (aBox ** 3.)), 'ro', alpha=0.9)
        for aBox, n1, n2 in zip(box_arr, N2_arr, N1_arr):
            myAx[0].plot(np.log10(n1 / (aBox ** 3.)), np.log10(n2 / (aBox ** 3.)), 'bo', alpha=0.9)
            myAx[1].plot((n1 / (aBox ** 3.)), (n2 / (aBox ** 3.)), 'bo', alpha=0.9)
        [anAx.plot(np.linspace(-10, 10, 10), np.linspace(-10, 10, 10), 'k--', lw=3, alpha=0.5) for anAx in myAx]
        myAx[1].set_xlim(lower_lims, upper_lims)
        myAx[1].set_ylim(lower_lims, upper_lims)
        myAx[0].set_xlim(np.log10(lower_lims), np.log10(upper_lims))
        myAx[0].set_ylim(np.log10(lower_lims), np.log10(upper_lims))

        [anAx.grid(lw=2, alpha=0.7) for anAx in myAx]
        plt.tight_layout()
        plt.show()

    def SetRunNameFor(self, sysName, runName):
        try:
            self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX'] = runName
        except KeyError:
            print("Failed! Did you type the correct system name?"
                  " Does this system exist?")

    def SetRunName_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.SetRunNameFor(aSys, aSys)

    def SetNumberOfRunsFor(self, sysName, run_num):
        try:
            self.SysInfo[sysName]['Runs'] = run_num
        except KeyError:
            print("Failed! Did you type the correct system name?"
                  " Does this system exist?")

    def SetNumberOfRuns_ForAll(self, run_num):
        for aSys in self.SysInfo.keys():
            self.SetNumberOfRunsFor(aSys, run_num)

    def MakeDirs_NoIntFor(self, sysName):

        dum_dir = self.SimulationPath
        dum_dir += 'OrthoLines/'
        MKDirCatch(dum_dir)
        dum_dir += sysName + '/'
        MKDirCatch(dum_dir)

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']

        for b_ID, a_tot_set in enumerate(zip(box_arr, num_arr)):
            a_box    = a_tot_set[0]
            a_con    = a_tot_set[1]
            this_dir = dum_dir + str(a_box) + '/'
            MKDirCatch(this_dir)
            # First Order
            con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/'
            MKDirCatch(con_dir)
            con_dir += 'NoInt/'
            MKDirCatch(con_dir)
            for a_run in range(num_of_runs):
                final_dir = con_dir + str(a_run + 1) + '/'
                MKDirCatch(final_dir)
            if a_con[0] != a_con[1]:
                # Reverse Order
                con_dir = this_dir + str(a_con[1]) + '_' + str(a_con[0]) + '/'
                MKDirCatch(con_dir)
                con_dir += 'NoInt/'
                MKDirCatch(con_dir)
                for a_run in range(num_of_runs):
                    final_dir = con_dir + str(a_run + 1) + '/'
                    MKDirCatch(final_dir)

    def MakeDirs_WIntFor(self, sysName):

        dum_dir = self.SimulationPath
        dum_dir += 'OrthoLines/'
        MKDirCatch(dum_dir)
        dum_dir += sysName + '/'
        MKDirCatch(dum_dir)

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']

        for b_ID, a_tot_set in enumerate(zip(box_arr, num_arr)):
            a_box = a_tot_set[0]
            a_con = a_tot_set[1]
            this_dir = dum_dir + str(a_box) + '/'
            MKDirCatch(this_dir)
            # First Order
            con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/'
            MKDirCatch(con_dir)
            con_dir += 'WInt/'
            MKDirCatch(con_dir)
            for a_run in range(num_of_runs):
                final_dir = con_dir + str(a_run + 1) + '/'
                MKDirCatch(final_dir)
            if a_con[0] != a_con[1]:
                # Reverse Order
                con_dir = this_dir + str(a_con[1]) + '_' + str(a_con[0]) + '/'
                MKDirCatch(con_dir)
                con_dir += 'WInt/'
                MKDirCatch(con_dir)
                for a_run in range(num_of_runs):
                    final_dir = con_dir + str(a_run + 1) + '/'
                    MKDirCatch(final_dir)

    def MakeDirs_For(self, sysName):
        self.MakeDirs_WIntFor(sysName)
        self.MakeDirs_NoIntFor(sysName)

    def MakeDirs_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.MakeDirs_For(aSys)

    def MakeDirs_ForAll_NoInt(self):
        for aSys in self.SysInfo.keys():
            self.MakeDirs_NoIntFor(aSys)

    def MakeDirs_ForAll_WInt(self):
        for aSys in self.SysInfo.keys():
            self.MakeDirs_WIntFor(aSys)

    def Write_ParamsWIntFor(self, sysName):
        dum_dir = self.SimulationPath+'OrthoLines/'+sysName+'/'

        param_copy = self.SysInfo[sysName]['Key File'][:]
        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']
        param_copy[1]['ENERGY_FILE'] = self.SysInfo[sysName]['Int Energy File']
        struc_folder = self.SysInfo[sysName]['Structure Folder']

        for b_ID, a_tot_set in enumerate(zip(box_arr, num_arr)):
            a_box = a_tot_set[0]
            a_con = a_tot_set[1]
            param_copy[1]['BOX_SIZE'] = a_box
            this_dir = dum_dir + str(a_box) + '/'
            # First Order
            con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
            param_copy[1]['STRUCT_FILE'] = struc_folder + str(a_con[0]) + '_' + str(a_con[1]) + '.prm'
            for a_run in range(num_of_runs):
                param_copy[1]['RANDOM_SEED'] = 0
                final_dir = con_dir + str(a_run + 1) + '/'
                with open(final_dir + 'param.key', 'w+') as pFile:
                    for a_key in param_copy[0]:
                        N_spcs = 25 - len(a_key)
                        pFile.write(a_key + ' ' * N_spcs + str(param_copy[1][a_key]) + '\n')
            if a_con[0] != a_con[1]:
                con_dir = this_dir + str(a_con[1]) + '_' + str(a_con[0]) + '/WInt/'
                param_copy[1]['STRUCT_FILE'] = struc_folder + str(a_con[1]) + '_' + str(a_con[0]) + '.prm'
                for a_run in range(num_of_runs):
                    param_copy[1]['RANDOM_SEED'] = 0
                    final_dir = con_dir + str(a_run + 1) + '/'
                    with open(final_dir + 'param.key', 'w+') as pFile:
                        for a_key in param_copy[0]:
                            N_spcs = 25 - len(a_key)
                            pFile.write(a_key + ' ' * N_spcs + str(param_copy[1][a_key]) + '\n')

    def Write_ParamsNoIntFor(self, sysName):
        dum_dir = self.SimulationPath+'OrthoLines/'+sysName+'/'

        param_copy = self.SysInfo[sysName]['Key File'][:]

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']
        param_copy[1]['ENERGY_FILE'] = self.SysInfo[sysName]['NoInt Energy File']
        param_copy[1]['MV_STROT_FREQ'] = '0.0';
        param_copy[1]['MV_COLOCAL_FREQ'] = '0.0';
        param_copy[1]['MV_SMCLSTR_FREQ'] = '0.0';
        param_copy[1]['MV_CLSTR_FREQ'] = '0.0';
        param_copy[1]['REPORT_NETWORK_FREQ'] = '0';
        struc_folder = self.SysInfo[sysName]['Structure Folder']

        for b_ID, a_tot_set in enumerate(zip(box_arr, num_arr)):
            a_box = a_tot_set[0]
            a_con = a_tot_set[1]
            param_copy[1]['BOX_SIZE'] = a_box
            this_dir = dum_dir + str(a_box) + '/'
            # First Order
            con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/NoInt/'
            param_copy[1]['STRUCT_FILE'] = struc_folder + str(a_con[0]) + '_' + str(a_con[1]) + '.prm'
            for a_run in range(num_of_runs):
                param_copy[1]['RANDOM_SEED'] = 0
                final_dir = con_dir + str(a_run + 1) + '/'
                with open(final_dir + 'param.key', 'w+') as pFile:
                    for a_key in param_copy[0]:
                        N_spcs = 25 - len(a_key)
                        pFile.write(a_key + ' ' * N_spcs + str(param_copy[1][a_key]) + '\n')
            # Reverse Order
            if a_con[0] != a_con[1]:
                con_dir = this_dir + str(a_con[1]) + '_' + str(a_con[0]) + '/NoInt/'
                param_copy[1]['STRUCT_FILE'] = struc_folder + str(a_con[1]) + '_' + str(a_con[0]) + '.prm'
                for a_run in range(num_of_runs):
                    param_copy[1]['RANDOM_SEED'] = 0
                    final_dir = con_dir + str(a_run + 1) + '/'
                    with open(final_dir + 'param.key', 'w+') as pFile:
                        for a_key in param_copy[0]:
                            N_spcs = 25 - len(a_key)
                            pFile.write(a_key + ' ' * N_spcs + str(param_copy[1][a_key]) + '\n')

    def Write_ParamsFor(self, sysName):
        self.Write_ParamsWIntFor(sysName)
        self.Write_ParamsNoIntFor(sysName)

    def Write_ParamsWInt_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.Write_ParamsWIntFor(aSys)

    def Write_ParamsNoInt_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.Write_ParamsNoIntFor(aSys)

    def Write_Params_ForAll(self):
        for aSys in self.SysInfo.keys():
            self.Write_ParamsWIntFor(aSys)
            self.Write_ParamsNoIntFor(aSys)

    def SubmitWIntJobs_ToQueueFor(self, sysName, WIntQueue):
        dum_dir = self.SimulationPath+'OrthoLines/'+sysName+'/'

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']

        for b_ID, a_tot_set in enumerate(zip(box_arr, num_arr)):
            a_box = a_tot_set[0]
            a_con = a_tot_set[1]
            this_dir = dum_dir + str(a_box) + '/'

            #First Order
            con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
            for a_run in range(num_of_runs):
                final_dir = con_dir + str(a_run + 1) + '/'
                run_dir = final_dir
                run_num = ['W' + sysName[:2], b_ID, self.QSubIter]
                run_num = [str(a_val) for a_val in run_num]
                run_num = "_".join(run_num)
                run_command = self.QSubCommand + ' ' + WIntQueue + ' ' + run_num + ' /.' + self.PathToLASSI
                try:
                    os.chdir(run_dir)
                except IOError:
                    print("Did you make all the directories?")
                    raise
                try:
                    ret_code = sproc.check_output(run_command, shell=True, stderr=sproc.STDOUT)
                    if ret_code < 0:
                        print("The submission failed! Signal: ", -retcode)
                    else:
                        print(ret_code)
                except OSError as myErr:
                    print("Couldn't submit because ", myErr)
                    raise
                time.sleep(2)
                self.QSubIter += 1
            #Reverse Order
            if a_con[0] != a_con[1]:
                con_dir = this_dir + str(a_con[1]) + '_' + str(a_con[0]) + '/WInt/'
                for a_run in range(num_of_runs):
                    final_dir = con_dir + str(a_run + 1) + '/'
                    run_dir = final_dir
                    run_num = ['W' + sysName[:2], b_ID, self.QSubIter]
                    run_num = [str(a_val) for a_val in run_num]
                    run_num = "_".join(run_num)
                    run_command = self.QSubCommand + ' ' + WIntQueue + ' ' + run_num + ' /.' + self.PathToLASSI
                    try:
                        os.chdir(run_dir)
                    except IOError:
                        print("Did you make all the directories?")
                        raise
                    try:
                        ret_code = sproc.check_output(run_command, shell=True, stderr=sproc.STDOUT)
                        if ret_code < 0:
                            print("The submission failed! Signal: ", -retcode)
                        else:
                            print(ret_code)
                    except OSError as myErr:
                        print("Couldn't submit because ", myErr)
                        raise
                    time.sleep(2)
                    self.QSubIter += 1
        os.chdir(self.CurrentDir)

    def SubmitNoIntJobs_ToQueueFor(self, sysName, NoIntQueue):
        dum_dir = self.SimulationPath+'OrthoLines/'+sysName+'/'

        box_arr = self.SysInfo[sysName]['Boxes'];
        num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        num_of_runs = self.SysInfo[sysName]['Runs']

        for b_ID, a_tot_set in enumerate(zip(box_arr, num_arr)):
            a_box = a_tot_set[0]
            a_con = a_tot_set[1]
            this_dir = dum_dir + str(a_box) + '/'

            #First Order
            con_dir = this_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/NoInt/'
            for a_run in range(num_of_runs):
                final_dir = con_dir + str(a_run + 1) + '/'
                run_dir = final_dir
                run_num = ['N' + sysName[:2], b_ID, self.QSubIter]
                run_num = [str(a_val) for a_val in run_num]
                run_num = "_".join(run_num)
                run_command = self.QSubCommand + ' ' + NoIntQueue + ' ' + run_num + ' /.' + self.PathToLASSI
                try:
                    os.chdir(run_dir)
                except IOError:
                    print("Did you make all the directories?")
                    raise
                try:
                    ret_code = sproc.check_output(run_command, shell=True, stderr=sproc.STDOUT)
                    if ret_code < 0:
                        print("The submission failed! Signal: ", -retcode)
                    else:
                        print(ret_code)
                except OSError as myErr:
                    print("Couldn't submit because ", myErr)
                    raise
                time.sleep(2)
                self.QSubIter += 1

            #Reverse Order
            if a_con[0] != a_con[1]:
                con_dir = this_dir + str(a_con[1]) + '_' + str(a_con[0]) + '/NoInt/'
                for a_run in range(num_of_runs):
                    final_dir = con_dir + str(a_run + 1) + '/'
                    run_dir = final_dir
                    run_num = ['N' + sysName[:2], b_ID, self.QSubIter]
                    run_num = [str(a_val) for a_val in run_num]
                    run_num = "_".join(run_num)
                    run_command = self.QSubCommand + ' ' + NoIntQueue + ' ' + run_num + ' /.' + self.PathToLASSI
                    try:
                        os.chdir(run_dir)
                    except IOError:
                        print("Did you make all the directories?")
                        raise
                    try:
                        ret_code = sproc.check_output(run_command, shell=True, stderr=sproc.STDOUT)
                        if ret_code < 0:
                            print("The submission failed! Signal: ", -retcode)
                        else:
                            print(ret_code)
                    except OSError as myErr:
                        print("Couldn't submit because ", myErr)
                        raise
                    time.sleep(2)
                    self.QSubIter += 1

        os.chdir(self.CurrentDir)

    def SubmitWIntJobs_ForAll(self, WIntQueue):
        for aSys in self.SysInfo.keys():
            self.SubmitWIntJobs_ToQueueFor(aSys, WIntQueue)

    def SubmitNoIntJobs_ForAll(self, NoIntQueue):
        for aSys in self.SysInfo.keys():
            self.SubmitWIntJobs_ToQueueFor(aSys, NoIntQueue)

    def SubmitJobs_ForAll(self):
        try:
            self.QSubCommand
            self.QSUB_WIntQ
            self.QSUB_NoIntQ
        except AttributeError:
            print("You need to set up the qsub command, and which queues to use!")
            return
        self.Reset_QSUB_Iter()
        for aSys in self.SysInfo.keys():
            self.SubmitNoIntJobs_ToQueueFor(aSys, self.QSUB_NoIntQ)
            self.SubmitWIntJobs_ToQueueFor(aSys, self.QSUB_WIntQ)

    def SetAndRun_For(self, sysName):
        self.MakeDirs_For(sysName)
        self.Write_ParamsFor(sysName)
        self.Reset_QSUB_Iter()
        self.SubmitWIntJobs_ToQueueFor(sysName, self.QSUB_WIntQ)
        self.SubmitNoIntJobs_ToQueueFor(sysName, self.QSUB_NoIntQ)

    def SetAndRunWInt_For(self, sysName):
        self.MakeDirs_WIntFor(sysName)
        self.Write_ParamsWIntFor(sysName)
        self.Reset_QSUB_Iter()
        self.SubmitWIntJobs_ToQueueFor(sysName, self.QSUB_WIntQ)

    def SetAndRunNoInt_For(self, sysName):
        self.MakeDirs_NoIntFor(sysName)
        self.Write_ParamsNoIntFor(sysName)
        self.Reset_QSUB_Iter()
        self.SubmitNoIntJobs_ToQueueFor(sysName, self.QSUB_NoIntQ)

class TwoComo_Ortho_Analysis(TwoComp_Ortho_Setup):
    def __init__(self, SimSetupInstance):
        try:
            self.SysInfo = SimSetupInstance.SysInfo
            self.GlobalParamFile = SimSetupInstance.GlobalParamFile
            self.SimulationPath = SimSetupInstance.SimulationPath
            self.CurrentDir = SimSetupInstance.CurrentDir
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise

    def Read_PDFFileFor(self, sysName, file_name, box_size, tot_temps):
        try:
            dum_dat = np.loadtxt(file_name)
            comp_list = self.SysInfo[sysName]['Comp List']
            tot_comps = self.SysInfo[sysName]['Tot Bead Types']
            tot_comps_pos = 1 + Index_RDF(tot_comps - 2, tot_comps - 1, tot_comps)
        except IOError:
            print("Does the PDF file exist? Did these runs finish and output data?")
            print("Failed for "+file_name)
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        actual_comps = len(comp_list)
        ret_dat = np.zeros((actual_comps, tot_temps, box_size * 4))
        for a_temp in range(tot_temps):
            for compID, compPair in enumerate(comp_list):
                compA = int(compPair[0]);
                compB = int(compPair[1])
                ret_dat[compID][a_temp] = dum_dat[tot_comps_pos * a_temp + Index_RDF(compA, compB, tot_comps)]
        return ret_dat

    def Read_CLUSFileFor(self, file_name):
        try:
            dum_dat = np.loadtxt(file_name)
        except IOError:
            print("Does the CLUS file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Read_COMDenFileFor(self, sysName, file_name, box_size, tot_temps):
        try:
            dum_dat   = np.loadtxt(file_name)
            comp_list = self.SysInfo[sysName]['COMDen List']
            num_mols = len(np.unique(comp_list.T[0]))-1
            tot_comps = num_mols*(num_mols+1)
        except IOError:
            print("Does the PDF file exist? Did these runs finish and output data?")
            print("Failed for " + file_name)
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        ret_dat = np.zeros((tot_comps, tot_temps, box_size * 4))
        for a_temp in range(tot_temps):
            for compID, compPair in enumerate(comp_list):
                compA = int(compPair[0]);
                compB = int(compPair[1])
                ret_dat[compID][a_temp] = dum_dat[tot_comps * a_temp + Index_COMDen(compA, compB, num_mols)]
        return ret_dat

    def Read_MolClusFileFor(self, file_name):
        try:
            dum_dat = np.loadtxt(file_name)
        except IOError:
            print("Does the MolClus file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Read_GYRADFileFor(self, file_name):
        try:
            dum_dat = np.loadtxt(file_name)
        except IOError:
            print("Does the GYRRAD file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Save_NoIntPDF_For(self, sysName):
        try:
            dum_dat = self.Collect_NoInt_PDFs_For(sysName)
            file_name = sysName + 'ORTH_N_PDF.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_NoInt_PDF'] = file_name
            pass
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntPDF_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_PDFs_For(sysName)
            file_name = sysName + 'ORTH_W_PDF.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_WInt_PDF'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntCOMDen_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_COMDens_For(sysName)
            file_name = sysName + 'ORTH_W_COMDen.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_WInt_COMDen'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_NoIntCOMDen_For(self, sysName):
        try:
            dum_dat = self.Collect_NoInt_COMDens_For(sysName)
            file_name = sysName + 'ORTH_N_COMDen.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_NoInt_COMDen'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntMolClus_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_MolClus_For(sysName)
            file_name = sysName + 'ORTH_W_MolClus.b'
            with open(file_name, "wb+") as save_file:
                pickle.dump(dum_dat, save_file)
            self.SysInfo[sysName]['Raw_WInt_MolClus'] = file_name
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Save_WIntClus_For(self, sysName):
        try:
            dum_dat = self.Collect_WInt_CLUS_For(sysName)
            file_name = sysName + 'ORTH_CLUS.b'
            np.save(file_name, dum_dat)
            self.SysInfo[sysName]['Raw_CLUS'] = file_name + '.npy'
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Set_Auto_PDFFileNames_For(self, sysName):
        self.SysInfo[sysName]['Raw_NoInt_PDF'] = sysName + 'ORTH_N_PDF.b'
        self.SysInfo[sysName]['Raw_WInt_PDF'] = sysName + 'ORTH_W_PDF.b'

    def Collect_NoInt_PDFs_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'PhaseDiags/'+sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        tot_raw_dat = []
        for b_ID, a_box in enumerate(box_arr):
            box_dir = over_dir + str(a_box)+'/'
            tot_raw_per_box = []
            for a_con in num_arr:
                con_dir = box_dir + str(a_con[0])+'_'+str(a_con[1])+'/NoInt/'
                tot_raw_dat_temp = []
                for a_run in range(run_num):
                    file_name = con_dir + str(a_run + 1) + '/' + run_name + '_RDF.dat'
                    tot_raw_dat_temp.append(self.Read_PDFFileFor(sysName, file_name, a_box, tot_temps))
                tot_raw_per_box.append(tot_raw_dat_temp)
            tot_raw_dat.append(tot_raw_per_box)
        return tot_raw_dat

    def Collect_WInt_PDFs_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'PhaseDiags/' + sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        tot_raw_dat = []
        for b_ID, a_box in enumerate(box_arr):
            box_dir = over_dir + str(a_box) + '/'
            tot_raw_per_box = []
            for a_con in num_arr:
                con_dir = box_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
                tot_raw_dat_temp = []
                for a_run in range(run_num):
                    file_name = con_dir + str(a_run + 1) + '/' + run_name + '_RDF.dat'
                    tot_raw_dat_temp.append(self.Read_PDFFileFor(sysName, file_name, a_box, tot_temps))
                tot_raw_per_box.append(tot_raw_dat_temp)
            tot_raw_dat.append(tot_raw_per_box)
        return tot_raw_dat

    def Collect_WInt_COMDens_For(self, sysName):
        try:
            box_arr   = self.SysInfo[sysName]['Boxes']
            key_dum   = self.SysInfo[sysName]['Key File'][1]
            over_dir  = self.SimulationPath + 'OrthoLines/' + sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name  = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num   = self.SysInfo[sysName]['Runs']
            num_arr   = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        # First component line
        tot_raw_dat = []
        tot_raw_per_box = []
        for b_ID, a_tot_set in enumerate(zip(box_arr, num_arr)):
            a_box = a_tot_set[0]
            a_con = a_tot_set[1]
            box_dir = over_dir + str(a_box) + '/'
            con_dir = box_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = con_dir + str(a_run + 1) + '/' + run_name + '_COMDen.dat'
                tot_raw_dat_temp.append(self.Read_COMDenFileFor(sysName, file_name, a_box, tot_temps))
            tot_raw_per_box.append(tot_raw_dat_temp)
        tot_raw_dat.append(tot_raw_per_box)

        # Second component line
        tot_raw_per_box = []
        for b_ID, a_tot_set in enumerate(zip(box_arr, num_arr)):
            a_box = a_tot_set[0]
            a_con = a_tot_set[1]
            box_dir = over_dir + str(a_box) + '/'
            con_dir = box_dir + str(a_con[1]) + '_' + str(a_con[0]) + '/WInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = con_dir + str(a_run + 1) + '/' + run_name + '_COMDen.dat'
                tot_raw_dat_temp.append(self.Read_COMDenFileFor(sysName, file_name, a_box, tot_temps))
            tot_raw_per_box.append(tot_raw_dat_temp)

        tot_raw_dat.append(tot_raw_per_box)

        return tot_raw_dat

    def Collect_NoInt_COMDens_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'OrthoLines/' + sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        # First component line
        tot_raw_dat = []
        tot_raw_per_box = []
        for b_ID, a_tot_set in enumerate(zip(box_arr, num_arr)):
            a_box = a_tot_set[0]
            a_con = a_tot_set[1]
            box_dir = over_dir + str(a_box) + '/'
            con_dir = box_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/NoInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = con_dir + str(a_run + 1) + '/' + run_name + '_COMDen.dat'
                tot_raw_dat_temp.append(self.Read_COMDenFileFor(sysName, file_name, a_box, tot_temps))
            tot_raw_per_box.append(tot_raw_dat_temp)

        tot_raw_dat.append(tot_raw_per_box)

        # Second component line
        tot_raw_per_box = []
        for b_ID, a_tot_set in enumerate(zip(box_arr, num_arr)):
            a_box = a_tot_set[0]
            a_con = a_tot_set[1]
            box_dir = over_dir + str(a_box) + '/'
            con_dir = box_dir + str(a_con[1]) + '_' + str(a_con[0]) + '/NoInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = con_dir + str(a_run + 1) + '/' + run_name + '_COMDen.dat'
                tot_raw_dat_temp.append(self.Read_COMDenFileFor(sysName, file_name, a_box, tot_temps))
            tot_raw_per_box.append(tot_raw_dat_temp)

        tot_raw_dat.append(tot_raw_per_box)

        return tot_raw_dat

    def Collect_WInt_MolClus_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'OrthoLines/' + sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        tot_raw_dat = []
        for b_ID, a_box in enumerate(box_arr):
            box_dir = over_dir + str(a_box) + '/'
            tot_raw_per_box = []
            for a_con in num_arr:
                con_dir = box_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
                tot_raw_dat_temp = []
                for a_run in range(run_num):
                    file_name = con_dir + str(a_run + 1) + '/' + run_name + '_MolClus.dat'
                    tot_raw_dat_temp.append(self.Read_MolClusFileFor(file_name))
                tot_raw_per_box.append(tot_raw_dat_temp)
            tot_raw_dat.append(tot_raw_per_box)
        return tot_raw_dat

    def Collect_WInt_CLUS_For(self, sysName):
        try:
            box_arr = self.SysInfo[sysName]['Boxes']
            key_dum = self.SysInfo[sysName]['Key File'][1]
            over_dir = self.SimulationPath + 'OrthoLines/' + sysName + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num = self.SysInfo[sysName]['Runs']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of TwoComp_Setup!")
            raise
        tot_raw_dat = []
        for b_ID, a_box in enumerate(box_arr):
            box_dir = over_dir + str(a_box) + '/'
            tot_raw_per_box = []
            for a_con in num_arr:
                con_dir = box_dir + str(a_con[0]) + '_' + str(a_con[1]) + '/WInt/'
                tot_raw_dat_temp = []
                for a_run in range(run_num):
                    file_name = con_dir + str(a_run + 1) + '/' + run_name + '_CLUS.dat'
                    tot_raw_dat_temp.append(self.Read_CLUSFileFor(file_name))
                tot_raw_per_box.append(tot_raw_dat_temp)
            tot_raw_dat.append(tot_raw_per_box)
        return tot_raw_dat

    def Gen_pRDFs_For(self, sysName, my_WInt=None, my_NoInt=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            self.Set_Auto_PDFFileNames_For(sysName)
        except KeyError:
            print("Is this the correct system?")
            return
        if my_WInt != None:
            this_sys['Raw_WInt_PDF'] = my_WInt
        if my_NoInt != None:
            this_sys['Raw_NoInt_PDF'] = my_NoInt
        nFile = this_sys['Raw_NoInt_PDF']
        iFile = this_sys['Raw_WInt_PDF']
        with open(iFile) as dumFile:
            i_pdfs = pickle.load(dumFile)
        with open(nFile) as dumFile:
            n_pdfs = pickle.load(dumFile)

        tot_rdf_arr = []
        for boxID, a_box in enumerate(box_ar):
            x_ar = np.arange(0, a_box, 0.25)
            n_ar = np.array(n_pdfs[boxID]);
            i_ar = np.array(i_pdfs[boxID]);
            n_run_ave, n_run_err = MeanAndError(n_ar, 3);
            n_temp_ave, n_temp_err = MeanAndError(n_run_ave, 1);
            n_temp_err = np.sqrt(n_temp_err**2. + np.sum(n_run_err**2.,axis=1))
            i_run_ave, i_run_err = MeanAndError(i_ar, 1);
            n_sup_ave_for_zeros = np.mean(np.mean(n_temp_ave, axis = 1), axis=0);
            zero_points = np.argwhere(n_sup_ave_for_zeros == 0.).T[0]
            x_ar = np.delete(x_ar, zero_points)
            rdf_per_num = []
            for numID, a_num in enumerate(zip(n_temp_ave, i_run_ave)):
                rdf_per_comp = []
                for compID, aComp in enumerate(zip(a_num[0], a_num[1])):
                    n_dum = np.delete(aComp[0], zero_points)
                    rdf_per_temp = []
                    for tempID, aTemp in enumerate(aComp[1]):
                        i_dum = np.delete(aTemp, zero_points)
                        p_rdf_val = i_dum / n_dum
                        rdf_per_temp.append(p_rdf_val)
                    rdf_per_comp.append(rdf_per_temp)
                rdf_per_num.append(rdf_per_comp)
            tot_rdf_arr.append([x_ar, rdf_per_num])
        rFile = sysName + '_CalcRDF.b'
        this_sys['RDF File'] = rFile
        with open(rFile, "wb+") as rdf_file:
            pickle.dump(tot_rdf_arr, rdf_file)

    def Gen_RhoBar_For(self, sysName, my_rdfFile=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            num_arr = self.SysInfo[sysName]['Molecule Numbers'];
            comp_list = this_sys['Comp List']
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
            this_sys['RDF File'] = sysName + '_CalcRDF.b'
        except KeyError:
            print("Is this the correct system?")
            return
        if my_rdfFile != None:
            this_sys['RDF File'] = my_rdfFile
        rdf_file = this_sys['RDF File']
        with open(rdf_file) as rFile:
            tot_rdfs = pickle.load(rFile)

        rho_bar_mat = np.zeros((len(box_ar), len(num_arr),  tot_temps, len(comp_list)))
        for boxID, a_box in enumerate(tot_rdfs):
            x_ar = a_box[0]
            pRDFs = a_box[1]
            box_size = box_ar[boxID]
            volume_norm = RDFVolumeElement(x_ar, box_size)
            for numID, aNum in enumerate(pRDFs):
                for compID, aComp in enumerate(aNum):
                    for tempID, aTemp in enumerate(aComp):
                        g_of_r = aTemp;
                        my_func = np.abs(g_of_r - 1.)
                        my_func_norm = my_func * volume_norm
                        rho_bar = np.trapz(my_func_norm, x_ar) / (box_size ** 3.)
                        # rho_bar = np.trapz(g_of_r*my_func_norm,x_ar)/np.trapz(g_of_r*volume_norm,x_ar)
                        rho_bar_mat[boxID][numID][tempID][compID] = rho_bar;
        rhoFile = sysName + '_rhobar.c'
        np.save(rhoFile, rho_bar_mat)
        rhoFile += '.npy'
        this_sys['RhoBar File'] = rhoFile

    def Gen_PhiC_For(self, sysName, my_clusFile=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            num_ar = self.SysInfo[sysName]['Molecule Numbers'];
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
            this_sys['Raw_CLUS'] = sysName + '_CLUS.b.npy'
        except KeyError:
            print("Is this the correct system?")
            return
        if my_clusFile != None:
            this_sys['Raw_CLUS'] = my_clusFile
        clus_file = this_sys['Raw_CLUS']
        tot_clus = np.load(clus_file)

        clus_run_avg = np.mean(tot_clus, axis=2)
        clus_run_err = np.std(tot_clus, axis=2)
        num_mols = float(num_ar[0,0]+num_ar[0,1])
        CMols    = np.linspace(1,num_mols,num_mols)
        phi_c_mat = np.zeros((len(box_ar), len(num_ar),  tot_temps, 5, 2))
        for boxID, aBox in enumerate(box_ar):
            for numID, aNum in enumerate(num_ar):
                for tempID in range(tot_temps):
                    clus_raw_m = clus_run_avg[boxID,numID,tempID]
                    clus_raw_e = clus_run_err[boxID,numID,tempID]
                    this_phi = phi_c_mat[boxID,numID,tempID]

                    this_phi[0,0] = clus_raw_m[0]/num_mols
                    this_phi[0,1] = clus_raw_e[0]/num_mols

                    CDist_Sum = np.sum(clus_raw_m[1:])
                    CDist_Norm = clus_raw_m[1:]/CDist_Sum
                    CErr_Norm = clus_raw_e[1:]/CDist_Sum
                    sigma_mean = np.dot(CMols,CDist_Norm)
                    sigma_m_err = np.sqrt(np.sum((CDist_Norm * CErr_Norm)**2.))
                    sigma_var = np.dot(CMols*CMols, CDist_Norm) - sigma_mean**2.
                    sigma_v_err1 = np.sqrt(np.sum((CMols*CMols*CErr_Norm)**2.))
                    sigma_v_err2 = 2.*sigma_m_err*sigma_mean
                    sigma_var_err = np.sqrt(sigma_v_err1**2. + sigma_v_err2**2.)
                    this_phi[1,0] = sigma_mean
                    this_phi[1,1] = sigma_m_err
                    this_phi[2,0] = sigma_var
                    this_phi[2,1] = sigma_var_err

                    CDist_Sum = np.sum(clus_raw_m[1:] * CMols)
                    CDist_Norm = clus_raw_m[1:] * CMols / CDist_Sum
                    CErr_Norm  = clus_raw_e[1:] * CMols / CDist_Sum
                    sigma_mean = np.dot(CMols, CDist_Norm)
                    sigma_m_err = np.sqrt(np.sum((CDist_Norm * CErr_Norm) ** 2.))
                    sigma_var = np.dot(CMols * CMols, CDist_Norm) - sigma_mean ** 2.
                    sigma_v_err1 = np.sqrt(np.sum((CMols * CMols * CErr_Norm) ** 2.))
                    sigma_v_err2 = 2. * sigma_m_err * sigma_mean
                    sigma_var_err = np.sqrt(sigma_v_err1 ** 2. + sigma_v_err2 ** 2.)
                    this_phi[3, 0] = sigma_mean
                    this_phi[3, 1] = sigma_m_err
                    this_phi[4, 0] = sigma_var
                    this_phi[4, 1] = sigma_var_err

        phiFile = sysName + '_phic.c'
        np.save(phiFile, phi_c_mat)
        phiFile += '.npy'
        this_sys['Perc File'] = phiFile

    def Gen_CorrectDen_For(self, sysName, my_WInt=None, my_NoInt=None):
        try:
            this_sys = self.SysInfo[sysName]
            box_ar = this_sys['Boxes']
            tot_mol_types = len(this_sys['Structure'])
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
        except KeyError:
            print("Is this the correct system?")
            return
        if my_WInt != None:
            this_sys['Raw_WInt_COMDen'] = my_WInt
        else:
            this_sys['Raw_WInt_COMDen'] = sysName + 'ORTH_W_COMDen.b'
        if my_NoInt != None:
            this_sys['Raw_NoInt_COMDen'] = my_NoInt
        else:
            this_sys['Raw_NoInt_COMDen'] = sysName + 'ORTH_N_COMDen.b'

        nFile = this_sys['Raw_NoInt_COMDen']
        iFile = this_sys['Raw_WInt_COMDen']

        with open(iFile) as dumFile:
            i_pdfs = pickle.load(dumFile)
        with open(nFile) as dumFile:
            n_pdfs = pickle.load(dumFile)

        tot_corr_den = []
        per_box_den  = []
        for boxID, a_box in enumerate(box_ar):
            x_ar = np.arange(0., a_box, 0.25)
            n_ar = np.array(n_pdfs[0][boxID])
            i_ar = np.array(i_pdfs[0][boxID])
            i_run_ave, i_run_err = MeanAndError(i_ar, 0);
            n_run_ave, n_run_err = MeanAndError(n_ar, 0);
            num_NAr = n_run_ave[:];
            num_WAr = i_run_ave[:];
            per_comp_den = []
            for compID, aComp in enumerate(zip(num_NAr, num_WAr)):
                comp_NAr = aComp[0]
                comp_IAr = aComp[1]
                n_temp_ave, n_temp_err = MeanAndError(comp_NAr, 0)
                zero_points = np.argwhere(n_temp_ave == 0).T[0]
                comp_xAr = np.delete(x_ar, zero_points)
                dum_nAr = np.delete(n_temp_ave, zero_points)
                per_temp_den = []
                for tempID, aTemp in enumerate(comp_IAr):
                    dum_iAr = np.delete(aTemp, zero_points)
                    per_temp_den.append(dum_iAr)
                per_comp_den.append([comp_xAr, dum_nAr, per_temp_den])
            per_box_den.append(per_comp_den)
        tot_corr_den.append(per_box_den)

        per_box_den = []
        for boxID, a_box in enumerate(box_ar):
            x_ar = np.arange(0., a_box, 0.25)
            n_ar = np.array(n_pdfs[1][boxID])
            i_ar = np.array(i_pdfs[1][boxID])
            i_run_ave, i_run_err = MeanAndError(i_ar, 0);
            n_run_ave, n_run_err = MeanAndError(n_ar, 0);
            num_NAr = n_run_ave[:];
            num_WAr = i_run_ave[:];
            per_comp_den = []
            for compID, aComp in enumerate(zip(num_NAr, num_WAr)):
                comp_NAr = aComp[0]
                comp_IAr = aComp[1]
                n_temp_ave, n_temp_err = MeanAndError(comp_NAr, 0)
                zero_points = np.argwhere(n_temp_ave == 0).T[0]
                comp_xAr = np.delete(x_ar, zero_points)
                dum_nAr = np.delete(n_temp_ave, zero_points)
                per_temp_den = []
                for tempID, aTemp in enumerate(comp_IAr):
                    dum_iAr = np.delete(aTemp, zero_points)
                    per_temp_den.append(dum_iAr)
                per_comp_den.append([comp_xAr, dum_nAr, per_temp_den])
            per_box_den.append(per_comp_den)
        tot_corr_den.append(per_box_den)



        rFile = sysName + 'ORTH_CorrectDen.b'
        this_sys['CorrectDen File'] = rFile
        with open(rFile, "wb+") as den_file:
            pickle.dump(tot_corr_den, den_file)

    def Gen_TieLineData_For(self, sysName, my_corrDenFile=None, idx_upto=13, idx_upfrom=-13, idx_uptill=-3):
        this_sys = self.SysInfo[sysName]
        if my_corrDenFile != None:
            this_sys['CorrectDen File'] = my_corrDenFile
        else:
            this_sys['CorrectDen File'] = sysName + 'ORTH_CorrectDen.b'
        dFile = this_sys['CorrectDen File'];
        with open(dFile) as dumFile:
            tot_den_data = pickle.load(dumFile)


        box_ar        = this_sys['Boxes']
        tot_struc     = this_sys['Structure']
        tot_mol_types = len(tot_struc)
        tot_poss_den  = tot_mol_types*(tot_mol_types+1)
        tot_temps     = this_sys['Key File'][1]['MC_CYCLE_NUM']
        temp_num_ar   = this_sys['Molecule Numbers'].T #This only has molecule numbers from the two sampled molecules
        #but ignores the rest, which we should manually add.
        tot_num_ar = []; tot_num_ar2 = []
        tot_num_ar.append(temp_num_ar[0]); tot_num_ar2.append(temp_num_ar[1]);
        tot_num_ar.append(temp_num_ar[1]); tot_num_ar2.append(temp_num_ar[0]);
        for additional_mol in tot_struc[2:]:
            this_num = additional_mol[0]
            this_num = this_num*np.ones(len(temp_num_ar[0]))
            tot_num_ar.append(this_num)
            tot_num_ar2.append(this_num);
        tot_num_ar  = np.array(tot_num_ar).T
        tot_num_ar2 = np.array(tot_num_ar2).T

        #The tie-line data are organized as follows:
        # Direction, Temperature, BoxSize And MoleculeNumber, Component, and then
        # [highCon, highCon_std, lowCon, lowCon_std, errSum/deltaCon]

        tie_line_data = np.zeros((2, tot_temps, tot_poss_den, 6, len(box_ar)))

        for bID, a_tot_set in enumerate(zip(box_ar, tot_num_ar)):
            aBox = a_tot_set[0]
            aNum = a_tot_set[1]
            inverse_vol = 1./(aBox**3.)
            for compID in range(tot_poss_den):
                this_comp_den = aNum[compID % tot_mol_types]*inverse_vol
                #this_comp_xAr = tot_den_data[0][bID][compID][0][1:]
                this_comp_nAr = tot_den_data[0][bID][compID][1][1:]
                for tempID in range(tot_temps):
                    this_comp_iAr = tot_den_data[0][bID][compID][2][tempID][1:]
                    this_comp_yAr = this_comp_den*this_comp_iAr/this_comp_nAr
                    highCon, highConEr = MeanAndError(this_comp_yAr[:idx_upto], 0)
                    lowCon, lowConEr   = MeanAndError(this_comp_yAr[idx_upfrom:idx_uptill], 0)
                    deltaCon           = np.abs(highCon-lowCon)
                    totErr             = np.sqrt(highConEr**2. + lowConEr**2.)

                    tie_line_data[0][tempID][compID][0][bID] = this_comp_den
                    tie_line_data[0][tempID][compID][1][bID] = highCon
                    tie_line_data[0][tempID][compID][2][bID] = highConEr
                    tie_line_data[0][tempID][compID][3][bID] = lowCon
                    tie_line_data[0][tempID][compID][4][bID] = lowConEr
                    tie_line_data[0][tempID][compID][5][bID] = totErr/deltaCon

        for bID, a_tot_set in enumerate(zip(box_ar, tot_num_ar2)):
            aBox = a_tot_set[0]
            aNum = a_tot_set[1]
            inverse_vol = 1./(aBox**3.)
            for compID in range(tot_poss_den):
                this_comp_den = aNum[compID % tot_mol_types]*inverse_vol
                #this_comp_xAr = tot_den_data[1][bID][compID][0][1:]
                this_comp_nAr = tot_den_data[1][bID][compID][1][1:]
                for tempID in range(tot_temps):
                    this_comp_iAr = tot_den_data[1][bID][compID][2][tempID][1:]
                    this_comp_yAr = this_comp_den*this_comp_iAr/this_comp_nAr
                    highCon, highConEr = MeanAndError(this_comp_yAr[:idx_upto], 0)
                    lowCon, lowConEr   = MeanAndError(this_comp_yAr[idx_upfrom:idx_uptill], 0)
                    deltaCon           = np.abs(highCon-lowCon)
                    totErr             = np.sqrt(highConEr**2. + lowConEr**2.)

                    tie_line_data[1][tempID][compID][0][bID] = this_comp_den
                    tie_line_data[1][tempID][compID][1][bID] = highCon
                    tie_line_data[1][tempID][compID][2][bID] = highConEr
                    tie_line_data[1][tempID][compID][3][bID] = lowCon
                    tie_line_data[1][tempID][compID][4][bID] = lowConEr
                    tie_line_data[1][tempID][compID][5][bID] = totErr/deltaCon
        tieFile = sysName + 'ORTH_TieLines.c'
        np.save(tieFile, tie_line_data)
        tieFile += '.npy'
        this_sys['TieLines File'] = tieFile

    def Gen_PlottingDataFor(self, sysName, comp1, comp2, my_TieLineFile=None):
        this_sys = self.SysInfo[sysName]
        if my_TieLineFile != None:
            this_sys['TieLines File'] = my_TieLineFile
        else:
            this_sys['TieLines File'] = sysName + 'ORTH_TieLines.c.npy'

        tie_file = this_sys['TieLines File']

        tot_tie_line_data = np.load(tie_file)

        tot_plot_data = []
        for a_dir in tot_tie_line_data:
            for tempID, aTemp in enumerate(a_dir):
                per_box_data = []
                for boxID in range(len(aTemp.T)):
                    cons_lin = ForPlot_Get_Ortho_TieLineConcs(a_dir, tempID, boxID, comp1, comp2)
                    cons_log = ForPlot_GenLogConcs(cons_lin[0], cons_lin[1])
                    cons_rat = ForPlot_GenRatioConcs(cons_lin[0], cons_lin[1])
                    per_box_data.append([cons_lin, cons_log, cons_rat])
            tot_plot_data.append(per_box_data)

        tot_plot_data = np.array(tot_plot_data)

        this_sys['Plotting Data'] = tot_plot_data