
import numpy as np
import scipy as sp
import os
import shutil
import subprocess as sproc
import time
import pickle
import matplotlib.pyplot as plt

__author__ = 'Furqan Dar'
__name__   = 'LASSI'

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
                                print(("The submission failed! Signal: ", -retcode))
                            else:
                                print(ret_code)
                        except OSError as myErr:
                            print(("Couldn't submit because ", myErr))
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


def sGen_Implicit_SymmBranched(Mol_Num, Branch_Num, Beads_PerBranch, Bead_Type, Hub_Type, Linker_Len):
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

            bond_list = [[i + 1 for i in range(Branch_Num)]]
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
    print(("Wrote structures to file: {:}".format(file_name)))


def Index_RDF(compA, compB, totComp):
    assert type(compA) == type(compB) == type(totComp) == int
    if compA > compB:
        return Index_RDF(compB, compA, totComp)
    elif compA == compB == -1:
        return 0
    else:
        return 1 + compA if compA == compB else totComp + compB - int((compA * (3 + compA - 2*totComp))/2)


def MeanAndError(theArray, theAxis):
    return np.mean(theArray, axis=theAxis), np.std(theArray, axis=theAxis)

"""
Functions to calculate the correct spherical volume element in a simple cubic lattice with periodic boundaries
"""
def RDFVol_F1(xVal, boxSize):
    return np.arctan(np.sqrt(-2. + (4.*(xVal**2.)/(boxSize**2.))))


def RDFVol_F2(xVal, boxSize):
    return (8.*xVal*np.arctan((2.*xVal*(-3. + (4.*(xVal**2.))/(boxSize**2.)))/
                              (boxSize*np.sqrt(-2. + (4.*(xVal**2.))/(boxSize**2.))*
                               (1. + (4.*(xVal**2.))/(boxSize**2.)))))/boxSize


def RDFVolumeElement(xAr, boxSize):
    '''
    Correctly calculates the spherical volume element in a cubic lattice with periodic boundary conditions
    '''
    new_x_Ar = np.ones(len(xAr));
    for xId, xVal in enumerate(xAr):
        if xVal <= boxSize/2.:
            new_x_Ar[xId] = 4.*np.pi*(xVal**2.);
        if xVal > boxSize/2. and xVal <= np.sqrt(2.)*boxSize/2.:
            new_x_Ar[xId] = 2.*np.pi*(3.*boxSize - 4.*xVal)*xVal
        if xVal > np.sqrt(2.)*boxSize/2. and xVal <= np.sqrt(3.)*boxSize/2.:
            new_x_Ar[xId] = 2.*boxSize*(-12.*RDFVol_F1(xVal, boxSize) + RDFVol_F2(xVal, boxSize) + 3.*np.pi)*xVal
    return new_x_Ar

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
        print("Initialized a LASSI setup with the following structures:")
        for aStr, aLin in zip(structures, linkers):
            self.SysInfo[aStr] = {}
            self.SysInfo[aStr]['Linker Length'] = aLin
            print((" " * 2 + "{:}:= Linker length {:} lattice units.".format(aStr, aLin)))
            self.SysInfo[aStr]['Structure'] = []
        self.AddParamFiles_ForAll(param_file)
        self.SetRunName_ForAll()
        self.SetNumberOfRuns_ForAll(2)

    def Set_SimulationPath(self, end_path_name):
        self.SimulationPath = end_path_name
        print(("Simulations shall be done in dir: {:}".format(self.SimulationPath)))

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
                    a_line = next(tot_file)
                    num_mol = int(a_line[:-1])
                    a_line = next(tot_file)
                    while (a_line[:4] != "}END"):
                        a_line = a_line[:-1].split()
                        this_line = [int(aVal) for aVal in a_line]
                        a_line = next(tot_file)
                        dum_struc_in.append(this_line)
                    dum_struc_tot.append([num_mol, np.array(dum_struc_in)])
        return dum_struc_tot

    def AddNewSystem(self, sysName, lin_len):
        if sysName in list(self.SysInfo.keys()):
            print("This structure name already exists! Doing nothing.")
        else:
            self.SysInfo[sysName] = {}
            self.SysInfo[sysName]['Linker Length'] = lin_len
            self.SysInfo[sysName]['Structure'] = []
            print(("Added {:}:= Linker length {:} lattice units.".format(sysName, lin_len)))
            self.AddParamFileTo(sysName, self.GlobalParamFile)
            self.SetRunNameFor(sysName, sysName)
            self.SetNumberOfRunsFor(sysName, 2)

    def AddStrucFileTo(self, sysName, file_name):
        try:
            self.SysInfo[sysName]['Structure File'] = file_name
            print(("{:} has structure file: {:}".format(sysName, file_name)))
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
            print(("{:} has energy files: Int:{:}  NoInt:{:}".format(sysName, self.SysInfo[sysName]['Int Energy File'],
                                                                    self.SysInfo[sysName]['NoInt Energy File'])))
            with open(self.CurrentDir + fileName_WInt) as WIntFile:
                WIntFile.readline() #First line is a comment
                tot_stickers = int(WIntFile.readline())
            self.SysInfo[sysName]['Tot Bead Types'] = tot_stickers
        except KeyError:
            print("Failed! Did you type the correct system name?")
            return

    def AddEnergyFiles_ForAll(self, fileName_NoInt, fileName_WInt):
        for aSys in list(self.SysInfo.keys()):
            self.AddEnergyFileTo(aSys, fileName_NoInt, fileName_WInt)

    def AddParamFileTo(self, sysName, fileName):
        self.SysInfo[sysName]['Key File'] = Read_ParamFile(fileName)

    def AddParamFiles_ForAll(self, fileName):
        for aSys in list(self.SysInfo.keys()):
            self.AddParamFileTo(aSys, fileName)

    def PrintParamsFor(self, sysName):
        try:
            dumKeys = self.SysInfo[sysName]['Key File']
            dumDict = dumKeys[1]
            dumKeys = dumKeys[0]
            for aKey in dumKeys:
                dum_spaces = 25 - len(aKey)
                print(('{:}{:}{:}'.format(aKey, ' ' * dum_spaces, dumDict[aKey])))
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

    def AddStruc_ExpLinear(self, sysName, mol_num, bead_num, bead_type, lin_len, lin_type):
        try:
            self.SysInfo[sysName]['Structure'].append(
                sGen_Explicit_Linear(mol_num, bead_num, bead_type, lin_len, lin_type))
        except KeyError:
            print("Failed! Did you type the correct system name? Or does this system even exist?")

    def AddStruc_ImpBranSymm(self, sysName, mol_num, branch_num, bead_per_branch, bead_type, hub_type, lin_len):
        try:
            self.SysInfo[sysName]['Structure'].append(
                sGen_Implicit_SymmBranched(mol_num, branch_num, bead_per_branch, bead_type, hub_type, lin_len))
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
        for aSys in list(self.SysInfo.keys()):
            sGen_WriteStructuresToFile(self.SysInfo[aSys]['Structure'], aSys + '_struc.prm')
            self.AddStrucFileTo(sysName, aSys + '_struc.prm')

    def SetTemperaturesFor(self, sysName, init_temp, final_temp, temp_steps, therm_temp):
        assert (type(temp_steps) == int) and (type(init_temp) == type(final_temp) == float)
        assert temp_steps > 0. and init_temp > 0. and final_temp > 0.
        delta_temp = (final_temp - init_temp) / (temp_steps-1)
        try:
            key_dum = self.SysInfo[sysName]['Key File'][1]
            key_dum['MC_TEMP'] = init_temp
            key_dum['MC_CYCLE_NUM'] = temp_steps
            key_dum['MC_DELTA_TEMP'] = delta_temp
            key_dum['PREEQ_TEMP'] = therm_temp
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetTemperatures_ForAll(self, init_temp, final_temp, temp_steps, therm_temp):
        for aSys in list(self.SysInfo.keys()):
            self.SetTemperaturesFor(aSys, init_temp, final_temp, temp_steps, therm_temp)

    def SetMCStepsFor(self, sysName, therm_steps, run_steps):
        try:
            key_dum = self.SysInfo[sysName]['Key File'][1]
            key_dum['N_STEPS'] = int(run_steps)
            key_dum['PREEQ_STEPS'] = int(therm_steps)
        except KeyError:
            print("Failed! Did you type the correct system name?")

    def SetMCSteps_ForAll(self, therm_steps, run_steps):
        for aSys in list(self.SysInfo.keys()):
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
        bead_types = np.array(bead_type_list)
        bead_types = bead_types.flatten()
        bead_types = np.unique(bead_types)
        bead_type_list = [[-1,-1]] #For the total
        for i in range(len(bead_types)):
            for j in range(i,len(bead_types)):
                bead_type_list.append([bead_types[i],bead_types[j]])
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
        for aSys in list(self.SysInfo.keys()):
            self.SetBoxSizesFor(aSys, low_con, high_con, tot_boxes)

    def SetBoxSizes_To(self, sysName):
        try:
            tot_ar = self.SysInfo[sysName]['Boxes']
        except KeyError:
            print("Failed! Did you type the correct system name? "
                  "Or have you setup the boxes for this system?")
            return
        for aSys in list(self.SysInfo.keys()):
            if aSys != sysName:
                self.SysInfo[aSys]['Boxes'] = tot_ar

    def SetRunNameFor(self, sysName, runName):
        try:
            self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX'] = runName
        except KeyError:
            print("Failed! Did you type the correct system name?"
                  " Does this system exist?")

    def SetRunName_ForAll(self):
        for aSys in list(self.SysInfo.keys()):
            self.SetRunNameFor(aSys, aSys)

    def SetNumberOfRunsFor(self, sysName, run_num):
        try:
            self.SysInfo[sysName]['Runs'] = run_num
        except KeyError:
            print("Failed! Did you type the correct system name?"
                  " Does this system exist?")

    def SetNumberOfRuns_ForAll(self, run_num):
        for aSys in list(self.SysInfo.keys()):
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
        for aSys in list(self.SysInfo.keys()):
            self.MakeDirs_For(aSys)

    def MakeDirs_ForAll_NoInt(self):
        for aSys in list(self.SysInfo.keys()):
            self.MakeDirs_NoIntFor(aSys)

    def MakeDirs_ForAll_WInt(self):
        for aSys in list(self.SysInfo.keys()):
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
                    param_copy[1]['RANDOM_SEED'] = a_run+1
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
                for a_run in range(self.SysInfo[sysName]['Runs']):
                    final_dir = this_dir + str(a_run + 1) + '/'
                    param_copy[1]['RANDOM_SEED'] = a_run
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
        for aSys in list(self.SysInfo.keys()):
            self.Write_ParamsWIntFor(aSys)

    def Write_ParamsNoInt_ForAll(self):
        for aSys in list(self.SysInfo.keys()):
            self.Write_ParamsNoIntFor(aSys)

    def Write_Params_ForAll(self):
        for aSys in list(self.SysInfo.keys()):
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
                            print(("The submission failed! Signal: ", -retcode))
                        else:
                            print(ret_code)
                    except OSError as myErr:
                        print(("Couldn't submit because ", myErr))
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
                            print(("The submission failed! Signal: ", -retcode))
                        else:
                            print(ret_code)
                    except OSError as myErr:
                        print(("Couldn't submit because ", myErr))
                        raise
                    time.sleep(2)
                    self.QSubIter += 1
        except KeyError:
            print("Have you setup the boxes yet?")
        os.chdir(self.CurrentDir)

    def SubmitWIntJobs_ForAll(self, WIntQueue):
        for aSys in list(self.SysInfo.keys()):
            self.SubmitWIntJobs_ToQueueFor(aSys, WIntQueue)

    def SubmitNoIntJobs_ForAll(self, NoIntQueue):
        for aSys in list(self.SysInfo.keys()):
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
        for aSys in list(self.SysInfo.keys()):
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
            self.SimulationPath  = SimSetupInstance.SimulationPath
            self.CurrentDir      = SimSetupInstance.CurrentDir
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            return

    def Read_PDFFileFor(self, sysName, file_name, box_size, tot_temps):
        try:
            dum_dat   = np.loadtxt(file_name)
            comp_list = self.SysInfo[sysName]['Comp List']
            tot_comps = self.SysInfo[sysName]['Tot Bead Types']
            tot_comps_pos = 1+Index_RDF(tot_comps-2, tot_comps-1, tot_comps)
        except IOError:
            print("Does the RDF file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        actual_comps = len(comp_list)
        ret_dat = np.zeros((actual_comps, tot_temps, box_size*4))
        for a_temp in range(tot_temps):
            for compID, compPair in enumerate(comp_list):
                compA = int(compPair[0]); compB = int(compPair[1])
                ret_dat[compID][a_temp] = dum_dat[tot_comps_pos*a_temp+Index_RDF(compA,compB,tot_comps)]
        return ret_dat

    def Read_CLUSFileFor(self, file_name):
        try:
            dum_dat   = np.loadtxt(file_name)
        except IOError:
            print("Does the RDF file exist? Did these runs finish and output data?")
            raise
        except KeyError:
            print("Did you type the system correctly?")
            raise
        return dum_dat

    def Read_GYRADFileFor(self, file_name):
        try:
            dum_dat   = np.loadtxt(file_name)
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
            file_name = sysName+'_N_PDF.b'
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
            np.save(file_name,dum_dat)
            self.SysInfo[sysName]['Raw_CLUS'] = file_name+'.npy'
            return
        except KeyError:
            print("Did you type the correct system name?")
            return

    def Set_PDFFileNames_For(self, sysName):
        self.SysInfo[sysName]['Raw_NoInt_PDF'] = sysName+'_N_PDF.b'
        self.SysInfo[sysName]['Raw_WInt_PDF']  = sysName+'_W_PDF.b'
    def Collect_NoInt_PDFs_For(self, sysName):
        try:
            box_arr   = self.SysInfo[sysName]['Boxes']
            key_dum   = self.SysInfo[sysName]['Key File'][1]
            over_dir  = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name  = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num   = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir+str(a_box)+'/NoInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir+str(a_run+1)+'/'+run_name+'_RDF.dat'
                tot_raw_dat_temp.append(self.Read_PDFFileFor(sysName, file_name, a_box, tot_temps))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat
    def Collect_WInt_PDFs_For(self, sysName):
        try:
            box_arr   = self.SysInfo[sysName]['Boxes']
            key_dum   = self.SysInfo[sysName]['Key File'][1]
            over_dir  = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name  = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num   = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir+str(a_box)+'/WInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir+str(a_run+1)+'/'+run_name+'_RDF.dat'
                tot_raw_dat_temp.append(self.Read_PDFFileFor(sysName, file_name, a_box, tot_temps))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat

    def Collect_WInt_CLUS_For(self, sysName):
        try:
            box_arr   = self.SysInfo[sysName]['Boxes']
            key_dum   = self.SysInfo[sysName]['Key File'][1]
            over_dir  = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name  = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num   = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir+str(a_box)+'/WInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir+str(a_run+1)+'/'+run_name+'_CLUS.dat'
                tot_raw_dat_temp.append(self.Read_CLUSFileFor(file_name))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat
    def Collect_GYRAD_For(self, sysName):
        try:
            box_arr   = self.SysInfo[sysName]['Boxes']
            key_dum   = self.SysInfo[sysName]['Key File'][1]
            over_dir  = self.SimulationPath + sysName + '/' + self.SysInfo[sysName]['Linker Length'] + '/'
            tot_temps = key_dum['MC_CYCLE_NUM']
            run_name  = self.SysInfo[sysName]['Key File'][1]['REPORT_PREFIX']
            run_num   = self.SysInfo[sysName]['Runs']
        except KeyError:
            print("Did you type the system name correctly?")
            raise
        except AttributeError:
            print("Analysis expects a fully functioning instance of Sim_Setup!")
            raise
        tot_raw_dat = []
        for a_box in box_arr:
            box_dir = over_dir+str(a_box)+'/WInt/'
            tot_raw_dat_temp = []
            for a_run in range(run_num):
                file_name = box_dir+str(a_run+1)+'/'+run_name+'_GR.dat'
                tot_raw_dat_temp.append(self.Read_GYRADFileFor(file_name))
            tot_raw_dat.append(tot_raw_dat_temp)
        return tot_raw_dat

    def Gen_pRDFs_For(self, sysName, my_WInt=None, my_NoInt=None):
        try:
            this_sys  = self.SysInfo[sysName]
            box_ar    = this_sys['Boxes']
            comp_list = this_sys['Comp List']
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
            self.Set_PDFFileNames_For(sysName)
        except KeyError:
            print("Is this the correct system?")
            return
        if my_WInt != None:
            this_sys['Raw_WInt_PDF'] = my_WInt
        if my_NoInt != None:
            this_sys['Raw_NoInt_PDF'] = my_WInt
        nFile = this_sys['Raw_NoInt_PDF']
        iFile = this_sys['Raw_WInt_PDF']
        with open(iFile) as dumFile:
            i_pdfs = pickle.load(dumFile)
        with open(nFile) as dumFile:
            n_pdfs = pickle.load(dumFile)

        tot_rdf_arr = []
        for boxID, a_box in enumerate(box_ar):
            x_ar = np.arange(0.,a_box,0.25)
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
                    p_rdf_ar = i_dum/n_dum
                    rdf_temp_ar.append(p_rdf_ar)
                rdf_box_ar.append(rdf_temp_ar)
            tot_rdf_arr.append([x_ar, rdf_box_ar])
        rFile = sysName+'_CalcRDF.b'
        this_sys['RDF File'] = rFile
        with open(rFile, "wb+") as rdf_file:
            pickle.dump(tot_rdf_arr, rdf_file)
    def Gen_RhoBar_For(self, sysName, my_rdfFile = None):
        try:
            this_sys  = self.SysInfo[sysName]
            box_ar    = this_sys['Boxes']
            comp_list = this_sys['Comp List']
            tot_temps = this_sys['Key File'][1]['MC_CYCLE_NUM']
            this_sys['RDF File'] = sysName+'_CalcRDF.b'
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
                    rho_bar_temp.append(rho_bar)
                rho_bar_comp.append(rho_bar_temp)
            rho_bar_mat.append(rho_bar_comp)
        rho_bar_mat = np.array(rho_bar_mat)
        rhoFile = sysName + '_rhobar.c'
        np.save(rhoFile, rho_bar_mat)
        rhoFile +='.npy'
        this_sys['RhoBar File'] = rhoFile
    def Gen_PhiC_For(self, sysName, my_clusFile = None):
        try:
            this_sys  = self.SysInfo[sysName]
            box_ar    = this_sys['Boxes']
            mol_num   = this_sys['Tot Molecules']
            this_sys['Raw_CLUS'] = sysName+'_CLUS.b.npy'
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
        phiFile +='.npy'
        this_sys['Perc File'] = phiFile

    def Gen_RhoBarInterp_For(self, sysName, of_comp = 0, temp_scale = 2., inter_points = 50):
        try:
            this_sys   = self.SysInfo[sysName]
            key_dum = this_sys['Key File'][1]
        except KeyError:
            print("Is this the correct system?")

        init_temp = key_dum['MC_TEMP']
        temp_steps = key_dum['MC_CYCLE_NUM']
        delta_temp = key_dum['MC_DELTA_TEMP']
        mol_num    = this_sys['Tot Molecules']
        box_ar = [[a_box for a_box in this_sys['Boxes']]] * temp_steps
        temp_ar    =[[init_temp + i * delta_temp for i in range(temp_steps)]] * len(box_ar[0])
        box_ar     = np.array(box_ar)
        temp_ar    = np.array(temp_ar)/temp_scale
        conc_ar    = np.log10(mol_num/(box_ar**3.))
        conc_ar    = conc_ar.flatten()
        temp_ar    = temp_ar.T.flatten()
        conc_grid = np.linspace(conc_ar.min(), conc_ar.max(), inter_points)
        temp_grid = np.linspace(temp_ar.min(), temp_ar.max(), inter_points)
        conc_grid, temp_grid = np.meshgrid(conc_grid, temp_grid)

        rho_barFile = this_sys['RhoBar File']
        rho_bar  = np.load(rho_barFile)
        this_rho = rho_bar[:,of_comp,:].T
        this_rho = this_rho.flatten()

        rho_func = sp.interpolate.Rbf(conc_ar, temp_ar, this_rho, function='linear')
        rho_grid = rho_func(conc_grid, temp_grid)
        return [conc_grid, temp_grid, rho_grid]
    def Gen_PhiCInterp_For(self, sysName,  temp_scale = 2., inter_points = 50):
        try:
            this_sys   = self.SysInfo[sysName]
            key_dum = this_sys['Key File'][1]
        except KeyError:
            print("Is this the correct system?")

        init_temp = key_dum['MC_TEMP']
        temp_steps = key_dum['MC_CYCLE_NUM']
        delta_temp = key_dum['MC_DELTA_TEMP']
        mol_num    = this_sys['Tot Molecules']
        box_ar = [[a_box for a_box in this_sys['Boxes']]] * temp_steps
        temp_ar    =[[init_temp + i * delta_temp for i in range(temp_steps)]] * len(box_ar[0])
        box_ar     = np.array(box_ar)
        temp_ar    = np.array(temp_ar)/temp_scale
        conc_ar    = np.log10(mol_num/(box_ar**3.))
        conc_ar    = conc_ar.flatten()
        temp_ar    = temp_ar.T.flatten()
        conc_grid = np.linspace(conc_ar.min(), conc_ar.max(), inter_points)
        temp_grid = np.linspace(temp_ar.min(), temp_ar.max(), inter_points)
        conc_grid, temp_grid = np.meshgrid(conc_grid, temp_grid)

        phi_cFile = this_sys['Perc File']
        phi_c  = np.load(phi_cFile)
        phi_c  = phi_c[:,:,0].T
        phi_c  = phi_c.flatten()/mol_num

        phi_func = sp.interpolate.Rbf(conc_ar, temp_ar, phi_c, function='linear')
        phi_grid = phi_func(conc_grid, temp_grid)
        return [conc_grid, temp_grid, phi_grid]
    def Gen_OrderParamInterp_For(self, sysName, of_comp = 0, temp_scale = 2., inter_points = 50):
        try:
            this_sys   = self.SysInfo[sysName]
            key_dum = this_sys['Key File'][1]
        except KeyError:
            print("Is this the correct system?")

        init_temp = key_dum['MC_TEMP']
        temp_steps = key_dum['MC_CYCLE_NUM']
        delta_temp = key_dum['MC_DELTA_TEMP']
        mol_num    = this_sys['Tot Molecules']
        box_ar = [[a_box for a_box in this_sys['Boxes']]] * temp_steps
        temp_ar    =[[init_temp + i * delta_temp for i in range(temp_steps)]] * len(box_ar[0])
        box_ar     = np.array(box_ar)
        temp_ar    = np.array(temp_ar)/temp_scale
        conc_ar    = np.log10(mol_num/(box_ar**3.))
        conc_ar    = conc_ar.flatten()
        temp_ar    = temp_ar.T.flatten()
        conc_grid  = np.linspace(conc_ar.min(), conc_ar.max(), inter_points)
        temp_grid  = np.linspace(temp_ar.min(), temp_ar.max(), inter_points)
        conc_grid, temp_grid = np.meshgrid(conc_grid, temp_grid)

        rho_barFile = this_sys['RhoBar File']
        rho_bar  = np.load(rho_barFile)
        this_rho = rho_bar[:,of_comp,:].T
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

    def Store_OrderParamInterp_For(self, sysName, of_comp = 0, temp_scale = 2., inter_points = 50):
        try:
            this_sys = self.SysInfo[sysName]
        except KeyError:
            print("Is this the correct system?")
        this_sys['Order Params'] = self.Gen_OrderParamInterp_For(sysName, of_comp, temp_scale, inter_points)
