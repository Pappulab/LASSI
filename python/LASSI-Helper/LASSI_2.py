"""
New LaSSI Helper modules. Contains the following modules:
    - SystemSetup
        Used for instantiating, defining, setting up, submitting jobs, and data collection.
    - Plotting
        Used to plot the various data that are generated.
        **Nothing here yet**
"""
__author__ = 'Furqan Dar'
__version__ = 0.1
import numpy as np
#import scipy as sp
import os
#import shutil
#import subprocess as sproc
#import time
import pickle
import errno
#import matplotlib.pyplot as plt
import copy
#import TrjProc


class _TopoUtils_Linear_Molecule_Back(object):
    """
    Helper function to generate linear molecules of arbitrary architectures
    """
    def __init__(self):
        self.Blocks = []
        self.Linkers = []
        self.InterBlockLinkers = []
        self.MolStructure = []

    def add_block(self, bead_list, linker_list):
        assert len(bead_list) == (len(linker_list) + 1), \
            "number of bead and linkers does not generate a correct molecule!"
        self.Blocks.append(bead_list)
        self.Linkers.append(linker_list)

    def add_inter_block_linker(self, linker_len):
        self.InterBlockLinkers.append(linker_len)

    def repeat_block(self, bead_list, linker_list, inter_linker, repeat_num=5):
        assert len(bead_list) == (len(linker_list) + 1), \
            "number of bead and linkers does not generate a correct molecule!"
        for a_block in range(repeat_num):
            self.Blocks.append(bead_list)
            self.Linkers.append(linker_list)
            if a_block < repeat_num - 1:
                self.add_inter_block_linker(inter_linker)

    def form_structure(self):
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
        tot_bead_id = 0
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

    def get_structure(self):
        self.form_structure()
        return self.MolStructure


class TopoUtils_Gen_Linear_Molecule(object):
    """
    General helper class to generate certain archetypal linear molecules
    """
    def __init__(self, num_of_mols=1000):
        """
        Initialize the molecule
        """
        self.MolNumber = num_of_mols
        self.MolStruc  = []

    @property
    def Structure(self):
        assert len(self.MolStruc) > 0, "The structure has not been formed yet."
        return [self.MolNumber, self.MolStruc]

    def set_mol_struc(self, new_struc):
        assert isinstance(new_struc, _TopoUtils_Linear_Molecule_Back), "The structure should be _TopoUtils_Linear_Molecule_Back Type object"
        self.MolStruc  = new_struc.get_structure()

    def gen_implicit_linear(self, bead_type=0, bead_num=7, lin_len=2):
        """
        Given the bead_type and bead_num, we generate a simple linear molecule with implicit linkers
        :param bead_type:
        :param bead_num:
        :param lin_len:
        :return: None
        """
        DumStruc = _TopoUtils_Linear_Molecule_Back()
        DumStruc.add_block([bead_type] * bead_num, [lin_len] * (bead_num - 1))
        self.MolStruc = DumStruc.get_structure()

    def gen_monomer(self, bead_type=0):
        """
        Given the bead_type and bead_num, we generate a simple linear molecule with implicit linkers
        :param bead_type:
        :return: None
        """
        DumStruc = _TopoUtils_Linear_Molecule_Back()
        DumStruc.add_block([bead_type] * 1, [])
        self.MolStruc = DumStruc.get_structure()

    def gen_dimer(self, bead_1=0, bead_2=1, lin_len=2):
        """
        Given the bead_type and bead_num, we generate a simple linear molecule with implicit linkers
        :param bead_1:
        :param bead_2:
        :param lin_len:
        :return: None
        """
        DumStruc = _TopoUtils_Linear_Molecule_Back()
        DumStruc.add_block([bead_1, bead_2], [lin_len] * 1)
        self.MolStruc = DumStruc.get_structure()


class TopoUtils_Gen_Linear_System(object):
    """
    Helpers to generate multiple molecules which can then be used as structures for the SystemSetup
    class below.
    """
    def __init__(self):
        self.MolStrucs = []
        self.MolNums   = 0
        self.StickNum  = 0

    def add_monomer(self, bead_type=0, mol_num=1000):
        self.MolNums += 1
        dum_struc = TopoUtils_Gen_Linear_Molecule(num_of_mols=mol_num)
        dum_struc.gen_monomer(bead_type)
        self.MolStrucs.append(dum_struc.Structure)

    def add_dimer(self, bead_1=0, bead_2=1, lin_len=2, mol_num=1000):
        self.MolNums += 1
        dum_struc = TopoUtils_Gen_Linear_Molecule(num_of_mols=mol_num)
        dum_struc.gen_dimer(bead_1, bead_2, lin_len=lin_len)
        self.MolStrucs.append(dum_struc.Structure)

    def add_implicit_linear(self, bead_type=0, bead_num=7, linker_len=2, mol_num=1000):
        self.MolNums += 1
        dum_struc = TopoUtils_Gen_Linear_Molecule(num_of_mols=mol_num)
        dum_struc.gen_implicit_linear(bead_type, bead_num, linker_len)
        self.MolStrucs.append(dum_struc.Structure)

    @property
    def SysStrucs(self):
        return self.MolStrucs

    @property
    def StickerNumber(self):
        dum_list = []
        for a_struc in self.SysStrucs:
            sticker_types = a_struc[1].T[1]
            dum_list.append(sticker_types)
        dum_list = np.array([a_type for a_set in dum_list for a_type in a_set])
        dum_counts, dum_vals = np.unique(dum_list, return_counts=True)
        self.StickNum = len(dum_counts)
        return self.StickNum

    def add_molecule(self, new_molecule):
        assert isinstance(new_molecule, TopoUtils_Gen_Linear_Molecule), "new_molecule must be TopoUtils_Gen_Linear_Molecule type"
        self.MolNums += 1
        self.MolStrucs.append(new_molecule.Structure)

    def __repr__(self):
        """
        Crude way to show the system. Gives each unique bead-type a letter. Linker lengths are represented as '-'
        :return: {str} that contains system representation
        """
        tot_stickers = self.StickerNumber
        bead_to_letter = {a_number:a_letter for a_letter, a_number in zip('ABCDEFGHIJKLMNOPQRSTUVWXYZ', range(26))}

        list_to_print = []
        for strucID, a_struc in enumerate(self.SysStrucs):
            dum_num = str(a_struc[0])
            dum_str = self._gen_mol_in_letters(a_struc[1], bead_to_letter)
            list_to_print.append(r''+dum_num+'\t'+dum_str)

        return "\n".join(list_to_print)

    def _gen_mol_in_letters(self, this_struc, bead_to_letter):
        dum_list = []
        bead_types = this_struc.T[1]
        linker_len = this_struc.T[2]
        start_bead, end_bead = bead_types[0], bead_types[-1]
        dum_list.append(bead_to_letter[start_bead])
        if len(this_struc) == 1:
            return "".join(dum_list)
        elif len(this_struc) == 2:
            dum_list.append('-'*linker_len[0])
            dum_list.append(bead_to_letter[end_bead])
            return "".join(dum_list)
        else:
            for beadID, (a_bead, a_lin) in enumerate(zip(bead_types[1:-1:2], linker_len[1:-1:2])):
                dum_list.append('-' * a_lin)
                dum_list.append(bead_to_letter[a_bead])
            dum_list.append('-' * linker_len[-1])
            dum_list.append(bead_to_letter[end_bead])
            return "".join(dum_list)


class EnergyUtils_Gen_File(object):
    def __init__(self, tot_stick_num):
        self.St_Num   = tot_stick_num
        self.Ovlp_En  = np.zeros((tot_stick_num, tot_stick_num))
        self.Cont_En  = np.zeros((tot_stick_num, tot_stick_num))
        self.Sti_En   = np.zeros((tot_stick_num, tot_stick_num))
        self.FSol_En  = np.zeros((tot_stick_num, tot_stick_num))
        self.TInd_En  = np.zeros((tot_stick_num, tot_stick_num))
        self.Tot_File = []
        self.FileQ    = False

    def Add_Ovlp_Btwn(self, energy=-0.5, a_pair=(0, 0)):
        comp1, comp2 = a_pair
        self.Ovlp_En[comp1, comp2] = energy
        self.Ovlp_En[comp2, comp1] = energy
        return None

    def Add_Cont_Btwn(self, energy=-0.5, a_pair=(0, 0)):
        comp1, comp2 = a_pair
        self.Cont_En[comp1, comp2] = energy
        self.Cont_En[comp2, comp1] = energy
        return None

    def Add_Sti_Btwn(self, energy=-0.5, a_pair=(0, 0)):
        comp1, comp2 = a_pair
        self.Sti_En[comp1, comp2] = energy
        self.Sti_En[comp2, comp1] = energy
        return None

    def Add_FSol_For(self, energy=-0.5, this_type=0):
        self.FSol_En[this_type, this_type] = energy
        return None

    def Add_TInd_For(self, energy=-0.5, this_type=0):
        self.TInd_En[this_type, this_type] = energy
        return None

    def Add_En_Btwn(self, energy=-0.5, a_pair=(0, 0), en_mode="Sti"):
        if en_mode == "Ovlp":
            self.Add_Ovlp_Btwn(energy, a_pair)
        if en_mode == "Cont":
            self.Add_Cont_Btwn(energy, a_pair)
        if en_mode == "Sti":
            self.Add_Sti_Btwn(energy, a_pair)
        else:
            raise NameError(en_mode + " is not a valid option. ['Ovlp', 'Cont', 'Sti'] are valid.")

    def Write_Matrix(self, this_mat):
        unique_ints = np.unique(this_mat)
        num_of_ints = len(unique_ints)
        write_mat = [];
        if num_of_ints == 1:
            write_mat.append(f"{this_mat[0, 0]:.3f}")
        else:
            for aRow in this_mat:
                this_line = []
                for aNum in aRow:
                    this_line.append(f"{aNum:.3f}")
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

        tot_file.append("#FSOL_POT")
        for aline in self.Write_Matrix(self.FSol_En):
            tot_file.append(aline)
        tot_file.append("")

        tot_file.append("#T_IND_POT")
        for aline in self.Write_Matrix(self.TInd_En):
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
        self.FileQ    = True

    def get_EnFile(self):
        """
        Getter for the energy method
        :return: Tot_File
        """
        assert self.FileQ, "File needs to have been formed using Form_Energy_File() method"
        return self.Tot_File

    def Print_File(self):
        """
        Prints the energy-file to the screen.
        :return:
        """
        assert self.FileQ, "File needs to have been formed using Form_Energy_File() method"
        for a_line in self.Tot_File:
            print(a_line)


class IO_Utils(object):
    """
    A collection of IO Utilities.
    """

    @staticmethod
    def MKDirCatch(this_dir):
        try:
            os.mkdir(this_dir)
        except OSError as myErr:
            if myErr.errno != errno.EEXIST:
                raise
            pass

    @staticmethod
    def gen_dir_str(boxSize=100, molSet=None, repNum=1, wInt=True):
        """
        Given a boxSize, a set of molecule-numbers, the replicate number, and the interaction state,
        we generate the str that corresponds to this particular simulation instance.
        """

        if molSet is None:
            molSet=[1000, ]

        dum_dir = []
        # dum_dir.append('')
        dum_dir.append(str(boxSize))
        molStr = "_".join([str(a_mol) for a_mol in molSet])  # Create underscored string from mol-nums
        dum_dir.append(molStr)
        if wInt:
            dum_dir.append("WInt")
        else:
            dum_dir.append("NoInt")

        dum_dir.append(str(repNum))
        if os.name == 'nt':
            return "\\".join(dum_dir) + "\\"
        else:
            return "/".join(dum_dir) + "/"

    def gen_dir_for_box(boxSize=100,
                        molList=None,
                        repNums=1,
                        wInt=True):
        """
        For a given box-size, and a list of mol_nums, we generate the strings corresponding to all the mol_nums
        :param molList:
        :param repNums:
        :param wInt:
        :return:
        """

        if molList is None:
            molList=[[100,100,], [200,200]]

        dum_dir_list = []
        for setID, a_set in enumerate(molList):
            for a_rep in range(1, repNums+1):
                dum_dir_list.append(IO_Utils.gen_dir_str(boxSize, a_set, a_rep, wInt))
        return dum_dir_list

    @staticmethod
    def gen_dir_for_box_nested(boxSize=100,
                        molList=None,
                        repNums=1,
                        wInt=True,
                        molCycles=1):
        """
        For a given box-size, and a list of mol_nums, we generate the strings corresponding to all the mol_nums
        :param molList:
        :param repNums:
        :param wInt:
        :return:
        """

        if molList is None:
            molList=[[100, 100,], [200, 200]]

        per_dir_list = []
        for cycID, a_cycle in enumerate(range(molCycles)):
            per_mol_list = []
            for setID, a_set in enumerate(molList):
                this_set = np.roll(a_set, a_cycle)
                per_rep_list = []
                for a_rep in range(1, repNums+1):
                    per_rep_list.append(IO_Utils.gen_dir_str(boxSize, this_set, a_rep, wInt))
                per_mol_list.append(per_rep_list)
            per_dir_list.append(per_mol_list)
        return per_dir_list
    @staticmethod
    def gen_dir_list(boxList=None,
                     molList=None,
                     repNums=3,
                     wInt=True):
        """
        Given the list of box-sizes (boxList), and the list of molecule numbers (molList), we generate a list of all the
        possible directories for the given interaction state (wInt) and replicate numbers (repNums))
        We check to make sure that the lengths of boxList and molList are the same
        """

        if boxList is None:
            boxList=[100, 110, ]

        if molList is None:
            molList=[[1000, 1000], [900, 900], ]

        assert len(boxList) == len(molList), "Number of boxes and number of molecule-tuples should be the same!"

        dum_dir_list = []
        for a_box, molSet in zip(boxList, molList):
            dum_dir_list.append(IO_Utils.gen_dir_for_box(a_box, molSet, repNums, wInt))

        return [a_dir for a_set in dum_dir_list for a_dir in a_set]

    @staticmethod
    def gen_dir_nested_list(boxList=None,
                            molList=None,
                            repNums=3,
                            wInt=True,
                            dir_prefix='/'):
        """
        Generate a nested list given a list of box-sizes, a list_of_mol_nums, replicate numbers, interaction state.
        molCycle is used to use np.roll(a_mol) to generate the permutations for each molecule number set.
        The idea being that it can be used to generate the different directions in the case of orthogonal sampling
        The outermost index is the direction, then box, then molecule number, then replicate.
        Again, remember that the boxList and molList are assumed to be from SystemSetup.get_independent_conditions()
        :param boxList List of box-sizes
        :param molList: List of molecule numbers
        :param repNums: Number of replicates
        :param wInt: Interaction state
        :return: nested_list_of_dirs = [a_dir, a_box, a_mol, a_rep]
        """

        if boxList is None:
            boxList=[100, 110, ]

        if molList is None:
            molList=[[1000, 1000], [900, 900], ]

        assert len(boxList) == len(molList), "Number of boxes and number of molecule-tuples should be the same!"
        per_box_list = []
        for totID, (a_box, a_mol_set) in enumerate(zip(boxList, molList)):
            per_mol_list = []
            for molID, a_mol in enumerate(a_mol_set):
                per_rep_list = []
                for repID, a_rep in enumerate(range(1, repNums+1)):
                    per_rep_list.append(dir_prefix + IO_Utils.gen_dir_str(boxSize=a_box, molSet=a_mol, repNum=a_rep,
                                                                          wInt=wInt))
                per_mol_list.append(per_rep_list)
            per_box_list.append(per_mol_list)

        return per_box_list

    @staticmethod
    def gen_run_conditions(boxList=None,
                           molList=None,
                           repNums=3):
        """
        Creates flattened arrays of boxes and molList which are 1-1 with the gen_dir_list list_of_directories.
        This is to be used in writing key-files and such
        :param molList:
        :param repNums:
        :param wInt:
        :return:
        """

        if boxList is None:
            boxList=[100, 110, ]

        if molList is None:
            molList=[[1000, 1000], [900, 900], ]

        dum_boxes = []
        dum_mols  = []

        for boxID, (a_box, a_mol) in enumerate(zip(boxList, molList)):
            for setID, a_set in enumerate(a_mol):
                for repID, a_rep in enumerate(range(repNums)):
                    dum_boxes.append(a_box)
                    dum_mols.append(a_set)

        return np.array(dum_boxes), np.array(dum_mols)
    @staticmethod
    def read_param_file(file_name):
        """
        Read a LASSI key-file.
        For each line in the key-file, we store the key-word and the value as a key-value pair in a large dictionary.
        Ignore all lines that are empty, or start with '#'
        :return: [list of keys, the full_dictionary]
        """
        dum_dict = {}
        dum_list = []
        with open(file_name) as pFile:
            totFile = pFile.readlines()
            for a_line in totFile:
                if a_line[0] == '#' or a_line[0] == '\n':
                    continue
                s_line = a_line.strip("\n").split(" ") #Skipping the \n character
                s_line = [a_key for a_key in s_line if a_key != '']
                s_line = s_line[:2]
                dum_list.append(s_line[0])
                dum_dict[s_line[0]] = s_line[1]
        return [dum_list, dum_dict]

    @staticmethod
    def write_param_file(param_obj,
                         file_path,
                         box_size=100,
                         run_name='DumR',
                         energy_file='DumE',
                         struc_file='DumS',
                         rng_seed=0):
        """
        Given a template key-file, we write the key-file for this particular simulation.

        :param file_path: Absolute path to where the key-file should be written.
        :param box_size: Size of simulation box
        :param run_name: Name of the run
        :param energy_file: Absolute path to the energy-file
        :param struc_file: Absolute path to the structure-file
        :param rng_seed: Seed for the RNG in C
        :return: None
        """
        dum_keys = param_obj[0][:]
        dum_vals = copy.deepcopy(param_obj[1])

        dum_vals['BOX_SIZE']       = str(box_size)
        dum_vals['REPORT_PREFIX']  = run_name
        dum_vals['STRUCT_FILE']    = struc_file
        dum_vals['ENERGY_FILE']    = energy_file
        dum_vals['RANDOM_SEED']    = str(rng_seed)
        with open(file_path, "w+") as pFile:
            for a_key in dum_keys:
                N_spcs = 25 - len(a_key)
                pFile.write(a_key + ' ' * N_spcs + dum_vals[a_key] + '\n')
        return None

    @staticmethod
    def write_struc_file(struc_obj,
                         file_path):
        """
        Write the structure object in the format readable by LASSI to file_path
        :param file_path:
        :return:
        """
        with open(file_path, "w+") as strucFile:
            for a_Mol in struc_obj:
                mol_num = a_Mol[0]
                mol_struc = a_Mol[1]
                beads_per_mol = len(np.unique(mol_struc.T[0]))
                strucFile.write("#New Molecule Type:= {:} beads per molecule\n".format(beads_per_mol))
                strucFile.write("NEW{\n")
                strucFile.write(str(mol_num) + "\n")
                for a_line in mol_struc:
                    this_line = [str(a_num) + "\t" for a_num in a_line]
                    this_line.append("\n")
                    strucFile.write("".join(this_line))
                strucFile.write("}END\n")
        return None

    @staticmethod
    def write_energy_file(energy_file_obj,
                          file_path):
        """
        Writes the energy file to the specified path
        :param energy_file_list: The output of
        :param file_path:
        :return: None
        """
        assert isinstance(energy_file_obj, EnergyUtils_Gen_File), "This method requires EnergyUtils_Gen_File objects."
        dum_file = energy_file_obj.get_EnFile()
        with open(file_path, "w") as myFile:
            for aline in dum_file:
                myFile.write(aline + "\n")
        return None

    @staticmethod
    def create_dirs_from_list(list_of_dirs):
        """
        Iterate over the given list to generate all the directories.
        :return:
        """
        for dirID, a_dir in enumerate(list_of_dirs):
            try:
                os.makedirs(a_dir)
            except OSError:
                pass

        return None

    @staticmethod
    def loop_over_dir_list(list_of_dirs, _passed_func, **kwargs):
        """
        Given the list of directories, we loop over each directory, change-directory into that directory,
        perform the function _passed_func(**kwargs) where **kwargs are assumed to be for _passed_func, store the
        return-values of that function, and continue.
        At the end, we change directory back to where we started evaluating this function.
        :param list_of_dirs
        :param _passed_func:
        :param kwargs:
        :return: list_of_return_vals
        """
        dum_start_dir = os.getcwd()
        ret_vals = []
        for dirID, a_dir in enumerate(list_of_dirs):
            os.chdir(a_dir)
            ret_vals.append(_passed_func(**kwargs))
        os.chdir(dum_start_dir)

        return ret_vals

    @staticmethod
    def loop_over_dirs_nested(nested_list, _passed_func, **kwargs):
        """
        Given a nested list of directories, we loop over each directory, and perform the function.
        The shape of the returned list will match that of the directory list.
        This is the nested version of loop_over_dirs_list()
        :param nested_list
        :param _passed_fun:
        :param kwargs:
        :return:
        """
        dum_start_dir = os.getcwd()
        per_box_vals = []
        for boxID, a_box in enumerate(nested_list):
            per_mol_vals = []
            for molID, a_mol in enumerate(a_box):
                per_rep_vals = []
                for repID, a_dir in enumerate(a_mol):
                    os.chdir(a_dir)
                    per_rep_vals.append(_passed_func(**kwargs))
                per_mol_vals.append(per_rep_vals)
            per_box_vals.append(per_mol_vals)

        os.chdir(dum_start_dir)

        return per_box_vals

    @staticmethod
    def read_nploadtxt(file_name='__CLUS.dat'):
        """
        To be used with directory looping wrapper functions, reads in file_name using np.loadtxt and returns the array.
        :return:
        """
        return np.loadtxt(file_name)


class SystemSetup(object):
    """
    General class that is used to store a particular system. A system is defined as having a unique interaction
    or topology set.
    """

    def __init__(self, comp_num=2, mol_min=100, mol_max=2000):
        self.CompNum = comp_num
        assert self.CompNum > 0, "Number of components must be more than 1!"

        self.MinMax = tuple([mol_min, mol_max])
        assert self.MinMax[0] <= self.MinMax[1], "MinMax should be (min, max), in that order where max>=min."

        self.Structure       = []
        self.StrucQ          = False
        self.WIntEnergyFile  = ''
        self.WIntEnQ         = False
        self.NoIntEnergyFile = ''
        self.NoIntEnQ        = False
        self.Boxes           = []
        self.BoxesQ          = False
        self.MolNums         = []
        self.MolNumsQ        = False
        self.OrthoQ          = False

    def reset_struc(self):
        """
        For the given SysName, we reinitialize the array that holds the topology/structure
        :return:
        """
        self.Structure = []
        self.StrucQ    = False
        return None

    def add_struc(self, sys_struc_ar):
        """
        Given a full structure array, like the ones produced by TopoUtils_Gen_Linear_System,
        we add the structure to this system.
        :param sys_struc_ar:
        :return:
        """
        assert len(sys_struc_ar) == self.CompNum, "Added structure does not match the number of components defined for this system."
        self.Structure = sys_struc_ar
        self.StrucQ    = True
        return None

    def gen_struc_from_mols(self, mol_list):
        """
        Given a set of molecule numbers, we renumber the molecules.
        :param mol_list:
        :return:
        """

        for compID, a_comp in enumerate(mol_list):
            self.Structure[compID][0] = a_comp
        return None

    def set_MinMax(self, mol_min=100, mol_max=2000):
        self.MinMax = tuple([mol_min, mol_max])
        return None

    def set_WInt_energy_file(self, file_name):
        """
        Set the absolute path for the WInt_energy_file for this system.
        :param file_name:
        :return: None
        """
        self.WIntEnergyFile = file_name
        self.WIntEnQ        = True

    def set_NoInt_energy_file(self, file_name):
        """
        Set the absolute path for the NoInt_energy_file for this system.
        :param file_name:
        :return: None
        """
        self.NoIntEnergyFile = file_name
        self.NoIntEnQ        = True

    def calc_mol_sample_linear(self, mol_min=100, mol_max=2000, num_of_bins=15):
        """
        Given the minimum and maximum molecule numbers, we generate a linear sampling.
        This is just a convenient wrapper for np.linspace
        """
        lin_mols = np.linspace(mol_min, mol_max - mol_min, num_of_bins + 1)
        lin_mols = np.array(lin_mols, dtype=int)
        return lin_mols

    def calc_mol_sample_log(self, mol_min=100, mol_max=2000, num_of_bins=15):
        """
        Given the minimum and maximum molecule numbers, we generate a log sampling.
        We log-sample until mol_max/2, then reverse the sampling to generate an even sampling around
        mol_max/2
        """
        actual_bins = int(np.floor((num_of_bins) / 2 + 1))
        log_mols = np.ceil(10.0 ** np.linspace(np.log10(mol_min), np.log10(mol_max / 2.), actual_bins))
        log_mols = np.array(log_mols, dtype=int)
        log_mols = np.sort(np.append(log_mols, mol_max - log_mols[-2::-1]))

        return log_mols

    def calc_mol_sample_hybrid(self, mol_min=100, mol_max=2000, num_of_bins=15):
        """
        Combining the two sampling techniques to generate the hybrid sampling.
        """
        lin_sam = self.calc_mol_sample_linear(mol_min, mol_max, num_of_bins)
        log_sam = self.calc_mol_sample_log(mol_min, mol_max, num_of_bins)

        dum_sam = np.append(lin_sam, log_sam)
        dum_sam = np.unique(dum_sam)

        return dum_sam

    def calc_mol_sample_N_comps(self, num_of_comps=2, mol_min=100, mol_max=2000, num_of_bins=15):
        """
        Used to generate a regular grid of hybrid sampled N-components.
        Suppose M is the mol_max.
        For 1-component, we just return back the hybrid sampling.
        For 2-components, we generate the hybrid sampling for 1-component. Then, since M=m1+m2, m2=M-m1.
        For 3-components and more, we have an iterative process:
            - We generate the hybrid sampling for 1-component, call it R.
            - Using R_{i=0}, we have now M' = M - R_{i=0}, which is used as a new max_mols for N-1 components.
                - This keeps iterating downwards till we reach N=2.
            - We then iterate over all R_i values
        """
        assert num_of_bins % 2 == 1, "We have to set maximum and minimum allowable molecule numbers!"
        hybrid_sam = self.calc_mol_sample_hybrid(mol_min, mol_max, num_of_bins)[::2]
        tot_sam = []
        if num_of_comps == 1:
            tot_sam.append(hybrid_sam)
            return np.array(tot_sam).T
        if num_of_comps == 2:
            tot_sam.append(hybrid_sam)
            tot_sam.append(mol_max - hybrid_sam)
            return np.array(tot_sam).T
        else:
            sup_dum_sam = []
            for numID, a_num in enumerate(hybrid_sam[:]):
                if mol_max - a_num < 2*mol_min:
                    continue
                dum_sam = self.calc_mol_sample_N_comps(num_of_comps=num_of_comps - 1,
                                                       mol_min=mol_min,
                                                       mol_max=mol_max - a_num,
                                                       num_of_bins=num_of_bins).T
                temp_sam = np.zeros((num_of_comps, dum_sam.shape[-1]), int)
                temp_sam[0] = a_num
                temp_sam[1:] = dum_sam
                sup_dum_sam.append(temp_sam.T)
            tot_sam = [a_trip for a_set in sup_dum_sam for a_trip in a_set]
            return np.array(tot_sam)

    def set_hybrid_molecule_numbers(self, num_of_bins=15):
        """
        Generate and set self.MolNums using calc_mol_sample_N_comps()
        :param num_of_bins: Estimated number of bins per component.
        :return: None
        """
        mol_min, mol_max = self.MinMax
        self.MolNums = self.calc_mol_sample_N_comps(num_of_comps=self.CompNum,
                                                    mol_min=mol_min,
                                                    mol_max=mol_max,
                                                    num_of_bins=num_of_bins)
        self.MolNumsQ = True
        return None

    def set_MolNums(self, mol_nums):
        """
        Manually set the list of list of molecule numbers. Remember that this ignores the preset molMin and molMax
        :param mol_nums: np.ndarray(n , self.CompNum)
        :return:
        """
        assert isinstance(mol_nums, np.ndarray),    "Mol Nums should be a NumPy array of NumPy arrays!"
        assert isinstance(mol_nums[0], np.ndarray), "Mol Nums should be a NumPy array of NumPy arrays!"
        assert mol_nums.shape[-1] == self.CompNum,  "Mol Nums does not have the correct number of components"
        self.MolNums = []

        for molID, a_mol in enumerate(mol_nums):
            self.MolNums.append(a_mol)

        self.MolNumsQ = True
        print("Remember that this method ignores mol_min and mol_max.")
        return None

    def calc_struc_beads_in_molecule(self, mol_struc):
        """
        Given a molecule structure, we calculate the total number of beads
        :param mol_struc:
        :return: number of beads in molecule (int)
        """
        this_struc = mol_struc[1][:] # First index is for molecule number
        bead_num = len(np.unique(this_struc.T[0]))
        return bead_num

    def calc_struc_beads_per_molecule(self):
        """
        Loop over each component in the structure to calculate the number of beads.
        :return: nd.array where each component is the number of beads for that component.
        """

        bead_nums = np.zeros(self.CompNum, int)
        for compID, a_struc in enumerate(self.Structure):
            bead_nums[compID] = self.calc_struc_beads_in_molecule(a_struc)
        return bead_nums

    def calc_struc_beads_max(self):
        """
        We loop over the structure for this system and calculate the total number of beads per molecule.
        :return: The maximum number of beads in this system.
        """
        assert self.StrucQ, "The structure has not been defined"
        mol_bead_nums = self.calc_struc_beads_per_molecule()
        bead_nums = np.zeros(len(self.MolNums), int)
        for setID, a_set in enumerate(self.MolNums):
            bead_nums[setID] = np.dot(a_set, mol_bead_nums)
        max_arg = bead_nums.argmax()
        return [max_arg, bead_nums[max_arg], self.MolNums[max_arg]]

    def calc_box_sample_log(self, low_con=1e-6, high_con=1e-2, tot_beads=2e5, tot_boxes=11):
        """
        Given the high and low concentrations, we generate box-sizes that are linearly spaced in log-space
        :param low_con:
        :param high_con:
        :param tot_beads:
        :param tot_boxes:
        :return: box-sizes
        """
        assert 0. < low_con < high_con, "Low conc should be lower than high conc!"
        dum_li = np.linspace(np.log10(low_con), np.log10(high_con), tot_boxes)
        dum_ar = 10. ** dum_li
        dum_ar = tot_beads / dum_ar
        dum_ar = np.array(dum_ar ** (1. / 3.), dtype=int)
        dum_ar_s = np.sort(dum_ar)
        return dum_ar_s

    def calc_box_sample_lin(self, low_con=1e-6, high_con=1e-2, tot_beads=2e5, tot_boxes=11):
        """
        Given the high and low concentrations, we generate box-sizes that are linearly spaced.
        :param low_con:
        :param high_con:
        :param tot_beads:
        :param tot_boxes:
        :return: box-sizes
        """
        assert 0. < low_con < high_con, "Low conc should be lower than high conc!"
        box_lo = int((tot_beads / low_con) ** (1. / 3.))
        box_hi = int((tot_beads / high_con) ** (1. / 3.))
        dum_li = np.linspace(box_lo, box_hi, tot_boxes, dtype=int)
        dum_ar = np.sort(dum_li)
        return dum_ar

    def calc_box_sample_hybrid(self, low_con=1e-6, high_con=1e-2, tot_beads=2e5, tot_boxes=11):
        """
        Use both box-sampling techniques to generate a hybrid box-size array
        :param low_con:
        :param high_con:
        :param tot_beads:
        :param tot_boxes:
        :return:
        """
        lin_box = self.calc_box_sample_lin(low_con, high_con, tot_beads, tot_boxes)
        log_box = self.calc_box_sample_log(low_con, high_con, tot_beads, tot_boxes)
        dum_box = np.append(lin_box, log_box)
        return np.unique(dum_box)[::2]

    def set_regular_boxes(self, low_con=1e-6, high_con=1e-2, tot_boxes=11):
        """
        We calculate the maximum number of beads for this system, and generate the hybrid_sampled box-sizes.
        :param low_con:
        :param high_con:
        :param tot_boxes:
        :return:
        """
        dum_beads   = self.calc_struc_beads_max()
        tot_beads   = dum_beads[1]
        self.Boxes  = self.calc_box_sample_hybrid(low_con, high_con, tot_beads, tot_boxes)
        self.BoxesQ = True
        print(r"Smallest box is {:}, and largest box is {:}".format(self.Boxes.min(), self.Boxes.max()))
        return None

    def set_boxes_to(self, list_of_boxes=np.array([100,110,120])):
        self.Boxes = list_of_boxes[:]
        self.BoxesQ = True
        print(r"Smallest box is {:}, and largest box is {:}".format(self.Boxes.min(), self.Boxes.max()))
        return None

    def calc_mol_boxes_for_conc(self, list_of_mol_nums=np.array([[100,200,300], [200,100,300]]),
                                this_comp=0,
                                target_con=1e-2):
        """
        Given a list of molecule numbers, we generate box sizes such that the concentration
        of 1 of the K-components is at the target_conc
        :param list_of_mol_nums:
        :param this_comp: Which component to look at in particular to set the target concentration
        :param target_con:
        :return: list_of_boxes
        """

        dum_boxes = list_of_mol_nums.T[this_comp] / target_con
        dum_boxes = dum_boxes ** (1./3.)
        dum_boxes = np.array(dum_boxes, int)

        return dum_boxes

    def calc_mol_ortho_boxes_for_conc(self, list_of_mol_nums=np.array([[100,200,300], [200,100,300]]),
                                      target_con=1e-2):
        """
        Given the target concentration and the list of mol_nums, we calculate what the box-size should be if each of the
        individual components' concentration was at target_con
        :param list_of_mol_nums:
        :param target_con:
        :return:
        """
        return np.array([np.array((a_comp / target_con) ** (1./3), dtype=int) for a_comp in list_of_mol_nums.T]).T

    def calc_mol_ortho_filter_boxes(self, list_of_boxes=np.array([[100, 100], [110, 110]]),
                                    this_comp=-1):
        """
        Given a generated list of box-sizes for a target concentration, we select box-sizes where only 1-component,
        this_comp, varies while all others are more-or-less constant. The slight differences are due to lattice-sizes
        being integer values.
        This function prunes the generated box-list to check
        :param list_of_boxes:
        :param this_comp: This component is allowed to vary
        :return: list_of_ids, list_of_boxes
        """

        dum_boxes    = np.delete(list_of_boxes.T, this_comp, axis=0).T
        list_of_ids  = []
        list_of_boxes= []
        for setID, a_set in enumerate(dum_boxes):
            dum_un = np.unique(a_set)
            if len(dum_un) == 1:
                list_of_ids.append(setID)
                list_of_boxes.append(a_set[0])
        return np.array(list_of_ids), np.array(list_of_boxes)


    def set_ortho_boxes(self, mol_num_boxes=11, target_con=1e-2,  comp_ig=-1):
        """
        Given the target concetration for all-components other than comp_ig, we generate the set of boxes
        and molecule numbers such that only 1 out of K-components varies. Simple cyclic
        permutations of the molecule-numbers per box-size produces different components to vary.
        :param mol_num_boxes: estimated number of bins -- usually a sever underestimate
        :param target_con: estimated target concentrations
        :param comp_ig: this component is allowed to vary
        :return:
        """
        self.OrthoQ = True
        self.set_hybrid_molecule_numbers(num_of_bins=mol_num_boxes)

        raw_molecule_ar = self.MolNums[:]
        raw_box_ar      = self.calc_mol_ortho_boxes_for_conc(raw_molecule_ar, target_con)

        pruned_ids, new_box_ar = self.calc_mol_ortho_filter_boxes(raw_box_ar, comp_ig)

        self.Boxes    = new_box_ar[:]
        self.BoxesQ   = True
        self.MolNums  = raw_molecule_ar[pruned_ids]
        print(r"Smallest box is {:}, and largest box is {:}".format(self.Boxes.min(), self.Boxes.max()))
        return None

    def get_number_of_conditions(self):
        """
        Calculate the total number of runs independent concentrations for this system. Let R be the number of
        independent concentrations
        If a regular grid, then we have N x M where N is the number of molecules sampled, and M is the number of
        boxes sampled
        If an orthogonal grid, then R = K x N where we have K components.
        :return:
        """
        assert self.BoxesQ, "The boxes for this system have not been defined!"
        assert self.MolNumsQ, "The molecule numbers for this system have not been calculated!"

        if not self.OrthoQ:
            return len(self.MolNums) * len(self.Boxes)
        else:
            return self.CompNum * len(self.MolNums)

    def get_independent_conditions(self):
        """
        For this system, we return the box-sizes and their corresponding molecule number sets.
        If a regular grid, then we have the same mol_num list for each box-size.
        If an ortho-grid, then for each box-size, we have the cycles of the molecule-numbers.
        :return: [box_list, list_of_list_of_mol_nums]
        """
        assert self.BoxesQ, "The boxes for this system have not been defined!"
        assert self.MolNumsQ, "The molecule numbers for this system have not been calculated!"

        list_of_mol_nums = []

        if not self.OrthoQ:
            list_of_mol_nums = [self.MolNums] * len(self.Boxes)
        else:
            for setID, a_set in enumerate(self.MolNums):
                dum_mol = [np.roll(a_set, a_comp) for a_comp in range(self.CompNum)]
                list_of_mol_nums.append(dum_mol)
        return self.Boxes, np.array(list_of_mol_nums)


class SimulationSetup(object):
    """
    Main class used to instantiate systems that will be simulated.
    """

    def __init__(self, key_file, sys_name_list=None, simulations_path=None, num_of_reps=2):
        """

        :param key_file: Absolute path to the key-file for this instance
        :param sys_name_list: The list of the unique system names.
        :param simulations_path: The absolute path where the simulations will take place.
        """

        if sys_name_list is None:
            sys_name_list = ["SysA"]

        self.KeyFileTemplate = IO_Utils.read_param_file(key_file)

        self.CurrentDirectory = self._gen_curDir()
        self.Sims_Path        = self._gen_simsPath(sims_path=simulations_path)
        self.Data_Path        = self._gen_dataPath(sims_path=simulations_path)
        for a_dir in [self.Sims_Path, self.Data_Path]:
            try:
                os.makedirs(a_dir)
            except OSError:
                pass
        self.PathToLASSI = '/project/fava/packages/bin/lassi'
        self.SysNames    = sys_name_list[:]
        self.Num_Reps    = num_of_reps
        self.Num_Temps   = int(self.KeyFileTemplate[1]['MC_CYCLE_NUM'])
        self.SysInfo     = {}
        for a_sys in sys_name_list:
            self.SysInfo[a_sys] = SystemSetup()
            for a_dir in [self._gen_struc_dir_prefix_for(a_sys),
                          self._gen_data_dir_prefix_for(a_sys)]:
                try:
                    os.makedirs(a_dir)
                except OSError:
                    pass

    def _gen_curDir(self) :
        """
        Generates the string for the current path. Wrapper for OS.
        :return: str for file path.
        """

        if os.name == 'nt':
            return os.getcwd() + '\\'
        else:
            return os.getcwd() + '/'

    def _gen_simsPath(self, sims_path=None):
        """
        Return string of directory where Runs shall be
        :return:
        """
        if os.name == 'nt':
            if sims_path is None:
                return self.CurrentDirectory + 'Runs\\'
            else:
                return sims_path + 'Runs\\'
        else:
            if sims_path is None:
                return self.CurrentDirectory + 'Runs/'
            else:
                return sims_path + 'Runs/'

    def _gen_dataPath(self, sims_path=None):
        """
        Return string of directory where Runs shall be
        :return:
        """
        if os.name == 'nt':
            if sims_path is None:
                return self.CurrentDirectory + 'Data\\'
            else:
                return sims_path + 'Data\\'
        else:
            if sims_path is None:
                return self.CurrentDirectory + 'Data/'
            else:
                return sims_path + 'Data/'

    def _gen_sims_dir_prefix_for(self, sys_name="SysA"):
        """
        Convenient wrapper to generate the directory prefix for the Runs directory
        :param sys_name:
        :return: {str}
        """
        if os.name == 'nt':
            return self.Sims_Path + sys_name + '\\'
        else:
            return self.Sims_Path + sys_name + '/'

    def _gen_data_dir_prefix_for(self, sys_name="SysA"):
        """
        Convenient wrapper to generate the directory prefix for the Data directory
        :param sys_name:
        :param wInt:
        :return: {str}
        """
        if os.name == 'nt':
            return self.Data_Path + sys_name + '\\'
        else:
            return self.Data_Path + sys_name + '/'

    def _gen_struc_dir_prefix_for(self, sys_name="SysA"):
        """
        Convenient wrapper to generate the directory prefix for the Data directory
        :param sys_name:
        :param wInt:
        :return: {str}
        """
        if os.name == 'nt':
            return self.CurrentDirectory + 'Structures\\' + sys_name + '\\'
        else:
            return self.CurrentDirectory + 'Structures/' + sys_name + '/'

    def _get_dirs_list_for(self, sys_name="SysA", wInt=True):
        """
        For this system, return all the directories for this interaction state.
        :param sys_name:
        :param wInt:
        :return:
        """
        dum_boxes, dum_mols = self.SysInfo[sys_name].get_independent_conditions()
        dum_dir_list = IO_Utils.gen_dir_list(dum_boxes,
                                             dum_mols,
                                             self.Num_Reps,
                                             wInt)
        dum_prefix = self._gen_sims_dir_prefix_for(sys_name=sys_name)

        dum_dir_list = [dum_prefix + a_dir for a_dir in dum_dir_list]

        return dum_dir_list

    def _get_dirs_nested_for(self, sys_name="SysA", wInt=True):
        """
        Generate the fully nested list of directories for every possible run condition for this system.
        :param sys_name:
        :param wInt:
        :return:
        """
        dum_boxes, dum_mols = self.SysInfo[sys_name].get_independent_conditions()
        dum_prefix = self._gen_sims_dir_prefix_for(sys_name)
        dum_dirs_nested_list = IO_Utils.gen_dir_nested_list(boxList=dum_boxes,
                                                            molList=dum_mols,
                                                            repNums=self.Num_Reps,
                                                            wInt=wInt,
                                                            dir_prefix=dum_prefix)

        return dum_dirs_nested_list

    def create_dirs_for(self, sys_name="SysA", wInt=True):
        """
        For this system, we create all the directories for this system with the given interaction state
        :param sys_name:
        :return:
        """

        dum_dir_list = self._get_dirs_list_for(sys_name=sys_name, wInt=wInt)
        IO_Utils.create_dirs_from_list(dum_dir_list)

        return None

    def write_strucs_for(self, sys_name):
        """
        We write all possible structures for this system to $CURDIR/Structures/sys_name/
        :param sys_name:
        :return: None
        """
        this_sys = self.SysInfo[sys_name]
        assert this_sys.StrucQ, "This system does not have a defined structure"
        assert this_sys.MolNumsQ, "This system does not have a calculated MolNums"

        dum_mols  = this_sys.MolNums[:]

        if not this_sys.OrthoQ:
            molCycles = 1;
        else:
            molCycles = this_sys.CompNum
        dum_prefix = self._gen_struc_dir_prefix_for(sys_name)
        for setID, a_set in enumerate(dum_mols):
            for a_comp in range(molCycles):
                this_set = np.roll(a_set, a_comp)
                this_sys.gen_struc_from_mols(this_set)
                dum_name   = dum_prefix + "_".join([str(a_num) for a_num in this_set])+".prm"
                IO_Utils.write_struc_file(this_sys.Structure, dum_name)

        return None

    def write_params_for(self, sys_name, wInt=True, fullRand=True):
        """
        Write key-files for LASSI in the appropriate directories
        :param sys_name:
        :param wInt:
        :return:
        """
        this_sys = self.SysInfo[sys_name]

        if wInt:
            assert this_sys.WIntEnQ, "Energy file for WInt has not been set yet!"
        else:
            assert this_sys.NoIntEnQ, "Energy file for NoInt has not been set yet!"

        dum_boxes, dum_mols = this_sys.get_independent_conditions()
        dum_nested_dirs     = self._get_dirs_nested_for(sys_name, wInt)

        dum_struc_dir_prefix = self._gen_struc_dir_prefix_for(sys_name)
        for boxID, (a_box_dir, a_box, a_mol_set) in enumerate(zip(dum_nested_dirs, dum_boxes, dum_mols)):
            for molID, (a_mol_dir, a_mol) in enumerate(zip(a_box_dir, a_mol_set)):
                for repID, (a_rep_dir, a_rep) in enumerate(zip(a_mol_dir, range(self.Num_Reps))):
                    dum_file_name  = a_rep_dir + 'param.key'
                    dum_struc_name = dum_struc_dir_prefix + "_".join([str(a_num) for a_num in a_mol]) + ".prm"
                    IO_Utils.write_param_file(param_obj=self.KeyFileTemplate,
                                              file_path=dum_file_name,
                                              box_size=a_box,
                                              run_name='_',
                                              energy_file=this_sys.WIntEnergyFile,
                                              struc_file=dum_struc_name,
                                              rng_seed=a_rep)

        return None

    def _read_raw_data_files_for(self, sys_name, wInt=True, file_name='__CLUS.dat'):
        """
        For this system, we generate the full nested_directory list, and then read the file_name.
        Remember that to generate the nested_list of directories, we do NOT use SystemSetup.get_independent_conditions()
        :param sys_name:
        :param wInt:
        :param file_name:
        :return:
        """

        dum_nest_list = self._get_dirs_nested_for(sys_name=sys_name, wInt=wInt)

        dum_ret_vals  = IO_Utils.loop_over_dirs_nested(nested_list=dum_nest_list,
                                                       _passed_func=IO_Utils.read_nploadtxt,
                                                       file_name=file_name)

        return dum_ret_vals

    def _save_raw_npy_data_for(self, sys_name, wInt=True, file_name='__CLUS.dat'):
        """
        Convenient wrapper to use np.save the particular data-file for this system.
        :param sys_name:
        :param wInt:
        :param file_name:
        :return:
        """

        dum_data = np.array(self._read_raw_data_files_for(sys_name, wInt, file_name))
        dum_file_name = self._gen_data_dir_prefix_for(sys_name) + file_name
        np.save(dum_file_name, dum_data)

    def _save_raw_pickle_data_for(self, sys_name, wInt=True, file_name='__COMDen.dat'):
        """
       Convenient wrapper to use pickle.dump the particular data-file for this system.
       This is for the types of data that cannot be converted to a numpy-array like the RDFs
       :param sys_name:
       :param wInt:
       :param file_name:
       :return:
       """

        dum_data = self._read_raw_data_files_for(sys_name, wInt, file_name)
        dum_file_name = self._gen_data_dir_prefix_for(sys_name) + file_name + '.b'

        with open(dum_file_name, 'wb') as dFile:
            pickle.dump(dum_data, dFile)

    def save_CorrDen_for(self, sys_name, path_to_norm_data='', norm_naming_func=None, COMDen_file_name='__COMDen.dat'):
        """
        Reads COMDen.dat for this system, and generate the CorrDen.dat and saves it to the data directory.
        We save the [xAr_list, corr_den] using pickle
        :param sys_name:
        :return:
        """

        this_sys = self.SysInfo[sys_name]

        dum_comden_name = self._gen_data_dir_prefix_for(sys_name) + COMDen_file_name + '.b'

        dum_comden_data = []
        with open(dum_comden_name, 'rb') as dFile:
            dum_comden_data = pickle.load(dFile)

        dum_comden_analysis = _COMDen_Utils(total_comden_data=dum_comden_data,
                                            temp_nums=self.Num_Temps,
                                            comp_nums=this_sys.CompNum,
                                            box_list=this_sys.Boxes)

        dum_comden_analysis.reshape_raw_data(wMode=1)
        dum_comden_analysis.gen_corr_den(path_to_norm_data=path_to_norm_data,
                                         naming_func=norm_naming_func)

        dum_corrden_name = self._gen_data_dir_prefix_for(sys_name) + '__CorrDen.dat.b'
        with open(dum_corrden_name, 'wb') as dFile:
            pickle.dump([dum_comden_analysis.xArr_list[:], dum_comden_analysis.corr_data[:]], dFile)

        return None

    def save_CorrDen_MolTypeCOM_for(self, sys_name, path_to_norm_data='', norm_naming_func=None, COMDen_file_name='__COMDen.dat'):
        """
        Reads COMDen.dat for this system, and generate the CorrDen.dat and saves it to the data directory.
        We save the [xAr_list, corr_den] using pickle
        :param sys_name:
        :return:
        """

        this_sys = self.SysInfo[sys_name]

        dum_comden_name = self._gen_data_dir_prefix_for(sys_name) + COMDen_file_name + '.b'

        dum_comden_data = []
        with open(dum_comden_name, 'rb') as dFile:
            dum_comden_data = pickle.load(dFile)

        dum_comden_analysis = _COMDen_Utils(total_comden_data=dum_comden_data,
                                            temp_nums=self.Num_Temps,
                                            comp_nums=this_sys.CompNum,
                                            box_list=this_sys.Boxes)

        dum_comden_analysis.reshape_raw_data(wMode=0)
        dum_comden_analysis.gen_corr_den(path_to_norm_data=path_to_norm_data,
                                         naming_func=norm_naming_func)

        dum_corrden_name = self._gen_data_dir_prefix_for(sys_name) + '__CorrDen_MolType.dat.b'
        with open(dum_corrden_name, 'wb') as dFile:
            pickle.dump([dum_comden_analysis.xArr_list[:], dum_comden_analysis.corr_data[:]], dFile)

        return None


class _COMDen_Utils(object):
    """
    Collection of funcitons to analyze and manipulate COMDen data files. It is assumed that the data are generated by
    SimulationSetup._save_raw_pickle_data_for() for the COMDen data.
    Therefore, the indexing structure of COMDen_data should be COMDen_Data[boxID][numID][repID][compID & tempIDs mixed].
    Since the last index is a convolution, the first thing we do is reshape the data into a more tractable form.

    reshaped_data[boxID][tempID][numID][compID][repID][:] so then each box-size can be stored as a numpy array
    """

    def __init__(self, total_comden_data, temp_nums=3, comp_nums=3, box_list=np.array([100, 110, 120])):
        """
        Given the total_comden_data, we initialize this object by storing the relavent dimensions.
        :param total_comden_data:
        :param temp_nums: The total number of different temperatures probed in this data.
        :param comp_nums: The total number of different components in the system.
        """
        self.raw_data         = total_comden_data[:] #Store a copy to be careful
        self.ReshapeQ         = False
        self.CorrDenQ         = False
        self.num_boxes        = len(self.raw_data)
        self.num_mols         = len(self.raw_data[0])
        self.num_reps         = len(self.raw_data[0][0])
        self.num_temps        = temp_nums
        self.num_comps        = comp_nums
        self.tot_comden_comps = self.num_comps * (self.num_comps + 1) # Total number of data lines per temperature per mode
        self.Boxes            = box_list[:]
        self.corr_data        = []
        self.xArr_list        = []

    def reshape_raw_data(self, wMode=1):
        """
        :param wMode = 0 for COM from MolType largest cluster, 1 for all molecules of MolType
        This reshapes the data to have the following structure
        :return:
        """
        dum_comps_list = self.gen_pairs_list(totMols=self.num_comps)
        per_box_data = []
        for boxID, com_den_of_box in enumerate(self.raw_data[:]):
            dum_box_data  = np.array(com_den_of_box)
            per_temp_data = []
            for tempID in range(self.num_temps):
                per_num_data = []
                for numID in range(self.num_mols):
                    per_comp_data = []
                    for compID, (compA, compB) in enumerate(dum_comps_list):
                        per_rep_data = []
                        for repID in range(self.num_reps):
                            this_idx = self.index_comp_to_comp_wMode_of_temp(compA=compA, compB=compB, temp_id=tempID, mode=wMode)
                            # print(this_idx)
                            per_rep_data.append(dum_box_data[numID, repID, this_idx])
                        per_comp_data.append(per_rep_data)
                    per_num_data.append(per_comp_data)
                per_temp_data.append(per_num_data)
            per_box_data.append(np.array(per_temp_data))
        self.raw_data = per_box_data
        self.ReshapeQ = True

        return None

    def gen_pairs_list(self, totMols=None):
        """
        Generates the comp-to-comp list given this many total_molecules.
        The list is like [[-1,0], [-1,1], ...,
                          [0,0],  [0,1], ...,
                          ]
        Every possible pair, and in order. With the addition of the -1,n pairs which correspond to the system's COM
        :param totMols:
        :return:
        """
        if totMols is not None:
            totMols = self.num_comps
        dum_pairs_list = []

        for i in range(totMols):
            dum_pairs_list.append([-1, i])

        for i in range(totMols):
            for j in range(totMols):
                dum_pairs_list.append([i,j])
        return np.array(dum_pairs_list)

    def index_comp_to_comp(self, compA=0, compB=0, tot_mols=None):
        """
        From the unshaped raw data, we return the index that corresponds the the density profile of compB with respect
        to the COM of compA
        :param compA:
        :param compB:
        :param tot_mols:
        :return:
        """

        if tot_mols is None:
            tot_mols = self.num_comps

        if compA < 0:
            return compB
        else:
            return tot_mols + compB + tot_mols * compA

    def index_comp_to_comp_wMode(self, compA=0, compB=0, tot_comden_comps=None, mode=0):
        """
        The raw index adjusted for the internal shift in LaSSI where if mode=0, we only look at the beads in the largest
        cluster for compA, and with mode=1 we look at all beads.
        :param compA:
        :param compB:
        :param tot_comden_comps:
        :param mode:
        :return:
        """
        assert 1 >= mode >= 0, "Mode can only be 0 or 1"

        if tot_comden_comps is None:
            tot_comden_comps = self.tot_comden_comps

        return tot_comden_comps * mode + self.index_comp_to_comp(compA, compB)

    def index_comp_to_comp_wMode_of_temp(self, compA=0, compB=0, tot_comden_comps=None, mode=0, temp_id=0):
        """
        Convenient wrapper that includes the shift due to tempID
        :param compA:
        :param compB:
        :param tot_comden_comps:
        :param mode:
        :param temp_id
        :return:
        """
        if tot_comden_comps is None:
            tot_comden_comps = self.tot_comden_comps
        return temp_id * 2 * tot_comden_comps + self.index_comp_to_comp_wMode(compA, compB, tot_comden_comps, mode)

    def gen_corr_den(self, path_to_norm_data, naming_func=None):
        """
        Given the total path to the normalization data, we loop over every box and normalize the number distributions
        to generate the density distributions. Particularly, we have N(r_i) / N_0(r_i) and we only perform the divisions
        where N_0(r_i) != 0.
        One can provide the naming function for the normalization files where the naming_func only takes box_size
        as the argument.
        :param path_to_norm_data:
        :return:
        """
        assert self.ReshapeQ, "The raw data needs to be reshaped first. Use the reshape_raw_data() method!"

        self.xArr_list = []

        per_box_data   = []
        for boxID, a_box in enumerate(self.Boxes):
            dum_norm_data = self._get_norm_data(path_to_norm_data=path_to_norm_data,
                                                box_size=a_box,
                                                naming_func=naming_func)
            good_pts      = np.where(dum_norm_data != 0)
            norm_den_cor  = dum_norm_data[good_pts]
            dum_xAr       = np.arange(0, a_box, 0.25)
            xAr_cor       = dum_xAr[good_pts]
            self.xArr_list.append(xAr_cor)

            per_temp_data = []
            for tempID in range(self.num_temps):
                per_num_data = []
                for numID in range(self.num_mols):
                    per_comp_data = []
                    for compID in range(self.tot_comden_comps):
                        per_rep_data = []
                        for repID in range(self.num_reps):
                            dum_yAr = self.raw_data[boxID][tempID, numID, compID, repID]
                            per_rep_data.append(dum_yAr[good_pts] / norm_den_cor)
                        per_comp_data.append(per_rep_data)
                    per_num_data.append(per_comp_data)
                per_temp_data.append(per_num_data)
            per_box_data.append(np.array(per_temp_data))

        self.corr_data = per_box_data[:]
        self.CorrDenQ  = True

        return None

    def _get_norm_data(self, path_to_norm_data, box_size, naming_func=None):
        """
        Fetch the normalization data given the box-size.
        :param path_to_norm_data:
        :param box_size:
        :param naming_func:
        :return:
        """
        assert box_size > 1, "Box-size should be bigger than 1!"

        if naming_func is None:
            naming_func = self._gen_norm_filename

        dum_file_name = path_to_norm_data + naming_func(box_size)

        return np.loadtxt(dum_file_name)

    @staticmethod
    def _gen_norm_filename(box_size):
        """
        Given the box-size, we generate the file-name for the normalization. This internal function assumes the file-names
        are P0_S_{box_size}.dat
        :param box_size:
        :return:
        """
        return f"P0_S_{box_size}.dat"

    @staticmethod
    def _gen_pairs_list(totMols=1):
        """
        Generates the comp-to-comp list given this many total_molecules.
        The list is like [[-1,0], [-1,1], ...,
                          [0,0],  [0,1], ...,
                          ]
        Every possible pair, and in order. With the addition of the -1,n pairs which correspond to the system's COM
        :param totMols:
        :return:
        """
        assert isinstance(totMols, int) and totMols > 0, "Positive integers only."

        dum_pairs_list = []
        for i in range(totMols):
            dum_pairs_list.append([-1, i])

        for i in range(totMols):
            for j in range(totMols):
                dum_pairs_list.append([i,j])
        return np.array(dum_pairs_list)


class _COMDen_Analysis(object):
    """
    Collection of functions to help with the analysis of COMDen data to calculate the densities of coexisting phases.
    It is assumed that initializer is given the [xAr_list, corr_data] generated from _COMDen_Utils() either directly from
    the class, or read from a previously saved file.
    """

    def __init__(self, tot_corr_den_obj, num_of_comps=3):

        self.xArList   = tot_corr_den_obj[0][:]
        self.CorrData  = tot_corr_den_obj[1][:]
        self.Num_Comps = num_of_comps
        self.PairsList = _COMDen_Utils.gen_pairs_list(self.Num_Comps)


class _Radial_Func_Utils(object):
    """
    Collection of functions to manipulate data that is sampled in radial distance like density distributions and
    pair-correlations.
    """

    def _gen_smooth_sum(dAr, step_size=4):
        """
        Extremely simple implementation of smoothing using sums.
        :param dAr: The data-array to be smoothed
        :param nAr:
        :param pAr:
        :param step_size:
        :return:
        """

        assert len(dAr) > step_size, "Step-size should be smaller than the input"
        dum_l = len(dAr)
        new_d = []

        for i in range(0, int(dum_l/step_size)):
            new_d.append(dAr[i * step_size : (i+1) * step_size].sum())

        return np.array(new_d)


class _DataWithError(object):
    """
    Since I am tired of implementing error propoagation for data that are like $(x, f(x), \delta f(x))$, I
    have decided to write a general-ish object that can be used to manipulate such data.
    In particular, simple operations like adding, subtracting, multiplying, dividing, and taking the logarithm will
    be implemented.

    - For addition and subtraction, the error is quadrature summed: $\sqrt{ \delta a^2 + \delta b^2 + ...}$
    - For multiplication, if $g = a \times b$, we have $ g \times \sqrt{ (\delta a / a)^2 + (\delta b / b)^2 }$.
    - For division, if $ g = a / b$, then we have to make sure that we only calculate the ratio where $b \neq 0$.
      The error is the same as multiplication. Furthermore, we filter out all the non-calculated values from x as well.
    - For logarithms, again we make sure to only calculate for non-zero quantities, and filter. If $ g = \log(a)$, then
      $\delta g = \delta a / a$, as the relative error.
    """

    def __init__(self, x_vals, y_vals, e_vals):
        """
        Initializer must get x-coordinates, y-values and their corresponding error-values.
        """
        assert len(x_vals) == len(y_vals) == len(e_vals), "All arrays must be the same length!"
        assert isinstance(x_vals, np.ndarray) and isinstance(y_vals, np.ndarray) and isinstance(e_vals, np.ndarray)
        self._x   = x_vals
        self._y   = y_vals
        self._e   = e_vals
    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, this_val):
        self._x = this_val
    @x.deleter
    def x(self):
        del self._x

    @property
    def y(self):
        return self._y
    @y.setter
    def y(self, this_val):
        self._y = this_val
    @y.deleter
    def y(self):
        del self._y

    @property
    def e(self):
        return self._e
    @e.setter
    def e(self, this_val):
        self._e = this_val
    @e.deleter
    def e(self):
        del self._e

    def __len__(self):
        return len(self.x)

    def __add__(self, other):
        """
        Addition of two data with error. The y-values are added, while the errors are added
        in quadrature.
        """
        assert isinstance(other, _DataWithError), "Both objects need to be of DataWithError type"
        assert len(self) == len(other), "Both objects have to be the same length."

        new_vals = self.y + other.y
        new_err  = np.sqrt(self.e**2. + other.e**2.)

        return _DataWithError(self.x, new_vals, new_err)

    def __sub__(self, other):
        """
        Subtraction of two data with error. The y-values are subtracted, while the errors are added
        in quadrature.
        """
        assert isinstance(other, _DataWithError), "Both objects need to be of DataWithError type"
        assert len(self) == len(other), "Both objects have to be the same length."

        new_vals = self.y - other.y
        new_err  = np.sqrt(self.e**2. + other.e**2.)

        return _DataWithError(self.x, new_vals, new_err)

    def __mul__(self, other):
        """
        Multiplication of two data with error. The y-values are multiplied, while the relative errors are added
        in quadrature.
        """
        assert isinstance(other, _DataWithError), "Both objects need to be of DataWithError type"
        assert len(self) == len(other), "Both objects have to be the same length."

        new_vals = self.y * other.y
        new_err  = new_vals * np.sqrt((self.e / self.y)**2. + (other.e / other.y)**2.)

        return _DataWithError(self.x, new_vals, new_err)

    def __truediv__(self, other):
        """
        Division of two data with error. Firstly, we find indecies where other.y is non-zero.
        The y-values are multiplied, while the relative errors are added
        in quadrature.
        """
        assert isinstance(other, _DataWithError), "Both objects need to be of DataWithError type"
        assert len(self) == len(other), "Both objects have to be the same length."

        gd_pts   = other.y != 0

        new_vals = self.y[gd_pts] / other.y[gd_pts]
        new_err  = new_vals * np.sqrt((self.e[gd_pts] / self.y[gd_pts])**2. + (other.e[gd_pts] / other.y[gd_pts])**2.)

        return _DataWithError(self.x[gd_pts], new_vals, new_err)

    def __pow__(self, power=2):
        """
        Raise the data to the power. For error, if $f = a^n$, then $\delta f = n * ( a ^ {n-1}) \delta a$
        """
        assert power is not 0, "Does not make sense for exponent of 0!"
        if power >= 1:
            new_x = self.x
            new_y = self.y ** power
            new_e = power * (self.y ** (power-1)) * self.e
        if power <= 1:
            gd_pts= self.y != 0

            new_x = self.x[gd_pts]
            new_y = self.y[gd_pts] ** power
            new_e = power * (self.y[gd_pts] ** (power-1)) * self.e[gd_pts]
        else:
            raise ValueError("Exponent of 0 is non-sensical here.")

        return _DataWithError(new_x, new_y, np.abs(new_e))

    def __abs__(self):
        """
        Return the absolute of the y-values. The error is left unchanged.
        """
        return _DataWithError(self.x, np.abs(self.y), self.e)

    def _non_negative(self):
        """
        Look for indecies where self.y is non-zero. Then we return a new object that has only the non-zero values.
        """
        gd_pts   = self.y != 0

        return _DataWithError(self.x[gd_pts], self.y[gd_pts], self.e[gd_pts])

    def _log10(self):
        """
        Calculate the log10 of the data and propagate the error through, which is relative error.
        """

        pos_data = self._non_negative()

        new_y    = np.log10(pos_data.y)
        new_e    = pos_data.e / pos_data.y

        return _DataWithError(pos_data.x, new_y, new_e)












