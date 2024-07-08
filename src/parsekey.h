#ifndef _PARSEKEY_H_ // include guard
#define _PARSEKEY_H_

#include "global.h"

//#include <stdio.h>
//#include <string.h>

int Parse_Keyfile(char* filename);

int Parse_EnergyFile(char* strEnFile);

void Parse_StructureFile(char* filename);

void Parse_StructureFile_CalcBeadsAndChains(char* filename, size_t* n_bead_num, size_t* n_chain_types,
                                            size_t* n_chain_num);

void CreateBeadsAndChains(size_t n_bead_num, size_t n_chain_num);
#endif // _PARSEKEY_H_
