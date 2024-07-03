Single Molecule Example - DiBlock Copolymer
-------------------------------------------

This example illustrates the basics of making key files, energy files and structure files.
Go through `structure_SM.prm` `energy_SM.prm` and `param_SM.key` in that order.
Then, assuming that the code is compiled, run using `./../../src/lassi -k param_SM.key`

After the simulation/run has finished, you can use VMD to visualize the trajectory using `Visualize_SM.vmd`

For the format of the energy file (the energy file only accepts
a specific set of keywords, and no comments), we have:

```
#STICKERS
2

## '#STICKERS' must be the first line in an energy file. (It tells LASSI how to build build and store energy matrices)
## The line following #STICKERS defines the total number of different bead types in the system. Since we have a
## di-block co-polymer, we only have 2 different bead types.

## The following keyword defines the isotropic interactions between the different beads, and is in the form of a symmetric
## matrix. This interaction is the overlap cost in a BFM. The matrix has the form:
##          |BeadType
## Bead_Type|0|1|
##          |1|1|
## where E_OVLP[0][0] defines the energy between 0 and 0 bead types, and so on.

#OVERLAP_POT
0.2 -0.2
-0.2 0.2

## The following two keywords are not yet fully implemented, but should add more isotropic interactions between beads
## where the _POT defines the energies, and _RAD defines the cut-off radii for those interactions
## NOTE: LEAVE THESE AS 0 for the time being.

#CONTACT_POT
0.0 0.0
0.0 0.0

#CONTACT_RAD
0.0 0.0
0.0 0.0

## The following keyword defines the anisotropic interaction energies between beads, or stickers in this case.

#SC_SC_POT
0.0 -2.0
-2.0 0.0

## The following keyword defines an internal linker length where the linker length between beads as defined in the
## structure file is LINKER_LENGTH*(length from structure).
## I usually just leave this at 1, and change the structure file, but if I want all linker lengths to be longer, changing
## this is faster.

#LINKER_LENGTH
1.0


## The following two keywords are not yet implemented, but add a harmonic bond interaction between beads, given a stiffness
## and an equilibrium length.

#LINKER_SPRCON
0.0

#LINKER_EQLEN
1.0
```