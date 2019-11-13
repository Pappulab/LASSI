LASSI - Examples
----------------

Assuming that you have compiled the code, for the following 4 examples, simply go to the appropriate directory,
and read the README, and structure file. Then, run with `./../../src/lassi -k {name of key file}`

##### A - Single Molecule



##### B - Di-block Copolymer System


##### C - Two Component System


##### D - Branched System



#### Miscellaneous
The general workflow for running simulations in general is:
1. Prepare a keyfile - call it keyFileName
2. Prepare an energy file.
3. Prepare a structure file.
4. Run using `./${PATH_TO_LASSI_EXEC} -k ${keyFileName}`

In order to generate phase diagrams in the concentration vs. temperature
plane, one has to run independent simulations for each concentration,
or the size of the simulation box within the keyfile. For each box size,
one then defines the initial and final temperatures, and the number of
temperature samples to take in a linear scale.

The LASSI module in python provided in `python/LASSI-Helper` can be used
to make some of the process of generating a binodal easier. For example,
making simple structures, calculating the correct box sizes for a
concentration range, and generating the appropriate directory tree and 
writing out independent keyfiles for each independent simulation. Then,
all that is required is to run the simulations whichever way your
architecture deems right. Then further, the LASSI module can be used
to collect and process the data, and generate the appropriate measurables
to detect the phase transitions.

