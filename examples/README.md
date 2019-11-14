LASSI - Examples
----------------

Assuming that you have compiled the code, for the following 4 examples, simply go to the appropriate directory,
and read the README, and structure file. Then, run with `./../../src/lassi -k {name of key file}`.

It is recommended that the examples are done in order because they build upon the previous
ones.

#### A - Single Molecule (Di-block Copolymer)

In order to illustrate the basics of making key files, energy files, and structure files,
we have a single di-block copolymer of the form `AAAAAAAAAABBBBBBBBBB`. We have a chain
that contains two different *stickers*, where the only anisotropic interactions are between
`A` and `B`.

#### B - Di-block Copolymer System

We then simulate 500 molecules from example 5, where we also have an annealing cycle to take
the system from a lower temperature where the system phase separates to a higher temperature
where the condensate dissolves.

#### C - Two Component System

Now, we take the previous system and break it into two co-polymers of `AAAAAAAAAA` and `BBBBBBBBBB`,
with a similar annealing scheme as example B. Thus, we have a multi-component system containing 
two types of molecules.

#### D - Branched System

Lastly, this examples  illustrates how to simulate systems where the molecules
are branched. The system contains H-shaped branched polymers, where the central bead
has a high excluded volume cost -- mimicking a folded domain. In this example, the 
central bead is a pure *spacer* because it does not interact anisotropically.

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

The `LASSI` module in python provided in `python/LASSI-Helper` can be used
to make some of the process of generating a binodal easier. For example,
making simple structures, calculating the correct box sizes for a
concentration range, and generating the appropriate directory tree and 
writing out independent keyfiles for each independent simulation. Then,
all that is required is to run the simulations whichever way your
architecture deems right. Then further, the `LASSI` module can be used
to collect and process the data, and generate the appropriate order parameters
to detect the phase transitions.

