==========================================================================
LASSI - LAttice Simulation Engine For Stickers And Spacer Interactions
==========================================================================

`version 0.1.0 - November 2019`

------------
Installation
------------

To install, run src/compile.sh in bash, and add the src/ folder to your `$PATH`

To test if the code compiled, go to `/runs/Test_Simulation` and run `./lassi`

------------
Introduction
------------

LASSI is a versatile simulation engine designed to perform Monte-Carlo (MC) simulations for phase transitions of explicit polymer instantiations of the *Stickers-and-Spacers* formalism. Furthermore, LASSI is designed, from the ground up, around  anisotropic interactions between individual monomers. LASSI can be used to generate full phase diagrams for a given polymer architecture (although rings aren't yet implemented).

-------------------
Tutorial & Examples
-------------------

There are 4 examples in the `examples/` directory that go over the basic ways to use LASSI with different types of systems. It is highly recommended that new users go through them in order.

Furthermore, there is a short tutorial in `/python/example` that illustrates using some python scripts to help generate phase diagrams. In effect, all that is required to generate a phase diagram is a structure file for the system of interest, a parameter file, an energy file with interactions and one without. Then, we pick a set of concentrations and run independent runs for each energy file, for each concentration (or box-size). The python script is my way of avoiding bash scripting, and also has convenient functions that will generate the correct directory trees, generate keyfiles for each condition, run the simulations (or submit them to a queue if you want), collect and process all the data from the runs. Furthermore, there are also functions to generate structure files that can contain simple linear polymers, linear polymers with explicit linkers, and symmetric branched polymers. Note that the tutorial is specific to the computing architecture at PappuLab and changes will be necessary for it to work on your machine.

-------------
Documentation
-------------

For detailed information about the underlying physics and applied examples, please see [LASSI: A lattice model for simulating phase transitions of multivalent proteins](https://doi.org/10.1371/journal.pcbi.1007028).
