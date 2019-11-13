Branched Molecule Example
-------------------------

To illustrate the use of LASSI to simulate branched polymers, we have
a system that contains 2 different H-shaped molecules. In keeping with
the previous examples, the two molecule types have only one type of
sticker each. Thus, in-cis anisotropic interactions aren't possible.
Furthermore, in order to mimic a folded domain like structure, there
is an additional sticker type that acts as the branching center for
the H and has a very large overlap energy.

As stated in `structure_BR.prm`, it is usually a good idea to *draw*
the desired structure on paper, for example, before attempting to
write the structure. 

The python helper scripts can be used to generate branched structures,
but for now only simple symmetric branched polymers can be generated.

Run using `./../../src/lassi -k param_BR.key`.



