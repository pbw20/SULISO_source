# SULISO_sample_files
Test files for the SULISO suite, showing various capabilities of 
the programs

The first three letters of the extensions refer to the program 
they are input for. These same files with the "out" suffix, show 
what can be expected as output for each program.

These examples illustrates the different roles of the LIPFR and 
UJISO programs. The same structure (a molecular cluster of a methyl 
group with five water molecules (arranged in a highly constrained 
symmetrical manner) is used with each of CAMVIB, LIPFR, CUTOFF and 
UJISO. This system yields 5 imaginary frequencies associated with 
genuine vibrational modes (as well as 6 zero-frequencies for 
translational and rotational motions of the whole molecular cluster) 
as a consequence of the symmetry constraints deliberately imposed 
on this system. The point is that the isotopic partition function 
ratios may be computed even in the presence of imaginary frequencies 
(which are treated as if they were real for this purpose); we have 
shown elsewhere (Ruiz-Pernía JJ, Williams IH. Chem Eur J 2012;18:
9405-9414) that this procedure leads to better behaviour in regard 
to their isotopic sensitivity than if the imaginary frequencies are 
ignored. Other codes simply remove or ignore imaginaries and the 
6 translations and rotations, whereas LIPFR and UJISO allow the user 
to decide how they should be treated. Application of CUTOFF followed 
by UJISO to this same example excludes the atoms and degrees of 
freedom responsible for the imaginary frequencies while retaining 
the smaller subset of atoms and degrees of freedom responsible for 
the isotope effects themselves. 

