cgHeliParm.py for DNA CG Martini trajectories
========
cgHeliParm.py has been developed for the analysis of the DNA double stranded structures in molecular dynamics (MD) simulations. Like 3DNA or Curves+, it defines a geometric reference frame for each base/base pair to calculate several structural descriptors of DNA from the GROMACS CG Martini MD trajectory. As output the helical descriptors are saved in individual files as a function of time.

<strong> Note: </strong> MDAnalysis must be installed together with the NumPy module. PDB files can also be analyzed using the script cgHeliParmPDB.py. When running cgHeliParm.py, look for the library_path that links to the 'data' folder in your system and modify it accordingly.

**Last Update: 11 Ago. 2015**

For a tutorial on how to use cgHeliParm.py, please visit [Martini's group home-page](http://md.chem.rug.nl/cgmartini/index.php/tutorial-martini-dna?a_id=358).

<strong> Please cite the following publications:</strong>

Lu, X.-J. and Olson, W.K. (2003)
[3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures.](http://nar.oxfordjournals.org/content/31/17/5108.short)
_Nucleic Acids Res_. 31(17), 5108-21.

Michaud-Agrawal, N. et al., (2011) 
[MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.] (http://onlinelibrary.wiley.com/doi/10.1002/jcc.21787/full)
_J Comp Chem_, 32(10), 2319–2327.

van der Walt, S. et al. (2011) 
[The NumPy Array: A Structure for Efficient Numerical Computation.] (http://ieeexplore.ieee.org/document/5725236/)
_Comput Sci Eng_, 13, 22–30.
