cgHeliParm for DNA CG Martini trajectories
========
cgHeliParm has been developed for the analysis of the DNA double stranded structures in molecular dynamics (MD) simulations. Like 3DNA, it defines a reference frame for each base/base pair to calculate several structural descriptors of DNA from the GROMACS CG Martini MD trajectory. As output, the helical descriptors are saved in individual files as a function of time, either in txt or json formats.

<strong> Note: </strong> cgHeliParm uses MDAnalysis and this must be installed together with the NumPy module with python2.7. Currently MDAnalysis does not support python3, which hinders cgHeliParm supporting python3 for the moment. Once there is a version compatible with python3, we will add it. It also uses mdreader, a subclass in MDAnalysis. You can download the last version of mdreader.py from [here](https://github.com/mnmelo/MDreader). PDB files can be analyzed using the script cgHeliParmPDB.py. When running cgHeliParm.py, look for the library_path that links to the 'data' folder in your system and modify it accordingly. Finally, there is no need to download Gromacs to use cgHeliParm.py (or cgHeliParmPDB.py).

**Last Update: 11 Ago. 2017**

For a tutorial on how to use cgHeliParm.py, please visit [Martini's group home-page](http://md.chem.rug.nl/cgmartini/index.php/tutorial-martini-dna?a_id=358).

<strong> Please cite the following publications:</strong>

Faustino, I. and Marrink, S.J. (2017)
[cgHeliParm: analysis of dsDNA helical parameters for coarse-grained MARTINI molecular dynamics simulations.](https://doi.org/10.1093/bioinformatics/btx444)
_Bioinformatics_. btx444.

Lu, X.-J. and Olson, W.K. (2003)
[3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures.](http://nar.oxfordjournals.org/content/31/17/5108.short)
_Nucleic Acids Res_. 31(17), 5108-21.

Michaud-Agrawal, N. et al., (2011) 
[MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.](http://onlinelibrary.wiley.com/doi/10.1002/jcc.21787/full)
_J Comp Chem_, 32(10), 2319–2327.

van der Walt, S. et al. (2011) 
[The NumPy Array: A Structure for Efficient Numerical Computation.](http://ieeexplore.ieee.org/document/5725236/)
_Comput Sci Eng_, 13, 22–30.
