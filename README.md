cgHeliParm.py for DNA/RNA CG Martini trajectories
========
cgHeliParm.py has been developed for the analysis of the DNA/RNA double stranded structures in molecular dynamics (MD) simulations. Like 3DNA or Curves+, it defines a geometric reference frame for each base/base pair to calculate several structural descriptors of DNA/RNA from the GROMACS CG Martini MD trajectory. As output the helical descriptors are saved in individual files as a function of time.

<strong> Note: </strong> since cgHeliParm.py uses MDAnalysis, it can be used with trajectory files from other MD packages such as NAMD or AMBER (only NETDCF files are accepted).  PDB files can also be analyzed.

To execute cgHeliParm.py, MDAnalysis python package should be installed.

**Last Update: 11 Ago. 2015**

For a tutorial on how to use cgHeliParm.py, please visit [Martini's group home-page](http://md.chem.rug.nl/cgmartini/index.php/tutorial-martini-dna?a_id=358).

<strong> Please cite the following publications:</strong>
Xiang-Jun Lu & Wilma K. Olson (2003)
[3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures.](http://nar.oxfordjournals.org/content/31/17/5108.short)
_Nucleic Acids Res_. 31(17), 5108-21.

Michaud-Agrawal, N. et al., (2011) 
[MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.] (http://onlinelibrary.wiley.com/doi/10.1002/jcc.21787/full)
_J Comp Chem_, 32(10), 2319–2327.
