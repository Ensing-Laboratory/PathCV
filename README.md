![figure](https://aip.scitation.org/na101/home/literatum/publisher/aip/journals/content/jcp/2018/jcp.2018.149.issue-7/1.5027392/20180613/images/large/1.5027392.figures.online.f1.jpeg)


# Path collective variables with Plumed 2.5

The `PathCV.cpp` source file found in this repository is the latest stable version of our path collective variable (Path-CV) code. It can be downloaded 
and then be compiled together with the [Plumed](https://www.plumed.org/download) program as described below.

Plumed is a software package that can be linked to various molecular simulation programs to carry out enhanced sampling simulations 
and study activated molecular transitions, such as chemical reactions, phase transitions, protein folding and conformational transitions, 
and so forth. 

By using a Path-CV, complex transitions that are difficult to capture with just a few collective variables can be described and biased along using
a biasing method, such as metadynamics, umbrella sampling, constrained MD, ABF, OPES, etc. A Path-CV is an adaptive function of other 
descriptive collective variables, such as bond distances, angles, coordination numbers and so on, which can be used to approach and describe as good
as possible the unknown intrinsic reaction coordinate of a molecular transition, by optimizing the function towards the mean transition path
or the minimum free energy path. The computational cost of biasing along the one-dimensional path-CV scales almost trivially with the number
of CVs used in the path-CV.

### Installation

To carry out path-metadynamics or other enhanced sampling simulations using our path collective variable as described in Refs. [1, 2], 
download the `PathCV.cpp` source file.

Save this file in your plumed-2 distribution in the directory `src/function/PathCV.cpp`. 
We are using [Plumed version 2.5](https://www.plumed.org/download) at the moment.

After recompiling the code (`make, make install, make doc`) the documentation on how to use the path-CV is found in `user-doc/html/index.html` 
under `Collective variables / Functions / PATHCV`.

Notice that [Plumed](https://www.plumed.org/download) already comes with a version of our Path-CV implementation. 
See for example [here](https://www.plumed.org/doc-v2.8/user-doc/html/_p_a_t_h.html). However, at the moment (January 2023), 
that version does not contain several useful features, such as `genpath` (to generate an initial linear path) and multiple walkers.

### Plumed eggs with example input files 

For a [Plumed Nest](https://www.plumed-nest.org) archive with example input files for a path-metadynamics simulation of DNA base-pairing 
transitions using Plumed and Gromacs 
click on [this link](https://www.plumed-nest.org/eggs/21/033/). Further instructions: This allows to run multiple-path-metadynamics on 
the Watson-Crick-Franklin to Hoogsteen base-pairing transition of an A-T pair in DNA. It requires PLUMED compiled with MPI and with the 
Path-CV code provided within. It also requires an MD engine that can run parallel replicas. We use GROMACS 2018.8 compiled with MPI. 
Notice that in the PLUMED input files "WALKERS_ID" must be adjusted for the different walkers. 
This example is taken from the work published in Ref. [3].

A second Plumed Nest example is found [here](https://www.plumed-nest.org/eggs/21/049/).


### References
[1] Path finding on high-dimensional free energy landscapes. Grisell Díaz Leines and Bernd Ensing, _Phys. Rev. Lett._ __109__ (2012), 020601,  [DOI: 10.1103/PhysRevLett.109.020601](https://doi.org/10.1103/PhysRevLett.109.020601).

[2] Advances in enhanced sampling along adaptive paths of collective variables. Alberto Pérez de Alba Ortíz , Ambuj Tiwari, Rakesh C. Puthenkalathil, and Bernd Ensing, _J. Chem. Phys._ __149__ (2018), 072320, [DOI: 10.1063/1.5027392](https://doi.org/10.1063/1.5027392).

[3] Sequence dependence of transient Hoogsteen base pairing in DNA. Alberto Pérez de Alba Ortíz, Jocelyne Vreede, and Bernd Ensing,  _PLOS Computational Biology_ __18__ (2022), e1010113, [DOI: 10.1371/journal.pcbi.1010113](https://doi.org/10.1371/journal.pcbi.1010113).
