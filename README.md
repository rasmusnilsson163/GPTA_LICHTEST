# How to install:
## Normal installation
Download or clone the repository. After entering the repository, do
```
mkdir build
cd build
cmake .. -D OMP=ON -D XDR=ON
make -j
cd ../
``` 
The gpta executable is then found in `/build/gpta`. 
## Debug mode installation
If the user wishes to install GPTA in debug mode, do
```
mkdir build_debug
cd build_debug
cmake -D CMAKE_BUILD_TYPE=DEBUG ..
make -j
cd ../
```

# What is LICH-TEST?
LICH-TEST is an algorithm used to identify and categorize ice-like atoms. The algorithm was first published [here](https://pubs.acs.org/doi/10.1021/acs.jpcb.1c01926) (10.1021/acs.jpcb.1c01926), and has previously been implemented in MATLAB and Python, which have their own git-hub page [here](https://github.com/opakarin/lich-test) (https://github.com/opakarin/lich-test).

# How to use LICH-TEST in GPTA
The LICH-TEST action is called by using the `--lichtest` command. With the command one has to specify which atoms should be selected using any of the selection tools, such as `+s`. If no selection is specified the program will crash. 

The user can also specify the cutoff radius and minimal score (referred to *S*<sub>min</sub> in the [article](https://pubs.acs.org/doi/10.1021/acs.jpcb.1c01926)) using the `+rcut` and `+minscore` arguments. `+rcut` is given in Ångström and has a default value of 3.5 Å. `+minscore` is unitless and has a default value of 0.5. 

Here are some examples of using LICH-TEST in GPTA:
```
gpta.x --i coord.pdb traj.dcd --lichtest +s Ow,o +out lichtest.out
gpta.x --i coord.pdb --lichtest +s mW --o labeled_ice.pdb
gpta.x --i coord.pdb --lichtest +s mW +minscore 0.6 --o labeled_ice.pdb
gpta.x --i coord.pdb --lichtest +s o,oh,ohs +rcut 2.5 --o labeled_ice.pdb
```

# Bugs and Notes
## Testing has been limited
While the implementation has been tested for usage in systems we have used, it has not been extensively tested on a wide set of systems. Thus, the user should check if the results seem reasonable. 

## Makefile is redundant
LICH-TEST was not added to the regular Makefile, only the cMakefile. Thus the user should only use cmake to compile this version. Part of this is due to some issues getting the Makefile to compile GPTA.

## No MPI parallelization for LICH-TEST
While the LICH-TEST action is parallelized using OpenMP, no MPI parallelization was implemented. Merge requests that optimize the parallelization are welcome. 

# Citation
Please cite the [LICH-TEST article](https://pubs.acs.org/doi/10.1021/acs.jpcb.1c01926) if it has been used in your project. You can do it by e.g. adding the following to your .bib file
```
@article{Roudsari2021,
    author = {Roudsari, Golnaz and Veshki, Farshad G. and Reischl, Bernhard and Pakarinen, Olli H.},
    doi = {10.1021/acs.jpcb.1c01926},
    issn = {1520-6106},
    journal = {The Journal of Physical Chemistry B},
    month = {apr},
    number = {15},
    pages = {3909--3917},
    title = {{Liquid Water and Interfacial, Cubic, and Hexagonal Ice Classification through Eclipsed and Staggered Conformation Template Matching}},
    url = {https://pubs.acs.org/doi/10.1021/acs.jpcb.1c01926},
    volume = {125},
    year = {2021}
}
```

