# What GPTA_LCIH-TEST?
This project is an implementation of the LICH-TEST algorithm in the program GPTA. 

LICH-TEST is an algorithm used to identify and categorize ice-like atoms. The algorithm was first published [here](https://pubs.acs.org/doi/10.1021/acs.jpcb.1c01926) (10.1021/acs.jpcb.1c01926), and has previously been implemented in MATLAB and Python, which have their own git-hub page [here](https://github.com/opakarin/lich-test).

GPTA is shorthand for General Purpose Trajectory Analyzer. It is developed by Paolo Raiteri and has its own github page [here](https://github.com/praiteri/GPTA/tree/main). The program is used to make analysis and manipulation of structure files easy. 

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

# How to use LICH-TEST in GPTA
## Basics
The LICH-TEST action is called by using the `--lichtest` command. With the command one has to specify which atoms should be selected using any of the selection tools, such as `+s`. If no selection is specified the program will crash. 

The user can also specify the cutoff radius and minimal score (referred to *S*<sub>min</sub> in the [article](https://pubs.acs.org/doi/10.1021/acs.jpcb.1c01926)) using the `+rcut` and `+minscore` arguments. `+rcut` is given in Ångström and has a default value of 3.5 Å. `+minscore` is unitless and has a default value of 0.5. 

## Output
### Outfile
The LICH-TEST action will always write an output file with the ice counts. The default name is `lichtest.out`, but this can be changed by using the `+out` flag followed by the name of the output file. The columns and their description can be found in the following table.
| Column number | Column name   |            Description                        |
| ------------- | ------------- | --------------------------------------------  |
|            1  | Frame Number  | Frame which was analyzed                      | 
|            2  | Liquid        | Number of liquid atoms                        | 
|            3  | Cubic         | Number of cubic ice (I<sub>c</sub>) atoms     |
|            4  | Hexagonal     | Number of hexagonal ice (I<sub>h</sub>) atoms |
|            5  | Mixed         | Number of mixed-interfacial atoms             |
|            6  | Cubic-if.     | Number of cubic-interfacial atoms             |
|            7  | Hexagonal-if. | Number of hexagonal-interfacial atoms         |
|            8  | Clath. hydr.  | Number of clathrate hydrate atoms             | 
|            9  | Interfacial   | Number of other interfacial atoms             |
|           10  | Clathrate-if. | Number of clathrate-interfacial atoms         |

A more thorough definition for each types can be found in [the paper](https://pubs.acs.org/doi/10.1021/acs.jpcb.1c01926).

### Labeled trajectory
By using the `--o` action in combination with the `--lichtest` action, the user can output a labeled strucutre file or trajecotry. In the labeled strucutre file the names of the atoms are replaced with abbreviation of the classification types. Below a table of which name corresponds to which type. The output formats the user can chose from are somewhat limited, as not all trajecotry formats store the name of the atoms at in each frame. It is thus recomended to use the `.pdb` format as this has been tested.
| Abbreviated atom name | Classification type | 
| ------------------- | ------------------- |
| OL | Liquid |
| OC | Cibic ice |
| OH | Hexagonal ice |
| OM | Mixed-interfacial |
| OIC | Cubic-interfacial | 
| OIH | Hexagonal-interfacial | 
| OCH | Clathrate hydrate |
| OI  | Other interfacial |
| OICH | Clathrate-interfacial |

## Examples
Here are some examples of using LICH-TEST in GPTA:
```
gpta.x --i coord.pdb traj.dcd --lichtest +s Ow,o +out lichtest.out
gpta.x --i coord.pdb --lichtest +s mW --o labeled_ice.pdb
gpta.x --i coord.pdb --lichtest +s mW +minscore 0.6 --o labeled_ice.pdb
gpta.x --i coord.pdb --lichtest +s o,oh,ohs +rcut 2.5 --o labeled_ice.pdb
```

# Bugs and Notes
## Limited testing was done
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

