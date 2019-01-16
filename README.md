## QCP IMPLEMENTATION AND DEMO IN FORTRAN 95/2003

### Subroutines and a sample program for implementing QCP for the sake of computing RMSD and superimpose molecules

#### Author: Veselin Kolev <vesso.kolev@gmail.com>
#### Released under BSD licence.
#### Version: 2019011500

#### Content:

#### 1. Introduction
#### 2. How to download the source code of the project
#### 3. Compiling of the source code
#### 4. Executing the demo
#### 5. Visualizing the result (the proximity)


_1. Introduction_

The included subroutines are part of the implementation of Quaternion Characteristic Polynomial (QCP) method, published by Naoto Hori at:

```
https://github.com/naotohori/fQCP
```

under BSD license. His code, in turn, is based on the C code written by Pu Liu and Douglas L. Theobald.

The Fortran code uploaded here is improved and partially revised by Veselin Kolev <vesso.kolev@gmail.com>.

To cite QCP:

 1. Douglas L. Theobald (2005) "Rapid calculation of RMSD using a quaternion-based characteristic
    polynomial." Acta Crystallographica A 61(4):478-480.

 2. Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2010) "Fast determination of the optimal
    rotational matrix for do_counted  superpositions." Journal of Computational Chemistry 31(7):1561-1563



_2. How to download the source code of the project_

The preferable method for obtaining the code is to use the tool ``git`` and clone the source tree locally to your file system:

```
git clone https://github.com/vessokolev/rmsd_f2003.git
```

You may download the source as a ZIP-archive by pressing the button "Clone or download" in the web-interface on GitHub, or by using wget:

```
wget https://github.com/vessokolev/rmsd_f2003/archive/master.zip
```


_3. Compiling of the source code_

To compile the code one need to have ``gfortran`` compiler. If you need to compile the code by using different compiler, set the environmental variable FC to point to the compiler executable file name. Then execute:

```
$ make
```

inside the folder containing the code.


_4. Executing the demo_

After successfully performing the compilation (see above) run:

```
$ ./demo
```

The computed RMSD (in Angstroms), numerically respresenting the proximity between the coordinates inside ``molecule_1.pdb`` and ``molecule_2.pdb``, will appear on the standard output, as a floating point number. Also, the file ``molecule_2_rotated.pdb`` will be created.


_5. Visualizing the result (the proximity)_

The input coordinates are stored inside the files ``molecule_1.pdb`` and ``molecule_2.pdb`` (following the descriptions recommended by the PDB file records specifications). During the execution of the demo program, the coordinates stored in ``molecule_2.pdb`` become a subject of sequential translations and rotations, in attempt to find the best match (proximity) between them and the coordinates in ``molecule_1.pdb``. The result (corresponding to the computed RSMD) is stored inside the file ``molecule_2_rotated.pdb``.

One possible way to visualize the proximity between the atoms in ``molecule_1.pdb`` and those in ``molecule_2_rotated.pdb`` is to load those two PDB files in VMD or UCSF Chimera.

For example:

```
$ chimera molecule_1.pdb molecule_2_rotated.pdb
```



