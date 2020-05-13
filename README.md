# Slave Spins

A flexbile slave spins mean-field code, with support for multi-band structure and inequivalent sites. The code structure is divided in two parts 

* A slave spin *library*, i.e. a set of tree-like connected modules rooted by a single module `SLAVE_SPINS` exposing to the user the necessary procedures. 
* A model dependent *driver*, to be placed in the directory `drivers`. The driver is the main program to be compiled and imports from  `SLAVE_SPINS` the necessary methods. An example, solving the Hubbard model the multi-orbital Bethe lattice case, is contained in the file `drivers/ss_hm_bethe.f90`. 

**The default driver is:** `ss_DFT`, which interfaces to conventional Wannier90 output and requires a number of input files. 

* `hij.conf`: the Wannier90 real-space Hamiltonian file in its original form
* `lat.conf`: a file containing the real-space lattice basis vectors $[e_1,e_2,e_3]$ in coordinates with respect to default orthonormal basis, i.e. $e_i=a_x\hat{x}+a_y\hat{y}+a_z\hat{z}$. 
* `kpath.conf`: a file containing a list of k-space points for band-structure visualisation. Each entry of the list has the form: `"Point_Symbol",[P_x,P_y,P_z]`, where `Point_Symbol` is a canonical high symmetry point denomination, e.g. $\Gamma$, $M$, etc.. and `P_x,P_y,P_z` are the coordinates of the points in a $[0:1]$ basis, i.e. $M=[0.5,0.5,0.0]$. 



The code relies on the following two Fortran libraries, which must be installed and available on your system:  

* SciFortran [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran)  
* DMFT_Tools [https://github.com/aamaricci/DMFTtools](https://github.com/aamaricci/DMFTtools) 

The first provides the necessary support for the mathematical procedures used throughout the SS code. The second library is only used in some drivers to simplify some specific tasks, e.g. reading W90 input, get electronic bandstructure, get Fermi Surface, get Green's functions, etc. 



#### Installation

Clone the repository
`$ git clone https://github.com/aamaricci/slave_spins`

Compile using `cmake`

```bash
$ cd SlaveSpins
$ mkdir Build
$ cd Build
$ cmake .. 
$ make
```

At line 4 you can specify multiple options  (default in square brackets):

```cmake
-DEXE=<name of the driver without extension> [ss_hm_bethe]
-DPREFIX=<OS location where to put your compiled exe> [~/.bin]
-DUSE_MPI=<TRUE/FALSE> [TRUE]
-DVERBOSE=<TRUE/FALSE> [FALSE]
-DBUILD_TYPE=<RELEASE/DEBUG/TESTING> [RELEASE]
```

A suffix corresponding to your actualy GIT branch, if different from `master` , is appended to the generated executable. 

*A first run in a empty directory will generate a default input file*

Look into the default driver for more info.



#### Details on the implementation

*to do...*



--

***COPYRIGHT & LICENSING***  
Copyright 2019 -  (c), Adriano Amaricci, Laura Fanfarillo, Massimo Capone, Luca de'Medici  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

