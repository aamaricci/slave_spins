# Slave Spins

A flexbile slave spins mean-field code, with support for multi-band structure and inequivalent sites.   

The code is based on:  

* SciFortran [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran)  
* DMFT_Tools [https://github.com/aamaricci/DMFTtools](https://github.com/aamaricci/DMFTtools) [For the driver part, see below]

The code structure is as follow:  

* The set of modules compile into a top layer named `SLAVE_SPINS.f90`  
* The actual implementation of the model dependent  is case by case performed in a driver program, usually placed in the directory `drivers`, which must include the `SLAVE_SPINS` module. 

An example, solving the Hubbard model the multi-orbital Bethe lattice case, is contained in the file `drivers/ss_hm_bethe.f90`.

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

