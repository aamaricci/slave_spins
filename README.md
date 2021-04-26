# Slave Spins

A flexbile slave spins mean-field code, with support for multi-band structure and inequivalent sites. This code uses MPI for faster execution. 

### Dependencies

The code is written around the SciFortran library, which can be found here with installation notes.   

* SciFortran https://github.com/QcmPlab/SciFortran

### Installation

Installation is  available using CMake.    

Clone the repo:

`git clone https://github.com/QcmPlab/Slave-Spins Slave_Spins`

and from the just created directory make a standard out-of-source CMake compilation:

`mkdir build`  
 `cd build`  
`cmake ..`     
`make`     
`make install`   

The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/slave_spins/PLAT/[VERSION]>` 
* `-DUSE_MPI=<yes>/no`  
* `-DVERBOSE=yes/<no> `  
* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  

### System loading:

The library can be loaded into the operative system using one of the following, automatically generated, methods:    

* environment module file `~/.modules.d/slave_spins/<PLAT>`  
* homebrew `bash` script `<PREFIX>/bin/configvars.sh`
* pkg-config file in `~/.pkg-config.d/slave_spins.pc`

### Python binding

**The Python binding is still in beta version**
Python binding (API) through module `sspy` can be  installed, once the library is successfully loaded in the OS, using the conventional toolchain:

`export F90=mpif90` (required if library has been compiled and installed with MPI support)  

1. `python setup.py install`
2. `pip install .`

Method 2. has the advantage of making `uninstall` operation feasible. 

### Uninstall

The library is removed with the command:

`make uninstall`

from the same building directory as for the installation part. 

### Driver program

**A default driver** named `ss_DFT` is available in the driver database:

`https://github.com/QcmPlab/Driver-Database`

This driver provides a generic  interface to conventional Wannier90 output and requires a number of additional input files. 

* `hij.conf`: the Wannier90 real-space Hamiltonian file in its original form
* `lat.conf`: a file containing the real-space lattice basis vectors $[e_1,e_2,e_3]$ in coordinates with respect to default orthonormal basis, i.e. $e_i=a_x\hat{x}+a_y\hat{y}+a_z\hat{z}$. 
* `kpath.conf`: a file containing a list of k-space points for band-structure visualisation. Each entry of the list has the form: `"Point_Symbol",[P_x,P_y,P_z]`, where `Point_Symbol` is a canonical high symmetry point denomination, e.g. $\Gamma$, $M$, etc.. and `P_x,P_y,P_z` are the coordinates of the points in a $[0:1]$ basis, i.e. $M=[0.5,0.5,0.0]$. 

*A first run in a empty directory will generate a default input file*

Look into the default driver for more info.



***COPYRIGHT & LICENSING***  
Copyright 2019 -  (c), Adriano Amaricci, Laura Fanfarillo, Massimo Capone, Luca de'Medici  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

<!--This program is free software: you can redistribute it and/or modify-->
<!--it under the terms of the GNU General Public License as published by-->
<!--the Free Software Foundation, either version 3 of the License, or-->
<!--(at your option) any later version.-->

<!--You should have received a copy of the GNU General Public License-->
<!--along with this program.  If not, see <http://www.gnu.org/licenses/>.-->

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. 

