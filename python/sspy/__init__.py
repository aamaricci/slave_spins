import numpy as np
import ss2py                         #import the Fortran module
import sspy                          #import the PYTHON library itself
from ss2py import *
from sspy import *

global Nlat         #Nlat =# of ineq sites
global Norb         #Norb =# of impurity orbitals
global Nspin        #Nspin=# spin degeneracy (max 2)
global filling      #filling
global Nloop        #max dmft loop variables
global xmu          #chemical potential
global beta         #inverse temperature
global eps          #broadening
global wini,wfin    #frequency range
global Nsuccess     #Number of repeated success to fall below convergence threshold  
global Lmats        #Number of Matsuabra Freq.
global Lreal        #Number of real-axis freq.


#FROM SS_INPUT_VARS
Nlat       = ss2py.ss_input_vars.nlat
Norb       = ss2py.ss_input_vars.norb
Nspin      = ss2py.ss_input_vars.nspin
filling    = ss2py.ss_input_vars.filling
Nloop      = ss2py.ss_input_vars.nloop
xmu        = ss2py.ss_input_vars.xmu
beta       = ss2py.ss_input_vars.beta
eps        = ss2py.ss_input_vars.eps
wini       = ss2py.ss_input_vars.wini
wfin       = ss2py.ss_input_vars.wfin
Nsuccess   = ss2py.ss_input_vars.nsuccess
Lmats      = ss2py.ss_input_vars.lmats
Lreal      = ss2py.ss_input_vars.lreal


#FROM SS_MAIN:
def solve(hk,dbands,hloc,ineq_sites):
    DimArg=len(np.shape(hk))
    if(DimArg==3):
        ss2py.solve_hk(hk,ineq_sites)
    elif(DimArg==2):
        ss2py.solve_dos(hk,dbands,hloc,ineq_sites)
    else:
        raise ValueError('Dim(Hk) /= {[Nlso,Nlso,Nk],[Nlso,Ne]}')
    return ;






