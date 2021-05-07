MODULE SS_VARS_GLOBAL
  USE SS_INPUT_VARS
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only: save_array,read_array,free_unit,free_units,file_length
  USE SF_LINALG
  USE SF_PAULI
  USE SF_SP_LINALG, only: sp_eigh
  USE SF_SPECIAL, only: fermi,heaviside,step
  USE SF_MISC,    only: assert_shape
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none

  integer,save                            :: Ns,Nss       !# of levels total = 2*Nlat*Norb
  integer                                 :: Nk
  integer                                 :: Nlso,Niso,Nso
  integer                                 :: Nineq
  integer,dimension(3)                    :: nDefOrder
  character(len=5),dimension(3)           :: DefOrder
  !
  complex(8),dimension(:,:,:),allocatable :: ss_Hk    !Inter-cell Hamiltionian, k-dependent
  complex(8),dimension(:,:),allocatable   :: ss_Hloc  !non-local, intra-cell hopping terms 
  complex(8),dimension(:,:),allocatable   :: ss_Hhyb  !local, intra-cell, hybridazation 
  real(8),dimension(:),allocatable        :: ss_Hdiag !local, intra-cell, diagonal terms 
  real(8),dimension(:,:),allocatable      :: ss_Wtk

  !
  real(8),dimension(:,:),allocatable      :: ss_c      ,ss_c_ineq       ![Nlat,Nso]
  real(8),dimension(:,:),allocatable      :: ss_dens   ,ss_dens_ineq    !..
  real(8),dimension(:,:),allocatable      :: ss_lambda0,ss_lambda0_ineq !..
  real(8),dimension(:,:),allocatable      :: ss_lambda ,ss_lambda_ineq  !..
  real(8),dimension(:,:),allocatable      :: ss_Sz     ,ss_Sz_ineq      !..
  real(8),dimension(:,:),allocatable      :: ss_Op     ,ss_Op_ineq      !..
  real(8),dimension(:,:),allocatable      :: ss_Heff   ,ss_Heff_ineq    !..
  real(8),dimension(:,:,:),allocatable    :: ss_OdgOp  ,ss_OdgOp_ineq   ![Nlat,Nso,Nso]
  real(8),dimension(:,:,:),allocatable    :: ss_Jhybr  ,ss_Jhybr_ineq   ![Nlat,Nso,Nso]
  real(8),dimension(:,:,:,:),allocatable  :: ss_SzSz   ,ss_SzSz_ineq    ![Nlat,4,Nso,Nso]
  !
  integer,dimension(:),allocatable        :: ss_ilat2ineq
  integer,dimension(:),allocatable        :: ss_ineq2ilat
  !
  real(8)                                 :: zeta_function
  logical                                 :: is_dos=.false.
  !
  character(len=32)                       :: ss_file_suffix=""
  !
  logical                                 :: master=.true.
  integer                                 :: mpi_rank=0
  integer                                 :: mpi_size=1


  interface ss_findloc
     module procedure :: ss_findloc_char
     module procedure :: ss_findloc_int
  end interface ss_findloc

contains


  function ss_findloc_char(array,val) result(pos)
    character(len=*),dimension(:) :: array
    character(len=*)              :: val
    integer                       :: pos,i
    pos=0
    do i=1,size(array)
       if(array(i)==val)then
          pos = i
          exit
       endif
    enddo
    return
  end function ss_findloc_char

  function ss_findloc_int(array,val) result(pos)
    integer,dimension(:) :: array
    integer              :: val
    integer              :: pos,i
    pos=0
    do i=1,size(array)
       if(array(i)==val)then
          pos = i
          exit
       endif
    enddo
    return
  end function ss_findloc_int

END MODULE SS_VARS_GLOBAL
