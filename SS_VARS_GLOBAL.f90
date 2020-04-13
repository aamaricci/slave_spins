MODULE SS_VARS_GLOBAL
  USE SS_INPUT_VARS
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only: save_array,read_array,free_unit,free_units,file_length
  USE SF_LINALG,  only: kron,eigh,diag,diagonal,operator(.x.),ones,eye,outerprod
  USE SF_SP_LINALG, only: sp_eigh
  USE SF_SPECIAL, only: fermi,heaviside,step
  USE SF_MISC,    only: assert_shape
  implicit none


  integer,save                            :: Ns,Nss       !# of levels total = 2*Nlat*Norb
  integer                                 :: Nk
  integer                                 :: Nlso,Niso
  integer                                 :: Nineq
  integer,dimension(3)                    :: nDefOrder
  character(len=5),dimension(3)           :: DefOrder
  !
  complex(8),dimension(:,:,:),allocatable :: ss_Hk
  real(8),dimension(:,:,:),allocatable    :: ss_Wtk
  complex(8),dimension(:,:),allocatable   :: ss_Hloc           !local hamiltonian
  !
  real(8),dimension(:,:),allocatable      :: ss_c      ,ss_c_ineq
  real(8),dimension(:,:),allocatable      :: ss_dens   ,ss_dens_ineq
  real(8),dimension(:,:),allocatable      :: ss_lambda0,ss_lambda0_ineq 
  real(8),dimension(:,:),allocatable      :: ss_lambda ,ss_lambda_ineq 
  real(8),dimension(:,:),allocatable      :: ss_Sz     ,ss_Sz_ineq
  real(8),dimension(:,:,:,:),allocatable  :: ss_SzSz   ,ss_SzSz_ineq
  real(8),dimension(:,:),allocatable      :: ss_Op     ,ss_Op_ineq  
  complex(8),dimension(:,:),allocatable   :: ss_weiss  ,ss_weiss_ineq   !4<S^x_{m\sigma}><E_el>  

  !
  integer,dimension(:),allocatable        :: ss_ilat2ineq
  integer,dimension(:),allocatable        :: ss_ineq2ilat
  !
  real(8)                                 :: zeta_function
  logical                                 :: is_dos=.false.
  !
  character(len=32)                       :: ss_file_suffix=""







END MODULE SS_VARS_GLOBAL
