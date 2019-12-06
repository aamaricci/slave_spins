MODULE SS_VARS_GLOBAL
  USE SS_INPUT_VARS
  USE SS_SPARSE_MATRIX
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only: save_array,read_array,free_unit,free_units,file_length
  USE SF_LINALG,  only: kron,eigh,diag,diagonal,operator(.x.),ones,eye,outerprod
  USE SF_SPECIAL, only: fermi,heaviside,step
  USE SF_MISC,    only: assert_shape
  implicit none

  integer,save                              :: Ns       !Number of levels = 2*Norb
  integer                                   :: Nk
  integer                                   :: Nso
  integer                                   :: Ndim
  !
  real(8),dimension(:),allocatable          :: ss_lambda    !lambda_{m\sigma}
  real(8),dimension(:),allocatable          :: ss_lambda0   !lambda^0_{m\sigma}
  real(8),dimension(:),allocatable          :: ss_weiss     !4<S^x_{m\sigma}><E_fermion>
  real(8),dimension(:),allocatable          :: ss_c
  complex(8),dimension(:,:,:),allocatable   :: ss_Hk
  real(8),dimension(:,:,:),allocatable      :: ss_Wtk
  real(8),dimension(:,:),allocatable        :: ss_Hloc           !local hamiltonian
  !
  real(8),dimension(:),allocatable          :: ss_dens      ![Ns: 1:Norb_up, 1:Norb_dw]
  real(8),dimension(:),allocatable          :: ss_zeta      !zeta_{m\sigma}
  real(8),dimension(:),allocatable          :: ss_Sz  
  !
  real(8),dimension(:,:,:),allocatable      :: ss_SzSz
  real(8),dimension(:),allocatable          :: ss_Op  
  !
  real(8)                                   :: ss_Ef
  real(8)                                   :: zeta_function
  logical                                   :: is_dos=.false.



  !Spin problems sparse Hamiltonian and dense Eigensolution
  !so far we are not using quantum numbers (how bad...)
  !=========================================================  
  type(sparse_matrix_csr)                   :: spHs
  !
  real(8),allocatable,dimension(:,:),target :: ss_Evecs
  real(8),allocatable,dimension(:)          :: ss_Evals
  integer                                   :: ss_Ndegen


  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                         :: ss_file_suffix=""       !suffix string attached to the output files.


END MODULE SS_VARS_GLOBAL
