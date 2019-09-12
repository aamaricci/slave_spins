MODULE SS_VARS_GLOBAL
  USE SF_CONSTANTS
  implicit none






  integer,save                            :: Ns       !Number of levels = 2*Norb
  integer                                 :: Nk
  !
  real(8),dimension(:),allocatable        :: ss_lambda    !lambda_{m\sigma}
  real(8),dimension(:),allocatable        :: ss_lambda0   !lambda^0_{m\sigma}
  real(8),dimension(:),allocatable        :: ss_zeta      !zeta_{m\sigma}
  real(8),dimension(:),allocatable        :: ss_weiss     !4<S^x_{m\sigma}><E_fermion>
  real(8),dimension(:),allocatable        :: ss_c
  real(8),dimension(:),allocatable        :: ss_Sz
  complex(8),dimension(:,:,:),allocatable :: ss_Hk
  real(8),dimension(:),allocatable        :: ss_Wtk
  real(8),dimension(:,:),allocatable      :: ss_Hloc           !local hamiltonian
  real(8),dimension(:),allocatable        :: ss_dens ![Ns: 1:Norb_up, 1:Norb_dw]

  logical                                 :: ss_Hdiag !
  !
  real(8)                                 :: zeta_function
  logical                                 :: is_bethe=.false.

  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                      :: ss_file_suffix=""       !suffix string attached to the output files.
  character(len=200)                     :: ss_input_file=""


END MODULE SS_VARS_GLOBAL
