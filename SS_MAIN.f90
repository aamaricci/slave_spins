MODULE SS_MAIN
  USE SS_INPUT_VARS
  USE SS_SETUP
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_SPIN
  !
  USE SF_LINALG, only: diag,diagonal
  !
  USE DMFT_TOOLS
  implicit none
  public :: ss_init
  public :: ss_solve

contains




  !< massive allocation of all parameters
  subroutine ss_init(hk_user,wtk_user,hloc)
    logical,save                                         :: isetup=.true.
    complex(8),dimension(:,:,:)                          :: hk_user  ![Nlso,Nlso,Nk]
    real(8),dimension(size(hk_user,3))                   :: wtk_user ![Nk]
    complex(8),dimension(Nspin*Norb,Nspin*Norb)          :: Htmp
    complex(8),dimension(Nspin*Norb,Nspin*Norb),optional :: Hloc
    integer :: ik
    !
    !< Init the SS structure + memory allocation
    Nk = size(hk_user,3)
    if(isetup)call ss_setup_structure()
    !
    !< Guess/Read the lambda/zeta input
    call ss_init_params()
    isetup=.false.
    !
    Htmp    = sum(Hk_user,dim=3);where(abs(Htmp)<1d-6)Htmp=zero
    forall(ik=1:Nk)&
         ss_Hk(:,:,ik) = Hk_user(:,:,ik) - Htmp    
    ss_Wtk  = Wtk_user
    !
    if(present(Hloc))then
       if(sum(abs(Htmp))>1d-12)stop "ss_init: Hloc seems to be present twice: in _Hloc and _Hk"
       ss_Hloc = Hloc           !local part
    else
       ss_Hloc = Htmp
    endif
    !
    ss_Hdiag=.true.
    do ik=1,Nk
       Htmp = ss_Hk(:,:,ik) + ss_Hloc
       if(sum(abs(Htmp - diag(diagonal(Htmp)))) > 1d-6)ss_Hdiag=.false.
    enddo
    !
  end subroutine ss_init




  subroutine ss_solve()
    !<< optimize lambda + constraint
    !
    ! < call ss_fermion
    call ss_solve_fermions
    !
    !< call ss_pspins
    call ss_solve_spins
    !
  end subroutine ss_solve




END MODULE SS_MAIN
