MODULE SS_MAIN
  USE SS_INPUT_VARS
  USE SS_SETUP
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_SPIN
  
  !
  USE SF_LINALG, only: diag,diagonal
  USE SF_OPTIMIZE, only: broyden1
  USE SF_IOTOOLS,only: save_array
  !
  USE DMFT_TOOLS
  implicit none
  public :: ss_init
  public :: ss_solve


  real(8),dimension(:),allocatable :: ss_lambda_init
  real(8),dimension(:),allocatable :: ss_zeta_init
  
  
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
    allocate(ss_lambda_init(Ns));ss_lambda_init=ss_lambda
    allocate(ss_zeta_init(Ns))  ;ss_zeta_init  =ss_zeta
  end subroutine ss_init




  subroutine ss_solve()


    call broyden1(broyden_ss_solve,ss_lambda)

    open(100,file="ss_zeta.dat")
    write(100,*)uloc(1),ss_zeta
    close(100)

    open(100,file="ss_lambda.dat")
    write(100,*)uloc(1),ss_lambda
    close(100)

    call save_array(trim(Pfile)//trim(ss_file_suffix)//".restart",[ss_lambda,ss_zeta])

  contains
    !
    !< solve the SS problem and optmize ss_lambda using Broyden
    function broyden_ss_solve(lambda) result(fss)
      real(8),dimension(:),intent(in) :: lambda
      real(8),dimension(size(lambda)) :: fss
      !
      integer                         :: iter
      logical                         :: z_converged
      integer                         :: Nitermax=100,Nsuccess=0
      real(8) :: tol=1d-6

      if(size(lambda)/=Ns)stop "broyden_ss_solver ERROR: size(lambda)!=Ns"

      ss_lambda=lambda
      !ss_zeta = ss_zeta_init
      z_converged=.false. ; iter=0

      do while(.not.z_converged.AND.iter<=Nitermax)
         iter=iter+1
         call start_loop(iter,Nitermax,"Z-iter")

         call ss_solve_fermions

         call ss_solve_spins

         z_converged = check_convergence_local(ss_zeta,tol,Nsuccess,Nitermax)

         call end_loop()
      end do

      !<constraint:
      fss = ss_Dens - (ss_Sz + 0.5d0)
      print*,ss_dens,ss_Sz+0.5d0
      print*,lambda,fss
    end function broyden_ss_solve


  end subroutine ss_solve




END MODULE SS_MAIN
