MODULE SS_MAIN
  USE SS_INPUT_VARS
  USE SS_SETUP
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_SPIN

  !
  USE SF_TIMER, only: start_timer,stop_timer
  USE SF_LINALG, only: diag,diagonal,kron
  USE SF_OPTIMIZE, only: broyden1,fsolve,broyden_mix
  USE SF_IOTOOLS,only: save_array
  USE SF_MISC,only: assert_shape
  !
  USE DMFT_TOOLS
  implicit none


  interface ss_init
     module procedure :: ss_init_hk
     module procedure :: ss_init_dos
  end interface ss_init

  public :: ss_init
  public :: ss_init_dos
  public :: ss_solve


  real(8),dimension(:),allocatable :: ss_lambda_init
  real(8),dimension(:),allocatable :: ss_zeta_init


contains




  !< massive allocation of all parameters
  subroutine ss_init_dos(Ebands,Dbands,Hloc)
    real(8),dimension(:,:)                      :: Ebands  ![Nlso,Ne]
    real(8),dimension(:,:)                      :: Dbands ![Nlso,Ne]
    real(8),dimension(:)                        :: Hloc
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Htmp
    real(8),dimension(Nspin*Norb,Nspin*Norb)    :: Wtmp
    integer                                     :: ie,io
    logical,save                                :: isetup=.true.
    !
    !< Init the SS structure + memory allocation
    Nk = size(Ebands,2)
    call assert_shape(Ebands,[Nspin*Norb,Nk],"ss_init_dos","Ebands")
    call assert_shape(Hloc,[Nspin*Norb],"ss_init_dos","Hloc")
    !
    if(isetup)call ss_setup_structure()
    !
    !< Guess/Read the lambda/zeta input
    call ss_init_params()
    isetup=.false.
    !
    ss_Hk = zero
    ss_Wtk= 0d0
    do ie=1,Nk
       Htmp=zero
       Wtmp=0d0
       do io=1,Nspin*Norb
          Htmp(io,io)  = one*Ebands(io,ie) - one*Hloc(io)
          Wtmp(io,io)  = Dbands(io,ie)
       end do
       select case(Nspin)
       case default
          ss_Hk(:,:,ie)  = kron(pauli_0,Htmp)
          ss_Wtk(:,:,ie) = kron(pauli_0,one*Wtmp)
       case (2)
          ss_Hk(:,:,ie)  = Htmp
          ss_Wtk(:,:,ie) = Wtmp
       end select
    end do
    !
    select case(Nspin)
    case default
       ss_Hloc = kron(pauli_0,one*diag(Hloc))
    case (2)
       ss_Hloc = diag(Hloc)
    end select
    !
    ss_Hdiag=.true.
    !
    allocate(ss_lambda_init(Ns));ss_lambda_init=ss_lambda
    allocate(ss_zeta_init(Ns))  ;ss_zeta_init  =ss_zeta
    !
    is_bethe=.true.
    !
  end subroutine ss_init_dos



  subroutine ss_init_hk(hk_user,wtk_user,hloc)
    complex(8),dimension(:,:,:)                          :: hk_user  ![Nlso,Nlso,Nk]
    real(8),dimension(size(hk_user,3))                   :: wtk_user ![Nk]
    complex(8),dimension(Nspin*Norb,Nspin*Norb),optional :: Hloc
    complex(8),dimension(Nspin*Norb,Nspin*Norb)          :: Htmp
    complex(8),dimension(:,:),allocatable                :: Hcheck
    logical,save                                         :: isetup=.true.
    integer :: ik,io
    !
    !< Init the SS structure + memory allocation
    Nk = size(hk_user,3)
    call assert_shape(hk_user,[Nspin*Norb,Nspin*Norb,Nk],"ss_init_hk","hk_user")
    if(isetup)call ss_setup_structure()
    !
    !< Guess/Read the lambda/zeta input
    call ss_init_params()
    isetup=.false.
    !
    Htmp    = sum(Hk_user,dim=3);where(abs(Htmp)<1d-6)Htmp=zero
    do ik=1,Nk
       select case(Nspin)
       case default
          ss_Hk(:,:,ik)  = kron(pauli_0,Hk_user(:,:,ik) - Htmp)
       case (2)
          ss_Hk(:,:,ik)  = Hk_user(:,:,ik) - Htmp
       end select
       ss_Wtk(:,:,ik) = Wtk_user(ik)
    end do
    !
    if(present(Hloc))then
       if(sum(abs(Htmp))>1d-12)stop "ss_init: Hloc seems to be present twice: in _Hloc and _Hk"
       select case(Nspin)
       case default
          ss_Hloc = kron(pauli_0,Hloc)
       case (2)
          ss_Hloc = Hloc           !local part
       end select
    else
       select case(Nspin)
       case default
          ss_Hloc = kron(pauli_0,Htmp)
       case (2)
          ss_Hloc = Htmp
       end select
    endif
    !
    ! ss_Hdiag=.true.
    ! allocate(Hcheck(Ns,Ns))
    ! do ik=1,Nk
    !    Hcheck = ss_Hk(:,:,ik) + ss_Hloc
    !    if(sum(abs(Hcheck - diag(diagonal(Hcheck)))) > 1d-6)ss_Hdiag=.false.
    ! enddo
    ! deallocate(Hcheck)
    !
    ! print*,ss_Hdiag
    call ss_get_lambda0()
    !
    allocate(ss_lambda_init(Ns));ss_lambda_init=ss_lambda
    allocate(ss_zeta_init(Ns))  ;ss_zeta_init  =ss_zeta
    !
  end subroutine ss_init_hk







  subroutine ss_solve()
    real(8),dimension(Nspin*Norb) :: lambda
    !
    select case(Nspin)
    case (1)
       lambda = ss_lambda(:Norb)
    case (2)
       lambda = ss_lambda
    case default
       stop "ss_solve ERROR: Nspin>2"
    end select
    !
    select case(solve_method)
    case ("broyden")       
       call broyden1(ss_solve_function,lambda)
    case ("hybrd")
       call fsolve(ss_solve_function,lambda,tol=1d-6)
    case default
       stop "ss_solve ERROR: solve_method not supported"
    end select
    !
    open(100,file="ss_zeta.dat")
    write(100,*)uloc(1),ss_zeta
    close(100)
    !
    open(100,file="ss_lambda.dat")
    write(100,*)uloc(1),ss_lambda
    close(100)
    !
    call save_array(trim(Pfile)//trim(ss_file_suffix)//".restart",[ss_lambda,ss_zeta])
    !
  contains
    !
    !< solve the SS problem and optmize ss_lambda using Broyden
    function ss_solve_function(lambda) result(fss)
      real(8),dimension(:),intent(in) :: lambda
      real(8),dimension(size(lambda)) :: fss
      !
      integer                         :: iter,Nsuccess=0
      logical                         :: z_converged
      real(8)                         :: zeta(Ns),Fzeta(Ns)
      !
      select case(Nspin)
      case (1)
         ss_lambda(:Norb) = lambda 
         call ss_spin_symmetry(ss_lambda)
      case (2)
         ss_lambda = lambda
      case default
         stop "ss_solve_function ERROR: Nspin>2"
      end select
      !
      z_converged=.false. ; iter=0
      !
      ss_zeta = ss_zeta_init
      !
      if(verbose>2)call start_timer()
      do while(.not.z_converged.AND.iter<=zeta_Nitermax)
         iter=iter+1
         call start_loop(iter,zeta_Nitermax,"Z-iter")
         zeta = ss_zeta
         call ss_solve_fermions
         call ss_solve_spins
         Fzeta= ss_zeta
         !
         call broyden_mix(zeta,Fzeta,zeta_Wmix,5,iter)
         !
         z_converged = check_convergence(ss_zeta,zeta_tolerance,Nsuccess,zeta_Nitermax)
         call end_loop()
      end do
      !
      !<constraint:
      fss = ss_Dens - (ss_Sz + 0.5d0)
      if(verbose>1)then
         write(*,"(A7,12G18.9)")"Dens  =",ss_dens
         write(*,"(A7,12G18.9)")"Sz+1/2=",ss_Sz+0.5d0
         write(*,*)""
         write(*,"(A7,12G18.9)")"Lambda=",lambda
         write(*,"(A7,12G18.9)")"F_ss  =",fss
      endif
      if(verbose>2)call stop_timer()
      if(verbose>1)then
         write(*,*)""
         write(*,*)""
      endif
    end function ss_solve_function
    !
  end subroutine ss_solve




END MODULE SS_MAIN
