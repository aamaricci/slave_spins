MODULE SS_MAIN
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_SPIN
  !
  USE SF_TIMER, only: start_timer,stop_timer
  USE SF_LINALG, only: diag,diagonal,kron
  USE SF_OPTIMIZE, only: broyden1,fsolve,broyden_mix,linear_mix,adaptive_mix
  USE SF_IOTOOLS,only: save_array
  USE SF_MISC,only: assert_shape
  !
  USE DMFT_TOOLS
  implicit none
  private

  interface ss_init
     module procedure :: ss_init_hk
     ! module procedure :: ss_init_dos
  end interface ss_init

  public :: ss_init
  public :: ss_solve


  real(8),dimension(:),allocatable :: ss_lambda_init
  real(8),dimension(:),allocatable :: ss_zeta_init

  integer,save             :: siter=0
  integer                  :: fiter,info
  logical                  :: fconverged

contains






  subroutine ss_init_hk(hk_user,wtk_user,hloc)
    complex(8),dimension(:,:,:)            :: hk_user  ![Nlso,Nlso,Nk]
    real(8),dimension(size(hk_user,3))     :: wtk_user ![Nk]
    complex(8),dimension(Nspin*Norb,Nspin*Norb),optional :: Hloc
    complex(8),dimension(Nspin*Norb,Nspin*Norb)          :: Htmp
    complex(8),dimension(:,:),allocatable  :: Hcheck
    logical,save                           :: isetup=.true.
    integer                                :: ik,io
    !
    !< Init the SS structure + memory allocation
    Nk = size(hk_user,3)
    call assert_shape(hk_user,[Nspin*Norb,Nspin*Norb,Nk],"ss_init_hk","hk_user")
    if(isetup)call ss_setup_structure()
    !
    !< Init the Hk structures
    Htmp    = sum(Hk_user,dim=3);where(abs(Htmp)<1d-6)Htmp=zero
    if(present(Hloc).AND.(sum(abs(Htmp))>1d-12))stop "ss_init: Hloc seems to be present twice: in _Hloc and _Hk"
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
    !< Init local non-interacting part 
    if(present(Hloc))Htmp = Hloc !Htmp should necessarily be zero (if condition above)
    select case(Nspin)
    case default
       ss_Hloc = kron(pauli_0,Htmp)
    case (2)
       ss_Hloc = Htmp
    end select
    !
    ss_Hdiag=.true.
    allocate(Hcheck(Ns,Ns))
    do ik=1,Nk
       Hcheck = ss_Hk(:,:,ik) + ss_Hloc
       if(sum(abs(Hcheck - diag(diagonal(Hcheck)))) > 1d-6)ss_Hdiag=.false.
    enddo
    deallocate(Hcheck)
    !
    !< Init/Read the lambda/zeta input
    if(filling/=dble(Norb))call ss_get_lambda0()
    call ss_init_params()
    if(verbose>2)then
       write(*,"(A7,12G18.9)")"Lam0  =",ss_lambda0
       write(*,"(A7,12G18.9)")"Lam   =",ss_lambda
       write(*,"(A7,12G18.9)")"Z     =",ss_zeta
    endif
    !
    !< Internal use:
    allocate(ss_lambda_init(Ns));ss_lambda_init=ss_lambda
    allocate(ss_zeta_init(Ns))  ;ss_zeta_init  =ss_zeta
    !
    isetup=.false.
    return
  end subroutine ss_init_hk




  subroutine ss_solve()
    real(8),dimension(Nso)   :: lambda
    real(8),dimension(2*Nso) :: params
    real(8),dimension(Nso)   :: Xvec,Fvec
    !
    select case(solve_method)
    case ("broyden")
       params =  [ss_lambda(:Nso),ss_zeta(:Nso)]
       call broyden1(ss_solve_full,params,tolf=solve_tolerance)
       !
    case ("fsolve")
       params =  [ss_lambda(:Nso),ss_zeta(:Nso)]
       call fsolve(ss_solve_full,params,tol=solve_tolerance)
       !
    case ("gg_broyden")
       lambda = ss_lambda(:Nso)
       call broyden1(ss_solve_lambda,lambda,tolf=solve_tolerance)
       !
    case ("gg_fsolve")
       lambda = ss_lambda(:Nso)
       call fsolve(ss_solve_lambda,lambda,tol=solve_tolerance)
       !
    case ("lf_solve")       
       include "ss_main_lf_solve.h90"
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
    open(100,file="ss_dens.dat")
    write(100,*)uloc(1),ss_dens
    close(100)
    !
    call save_array(trim(Pfile)//trim(ss_file_suffix)//".restart",[ss_lambda,ss_zeta])
    !
  end subroutine ss_solve




  function ss_solve_full(aparams) result(fss)
    real(8),dimension(:),intent(in)  :: aparams !2*Nso
    real(8),dimension(size(aparams)) :: fss
    real(8)                          :: zeta(Nso)
    !
    ss_lambda(1:Nso) = aparams(1:Nso)
    ss_zeta(1:Nso)   = aparams(Nso+1:2*Nso)  
    !
    zeta = ss_zeta(1:Nso)
    call ss_solve_fermions
    call ss_solve_spins
    if(verbose>3)write(*,"(A6,12G18.9)")"Ef   =",xmu
    if(verbose>3)write(*,"(A6,12G18.9)")"C    =",ss_c
    !
    !<constraint:
    fss(1:Nso) = ss_Dens(1:Nso) - (ss_Sz(1:Nso) + 0.5d0)
    fss(Nso+1:)= ss_zeta(1:Nso) - zeta
    !
    if(verbose>1)then
       write(*,"(A7,12G18.9)")"N     =",ss_dens(:Nso),sum(ss_dens),filling
       write(*,"(A7,12G18.9)")"Lambda=",ss_lambda(:Nso)
       write(*,"(A7,12G18.9)")"Z_ss  =",ss_zeta(:Nso)      
       write(*,"(A7,12G18.9)")"F_ss  =",fss
       write(*,*)""
    endif
    open(100,file="ss_zeta.dat")
    write(100,*)uloc(1),ss_zeta
    close(100)
    !
    open(100,file="ss_lambda.dat")
    write(100,*)uloc(1),ss_lambda
    close(100)
    !
    open(100,file="ss_dens.dat")
    write(100,*)uloc(1),ss_dens
    close(100)
  end function ss_solve_full


  !< solve the SS problem and optmize ss_lambda using Broyden/Hybrd optimization
  function ss_solve_lambda(lambda) result(fss)
    real(8),dimension(:),intent(in) :: lambda
    real(8),dimension(size(lambda)) :: fss
    integer                         :: iter,Nsuccess=0
    logical                         :: z_converged
    real(8)                         :: zeta(Ns),Fzeta(Ns)
    !
    ss_lambda(:Nso) = lambda 
    if(Nspin==1)call ss_spin_symmetry(ss_lambda)
    !
    if(zeta_restart_init)ss_zeta = ss_zeta_init
    !
    if(verbose>2)call start_timer()
    z_converged=.false. ; iter=0
    do while(.not.z_converged.AND.iter<=loop_Nitermax)
       iter=iter+1
       call start_loop(iter,loop_Nitermax,"Z-iter")
       zeta = ss_zeta
       call ss_solve_fermions
       if(verbose>3)write(*,"(A6,12G18.9)")"Ef   =",xmu
       if(verbose>3)write(*,"(A6,12G18.9)")"C    =",ss_c
       if(verbose>2)write(*,"(A6,12G18.9)")"N    =",ss_dens(:Nso),sum(ss_dens),filling
       call ss_solve_spins
       if(verbose>2)write(*,"(A6,12G18.9)")"Z_ss =",ss_zeta(:Nso)
       Fzeta= ss_zeta - zeta
       !
       select case(loop_MixType)
       case ("linear")            
          call linear_mix(ss_zeta,Fzeta,loop_Wmix)
       case ("adaptive")
          call adaptive_mix(ss_zeta,Fzeta,loop_Wmix,iter)
       case default
          call broyden_mix(zeta,Fzeta,loop_Wmix,loop_Nmix,iter)
       end select
       !
       z_converged = check_convergence_local(ss_zeta-zeta,loop_tolerance,Nsuccess,loop_Nitermax)
       call end_loop() 
    end do
    !
    !<constraint:
    fss = ss_Dens - (ss_Sz + 0.5d0)
    !
    open(100,file="ss_zeta.dat")
    write(100,*)uloc(1),ss_zeta
    close(100)
    !
    open(100,file="ss_lambda.dat")
    write(100,*)uloc(1),ss_lambda
    close(100)
    !
    open(100,file="ss_dens.dat")
    write(100,*)uloc(1),ss_dens
    close(100)
    !
    if(verbose>1)then
       write(*,"(A7,12G18.9)")"Dens  =",ss_dens(:Nso)
       write(*,"(A7,12G18.9)")"Sz+1/2=",ss_Sz(:Nso)+0.5d0
       write(*,*)""
       write(*,"(A7,12G18.9)")"Lambda=",lambda(:Nso)
       write(*,"(A7,12G18.9)")"F_ss  =",fss(:Nso)
       if(verbose>2)call stop_timer()
       write(*,*)""
       write(*,*)""
    endif
  end function ss_solve_lambda


  function solve_Hps(lambda) result(Fss)
    real(8),dimension(:),intent(in) :: lambda
    real(8),dimension(size(lambda)) :: fss
    !
    siter=siter+1
    ss_lambda(1:Nso) = lambda
    if(Nspin==1)call ss_spin_symmetry(ss_lambda)
    !
    call ss_solve_spins
    !
    fss = (ss_Dens(:Nso) - (ss_Sz(:Nso) + 0.5d0))
    if(verbose>1)write(*,"(I18,20G16.7)")siter,Fss(:Nso),ss_lambda(:Nso)
  end function solve_Hps



END MODULE SS_MAIN




! !< massive allocation of all parameters
! subroutine ss_init_dos(Ebands,Dbands,Hloc)
!   real(8),dimension(:,:)                      :: Ebands  ![Nlso,Ne]
!   real(8),dimension(:,:)                      :: Dbands ![Nlso,Ne]
!   real(8),dimension(:)                        :: Hloc
!   complex(8),dimension(Nso,Nso) :: Htmp
!   real(8),dimension(Nso,Nso)    :: Wtmp
!   integer                                     :: ie,io
!   logical,save                                :: isetup=.true.
!   !
!   !< Init the SS structure + memory allocation
!   Nk = size(Ebands,2)
!   call assert_shape(Ebands,[Nso,Nk],"ss_init_dos","Ebands")
!   call assert_shape(Hloc,[Nso],"ss_init_dos","Hloc")
!   !
!   if(isetup)call ss_setup_structure()
!   !
!   !< Guess/Read the lambda/zeta input
!   call ss_init_params()
!   isetup=.false.
!   !
!   ss_Hk = zero
!   ss_Wtk= 0d0
!   do ie=1,Nk
!      Htmp=zero
!      Wtmp=0d0
!      do io=1,Nso
!         Htmp(io,io)  = one*Ebands(io,ie) - one*Hloc(io)
!         Wtmp(io,io)  = Dbands(io,ie)
!      end do
!      select case(Nspin)
!      case default
!         ss_Hk(:,:,ie)  = kron(pauli_0,Htmp)
!         ss_Wtk(:,:,ie) = kron(pauli_0,one*Wtmp)
!      case (2)
!         ss_Hk(:,:,ie)  = Htmp
!         ss_Wtk(:,:,ie) = Wtmp
!      end select
!   end do
!   !
!   select case(Nspin)
!   case default
!      ss_Hloc = kron(pauli_0,one*diag(Hloc))
!   case (2)
!      ss_Hloc = diag(Hloc)
!   end select
!   !
!   ss_Hdiag=.true.
!   !
!   allocate(ss_lambda_init(Ns));ss_lambda_init=ss_lambda
!   allocate(ss_zeta_init(Ns))  ;ss_zeta_init  =ss_zeta
!   !
!   is_bethe=.true.
!   !
! end subroutine ss_init_dos



! subroutine ss_solve_routine() 
!   real(8),dimension(Nso) :: lambda
!   real(8),dimension(Nso) :: Xvec,Fvec
!   integer                       :: fiter,info
!   logical                       :: fconverged
!   !
!   !< Setup:
!   if(Nspin==1)call ss_spin_symmetry(ss_lambda)
!   !
!   if(verbose>2)call start_timer()
!   fconverged = .false.; fiter=0 
!   do while(.not.fconverged.AND.fiter<=loop_Nitermax)
!      fiter=fiter+1
!      call start_loop(fiter,loop_Nitermax,"Fermion-iter")
!      !
!      !
!      lambda = ss_lambda(:Nso)    !get a guess to start the Css: n=Sz+1/2 optimization
!      call ss_solve_fermions
!      if(verbose>3)write(*,"(A6,12G18.9)")"Ef   =",xmu
!      if(verbose>3)write(*,"(A6,12G18.9)")"C    =",ss_c
!      if(verbose>2)write(*,"(A6,12G18.9)")"N    =",ss_dens(:Nso),sum(ss_dens),filling
!      if(verbose>1)write(*,"(3A18)")"Iter","Fss","Lambda"
!      !
!      siter  = 0
!      Xvec   = ss_zeta(:Nso)
!      call broyden1(solve_Hps,lambda,tolf=solve_tolerance)
!      !call fsolve(solve_Hps,lambda,tol=solve_tolerance,info=info)
!      Fvec   = ss_zeta(:Nso)! - Xvec
!      !
!      ! call linear_mix(Xvec,Fvec,loop_Wmix)
!      !call adaptive_mix(Xvec,Fvec,loop_Wmix,fiter)
!      !call broyden_mix(Xvec,Fvec,loop_Wmix,loop_Nmix,fiter)
!      ! ss_zeta(:Nso)     = Xvec
!      ss_zeta(:Nso) = loop_Wmix*ss_zeta(:Nso) + (1d0-loop_Wmix)*Fvec
!      if(verbose>2)write(*,"(A6,12G18.9,I4)")"Z_ss =",ss_zeta(:Nso),info
!      !
!      fconverged = check_convergence(ss_zeta(:Nso),loop_tolerance,Nsuccess,loop_Nitermax)
!      call end_loop()
!   end do
!   !
!   if(verbose>1)then
!      write(*,"(A7,12G18.9)")"Dens  =",ss_dens(:Nso)
!      write(*,"(A7,12G18.9)")"Sz+1/2=",ss_Sz(:Nso)+0.5d0
!      write(*,*)""
!      write(*,"(A7,12G18.9)")"Lambda=",lambda(:Nso)
!      write(*,"(A7,12G18.9)")"F_ss  =",ss_Dens(:Nso) - (ss_Sz(:Nso) + 0.5d0)
!      if(verbose>2)call stop_timer()
!      write(*,*)""
!      write(*,*)""
!   endif
! end subroutine ss_solve_routine
