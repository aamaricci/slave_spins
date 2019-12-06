MODULE SS_MAIN
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_SPIN
  !
  USE SF_TIMER, only: start_timer,stop_timer
  USE SF_OPTIMIZE, only: broyden1,fsolve,broyden_mix,linear_mix,adaptive_mix
  !
  USE DMFT_TOOLS
  implicit none
  private

  interface ss_init
     module procedure :: ss_init_hk
     module procedure :: ss_init_dos
  end interface ss_init

  public :: ss_init
  public :: ss_solve
  public :: ss_get_Hf
  public :: ss_get_dens
  public :: ss_get_zeta
  public :: ss_get_lambda
  public :: ss_get_sz


  real(8),dimension(:),allocatable     :: ss_lambda_init
  real(8),dimension(:),allocatable     :: ss_zeta_init

  integer,save                         :: siter=0
  integer                              :: fiter,info
  logical                              :: fconverged

contains


  !< Init SS calculation by passing the Hamiltonian H(k), the k-point weight W(k) and the local
  !  part of the Hamiltonian H_loc = sum_k H(k)
  subroutine ss_init_hk(hk_user,wtk_user,hloc)
    complex(8),dimension(:,:,:)                          :: hk_user  ![Nlso,Nlso,Nk]
    real(8),dimension(size(hk_user,3))                   :: wtk_user ![Nk]
    complex(8),dimension(Nspin*Norb,Nspin*Norb),optional :: Hloc
    complex(8),dimension(Nspin*Norb,Nspin*Norb)          :: Htmp
    complex(8),dimension(:,:),allocatable                :: Hcheck
    logical,save                                         :: isetup=.true.
    integer                                              :: ik,io
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
    !< Init/Read the lambda/zeta input
    if(filling/=dble(Norb))call ss_get_lambda0()
    call ss_init_params()
    if(verbose>2)then
       write(*,"(A7,12G18.9)")"Lam0  =",ss_lambda0
       write(*,"(A7,12G18.9)")"Lam   =",ss_lambda
       write(*,"(A7,12G18.9)")"Z     =",ss_zeta
       write(*,"(A7,12G18.9)")"Ef    =",ss_Ef
       write(*,*)" "
    endif
    !
    !< Internal use:
    allocate(ss_lambda_init(Ns));ss_lambda_init=ss_lambda
    allocate(ss_zeta_init(Ns))  ;ss_zeta_init  =ss_zeta
    !
    isetup=.false.
    return
  end subroutine ss_init_hk

  !< Init SS calculation by passing the Density of states D(e), the dispersions E(e) and the local
  !  part of the Hamiltonian H_loc = sum_k H(k)
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
    call assert_shape(Dbands,[Nspin*Norb,Nk],"ss_init_dos","Dbands")
    call assert_shape(Hloc,[Nspin*Norb],"ss_init_dos","Hloc")
    !
    if(isetup)call ss_setup_structure()
    is_dos=.true.
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
       !
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
    !< Init/Read the lambda/zeta input
    if(filling/=dble(Norb))call ss_get_lambda0()
    call ss_init_params()
    if(verbose>2)then
       write(*,"(A7,12G18.9)")"Lam0  =",ss_lambda0
       write(*,"(A7,12G18.9)")"Lam   =",ss_lambda
       write(*,"(A7,12G18.9)")"Z     =",ss_zeta
       write(*,"(A7,12G18.9)")"Ef    =",ss_Ef
       write(*,*)" "
    endif
    !
    !< Internal use:
    allocate(ss_lambda_init(Ns));ss_lambda_init=ss_lambda
    allocate(ss_zeta_init(Ns))  ;ss_zeta_init  =ss_zeta
    !
    isetup=.false.
    !
    return
  end subroutine ss_init_dos




  !< Execute the SS calculation using different algorithmes as from the input file
  subroutine ss_solve()
    real(8),dimension(Nso)     :: lambda
    real(8),dimension(2*Nso)   :: params
    real(8),dimension(Nso+1)   :: lambda1
    real(8),dimension(2*Nso+1) :: params1
    real(8),dimension(Nso)     :: Xvec,Fvec
    !
    call save_array(trim(Pfile)//trim(ss_file_suffix)//".used",[ss_lambda,ss_zeta,xmu])
    !
    call start_timer
    select case(solve_method)
    case ("fsolve")
       params  = [ss_lambda(:Nso),ss_zeta(:Nso)]
       params1 = [params,xmu]   !we use xmu here to keep it generic with respect to mu, fixed or varied.
       if(filling==0d0)then
          call fsolve(ss_solve_full,params,tol=solve_tolerance)
       else
          call fsolve(ss_solve_full,params1,tol=solve_tolerance)
       end if
       !
    case ("broyden")
       params  = [ss_lambda(:Nso),ss_zeta(:Nso)]
       params1 = [params,xmu]   !we use xmu here to keep it generic with respect to mu, fixed or varied.
       if(filling==0d0)then
          call broyden1(ss_solve_full,params,tolf=solve_tolerance)
       else
          call broyden1(ss_solve_full,params1,tolf=solve_tolerance)
       endif
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
       call ss_lf_solve()
       !
    case default
       write(*,*)"ERROR in ss_solve(): no solve_method named as input *"//str(solve_method)//"*"
       stop
    end select
    !
    call stop_timer
    !
    call save_array(trim(Pfile)//trim(ss_file_suffix)//".restart",[ss_lambda,ss_zeta,xmu])
    !
    call ss_spin_corr()
    call ss_write_last()
    !
  end subroutine ss_solve



  !< solve the SS problem by optimizing the parameter vector X=[lambda,Z{,Ef}]
  !  with the conditions:
  !  n - Sz -1/2 =0
  !  z_i - z_{i-1}=0
  !  {n - filling=0}
  include "ss_main_solve_full.h90"


  !< solve the SS problem by optimizing separately lambda or [lambda,Ef] and Z
  !GG method: optimize broyden/fsolve in lambda and iterate over Z for any fixed lambda 
  include "ss_main_solve_gg.h90"
  !LF method: iterate over lambda, solve fermion at fixed lambda + broyden/fsolve for spins changing lambda, fix Z
  include "ss_main_solve_lf.h90"




  subroutine ss_write_all()
    integer :: unit
    open(free_unit(unit),file="lambda_all.ss",position='append')
    write(unit,*)ss_lambda
    close(unit)
    !
    open(free_unit(unit),file="zeta_all.ss",position='append')
    write(unit,*)ss_zeta
    close(unit)
    !
    open(free_unit(unit),file="Ef_all.ss",position='append')
    write(unit,*)ss_Ef
    close(unit)
    !
    open(free_unit(unit),file="dens_all.ss",position='append')
    write(unit,*)ss_dens
    close(unit)
    !
    open(free_unit(unit),file="sz_all.ss",position='append')
    write(unit,*)ss_Sz
    close(unit)
  end subroutine ss_write_all



  subroutine ss_write_last()
    integer :: unit,units(4),i,iorb,jorb
    real(8) :: SzSz(4)
    open(free_unit(unit),file="hubbards.ss")
    write(unit,"(90F15.9)")(uloc(iorb),iorb=1,Norb),Ust,Jh
    close(unit)
    !
    open(free_unit(unit),file="lambda_last.ss")
    write(unit,*)ss_lambda
    close(unit)
    !
    open(free_unit(unit),file="zeta_last.ss")
    write(unit,*)ss_zeta
    close(unit)
    !
    open(free_unit(unit),file="Ef_last.ss")
    write(unit,*)ss_Ef
    close(unit)
    !
    open(free_unit(unit),file="dens_last.ss")
    write(unit,*)ss_dens
    close(unit)
    !
    open(free_unit(unit),file="sz_last.ss")
    write(unit,*)ss_Sz
    close(unit)
    !
    open(free_unit(unit),file="Op_last.ss")
    write(unit,*)ss_Op
    close(unit)
    !
    units = free_units(4)
    open(units(1),file="SzSz_uu.ss")
    open(units(2),file="SzSz_dd.ss")
    open(units(3),file="SzSz_ud.ss")
    open(units(4),file="SzSz_du.ss")
    do iorb=1,Norb
       do jorb=1,Norb
          do i=1,4
             write(units(i),*)iorb,jorb,ss_SzSz(i,iorb,jorb)
          enddo
       enddo
    enddo
    do i=1,4
       close(units(i))
    enddo
  end subroutine ss_write_last




  subroutine ss_get_Hf(Hk)
    complex(8),dimension(:,:,:) :: Hk
    complex(8),dimension(Ns,Ns) :: diagZ,Hk_f
    real(8),dimension(Ns)       :: sq_zeta
    integer                     :: ik
    call assert_shape(Hk,[Nspin*Norb,Nspin*Norb,Nk],"ss_get_Hf","Hk")
    sq_zeta = sqrt(ss_zeta)
    diagZ   = diag(sq_zeta)
    do ik=1,Nk
       Hk_f       = (diagZ .x. ss_Hk(:,:,ik)) .x. diagZ
       Hk(:,:,ik) = Hk_f(:Nso,:Nso)
    enddo
  end subroutine ss_get_Hf

  subroutine ss_get_dens(dens)
    real(8),dimension(Nspin,Norb) :: dens
    integer :: iorb,ispin
    do ispin=1,Nspin
       do iorb=1,Norb
          dens(ispin,iorb) = ss_dens(iorb+(ispin-1)*Norb)
       enddo
    enddo
  end subroutine ss_get_dens

  subroutine ss_get_zeta(zeta)
    real(8),dimension(Nspin,Norb) :: zeta
    integer :: iorb,ispin
    do ispin=1,Nspin
       do iorb=1,Norb
          zeta(ispin,iorb) = ss_zeta(iorb+(ispin-1)*Norb)
       enddo
    enddo
  end subroutine ss_get_zeta

  subroutine ss_get_sz(sz)
    real(8),dimension(Nspin,Norb) :: sz
    integer :: iorb,ispin
    do ispin=1,Nspin
       do iorb=1,Norb
          sz(ispin,iorb) = ss_sz(iorb+(ispin-1)*Norb)
       enddo
    enddo
  end subroutine ss_get_sz

  subroutine ss_get_lambda(lambda)
    real(8),dimension(Nspin,Norb) :: lambda
    integer :: iorb,ispin
    do ispin=1,Nspin
       do iorb=1,Norb
          lambda(ispin,iorb) = ss_lambda(iorb+(ispin-1)*Norb)
       enddo
    enddo
  end subroutine ss_get_lambda

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




! ss_Hdiag=.true.
! allocate(Hcheck(Ns,Ns))
! do ik=1,Nk
!    Hcheck = ss_Hk(:,:,ik) + ss_Hloc
!    if(sum(abs(Hcheck - diag(diagonal(Hcheck)))) > 1d-6)ss_Hdiag=.false.
! enddo
! deallocate(Hcheck)
