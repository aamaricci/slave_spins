MODULE SS_MAIN
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_SPIN
  !
  USE SF_PAULI
  USE SF_LINALG, only: kron
  USE SF_TIMER, only: start_timer,stop_timer
  USE SF_OPTIMIZE, only: broyden1,fsolve!,broyden_mix,linear_mix,adaptive_mix
  !
  implicit none
  private

  interface ss_init
     module procedure :: ss_init_hk
     module procedure :: ss_init_dos
  end interface ss_init


  public :: ss_init
  public :: ss_solve

  real(8),dimension(:),allocatable     :: ss_lambda_init
  real(8),dimension(:),allocatable     :: ss_zeta_init
  integer,save                         :: siter=0
  integer                              :: fiter,info
  logical                              :: fconverged

contains



  !< Init SS calculation by passing the Hamiltonian H(k)
  subroutine ss_init_hk(hk_user,UserOrder,Hloc,ineq_sites)
    complex(8),dimension(:,:,:)                                    :: hk_user  ![Nlso,Nlso,Nk]
    character(len=*),dimension(3),optional                         :: UserOrder
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb),optional :: Hloc
    integer,dimension(Nlat),optional                               :: ineq_sites
    !
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb)          :: Htmp,Hk
    complex(8),dimension(:,:),allocatable                          :: Hcheck
    logical,save                                                   :: isetup=.true.
    integer                                                        :: ik,io,ilat
    character(len=5),dimension(3)                                  :: UserOrder_
    real(8),dimension(Nlat,2*Norb)                                 :: TmpLambda0,TmpZeta,TmpLambda
    !
    UserOrder_ = [character(len=5) :: "Norb","Nspin","Nlat"];
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Nk = size(hk_user,3)
    !
    !< Init SS parameters:
    if(isetup)then
       if(present(ineq_sites))then
          call ss_setup_structure(ineq_sites)
       else
          call ss_setup_structure()
       endif
    endif
    !
    !< Init the SS structure + memory allocation
    call assert_shape(hk_user,[Nlso,Nlso,Nk],"ss_init_hk","hk_user")
    !
    !< Init the Hk structures
    Htmp = sum(Hk_user,dim=3)/Nk;where(abs(Htmp)<1d-6)Htmp=zero
    !
    if(present(Hloc).AND.(sum(abs(Htmp))>1d-12))&
         stop "ss_init: Hloc seems to be present twice: in _Hloc and _Hk"
    !
    do ik=1,Nk
       !< if order of Hk_user is not correct set the SS_order function to actual reorder
       Hk = ss_user2ss(Hk_user(:,:,ik),UserOrder_)
       !
       select case(Nspin)
       case default
          ss_Hk(:,:,ik)  = kron(pauli_0,Hk - Htmp)
       case (2)
          ss_Hk(:,:,ik)  = Hk - Htmp
       end select
       ss_Wtk(:,:,ik) = 1d0/Nk  !Wtk_user(ik)
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
    !if(filling/=dble(Norb))call ss_solve_lambda0()
    if(filling/=0d0)call ss_solve_lambda0()
    call ss_init_params()
    !
    if(verbose>2)then
       call ss_spread_array(ss_lambda0,TmpLambda0)
       call ss_spread_array(ss_lambda,TmpLambda)
       call ss_spread_array(ss_zeta,TmpZeta)
       do ilat=1,Nlat
          write(*,"(A7,12G18.9)")"Lam0  =",TmpLambda0(ilat,:)
          write(*,"(A7,12G18.9)")"Lam   =",TmpLambda(ilat,:)
          write(*,"(A7,12G18.9)")"Z     =",TmpZeta(ilat,:)
       enddo
       write(*,"(A7,12G18.9)")"mu    =",xmu
       write(*,*)" "
    endif
    !
    isetup=.false.
    return
  end subroutine ss_init_hk



  !< Init SS calculation by passing the Density of states D(e), the dispersions E(e) and the local
  !  part of the Hamiltonian H_loc = sum_k H(k)
  subroutine ss_init_dos(Ebands,Dbands,UserOrder,Hloc,ineq_sites)
    real(8),dimension(:,:)                                :: Ebands  ![Nlso,Ne]
    real(8),dimension(:,:)                                :: Dbands ![Nlso,Ne]
    character(len=*),dimension(3),optional                :: UserOrder
    real(8),dimension(:),optional                         :: Hloc
    integer,dimension(Nlat),optional                      :: ineq_sites
    !
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb) :: Htmp
    real(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb)    :: Wtmp
    real(8),dimension(Nspin*Nlat*Norb)                    :: Eb,Db,Hloc_
    integer                                               :: ie,io,ilat
    logical,save                                          :: isetup=.true.
    character(len=5),dimension(3)                         :: UserOrder_
    real(8),dimension(Nlat,2*Norb)                        :: TmpLambda0,TmpZeta,TmpLambda
    !
    UserOrder_ = [character(len=5) :: "Norb","Nspin","Nlat"];
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Nk = size(Ebands,2)
    !
    if(isetup)then
       if(present(ineq_sites))then
          call ss_setup_structure(ineq_sites)
       else
          call ss_setup_structure()
       endif
    endif
    is_dos=.true.
    !
    !< Init the SS structure + memory allocation
    call assert_shape(Ebands,[Nspin*Nlat*Norb,Nk],"ss_init_dos","Ebands")
    call assert_shape(Dbands,[Nspin*Nlat*Norb,Nk],"ss_init_dos","Dbands")
    Hloc_      = zero
    if(present(Hloc))then
       call assert_shape(Hloc,[Nspin*Nlat*Norb],"ss_init_dos","Hloc")
       Hloc_ = Hloc
    endif
    !
    !
    ss_Hk = zero
    ss_Wtk= 0d0
    do ie=1,Nk
       Eb = ss_user2ss(Ebands(:,ie),UserOrder_)
       Db = ss_user2ss(Dbands(:,ie),UserOrder_)
       Htmp=zero
       Wtmp=0d0
       do io=1,Nspin*Nlat*Norb
          Htmp(io,io)  = one*Eb(io) - one*Hloc_(io)
          Wtmp(io,io)  = Db(io)
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
       ss_Hloc = kron(pauli_0,one*diag(Hloc_))
    case (2)
       ss_Hloc = diag(Hloc_)
    end select
    !
    !< Init/Read the lambda/zeta input
    ! if(filling/=dble(Norb))call ss_solve_lambda0()
    if(filling/=0d0)call ss_solve_lambda0()
    call ss_init_params()
    if(verbose>2)then
       call ss_spread_array(ss_lambda0,TmpLambda0)
       call ss_spread_array(ss_lambda,TmpLambda)
       call ss_spread_array(ss_zeta,TmpZeta)
       do ilat=1,Nlat
          write(*,"(A7,12G18.9)")"Lam0  =",TmpLambda0(ilat,:)
          write(*,"(A7,12G18.9)")"Lam   =",TmpLambda(ilat,:)
          write(*,"(A7,12G18.9)")"Z     =",TmpZeta(ilat,:)
       enddo
       write(*,"(A7,12G18.9)")"mu    =",xmu
       write(*,*)" "
    endif
    !
    isetup=.false.
    !
    return
  end subroutine ss_init_dos




  !< Execute the SS calculation using different algorithmes as from the input file
  subroutine ss_solve()
    real(8),dimension(Nlso)     :: lambda
    real(8),dimension(2*Nlso)   :: params
    real(8),dimension(Nlso+1)   :: lambda1
    real(8),dimension(2*Nlso+1) :: params1
    ! real(8),dimension(Nlso)     :: Xvec,Fvec
    !
    call save_array(trim(Pfile)//trim(ss_file_suffix)//".used",[ss_lambda,ss_zeta,xmu]) !2*Ns+1
    !
    call start_timer
    select case(solve_method)
    case ("fsolve")
       params  = [ss_lambda(:Nlso),ss_zeta(:Nlso)]
       params1 = [params,xmu] 
       if(filling==0d0)then
          call fsolve(ss_solve_full,params,tol=solve_tolerance)
       else
          call fsolve(ss_solve_full,params1,tol=solve_tolerance)
       end if
       !
    case ("broyden")
       params  = [ss_lambda(:Nlso),ss_zeta(:Nlso)]
       params1 = [params,xmu] 
       if(filling==0d0)then
          call broyden1(ss_solve_full,params,tol=solve_tolerance)
       else
          call broyden1(ss_solve_full,params1,tol=solve_tolerance)
       endif
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
    call ss_write_last()
    !
  end subroutine ss_solve



  !< solve the SS problem by optimizing the parameter vector X=[lambda,Z{,xmu}]
  !  with the conditions:
  !  n - Sz -1/2 =0
  !  z_i - z_{i-1}=0
  !  {n - filling=0}
  function ss_solve_full(aparams) result(fss)
    real(8),dimension(:),intent(in)  :: aparams !2*Nlso[+1]
    real(8),dimension(size(aparams)) :: fss
    real(8)                          :: zeta(Nlso)
    logical                          :: bool
    integer                          :: ilat,ineq,io,il,iorb,ispin
    real(8),dimension(Nlat,2*Norb)   :: TmpDens,TmpZeta,TmpLambda,TmpC
    !
    bool = (size(aparams)==2*Nlso+1)
    !
    siter = siter+1
    !
    if(verbose>2)call start_loop(siter,0,"FSOLVE-iter")
    !
    ss_lambda(1:Nlso) = aparams(1:Nlso)
    ss_zeta(1:Nlso)   = aparams(Nlso+1:2*Nlso)
    if(bool)xmu       = aparams(2*Nlso+1)
    !
    zeta = ss_zeta(1:Nlso)
    call ss_solve_fermions
    do ineq=1,Nineq
       call ss_solve_spins(ineq)
    enddo
    !< propagate spin-solution to all sites
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       do ispin=1,2
          do iorb=1,Norb
             io = ss_indices2i([iorb,ilat,ispin],[Norb,Nlat,2])
             il = ss_indices2i([iorb,ispin],[Norb,2])
             ss_Sz(io)   = ss_Sz_ineq(ineq,il)
             ss_Op(io)   = ss_Op_ineq(ineq,il)
          enddo
       enddo
       ss_SzSz(ilat,:,:,:) = ss_SzSz_ineq(ineq,:,:,:)
    enddo
    !
    !< Get Z:
    ss_Zeta = ss_Op**2
    if(Nspin==1)call ss_spin_symmetry(ss_zeta)
    !
    !<constraint:
    fss(1:Nlso)           = ss_Dens(1:Nlso) - (ss_Sz(1:Nlso) + 0.5d0)
    fss(Nlso+1:2*Nlso)    = ss_zeta(1:Nlso) - zeta
    if(bool)fss(2*Nlso+1) = sum(ss_dens) - filling
    !
    if(verbose>1)then
       call ss_spread_array(ss_c,TmpC)
       call ss_spread_array(ss_dens,TmpDens)
       call ss_spread_array(ss_lambda,TmpLambda)
       call ss_spread_array(ss_zeta,TmpZeta)
       do ilat=1,Nlat
          write(*,*)" SITE= "//str(ilat,4)
          if(verbose>3)write(*,"(A6,12G18.9)")"C    =",TmpC(ilat,:)
          write(*,"(A7,12G18.9)")"N     =",TmpDens(ilat,:Nspin*Norb),&
               sum(TmpDens(ilat,:Nspin*Norb)),filling
          write(*,"(A7,12G18.9)")"Lambda=",TmpLambda(ilat,:Nspin*Norb)
          write(*,"(A7,12G18.9)")"Z_ss  =",TmpZeta(ilat,:Nspin*Norb)
          write(*,*)" "
       enddo
       write(*,*)" - "
       write(*,"(A7,50G18.9)")"F_ss  =",fss
       write(*,*)""
    endif
    !
    if(verbose>2)call end_loop()
    !
    call ss_write_last()
    !
  end function ss_solve_full








  subroutine ss_write_last()
    integer                        :: unit,units(4),i,iorb,jorb,ilat
    real(8),dimension(Nlat,2*Norb) :: TmpDens,TmpZeta,TmpLambda,TmpOp,TmpSz
    open(free_unit(unit),file="hubbards.ss")
    write(unit,"(90F15.9)")(uloc(iorb),iorb=1,Norb),Ust,Jh
    close(unit)
    !
    call ss_spread_array(ss_Dens,TmpDens)
    call ss_spread_array(ss_Lambda,TmpLambda)
    call ss_spread_array(ss_Zeta,TmpZeta)
    call ss_spread_array(ss_Op,TmpOp)
    call ss_spread_array(ss_Sz,TmpSz)

    do ilat=1,Nlat
       open(free_unit(unit),file="lambda_site"//str(ilat,4)//".ss")
       write(unit,*)TmpLambda(ilat,:)
       close(unit)
       !
       open(free_unit(unit),file="zeta_site"//str(ilat,4)//".ss")
       write(unit,*)TmpZeta(ilat,:)
       close(unit)
       !
       open(free_unit(unit),file="dens_site"//str(ilat,4)//".ss")
       write(unit,*)TmpDens(ilat,:)
       close(unit)
       !
       open(free_unit(unit),file="sz_site"//str(ilat,4)//".ss")
       write(unit,*)TmpSz(ilat,:)
       close(unit)
       !
       open(free_unit(unit),file="Op_site"//str(ilat,4)//".ss")
       write(unit,*)TmpOp(ilat,:)
       close(unit)
       !
       units = free_units(4)
       open(units(1),file="SzSz_uu_site"//str(ilat,4)//".ss")
       open(units(2),file="SzSz_dd_site"//str(ilat,4)//".ss")
       open(units(3),file="SzSz_ud_site"//str(ilat,4)//".ss")
       open(units(4),file="SzSz_du_site"//str(ilat,4)//".ss")
       do iorb=1,Norb
          do jorb=1,Norb
             do i=1,4
                write(units(i),*)iorb,jorb,ss_SzSz(ilat,i,iorb,jorb)
             enddo
          enddo
       enddo
       do i=1,4
          close(units(i))
       enddo
       !
    enddo
    !
    open(free_unit(unit),file="mu.ss")
    write(unit,*)xmu
    close(unit)
    !
  end subroutine ss_write_last




  subroutine start_loop(loop,max,name,unit,id)
    integer                   :: loop
    integer,optional          :: max,unit,id
    character(len=*),optional :: name
    character(len=16)         :: loop_name
    integer                   :: unit_,id_
    loop_name="main-loop";if(present(name))loop_name=name
    unit_    =6          ;if(present(unit))unit_=unit
    id_      =0          ;if(present(id))id_=id
    write(unit_,*)
    if(.not.present(max))then
       write(unit_,"(A,I5)")"-----"//trim(adjustl(trim(loop_name))),loop,"-----"
    else
       write(unit_,"(A,I5,A,I5,A)")"-----"//trim(adjustl(trim(loop_name))),loop,&
            " (max:",max,")-----"
    endif
    call start_timer
  end subroutine start_loop


  subroutine end_loop(unit,id)
    integer,optional :: unit,id
    integer          :: unit_,id_
    unit_=6 ; if(present(unit))unit_=unit
    id_  =0 ; if(present(id))id_=id
    write(unit_,"(A)")"====================================="
    call stop_timer
    write(unit_,*)
    write(unit_,*)
  end subroutine end_loop


END MODULE SS_MAIN


!
! case ("gg_broyden")
!    lambda = ss_lambda(:Nlso)
!    call broyden1(ss_solve_lambda,lambda,tolf=solve_tolerance)
!    !
! case ("gg_fsolve")
!    lambda = ss_lambda(:Nlso)
!    call fsolve(ss_solve_lambda,lambda,tol=solve_tolerance)
!    !
! case ("lf_solve")
!    call ss_lf_solve()
!


! !< solve the SS problem by optimizing separately lambda or [lambda,xmu] and Z
! !GG method: optimize broyden/fsolve in lambda and iterate over Z for any fixed lambda 
! include "SS_MAIN/ss_main_solve_gg.h90"
! !LF method: iterate over lambda, solve fermion at fixed lambda + broyden/fsolve for spins changing lambda, fix Z
! include "SS_MAIN/ss_main_solve_lf.h90"


! !< Internal use:
! allocate(ss_lambda_init(Ns));ss_lambda_init=ss_lambda
! allocate(ss_zeta_init(Ns))  ;ss_zeta_init  =ss_zeta


! ss_Hdiag=.true.
! allocate(Hcheck(Ns,Ns))
! do ik=1,Nk
!    Hcheck = ss_Hk(:,:,ik) + ss_Hloc
!    if(sum(abs(Hcheck - diag(diagonal(Hcheck)))) > 1d-6)ss_Hdiag=.false.
! enddo
! deallocate(Hcheck)
