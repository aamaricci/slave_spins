MODULE SS_MAIN
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_SPIN
  !
  USE SF_PAULI
  USE SF_LINALG, only: kron
  USE SF_TIMER, only: start_timer,stop_timer
  USE SF_OPTIMIZE, only: broyden1,fsolve,leastsq!,broyden_mix,linear_mix,adaptive_mix
  !
  implicit none
  private

  interface ss_init
     module procedure :: ss_init_hk
     module procedure :: ss_init_dos
  end interface ss_init


  public :: ss_init
  public :: ss_solve

  ! real(8),dimension(:),allocatable :: ss_lambda_init
  real(8),dimension(:),allocatable :: ss_zeta_init
  integer,save                     :: siter=0
  integer                          :: fiter,info
  logical                          :: fconverged
  integer                          :: iorb,ispin,ilat,ineq,io,il


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
    integer                                                        :: ik
    character(len=5),dimension(3)                                  :: UserOrder_
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
    if(filling/=0d0)call ss_solve_lambda0()
    call ss_init_params()
    !
    if(verbose>2)then
       do ineq=1,Nineq
          write(*,"(A7,12G18.9)")"Lam0  =",ss_lambda0(ineq,:)
          write(*,"(A7,12G18.9)")"Lam   =",ss_lambda(ineq,:)
          write(*,"(A7,12G18.9)")"Z     =",ss_zeta(ineq,:)
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
    !
    if(verbose>2)then
       do ineq=1,Nineq
          write(*,"(A7,12G18.9)")"Lam0  =",ss_lambda0(ineq,:)
          write(*,"(A7,12G18.9)")"Lam   =",ss_lambda(ineq,:)
          write(*,"(A7,12G18.9)")"Z     =",ss_zeta(ineq,:)
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
    real(8),dimension(2*Niso)    :: params
    real(8),dimension(Niso)      :: apar
    real(8),dimension(2*Niso+1)  :: params1
    real(8),dimension(Niso+1)    :: apar1
    real(8),dimension(Nineq*Nss) :: lambda,zeta
    !
    call save_array(trim(Pfile)//trim(ss_file_suffix)//".used",[ss_lambda,ss_zeta,xmu]) !2*Ns+1
    !
    !< Reduce parameters to Inequivalent sites only:
    do ineq=1,Nineq
       ilat = ss_ineq2ilat(ineq)
       ss_lambda_ineq(ineq,:) = ss_lambda(ilat,:)
       ss_zeta_ineq(ineq,:)   = ss_zeta(ilat,:)
    enddo
    !< Pack array:
    lambda =  ss_pack_array(ss_lambda_ineq,Nineq)
    zeta   =  ss_pack_array(ss_zeta_ineq,Nineq)
    !
    call start_timer

    select case(solve_method)
    case ("fsolve")
       params  = [lambda(:Niso),zeta(:Niso)]
       params1 = [params,xmu] 
       if(filling==0d0)then
          call fsolve(ss_solve_full,params,tol=solve_tolerance)
       else
          call fsolve(ss_solve_full,params1,tol=solve_tolerance)
       end if
       !
    case ("broyden")
       params  = [lambda(:Niso),zeta(:Niso)]
       params1 = [params,xmu] 
       if(filling==0d0)then
          call broyden1(ss_solve_full,params,tol=solve_tolerance)
       else
          call broyden1(ss_solve_full,params1,tol=solve_tolerance)
       endif
       !
    case ("leastsq")
       apar  = [lambda(:Niso)]
       apar1 = [params,xmu] 
       if(filling==0d0)then
          call leastsq(ss_solve_least,apar,m=2*Niso,tol=solve_tolerance)
       else
          call leastsq(ss_solve_least,apar1,m=2*Niso+1,tol=solve_tolerance)
       end if
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
    logical                          :: bool
    integer                          :: ilat,ineq,io,il,iorb,ispin
    real(8),dimension(Nineq,Nss)     :: TmpZeta
    real(8),dimension(Nineq*Nss)     :: lambda,zeta,sz,Dens,zeta_prev
    !
    bool = (size(aparams)==2*Niso+1)
    !
    siter = siter+1
    if(verbose>2)call start_loop(siter,0,"FSOLVE-iter")
    !
    !< extract the input for ineq. sites
    lambda(:Niso) = aparams(1:Niso)
    zeta(:Niso)   = aparams(Niso+1:2*Niso)
    if(bool)xmu   = aparams(2*Niso+1)
    if(Nspin==1)then
       call ss_spin_symmetry(lambda,Nineq)
       call ss_spin_symmetry(zeta,Nineq)
    endif
    ss_lambda_ineq = ss_unpack_array(lambda,Nineq)
    ss_zeta_ineq   = ss_unpack_array(zeta,Nineq)
    !
    !< propagate input to all sites
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_lambda(ilat,:) = ss_lambda_ineq(ineq,:)
       ss_zeta(ilat,:)   = ss_zeta_ineq(ineq,:)
    enddo
    !
    !< save input zeta 
    zeta_prev = ss_pack_array(ss_zeta_ineq,Nineq)
    TmpZeta   = ss_zeta_ineq
    !
    !< Solve Fermions:
    call ss_solve_fermions    
    !
    !< Solve Spins:
    do ineq=1,Nineq
       call ss_solve_spins(ineq)
    enddo
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_Sz(ilat,:)       = ss_Sz_ineq(ineq,:)
       ss_Op(ilat,:)       = ss_Op_ineq(ineq,:)
       ss_Zeta(ilat,:)     = ss_Zeta_ineq(ineq,:)
       ss_SzSz(ilat,:,:,:) = ss_SzSz_ineq(ineq,:,:,:)       
    enddo
    if(Nspin==1)call ss_spin_symmetry(ss_zeta,Nlat)
    !
    ! < Constraint:
    do ineq=1,Nineq
       ilat = ss_ineq2ilat(ineq)
       ss_Dens_ineq(ineq,:) = ss_Dens(ilat,:)
    enddo
    Dens    = ss_pack_array(ss_Dens_ineq,Nineq)
    Sz      = ss_pack_array(ss_Sz_ineq,Nineq)
    Zeta    = ss_pack_array(ss_zeta_ineq,Nineq)
    !
    fss(1:Niso)           = Dens(1:Niso) - (Sz(1:Niso) + 0.5d0)
    fss(Niso+1:2*Niso)    = abs(zeta(1:Niso) - zeta_prev(1:Niso))
    if(bool)fss(2*Niso+1) = sum(ss_dens) - filling
    !
    if(verbose>1)then
       do ineq=1,Nineq
          ilat = ss_ineq2ilat(ineq)
          write(*,*)" SITE= "//str(ilat,4)
          !
          if(verbose>3)write(*,"(A7,12G18.9)")"C     =",ss_C_ineq(ineq,:Nspin*Norb)
          !
          write(*,"(A7,12G18.9)")"N     =",ss_Dens_ineq(ineq,:Nspin*Norb),&
               sum(ss_Dens_ineq(ineq,:Nspin*Norb))*(3-Nspin),filling/Nlat
          write(*,"(A7,12G18.9)")"Lambda=",ss_lambda_ineq(ineq,:Nspin*Norb)
          if(verbose>4)write(*,"(A7,12G18.9)")"Z_prev=",TmpZeta(ineq,:Nspin*Norb)
          write(*,"(A7,12G18.9)")"Z_ss  =",ss_zeta_ineq(ineq,:Nspin*Norb)
          write(*,*)" "
       enddo
       write(*,"(A11,50G18.9)")"F_lambda  =",fss(1:Niso)
       write(*,"(A11,50G18.9)")"F_zeta    =",fss(Niso+1:2*Niso)
       if(bool)write(*,"(A11,G18.9)")"F_filling =",fss(2*Niso+1)
       write(*,*)""
    endif
    !
    if(verbose>2)call end_loop()
    !
    call ss_write_last()
    !
    if(siter>=Niter)stop "SS_SOLVE_FULL warning: iter > Niter"
  end function ss_solve_full



  function ss_solve_least(aparams,m) result(fss)
    real(8),dimension(:)         :: aparams !2*Nlso[+1]
    integer                      :: m
    real(8),dimension(m)         :: fss
    logical                      :: bool
    integer                      :: ilat,ineq,io,il,iorb,ispin
    real(8),dimension(Nineq,Nss) :: TmpZeta
    real(8),dimension(Nineq*Nss) :: lambda,zeta,sz,Dens,zeta_prev
    !
    bool = (size(aparams)==Niso+1)
    !
    siter = siter+1
    if(verbose>2)call start_loop(siter,0,"FSOLVE-iter")
    !
    !< extract the input for ineq. sites
    lambda(:Niso) = aparams(1:Niso)
    if(bool)xmu   = aparams(Niso+1)
    if(Nspin==1)call ss_spin_symmetry(lambda,Nineq)
    ss_lambda_ineq = ss_unpack_array(lambda,Nineq)
    !
    !< propagate input to all sites
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_lambda(ilat,:) = ss_lambda_ineq(ineq,:)
       ss_zeta(ilat,:)   = ss_zeta_ineq(ineq,:) !from prev. loop
    enddo
    !
    !< save input zeta 
    zeta_prev = ss_pack_array(ss_zeta_ineq,Nineq)
    TmpZeta   = ss_zeta_ineq
    !
    !< Solve Fermions:
    call ss_solve_fermions    
    !
    !< Solve Spins:
    do ineq=1,Nineq
       call ss_solve_spins(ineq)
    enddo
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_Sz(ilat,:)       = ss_Sz_ineq(ineq,:)
       ss_Op(ilat,:)       = ss_Op_ineq(ineq,:)
       ss_Zeta(ilat,:)     = ss_Zeta_ineq(ineq,:)
       ss_SzSz(ilat,:,:,:) = ss_SzSz_ineq(ineq,:,:,:)       
    enddo
    if(Nspin==1)call ss_spin_symmetry(ss_zeta,Nlat)
    !
    ! < Constraint:
    do ineq=1,Nineq
       ilat = ss_ineq2ilat(ineq)
       ss_Dens_ineq(ineq,:) = ss_Dens(ilat,:)
    enddo
    Dens    = ss_pack_array(ss_Dens_ineq,Nineq)
    Sz      = ss_pack_array(ss_Sz_ineq,Nineq)
    Zeta    = ss_pack_array(ss_zeta_ineq,Nineq)
    !
    fss(1:Niso)           = Dens(1:Niso) - (Sz(1:Niso) + 0.5d0)
    fss(Niso+1:2*Niso)    = abs(zeta(1:Niso) - zeta_prev(1:Niso))
    if(bool)fss(2*Niso+1) = sum(ss_dens) - filling
    !
    if(verbose>1)then
       do ineq=1,Nineq
          ilat = ss_ineq2ilat(ineq)
          write(*,*)" SITE= "//str(ilat,4)
          !
          if(verbose>3)write(*,"(A7,12G18.9)")"C     =",ss_C_ineq(ineq,:Nspin*Norb)
          !
          write(*,"(A7,12G18.9)")"N     =",ss_Dens_ineq(ineq,:Nspin*Norb),&
               sum(ss_Dens_ineq(ineq,:Nspin*Norb))*(3-Nspin),filling/Nlat
          write(*,"(A7,12G18.9)")"Lambda=",ss_lambda_ineq(ineq,:Nspin*Norb)
          if(verbose>4)write(*,"(A7,12G18.9)")"Z_prev=",TmpZeta(ineq,:Nspin*Norb)
          write(*,"(A7,12G18.9)")"Z_ss  =",ss_zeta_ineq(ineq,:Nspin*Norb)
          write(*,*)" "
       enddo
       write(*,"(A11,50G18.9)")"F_lambda  =",fss(1:Niso)
       write(*,"(A11,50G18.9)")"F_zeta    =",fss(Niso+1:2*Niso)
       if(bool)write(*,"(A11,G18.9)")"F_filling =",fss(2*Niso+1)
       write(*,*)""
    endif
    !
    if(verbose>2)call end_loop()
    !
    call ss_write_last()
    !
    if(siter>=Niter)stop "SS_SOLVE_FULL warning: iter > Niter"
  end function ss_solve_least







  subroutine ss_write_last()
    integer                        :: unit,units(4),i,iorb,jorb,ilat
    real(8),dimension(Nlat,2*Norb) :: TmpDens
    open(free_unit(unit),file="hubbards.ss")
    write(unit,"(90F15.9)")(uloc(iorb),iorb=1,Norb),Ust,Jh
    close(unit)
    !
    TmpDens=ss_unpack_array(ss_Dens,Nlat)
    ! TmpLambda= ss_unpack_array(ss_Lambda)
    ! TmpZeta= ss_unpack_array(ss_Zeta)
    ! TmpOp= ss_unpack_array(ss_Op)
    ! TmpSz= ss_unpack_array(ss_Sz)

    do ilat=1,Nlat
       open(free_unit(unit),file="lambda_site"//str(ilat,4)//".ss")
       write(unit,*)ss_lambda(ilat,:)
       close(unit)
       !
       open(free_unit(unit),file="zeta_site"//str(ilat,4)//".ss")
       write(unit,*)ss_zeta(ilat,:)
       close(unit)
       !
       open(free_unit(unit),file="dens_site"//str(ilat,4)//".ss")
       write(unit,*)TmpDens(ilat,:)
       close(unit)
       !
       open(free_unit(unit),file="sz_site"//str(ilat,4)//".ss")
       write(unit,*)ss_Sz(ilat,:)
       close(unit)
       !
       open(free_unit(unit),file="Op_site"//str(ilat,4)//".ss")
       write(unit,*)ss_Op(ilat,:)
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
