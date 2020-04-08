MODULE SS_SOLVE_MAIN
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  USE SS_IO
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_SPIN
  !
  USE SF_PAULI
  USE SF_TIMER, only: start_timer,stop_timer
  USE SF_OPTIMIZE, only: broyden1,fsolve,leastsq
  !
  implicit none
  private


  public :: ss_solve_methods

  integer,save                     :: siter=0
  integer                          :: iorb,ispin,ilat,ineq,io,il


contains


  !< Execute the SS calculation using different algorithmes as from the input file
  subroutine ss_solve_methods()
    real(8),dimension(2*Niso)    :: params
    real(8),dimension(Niso)      :: apar
    real(8),dimension(2*Niso+1)  :: params1
    real(8),dimension(Niso+1)    :: apar1
    real(8),dimension(Nineq*Nss) :: lambda,zeta
    !
    ! call save_array(trim(Pfile)//trim(ss_file_suffix)//".used",[ss_lambda,ss_zeta,xmu])
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
    !
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
       apar1 = [apar,xmu] 
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
    !< store parameters:
    call save_array(trim(Pfile)//trim(ss_file_suffix)//".restart",&
         [ss_pack_array(ss_lambda,Nlat), ss_pack_array(ss_zeta,Nlat), xmu])
    !
    call ss_write_last()
    !
  end subroutine ss_solve_methods








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
    if(verbose>2)call start_loop(siter,Niter,"FSOLVE-iter")
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
    !< store parameters:
    call save_array(trim(Pfile)//trim(ss_file_suffix)//".used",&
         [ss_pack_array(ss_lambda,Nlat), ss_pack_array(ss_zeta,Nlat), xmu])
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
    call ss_print_screen
    if(verbose>1)then
       write(*,"(A11,50G18.9)")"F_cnstr   =",fss(1:Niso)
       write(*,"(A11,50G18.9)")"F_zeta    =",fss(Niso+1:2*Niso)
       if(bool)write(*,"(A11,G18.9)")"F_filling =",fss(2*Niso+1)
       write(*,*)""
    end if
    !
    if(verbose>2)call end_loop()
    !
    call ss_write_last()
    !
    if(siter>=Niter)stop "SS_SOLVE_FULL warning: iter > Niter"
  end function ss_solve_full





  function ss_solve_least(aparams,m) result(fss)
    real(8),dimension(:)         :: aparams
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
    if(verbose>2)call start_loop(siter,Niter,"LEASTSQ-iter")
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
    !< store parameters:
    call save_array(trim(Pfile)//trim(ss_file_suffix)//".used",&
         [ss_pack_array(ss_lambda,Nlat), ss_pack_array(ss_zeta,Nlat), xmu])
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
    call ss_print_screen
    if(verbose>1)then
       write(*,"(A11,50G18.9)")"F_cnstr   =",fss(1:Niso)
       write(*,"(A11,50G18.9)")"F_zeta    =",fss(Niso+1:2*Niso)
       if(bool)write(*,"(A11,G18.9)")"F_filling =",fss(2*Niso+1)
       write(*,*)""
    end if
    !
    if(verbose>2)call end_loop()
    !
    call ss_write_last()
    !
    if(siter>=Niter)stop "SS_SOLVE_FULL warning: iter > Niter"
  end function ss_solve_least




  subroutine ss_print_screen
    if(verbose>1)then
       do ineq=1,Nineq
          ilat = ss_ineq2ilat(ineq)
          write(*,*)" SITE= "//str(ilat,4)
          !
          if(verbose>3)write(*,"(A7,12G18.9)")"C     =",ss_C_ineq(ineq,:Nspin*Norb)
          !
          write(*,"(A7,12G18.9)")"N     =",ss_Dens_ineq(ineq,:Nspin*Norb),&
               sum(ss_Dens_ineq(ineq,:Nspin*Norb))*(3-Nspin),filling/Nlat
          write(*,"(A7,12G18.9)")"Sz    =",ss_Sz_ineq(ineq,:Nspin*Norb)
          write(*,"(A7,12G18.9)")"Lambda=",ss_lambda_ineq(ineq,:Nspin*Norb)
          if(verbose>3)write(*,"(A7,12G18.9)")"Op    =",ss_Op_ineq(ineq,:Nspin*Norb)
          ! if(verbose>4)write(*,"(A7,12G18.9)")"Z_prev=",TmpZeta(ineq,:Nspin*Norb)
          write(*,"(A7,12G18.9)")"Z_ss  =",ss_zeta_ineq(ineq,:Nspin*Norb)
          write(*,*)" "
       enddo
    endif
  end subroutine ss_print_screen

  
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


END MODULE SS_SOLVE_MAIN


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
