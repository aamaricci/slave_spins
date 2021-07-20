MODULE SS_SOLVE_MAIN
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  USE SS_IO
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_SPIN
  !
  USE SF_IOTOOLS, only: str
  USE SF_FONTS
  USE SF_PAULI
  USE SF_TIMER, only: start_timer,stop_timer
  USE SF_OPTIMIZE, only: broyden1,fsolve,leastsq,fmin_cg
  !
  implicit none
  private


  public :: ss_solve_methods

  integer,save                       :: siter=0
  integer                            :: iorb,ispin,ilat,ineq,io,il
  real(8),dimension(:,:),allocatable :: ss_lambda_init
  real(8),dimension(:,:),allocatable :: ss_Op_init

contains


  !< Execute the SS calculation using different algorithmes as from the input file
  subroutine ss_solve_methods()
    real(8),dimension(2*Niso)    :: params
    real(8),dimension(Niso)      :: apar
    real(8),dimension(2*Niso+1)  :: params1
    real(8),dimension(Niso+1)    :: apar1
    real(8),dimension(Niso)      :: lambda,Op
    integer                      :: iter
    real(8)                      :: chi
    !
    !< Reduce parameters to Inequivalent sites only:
    do ineq=1,Nineq
       ilat = ss_ineq2ilat(ineq)
       ss_lambda_ineq(ineq,:) = ss_lambda(ilat,:)
       ss_Op_ineq(ineq,:)     = ss_Op(ilat,:)
    enddo
    !< Pack array:
    lambda =  ss_pack_array(ss_lambda_ineq,Nineq)
    Op     =  ss_pack_array(ss_Op_ineq,Nineq)
    !
    allocate(ss_lambda_init(Nineq,Nso));ss_lambda_init=ss_lambda_ineq
    allocate(ss_Op_init(Nineq,Nso))    ;ss_Op_init  =ss_Op_ineq
    !
    params  = [lambda,Op]
    params1 = [params,xmu]
    !
    apar  = [lambda]
    apar1 = [apar,xmu]
    !
    if(master)call start_timer
    !
    select case(solve_method)
    case ("fsolve")
       if(filling==0d0)then
          call fsolve(ss_solve_full,params,tol=solve_tolerance,check=(verbose>3))
       else
          call fsolve(ss_solve_full,params1,tol=solve_tolerance,check=(verbose>3))
       end if
       !
    case ("broyden")
       if(filling==0d0)then
          call broyden1(ss_solve_full,params,tol=solve_tolerance)
       else
          call broyden1(ss_solve_full,params1,tol=solve_tolerance)
       endif
       !
    case ("leastsq")
       if(filling==0d0)then
          call leastsq(ss_solve_least,apar,m=2*Niso,tol=solve_tolerance)
       else
          call leastsq(ss_solve_least,apar1,m=2*Niso+1,tol=solve_tolerance)
       end if
       !
    case ("conjg")
       if(filling==0d0)then
          call fmin_cg(params, ss_solve_cg, iter, chi)
       else
          call fmin_cg(params1,ss_solve_cg, iter, chi)
       end if
       if(verbose>1)write(*,"(A,ES18.9,1x,I5)")"Chi^2_fss|iter:",chi,"|",iter
       !
    case ("gg_broyden")
       if(filling==0d0)then
          call broyden1(ss_solve_gg,apar,tol=solve_tolerance)
       else
          call broyden1(ss_solve_gg,apar1,tol=solve_tolerance)
       endif
       !
    case ("gg_fsolve")
       if(filling==0d0)then
          call fsolve(ss_solve_gg,apar,tol=solve_tolerance,check=(verbose>3))
       else
          call fsolve(ss_solve_gg,apar1,tol=solve_tolerance,check=(verbose>3))
       endif
       !
    case default
       write(*,*)"ERROR in ss_solve(): no solve_method named as input *"//str(solve_method)//"*"
       stop
    end select
    !
    if(master)call stop_timer
    !
    !< store parameters:
    if(master)call save_array(trim(Pfile)//trim(ss_file_suffix)//".restart",&
         [&
         ss_pack_array(ss_lambda,Nlat), &
         ss_pack_array(ss_Op,Nlat),     &
         pack(ss_pack_array(ss_OdgOp,Nlat),.true.),& !this creates a 2d array, packed into a 1d
         xmu])
    !
    call ss_write_last()
    !
  end subroutine ss_solve_methods







  !
  !< Solve F(x)=0, w/ F:\RRR^{2Nlso+1}--->\RRR^{2Nlso+1}.
  ! used in: FSOLVE, BROYDEN
  function ss_solve_full(aparams) result(fss)
    real(8),dimension(:),intent(in)  :: aparams !2*Nlso[+1]
    real(8),dimension(size(aparams)) :: fss
    logical                          :: bool
    integer                          :: ilat,ineq,io,il,iorb,ispin
    real(8),dimension(Nineq,Nso)     :: TmpZ
    real(8),dimension(Nineq*Nso)     :: lambda,Zeta,Op,sz,Dens,Zeta_prev
    !
    bool = (size(aparams)==2*Niso+1)
    !
    siter = siter+1
    call start_loop(siter,Niter,"FSOLVE-iter")
    !
    !< extract the input for ineq. sites
    lambda      = aparams(1:Niso)
    Op          = aparams(Niso+1:2*Niso)
    if(bool)xmu = aparams(2*Niso+1)
    !
    !< Map onto ineq. arrays
    ss_lambda_ineq = ss_unpack_array(lambda,Nineq)
    ss_Op_ineq     = ss_unpack_array(Op,Nineq)
    !
    !< propagate input to all sites
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_lambda(ilat,:) = ss_lambda_ineq(ineq,:)
       ss_Op(ilat,:)     = ss_Op_ineq(ineq,:)
    enddo
    !
    !< store parameters:
    if(master)call save_array(trim(Pfile)//trim(ss_file_suffix)//".used",&
         [ss_pack_array(ss_lambda,Nlat), ss_pack_array(ss_Op,Nlat), xmu])
    !
    !< save input
    Zeta_prev = ss_pack_array(ss_Op_ineq**2,Nineq)
    ! Zeta_prev = Zeta_Prev**2
    !
    !< Solve Fermions:
    if(master.AND.verbose>2)call start_timer
    call ss_solve_fermions
    if(master.AND.verbose>2)call stop_timer("solve fermions")
    !
    !< Solve Spins:
    if(master.AND.verbose>2)call start_timer
    do ineq=1,Nineq
       call ss_solve_spins(ineq)
    enddo
    if(master.AND.verbose>2)call stop_timer("solve spins")
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_Sz(ilat,:)       = ss_Sz_ineq(ineq,:)
       ss_Op(ilat,:)       = ss_Op_ineq(ineq,:)
       ss_SzSz(ilat,:,:,:) = ss_SzSz_ineq(ineq,:,:,:)
       ss_OdgOp(ilat,:,:)  = ss_OdgOp_ineq(ineq,:,:)
    enddo
    !
    !add here the formula updating lambda0
    if(self_lam0)call ss_update_lambda0()
    !
    ! < Constraint:
    Dens  = ss_pack_array(ss_Dens_ineq,Nineq)
    Sz    = ss_pack_array(ss_Sz_ineq,Nineq)
    Op    = ss_pack_array(ss_Op_ineq,Nineq)
    Zeta  = Op**2
    !
    fss(1:Niso)           = Dens - (Sz + 0.5d0)
    fss(Niso+1:2*Niso)    = Zeta - Zeta_prev
    if(bool)fss(2*Niso+1) = sum(ss_dens)*(3-Nspin) - filling
    !
    call ss_print_screen()
    if(master.AND.verbose>1)then
       write(*,"(A11,50G18.9)")"F_cnstr   =",fss(1:Niso)
       write(*,"(A11,50G18.9)")"F_z       =",fss(Niso+1:2*Niso)
       if(bool)write(*,"(A11,G18.9)")"F_filling =",fss(2*Niso+1)
       write(*,*)""
    end if
    !
    call end_loop()
    !
    call ss_write_last()
    !
    if(siter>=Niter)stop "SS_SOLVE_FULL warning: iter > Niter"
  end function ss_solve_full







  !
  !< Solve F(x)=0, w/ F:\RRR^{2Nlso+1}--->\RRR
  ! used in: CONJUGATE GRADIENT FMIN
  function ss_solve_cg(aparams) result(chi2)
    real(8),dimension(:)             :: aparams !2*Nlso[+1]
    real(8),dimension(size(aparams)) :: fss
    real(8)                          :: chi2
    logical                          :: bool
    integer                          :: ilat,ineq,io,il,iorb,ispin
    real(8),dimension(Nineq,Nso)     :: TmpZ
    real(8),dimension(Nineq*Nso)     :: lambda,Op,sz,Dens,Zeta,Zeta_prev
    !
    bool = (size(aparams)==2*Niso+1)
    !
    siter = siter+1
    call start_loop(siter,Niter,"FSOLVE-iter")
    !
    !< extract the input for ineq. sites
    lambda      = aparams(1:Niso)
    Op          = aparams(Niso+1:2*Niso)
    if(bool)xmu = aparams(2*Niso+1)
    !
    !< Map onto ineq. arrays
    ss_lambda_ineq = ss_unpack_array(lambda,Nineq)
    ss_Op_ineq     = ss_unpack_array(Op,Nineq)
    !
    !< propagate input to all sites
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_lambda(ilat,:) = ss_lambda_ineq(ineq,:)
       ss_Op(ilat,:)     = ss_Op_ineq(ineq,:)
    enddo
    !
    !< store parameters:
    if(master)call save_array(trim(Pfile)//trim(ss_file_suffix)//".used",&
         [ss_pack_array(ss_lambda,Nlat), ss_pack_array(ss_Op,Nlat), xmu])
    !
    !< save input  
    Zeta_prev = ss_pack_array(ss_Op_ineq,Nineq);Zeta_prev=Zeta_prev**2
    TmpZ      = ss_Op_ineq**2
    !
    !< Solve Fermions:
    if(master.AND.verbose>2)call start_timer
    call ss_solve_fermions
    if(master.AND.verbose>2)call stop_timer("solve fermions")
    !
    !< Solve Spins:
    if(master.AND.verbose>2)call start_timer
    do ineq=1,Nineq
       call ss_solve_spins(ineq)
    enddo
    if(master.AND.verbose>2)call stop_timer("solve spins")
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_Sz(ilat,:)       = ss_Sz_ineq(ineq,:)
       ss_Op(ilat,:)       = ss_Op_ineq(ineq,:)
       ss_SzSz(ilat,:,:,:) = ss_SzSz_ineq(ineq,:,:,:)
       ss_OdgOp(ilat,:,:)  = ss_OdgOp_ineq(ineq,:,:)
    enddo
    !
    !add here the formula updating lambda0
    if(self_lam0)call ss_update_lambda0()
    !
    ! < Constraint:
    do ineq=1,Nineq
       ilat = ss_ineq2ilat(ineq)
       ss_Dens_ineq(ineq,:) = ss_Dens(ilat,:)
    enddo
    Dens    = ss_pack_array(ss_Dens_ineq,Nineq)
    Sz      = ss_pack_array(ss_Sz_ineq,Nineq)
    Op      = ss_pack_array(ss_Op_ineq,Nineq)
    Zeta    = Op**2
    !
    Fss(1:Niso)           = Dens - (Sz + 0.5d0)
    Fss(Niso+1:2*Niso)    = Zeta - Zeta_prev
    if(bool)fss(2*Niso+1) = sum(ss_dens) - filling
    !
    chi2 = dot_product(Fss,Fss)/size(Fss)
    !
    call ss_print_screen(TmpZ)
    if(master.AND.verbose>1)then
       write(*,"(A11,50G18.9)")"F_cnstr   =",fss(1:Niso)
       write(*,"(A11,50G18.9)")"F_Z       =",fss(Niso+1:2*Niso)
       if(bool)write(*,"(A11,G18.9)")"F_filling =",fss(2*Niso+1)
       write(*,*)""
    endif
    !
    call end_loop()
    !
    call ss_write_last()
    !
    if(siter>=Niter)stop "SS_SOLVE_FULL warning: iter > Niter"
  end function ss_solve_cg






  !
  !< Solve F(x)=0, w/ F:\RRR^{Nlso+1}--->\RRR^{2Nlso+1}.
  ! used in: LEASTSQ
  function ss_solve_least(aparams,m) result(fss)
    real(8),dimension(:)         :: aparams
    integer                      :: m
    real(8),dimension(m)         :: fss
    logical                      :: bool
    integer                      :: ilat,ineq,io,il,iorb,ispin
    real(8),dimension(Nineq,Nso) :: TmpZ
    real(8),dimension(Nineq*Nso) :: lambda,Op,sz,Dens,Zeta,Zeta_prev
    !
    bool = (size(aparams)==Niso+1)
    !
    siter = siter+1
    call start_loop(siter,Niter,"LEASTSQ-iter")
    !
    !< extract the input for ineq. sites
    lambda      = aparams(1:Niso)
    if(bool)xmu = aparams(Niso+1)
    !
    !< Symmetrize:
    if(Nspin==1)call ss_spin_symmetry(ss_lambda,Nlat)
    !
    !< Map onto ineq. arrays
    ss_lambda_ineq = ss_unpack_array(lambda,Nineq)
    !
    !< propagate input to all sites
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_lambda(ilat,:) = ss_lambda_ineq(ineq,:)
       ss_Op(ilat,:)     = ss_Op_ineq(ineq,:) !From previous loop
    enddo
    !
    !< store parameters:
    if(master)call save_array(trim(Pfile)//trim(ss_file_suffix)//".used",&
         [ss_pack_array(ss_lambda,Nlat), ss_pack_array(ss_Op,Nlat), xmu])
    !
    !< save input  
    Zeta_prev = ss_pack_array(ss_Op_ineq**2,Nineq)
    TmpZ      = ss_Op_ineq**2
    !
    if(master.AND.verbose>2)call start_timer
    call ss_solve_fermions
    if(master.AND.verbose>2)call stop_timer("solve fermions")
    !
    !< Solve Spins:
    if(master.AND.verbose>2)call start_timer
    do ineq=1,Nineq
       call ss_solve_spins(ineq)
    enddo
    if(master.AND.verbose>2)call stop_timer("solve spins")
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_Sz(ilat,:)       = ss_Sz_ineq(ineq,:)
       ss_Op(ilat,:)       = ss_Op_ineq(ineq,:)
       ss_SzSz(ilat,:,:,:) = ss_SzSz_ineq(ineq,:,:,:)
       ss_OdgOp(ilat,:,:)  = ss_OdgOp_ineq(ineq,:,:)
    enddo
    !
    !add here the formula updating lambda0
    if(self_lam0)call ss_update_lambda0()
    !
    ! < Constraint:
    do ineq=1,Nineq
       ilat = ss_ineq2ilat(ineq)
       ss_Dens_ineq(ineq,:) = ss_Dens(ilat,:)
    enddo
    !
    Dens    = ss_pack_array(ss_Dens_ineq,Nineq)
    Sz      = ss_pack_array(ss_Sz_ineq,Nineq)
    Op      = ss_pack_array(ss_Op_ineq,Nineq)
    Zeta    = Op**2
    !
    fss(1:Niso)           = Dens - (Sz + 0.5d0)
    fss(Niso+1:2*Niso)    = Zeta - Zeta_prev
    if(bool)fss(2*Niso+1) = sum(ss_dens) - filling
    !
    call ss_print_screen(TmpZ)
    if(master.AND.verbose>1)then
       write(*,"(A11,50G18.9)")"F_cnstr   =",fss(1:Niso)
       write(*,"(A11,50G18.9)")"F_Z       =",fss(Niso+1:2*Niso)
       if(bool)write(*,"(A11,G18.9)")"F_filling =",fss(2*Niso+1)
       write(*,*)""
    endif
    !
    call end_loop()
    !
    call ss_write_last()
    !
    if(siter>=Niter)stop "SS_SOLVE_LEASTSQ warning: iter > Niter"
  end function ss_solve_least












  !
  !< Solve G(x)=0, w/ G:\RRR^{Nlso+1}--->\RRR^{Nlso+1}.
  ! used in: FSOLVE, BROYDEN, with G!=F and with an internal
  ! Z-loop
  function ss_solve_gg(aparams) result(fss)
    real(8),dimension(:),intent(in)  :: aparams
    real(8),dimension(size(aparams)) :: fss
    integer                          :: iter,Nsuccess=0
    logical                          :: z_converged,bool
    real(8),dimension(Nineq*Nso)     :: lambda,Op,sz,Dens,Zeta,Zeta_prev
    !
    bool = size(aparams)==Niso+1
    !
    ! > Retrieve ss_lambda:
    lambda      = aparams(1:Niso)
    if(bool)xmu = aparams(Niso+1)
    !
    !< Symmetrize:
    if(Nspin==1)call ss_spin_symmetry(ss_lambda,Nlat)
    !
    !< Map onto ineq. arrays
    ss_lambda_ineq = ss_unpack_array(lambda,Nineq)
    if(restart_init)ss_Op_ineq = ss_Op_init
    !
    !< propagate input to all sites
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_lambda(ilat,:) = ss_lambda_ineq(ineq,:)
       ss_Op(ilat,:)     = ss_Op_ineq(ineq,:)
    enddo
    !
    !
    !> Solve:
    if(master.AND.verbose>2)call start_timer()
    !
    z_converged=.false. ; iter=0
    do while(.not.z_converged.AND.iter<=Niter)
       iter=iter+1
       call start_loop(iter,Niter,"Z-iter")
       !
       Zeta_prev = ss_pack_array(ss_Op_ineq**2,Nineq)
       !
       if(master.AND.verbose>2)call start_timer
       call ss_solve_fermions
       if(master.AND.verbose>2)call stop_timer("solve fermions")
       !
       !< Solve Spins:
       if(master.AND.verbose>2)call start_timer
       do ineq=1,Nineq
          call ss_solve_spins(ineq)
       enddo
       if(master.AND.verbose>2)call stop_timer("solve spins")
       do ilat=1,Nlat
          ineq = ss_ilat2ineq(ilat)
          ss_Sz(ilat,:)       = ss_Sz_ineq(ineq,:)
          ss_Op(ilat,:)       = ss_Op_ineq(ineq,:)
          ss_SzSz(ilat,:,:,:) = ss_SzSz_ineq(ineq,:,:,:)
          ss_OdgOp(ilat,:,:)  = ss_OdgOp_ineq(ineq,:,:)
       enddo
       !
       !add here the formula updating lambda0
       if(self_lam0)call ss_update_lambda0()
       !
       do ineq=1,Nineq
          ilat = ss_ineq2ilat(ineq)
          ss_Dens_ineq(ineq,:) = ss_Dens(ilat,:)
       enddo
       Dens    = ss_pack_array(ss_Dens_ineq,Nineq)
       Sz      = ss_pack_array(ss_Sz_ineq,Nineq)
       Op      = ss_pack_array(ss_Op_ineq,Nineq)
       Zeta    = Op**2
       !
       call ss_print_screen
       !
       Zeta = loop_Wmix*Zeta + (1d0-loop_Wmix)*Zeta_prev
       !
       if(master)z_converged = check_convergence((Zeta-Zeta_prev),loop_tolerance,Nsuccess,Niter)
#ifdef _MPI
       if(check_MPI())call bcast_MPI(MPI_COMM_WORLD,z_converged)
#endif
       call end_loop() 
    end do
    !
    !<constraint:
    fss(1:Niso)         = Dens - (Sz + 0.5d0)
    if(bool)fss(Niso+1) = sum(ss_Dens) - filling
    !
    if(master.AND.verbose>1)then
       write(*,"(A11,50G18.9)")"F_cnstr   =",fss(1:Niso)
       if(bool)write(*,"(A11,G18.9)")"F_filling =",fss(Niso+1)
       write(*,*)""
    endif
    if(master.AND.verbose>2)call stop_timer()
    !
    call ss_write_last()
    !
  end function ss_solve_gg








  !##################################################################
  !
  !                     AUXILIARY FUNCTIONS
  !
  !##################################################################


  subroutine ss_print_screen(Tmp)
    real(8),dimension(Nineq,Nso),optional :: Tmp
    if(master.AND.verbose>1)then
       do ineq=1,Nineq
          ilat = ss_ineq2ilat(ineq)
          write(*,*)" SITE= "//str(ilat,4)
          !
          if(verbose>3)write(*,"(A7,12G18.9)")"C     =",ss_C_ineq(ineq,:)
          !
          write(*,"(A7,12G18.9)")"mu    =",xmu
          write(*,"(A7,12G18.9)")"N     =",ss_Dens_ineq(ineq,:),sum(ss_Dens)*(3-Nspin),filling
          write(*,"(A7,12G18.9)")"Sz+1/2=",ss_Sz_ineq(ineq,:)+0.5d0
          write(*,"(A7,12G18.9)")"Lambda=",ss_lambda_ineq(ineq,:)
          if(verbose>3)write(*,"(A7,12G18.9)")"Op    =",ss_Op_ineq(ineq,:)
          if(verbose>4.AND.present(Tmp))&
               write(*,"(A7,12G18.9)")"Z_prev=",Tmp(ineq,:)
          write(*,"(A7,12G18.9)")"Z_ss  =",ss_Op_ineq(ineq,:)**2
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
    if(master.AND.verbose>2)then
       write(unit_,*)
       if(.not.present(max))then
          write(unit_,"(A,I5)")"-----"//trim(adjustl(trim(loop_name))),loop,"-----"
       else
          write(unit_,"(A,I5,A,I5,A)")"-----"//trim(adjustl(trim(loop_name))),loop,&
               " (max:",max,")-----"
       endif
       call start_timer
    endif
  end subroutine start_loop


  subroutine end_loop(unit,id)
    integer,optional :: unit,id
    integer          :: unit_,id_
    unit_=6 ; if(present(unit))unit_=unit
    id_  =0 ; if(present(id))id_=id
    if(master.AND.verbose>2)then
       write(unit_,"(A)")"====================================="
       call stop_timer
       write(unit_,*)
       write(unit_,*)
    endif
  end subroutine end_loop



  function check_convergence(Xnew,eps,N1,N2) result(convergence)
    real(8),intent(in)            :: Xnew(:)
    real(8)                       :: eps
    integer                       :: N1,N2
    integer                       :: unit
    integer                       :: i,j,Msize1
    logical                       :: convergence  
    real(8)                       :: error(2),err
    real(8),dimension(size(Xnew)) :: Verror
    real(8),save,allocatable      :: Xold(:)
    integer,save                  :: success=0,check=1
    character(len=2)              :: label
    character(len=100)            :: file_
    file_='error.err'
    Msize1=size(Xnew)
    if(.not.allocated(Xold))then
       allocate(Xold(Msize1))
       Xold=0.d0
       open(free_unit(unit),file=str(file_),position="append")
    endif
    write(unit,*);close(unit)
    !
    Verror=abs(Xnew-Xold)
    if(check==1)Verror=1.d0
    err=sum(Verror)/dble(size(Verror))
    !
    Xold=Xnew
    !
    open(unit,file=str(file_),position="append")
    write(unit,"(I5,ES15.7)")check,err
    close(unit)
    !
    if(err < eps)then
       success=success+1
    else
       success=0
    endif
    convergence=.false.
    !
    if(success > N1)convergence=.true.
    if(check>=N2)then
       open(10,file="ERROR.README")
       write(10,*)""
       close(10)
       write(*,"(A,I4,A)")"Not converged after",N2," iterations."
       convergence=.true.
    endif
    if(convergence)then
       check=0
       deallocate(Xold)
    endif
    !
    if(convergence)then
       write(*,"(A,ES15.7)")bold_green("    error="),err
    else
       if(err < eps)then
          write(*,"(A,ES15.7)")bold_yellow("    error="),err
       else
          write(*,"(A,ES15.7)")bold_red("    error="),err
       endif
    endif
    !
    check=check+1
  end function check_convergence

END MODULE SS_SOLVE_MAIN



