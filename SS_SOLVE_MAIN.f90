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
  real(8),dimension(:,:),allocatable :: ss_zeta_init

contains


  !< Execute the SS calculation using different algorithmes as from the input file
  subroutine ss_solve_methods()
    real(8),dimension(2*Niso)    :: params
    real(8),dimension(Niso)      :: apar
    real(8),dimension(2*Niso+1)  :: params1
    real(8),dimension(Niso+1)    :: apar1
    real(8),dimension(Nineq*Nss) :: lambda,zeta
    integer :: iter
    real(8) :: chi
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
    allocate(ss_lambda_init(Nineq,Nss));ss_lambda_init=ss_lambda_ineq
    allocate(ss_zeta_init(Nineq,Nss))  ;ss_zeta_init  =ss_zeta_ineq
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
    case ("conjg")
       params  = [lambda(:Niso),zeta(:Niso)]
       params1 = [params,xmu] 
       if(filling==0d0)then
          call fmin_cg(params, ss_solve_cg, iter, chi)!, itmax=Niter, ftol=solve_tolerance)
       else
          call fmin_cg(params1,ss_solve_cg, iter, chi)!, itmax=Niter, ftol=solve_tolerance)
       end if
       if(verbose>1)write(*,"(A,ES18.9,1x,I5)")"Chi^2_fss|iter:",chi,"|",iter
       !
    case ("gg_broyden")
       call broyden1(ss_solve_gg,lambda(:Niso),tol=solve_tolerance)
       !
    case ("gg_fsolve")
       call fsolve(ss_solve_gg,lambda(:Niso),tol=solve_tolerance)
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


  function ss_solve_cg(aparams) result(chi2)
    real(8),dimension(:)             :: aparams !2*Nlso[+1]
    real(8),dimension(size(aparams)) :: fss
    real(8)                          :: chi2
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
    chi2 = sum(abs(Fss)**2)/size(Fss)
    print*,chi2
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
  end function ss_solve_cg



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







  function ss_solve_gg(aparams) result(fss)
    real(8),dimension(:),intent(in)  :: aparams
    real(8),dimension(size(aparams)) :: fss
    integer                          :: iter,Nsuccess=0
    logical                          :: z_converged,bool
    real(8),dimension(Nineq*Nss)     :: lambda,zeta,sz,Dens,zeta_prev,Fzeta
    !
    bool = size(lambda)==Niso+1
    !
    ! > Retrieve ss_lambda:
    lambda(:Niso) = aparams(1:Niso)
    if(bool)xmu   = aparams(Niso+1)
    if(Nspin==1)call ss_spin_symmetry(lambda,Nineq)
    ss_lambda_ineq = ss_unpack_array(lambda,Nineq)
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_lambda(ilat,:) = ss_lambda_ineq(ineq,:)
    enddo
    !
    !> Retrieve starting ss_zeta:
    if(zeta_restart_init)ss_zeta_ineq = ss_zeta_init
    do ilat=1,Nlat
       ineq = ss_ilat2ineq(ilat)
       ss_zeta(ilat,:)   = ss_zeta_ineq(ineq,:)
    enddo
    !
    !> Solve:
    if(verbose>2)call start_timer()
    z_converged=.false. ; iter=0
    do while(.not.z_converged.AND.iter<=Niter)
       iter=iter+1
       call start_loop(iter,Niter,"Z-iter")
       !
       zeta_prev = ss_pack_array(ss_zeta_ineq,Nineq)
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
       do ineq=1,Nineq
          ilat = ss_ineq2ilat(ineq)
          ss_Dens_ineq(ineq,:) = ss_Dens(ilat,:)
       enddo
       Dens    = ss_pack_array(ss_Dens_ineq,Nineq)
       Sz      = ss_pack_array(ss_Sz_ineq,Nineq)
       Zeta    = ss_pack_array(ss_zeta_ineq,Nineq)
       !
       call ss_print_screen
       !
       zeta = loop_Wmix*zeta + (1d0-loop_Wmix)*zeta_prev
       !
       ! Fzeta= zeta - zeta_prev
       !
       ! select case(loop_MixType)
       ! case ("linear")            
       !    call linear_mix(zeta,Fzeta,loop_Wmix)
       ! case ("adaptive")
       !    call adaptive_mix(zeta,Fzeta,loop_Wmix,iter)
       ! case default
       !    call broyden_mix(zeta,Fzeta,loop_Wmix,loop_Nmix,iter)
       ! end select
       !
       z_converged = check_convergence(zeta-zeta_prev,loop_tolerance,Nsuccess,Niter)
       call end_loop() 
    end do
    !
    !<constraint:
    fss(1:Niso)         = Dens(1:Niso) - (Sz(1:Niso) + 0.5d0)
    if(bool)fss(Niso+1) = sum(ss_Dens) - filling
    !
    if(verbose>1)then
       write(*,"(A11,50G18.9)")"F_cnstr   =",fss(1:Niso)
       if(bool)write(*,"(A11,G18.9)")"F_filling =",fss(Niso+1)
       write(*,*)""
    endif
    if(verbose>2)call stop_timer()
    !
    call ss_write_last()
    !
  end function ss_solve_gg








  subroutine ss_print_screen
    if(verbose>1)then
       do ineq=1,Nineq
          ilat = ss_ineq2ilat(ineq)
          write(*,*)" SITE= "//str(ilat,4)
          !
          if(verbose>3)write(*,"(A7,12G18.9)")"C     =",ss_C_ineq(ineq,:Nspin*Norb)
          !
          write(*,"(A7,12G18.9)")"mu    =",xmu
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
