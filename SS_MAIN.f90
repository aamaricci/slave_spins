MODULE SS_MAIN
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_MAIN
  !
  USE SF_PAULI
  USE SF_LINALG, only: kron
  USE SF_TIMER, only: start_timer,stop_timer
  !
  implicit none
  private

  interface ss_solve
     module procedure :: ss_solve_hk
     module procedure :: ss_solve_dos
  end interface ss_solve


  public :: ss_solve

  integer,save                     :: siter=0
  integer                          :: fiter,info
  logical                          :: fconverged
  integer                          :: iorb,ispin,ilat,ineq,io,il


contains



  !< Init SS calculation by passing the Hamiltonian H(k)
  subroutine ss_solve_hk(hk_user,UserOrder,Hloc,ineq_sites)
    complex(8),dimension(:,:,:)                           :: hk_user  ![Nlso,Nlso,Nk]
    character(len=*),dimension(3),optional                :: UserOrder
    real(8),dimension(Nspin*Nlat*Norb),optional           :: Hloc
    integer,dimension(Nlat),optional                      :: ineq_sites
    !
    real(8),dimension(Nspin*Nlat*Norb)                    :: Hdiag,Haux
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb) :: Htmp,Hk
    complex(8),dimension(:,:),allocatable                 :: Hcheck
    integer                                               :: ik
    character(len=5),dimension(3)                         :: UserOrder_
    !
#ifdef _MPI
    if(check_MPI())master = get_master_MPI()
#endif
    !
    UserOrder_ = [character(len=5) :: "Norb","Nspin","Nlat"];
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Nk = size(hk_user,3)
    !
    !< Init SS parameters:
    if(present(ineq_sites))then
       call ss_setup_structure(ineq_sites)
    else
       call ss_setup_structure()
    endif
    !
    !< Init the SS structure + memory allocation
    call assert_shape(hk_user,[Nlso,Nlso,Nk],"ss_init_hk","hk_user")
    !
    !< Init local non-interacting part and reorder it
    if(present(Hloc))then
       Hdiag = ss_user2ss(Hloc,UserOrder_) 
    else
       Haux  = diagonal( sum(Hk_user,dim=3)/Nk )
       Hdiag = ss_user2ss(Haux,UserOrder_ )
       where(abs(Hdiag)<1d-6)Hdiag=0d0
    endif
    !
    ss_Hdiag(:Nlso) = Hdiag
    if(Nspin==1)call ss_spin_symmetry(ss_Hdiag,Nlat)
    !
    !< Init the Hk structures, subtract the local, diagonal part, coupled to the density 
    do ik=1,Nk
       !< if order of Hk_user is not correct set the SS_order function to actual reorder
       Hk = ss_user2ss(Hk_user(:,:,ik),UserOrder_)
       !
       select case(Nspin)
       case default
          ss_Hk(:,:,ik)  = kron(pauli_0,Hk) - diag(ss_Hdiag)
       case (2)
          ss_Hk(:,:,ik)  = Hk - diag(ss_Hdiag)
       end select
    end do
    ss_Wtk(1,:) = 1d0/Nk
    !
    !
    !< Init/Read the lambda input
    if(filling/=0d0)call ss_solve_lambda0()
    call ss_init_params()
    !
    if(master.AND.verbose>2)then
       do ineq=1,Nineq
          write(*,"(A7,12G18.9)")"Lam0  =",ss_lambda0(ineq,:)
          write(*,"(A7,12G18.9)")"Lam   =",ss_lambda(ineq,:)
          write(*,"(A7,12G18.9)")"Op    =",ss_Op(ineq,:)
       enddo
       write(*,"(A7,12G18.9)")"mu    =",xmu
       write(*,*)" "
    endif
    !
    call ss_solve_methods()
    !
    return
  end subroutine ss_solve_hk



  !< Init SS calculation by passing the Density of states D(e), the dispersions E(e) and the local
  !  part of the Hamiltonian H_loc = sum_k H(k)
  subroutine ss_solve_dos(Ebands,Dbands,UserOrder,Hloc,ineq_sites)
    real(8),dimension(:,:)                                :: Ebands  ![Nlso,Ne]
    real(8),dimension(:,:)                                :: Dbands ![Nlso,Ne]
    character(len=*),dimension(3),optional                :: UserOrder
    real(8),dimension(:),optional                         :: Hloc
    integer,dimension(Nlat),optional                      :: ineq_sites
    !
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb) :: Htmp
    real(8),dimension(Nspin*Nlat*Norb)                    :: Wtmp
    real(8),dimension(Nspin*Nlat*Norb)                    :: Eb,Db,Hdiag
    integer                                               :: ie,io,ilat
    character(len=5),dimension(3)                         :: UserOrder_
    !
#ifdef _MPI
    if(check_MPI())master = get_master_MPI()
#endif
    !
    UserOrder_ = [character(len=5) :: "Norb","Nspin","Nlat"];
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Nk = size(Ebands,2)
    !
    is_dos=.true.
    !
    if(present(ineq_sites))then
       call ss_setup_structure(ineq_sites)
    else
       call ss_setup_structure()
    endif
    !
    !< Init the SS structure + memory allocation
    call assert_shape(Ebands,[Nspin*Nlat*Norb,Nk],"ss_init_dos","Ebands")
    call assert_shape(Dbands,[Nspin*Nlat*Norb,Nk],"ss_init_dos","Dbands")
    if(present(Hloc))call assert_shape(Hloc,[Nspin*Nlat*Norb],"ss_init_dos","Hloc")
    !
    !< Init local non-interacting part and reorder it
    Hdiag = 0d0
    if(present(Hloc))Hdiag = ss_user2ss(Hloc,UserOrder_)
    ss_Hdiag(:Nlso) = Hdiag
    if(Nspin==1)call ss_spin_symmetry(ss_Hdiag,Nlat)
    !
    ss_Hk = zero
    ss_Wtk= 0d0
    do ie=1,Nk
       Eb = ss_user2ss(Ebands(:,ie),UserOrder_)
       Db = ss_user2ss(Dbands(:,ie),UserOrder_)
       Htmp=zero
       do io=1,Nspin*Nlat*Norb
          Htmp(io,io)  = one*Eb(io)
       end do
       !
       select case(Nspin)
       case default
          ss_Hk(:,:,ie)  = kron(pauli_0,Htmp)
          ss_Wtk(:,ie)   = [Db,Db]
       case (2)
          ss_Hk(:,:,ie)  = Htmp
          ss_Wtk(:,ie)   = Db
       end select
    end do
    !
    !
    !< Init/Read the lambda input
    ! if(filling/=dble(Norb))call ss_solve_lambda0()
    if(filling/=0d0)call ss_solve_lambda0()
    call ss_init_params()
    !
    if(master.AND.verbose>2)then
       do ineq=1,Nineq
          write(*,"(A7,12G18.9)")"Lam0  =",ss_lambda0(ineq,:)
          write(*,"(A7,12G18.9)")"Lam   =",ss_lambda(ineq,:)
          write(*,"(A7,12G18.9)")"Op    =",ss_Op(ineq,:)
       enddo
       write(*,"(A7,12G18.9)")"mu    =",xmu
       write(*,*)" "
    endif
    !
    call ss_solve_methods()
    !
    return
  end subroutine ss_solve_dos





  subroutine ss_init_params()
    logical                          :: IOfile
    real(8),dimension(:),allocatable :: params
    integer                          :: Len
    inquire(file=trim(Pfile)//trim(ss_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       len = file_length(trim(Pfile)//trim(ss_file_suffix)//".restart")
       if( (len/=2*Ns) .AND. (len/=2*Ns+1) )stop "SS_INIT_PARAMS ERROR: len!=2*Ns OR 2*Ns+1"
       allocate(params(len))
       call read_array(trim(Pfile)//trim(ss_file_suffix)//".restart",params)
       ss_lambda = ss_unpack_array(params(1:Ns),Nlat)
       ss_op     = ss_unpack_array(params(Ns+1:2*Ns),Nlat)
       if(size(params)==2*Ns+1)xmu=params(2*Ns+1)
    else
       ss_lambda = -ss_lambda0
       ss_Op     = 1d0
    endif
  end subroutine ss_init_params


END MODULE SS_MAIN
