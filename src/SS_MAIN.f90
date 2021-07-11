MODULE SS_MAIN
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  USE SS_SOLVE_FERMION
  USE SS_SOLVE_MAIN
  !
  USE SF_PAULI
  USE SF_LINALG, only: kron
  USE SF_TIMER,  only: start_timer,stop_timer
  !
  implicit none
  private

  interface ss_solve
     module procedure :: ss_solve_hk
     module procedure :: ss_solve_dos
  end interface ss_solve


  public :: ss_solve

  integer,save :: siter=0
  integer      :: fiter,info
  logical      :: fconverged
  integer      :: iorb,jorb,ispin,jspin,ilat,ineq,io,il


contains



  !< Init SS calculation by passing the Hamiltonian H(k)
  subroutine ss_solve_hk(hk_user,ineq_sites,UserOrder)
    complex(8),dimension(:,:,:)                           :: hk_user  ![Nlso,Nlso,Nk]
    integer,dimension(Nlat),optional                      :: ineq_sites
    character(len=5),dimension(3),optional                :: UserOrder
    !
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb) :: Htmp
    integer                                               :: iorb,jorb,ispin,io,jo,ilat,jlat,ineq,ii,jj

    integer                                               :: ik
    character(len=5),dimension(3)                         :: UserOrder_
    !
#ifdef _MPI
    if(check_MPI())master = get_master_MPI()
#endif
    !
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

    !< Init the Hk structures, subtract the local, diagonal part, coupled to the density
    !< First: if order of Hk_user is not correct set the SS_order function to actual reorder
    if(present(UserOrder))then
       write(*,"(A)")"SS Reordering according to provided: UserOrder"
       do ik=1,Nk
          ss_Hk(:,:,ik) = ss_user2ss(Hk_user(:,:,ik),UserOrder_)
       end do
    endif
    where(abs(ss_Hk)<1d-6)ss_Hk=zero
    ss_Wtk(1,:) = 1d0/Nk
    !
    Htmp = sum(ss_Hk,dim=3)/Nk
    !
    !< get the diagonal part
    ss_Hdiag = diagonal(Htmp)   
    !
    do ispin=1,Nspin
       do ilat=1,Nlat
          do iorb=1,Norb
             io = ss_indices2i([iorb,ilat,ispin],[Norb,Nlat,Nspin])
             !
             do jspin=1,Nspin
                do jlat=1,Nlat
                   do jorb=1,Norb
                      jo = ss_indices2i([jorb,jlat,jspin],[Norb,Nlat,Nspin])
                      !
                      if(ilat/=jlat)then
                         !< get the inter-site part of the local Hamiltonian
                         ss_Hloc(io,jo) = Htmp(io,jo)
                      else
                         !< get the off-diagonal part local Hamiltonian
                         !io/=jo:
                         ![iorb/=jorb,ispin==jspin], [iorb==jorb,ispin/=jspin], [iorb/=jorb,ispin/=jspin]
                         if(io/=jo)ss_Hhyb(io,jo) = Htmp(io,jo)
                      endif
                      !
                   enddo
                enddo
             enddo
             !
          enddo
       enddo
    enddo
    !
    !< remove local parts from H(k)
    forall(ik=1:Nk)ss_Hk(:,:,ik) = ss_Hk(:,:,ik) - diag(ss_Hdiag) - ss_Hloc - ss_Hhyb 
    !
    !
    !< Init/Read the lambda input
    if(use_lam0)call ss_solve_lambda0()
    call ss_init_params()
    !
    if(master.AND.verbose>2)then
       do ineq=1,Nineq
          ilat = ss_ineq2ilat(ineq)
          write(*,"(A7,12G18.9)")"Lam   =",ss_lambda(ilat,:)
          write(*,"(A7,12G18.9)")"Op    =",ss_Op(ilat,:)
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
  subroutine ss_solve_dos(Ebands,Dbands,Hloc,ineq_sites,UserOrder)
    real(8),dimension(:,:)                                :: Ebands  ![Nlso,Ne]
    real(8),dimension(:,:)                                :: Dbands ![Nlso,Ne]
    real(8),dimension(Nspin*Nlat*Norb),optional           :: Hloc
    integer,dimension(Nlat),optional                      :: ineq_sites
    character(len=5),dimension(3),optional                :: UserOrder
    !
    real(8),dimension(Nspin*Nlat*Norb)                    :: Eb,Db
    integer                                               :: ie,io,ilat,ineq
    character(len=5),dimension(3)                         :: UserOrder_
    !
#ifdef _MPI
    if(check_MPI())master = get_master_MPI()
#endif
    !
    if(Nspin==2)then
       write(*,*)"WARNING: Nspin=2 with DOS interface. SOC is not allowed here yet. Ask Developers."
       call sleep(2)
    endif
    !
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
    ! if(present(Hloc))call assert_shape(Hloc,[Nspin*Nlat*Norb],"ss_init_dos","Hloc")
    !
    !< Init local non-interacting part and reorder it
    if(present(Hloc))ss_Hdiag = ss_user2ss(Hloc,UserOrder_)
    !
    !if required re-order the bands:
    if(present(UserOrder))write(*,"(A)")"SS warning: reordering according to provided: UserOrder"
    do ie=1,Nk
       if(present(UserOrder))then
          Eb = ss_user2ss(Ebands(:,ie),UserOrder_)
       else
          Eb = Ebands(:,ie)
       endif
       do io=1,Nspin*Nlat*Norb
          ss_Hk(io,io,ie)  = one*Eb(io)
       end do
       if(present(UserOrder))then
          ss_Wtk(:,ie) = ss_user2ss(Dbands(:,ie),UserOrder_)
       else
          ss_Wtk(:,ie) = Dbands(:,ie)
       endif
    end do
    !
    !
    !< Init/Read the lambda input
    ! if(filling/=dble(Norb))call ss_solve_lambda0()
    if(use_lam0)call ss_solve_lambda0()
    call ss_init_params()
    !
    if(master.AND.verbose>2)then
       do ineq=1,Nineq
          ilat = ss_ineq2ilat(ineq)
          write(*,"(A7,12G18.9)")"Lam   =",ss_lambda(ilat,:)
          write(*,"(A7,12G18.9)")"Op    =",ss_Op(ilat,:)
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
       if( (len/=2*Nlso) .AND. (len/=2*Nlso+1) )stop "SS_INIT_PARAMS ERROR: len!=2*Nso OR 2*Nso+1"
       allocate(params(len))
       call read_array(trim(Pfile)//trim(ss_file_suffix)//".restart",params)
       ss_lambda = ss_unpack_array(params(1:Nlso),Nlat)
       ss_op     = ss_unpack_array(params(Nlso+1:2*Nlso),Nlat)
       if(size(params)==2*Nlso+1)xmu=params(2*Nlso+1)
    else
       ss_lambda =  ss_lambda0
       ss_Op     =  1d0
    endif
  end subroutine ss_init_params


END MODULE SS_MAIN
