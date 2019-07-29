MODULE SS_VARS_GLOBAL
  USE SF_CONSTANTS
  implicit none





  ! !---------------- SECTOR-TO-FOCK SPACE STRUCTURE -------------------!
  ! type sector_map
  !    integer,dimension(:),allocatable :: map
  !    logical                          :: status=.false.
  ! end type sector_map

  ! interface map_allocate
  !    module procedure :: map_allocate_scalar
  !    module procedure :: map_allocate_vector
  ! end interface map_allocate

  ! interface map_deallocate
  !    module procedure :: map_deallocate_scalar
  !    module procedure :: map_deallocate_vector
  ! end interface map_deallocate



  !------------------ ABTRACT INTERFACES PROCEDURES ------------------!
  !SPARSE MATRIX-VECTOR PRODUCTS USED IN ED_MATVEC
  !dbleMat*dbleVec
  abstract interface
     subroutine dd_sparse_HxV(Nloc,v,Hv)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: v
       real(8),dimension(Nloc) :: Hv
     end subroutine dd_sparse_HxV
  end interface






  integer,save                            :: Ns       !Number of levels = 2*Norb
  integer                                 :: Nk
  !
  real(8),dimension(:),allocatable        :: ss_lambda    !lambda_{m\sigma}
  real(8),dimension(:),allocatable        :: ss_zeta      !zeta_{m\sigma}
  real(8),dimension(:),allocatable        :: ss_weiss     !4<S^x_{m\sigma}><E_fermion>
  real(8),dimension(:),allocatable        :: ss_c
  real(8),dimension(:),allocatable        :: ss_Sz
  complex(8),dimension(:,:,:),allocatable :: ss_Hk
  real(8),dimension(:),allocatable        :: ss_Wtk
  real(8),dimension(:,:),allocatable      :: ss_Hloc           !local hamiltonian
  real(8),dimension(:),allocatable        :: ss_dens ![Ns: 1:Norb_up, 1:Norb_dw]
  logical                                 :: ss_Hdiag !



  real(8) :: zeta_function

  
  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                      :: ss_file_suffix=""       !suffix string attached to the output files.
  character(len=200)                     :: ss_input_file=""


  ! contains


  !   !=========================================================
  !   subroutine map_allocate_scalar(H,N)
  !     type(sector_map) :: H
  !     integer          :: N
  !     if(H%status) call map_deallocate_scalar(H)
  !     allocate(H%map(N))
  !     H%status=.true.
  !   end subroutine map_allocate_scalar
  !   !
  !   subroutine map_allocate_vector(H,N)
  !     type(sector_map),dimension(:)       :: H
  !     integer,dimension(size(H))          :: N
  !     integer                             :: i
  !     do i=1,size(H)
  !        call map_allocate_scalar(H(i),N(i))
  !     enddo
  !   end subroutine map_allocate_vector


  !   !=========================================================
  !   subroutine map_deallocate_scalar(H)
  !     type(sector_map) :: H
  !     if(.not.H%status)then
  !        write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
  !        return
  !     endif
  !     if(allocated(H%map))deallocate(H%map)
  !     H%status=.false.
  !   end subroutine map_deallocate_scalar
  !   !
  !   subroutine map_deallocate_vector(H)
  !     type(sector_map),dimension(:) :: H
  !     integer                       :: i
  !     do i=1,size(H)
  !        call map_deallocate_scalar(H(i))
  !     enddo
  !   end subroutine map_deallocate_vector




  !   !=========================================================
  !   subroutine ed_set_MpiComm(comm)
  ! #ifdef _MPI
  !     integer :: comm,ierr
  !     ! call MPI_Comm_dup(Comm,MpiComm_Global,ierr)
  !     ! call MPI_Comm_dup(Comm,MpiComm,ierr)
  !     MpiComm_Global = comm
  !     MpiComm        = comm
  !     call Mpi_Comm_group(MpiComm_Global,MpiGroup_Global,ierr)
  !     MpiStatus      = .true.
  !     MpiSize        = get_Size_MPI(MpiComm_Global)
  !     MpiRank        = get_Rank_MPI(MpiComm_Global)
  !     MpiMaster      = get_Master_MPI(MpiComm_Global)
  ! #else
  !     integer,optional :: comm
  ! #endif
  !   end subroutine ed_set_MpiComm

  !   subroutine ed_del_MpiComm()
  ! #ifdef _MPI    
  !     MpiComm_Global = MPI_UNDEFINED
  !     MpiComm        = MPI_UNDEFINED
  !     MpiGroup_Global= MPI_GROUP_NULL
  !     MpiStatus      = .false.
  !     MpiSize        = 1
  !     MpiRank        = 0
  !     MpiMaster      = .true.
  ! #endif
  !   end subroutine ed_del_MpiComm



END MODULE SS_VARS_GLOBAL
