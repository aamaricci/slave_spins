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



  ! !------------------ ABTRACT INTERFACES PROCEDURES ------------------!
  ! !SPARSE MATRIX-VECTOR PRODUCTS USED IN ED_MATVEC
  ! !dbleMat*dbleVec
  ! abstract interface
  !    subroutine dd_sparse_HxV(Nloc,v,Hv)
  !      integer                 :: Nloc
  !      real(8),dimension(Nloc) :: v
  !      real(8),dimension(Nloc) :: Hv
  !    end subroutine dd_sparse_HxV
  ! end interface



  !SIZE OF THE PROBLEM
  !=========================================================
  integer,save                                       :: Ns       !Number of levels per spin
  integer,save                                       :: Nsectors !Number of sectors


  !local part of the Hamiltonian
  !INTERNAL USE (accessed thru functions)
  !=========================================================
  real(8),dimension(:,:,:,:),allocatable             :: impHloc           !local hamiltonian


  !Some maps between sectors and full Hilbert space (pointers)
  !PRIVATE:
  !=========================================================
  integer,allocatable,dimension(:)                   :: getDim             ! [Nsectors]



  ! !Variables for DIAGONALIZATION
  ! !PRIVATE
  ! !=========================================================  
  ! type(sparse_matrix_csr)                            :: spH0d !diagonal part
  ! type(sparse_matrix_csr),dimension(:),allocatable   :: spH0ups,spH0dws !reduced UP and DW parts
  ! !
  ! procedure(dd_sparse_HxV),pointer                   :: spHtimesV_p=>null()





  !Impurity Green's function and Self-Energies: (Nspin,Nspin,Norb,Norb,:)
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impSmats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impSreal
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impGmats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impGreal




  !Density and double occupancy
  !PRIVATE (now public but accessible thru routines)
  !=========================================================
  real(8),dimension(:),allocatable                   ::  ss_dens
  real(8),dimension(:),allocatable                   ::  ss_dens_up,ss_dens_dw
  real(8),dimension(:),allocatable                   ::  ss_docc



  !Local energies and generalized double occupancies
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  real(8)                                            :: ss_Ekin
  real(8)                                            :: ss_Epot
  real(8)                                            :: ss_Eint
  real(8)                                            :: ss_Ehartree
  real(8)                                            :: ss_Eknot





  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                                  :: ss_file_suffix=""       !suffix string attached to the output files.
  logical                                            :: Jhflag              !spin-exchange and pair-hopping flag.
  logical                                            :: offdiag_gf_flag=.false.
  character(len=200)                                 :: ss_input_file=""


  !   !This is the internal Mpi Communicator and variables.
  !   !=========================================================
  ! #ifdef _MPI
  !   integer                                            :: MpiComm_Global=MPI_COMM_NULL
  !   integer                                            :: MpiComm=MPI_COMM_NULL
  ! #endif
  !   integer                                            :: MpiGroup_Global=MPI_GROUP_NULL
  !   integer                                            :: MpiGroup=MPI_GROUP_NULL
  !   logical                                            :: MpiStatus=.false.
  !   logical                                            :: MpiMaster=.true.
  !   integer                                            :: MpiRank=0
  !   integer                                            :: MpiSize=1
  !   integer,allocatable,dimension(:)                   :: MpiMembers
  !   integer                                            :: mpiQup=0
  !   integer                                            :: mpiRup=0
  !   integer                                            :: mpiQdw=0
  !   integer                                            :: mpiRdw=0
  !   integer                                            :: mpiQ=0
  !   integer                                            :: mpiR=0
  !   integer                                            :: mpiIstart
  !   integer                                            :: mpiIend
  !   integer                                            :: mpiIshift
  !   logical                                            :: mpiAllThreads=.true.



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
