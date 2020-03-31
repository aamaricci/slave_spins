MODULE SS_SETUP
  USE SS_VARS_GLOBAL
  USE SCIFOR, only: uniq_array,arange
  !
  implicit none
  private

  interface ss_indices2i
     module procedure :: indices2i
  end interface ss_indices2i

  interface ss_i2indices
     module procedure :: i2indices
  end interface ss_i2indices

  interface ss_indx_reorder
     module procedure :: indx_reorder
  end interface ss_indx_reorder

  interface ss_spin_symmetry
     module procedure :: ss_spin_symmetry_nlso
     module procedure :: ss_spin_symmetry_nn
  end interface ss_spin_symmetry


  public :: ss_setup_structure
  public :: ss_init_params
  !
  public :: ss_reorder_hk
  public :: ss_reorder_bands
  !
  public :: ss_indices2i
  public :: ss_i2indices
  public :: ss_indx_reorder
  public :: ss_spread_array
  public :: ss_spin_symmetry


contains


  subroutine ss_setup_structure(ineq_sites)
    integer,dimension(Nlat),optional :: ineq_sites
    !
    if(Nspin>2)stop "ss_setup error: Nspin>2"
    if(Norb>5)stop "ss_setup error: Norb>5. Fix me at SS_SETUP.f90/line.25"
    !
    allocate(ss_ilat2ineq(Nlat))
    ss_ilat2ineq=arange(1,Nlat)                    !default list
    if(present(ineq_sites))&                       !
         ss_ilat2ineq=ineq_sites                   !if any use user provided list 
    call uniq_array(ss_ilat2ineq,ss_ineq2ilat)     !allocate and build ineq2ilat list\
    !
    call ss_setup_dimensions()
    !
    allocate(ss_lambda(Ns))
    allocate(ss_lambda0(Ns))
    allocate(ss_weiss(Ns))
    allocate(ss_c(Ns))
    allocate(ss_dens(Ns))
    ss_lambda = 0d0
    ss_lambda0= 0d0
    ss_weiss  = 0d0
    ss_c      = 0d0
    ss_dens   = 0d0
    !
    allocate(ss_zeta(Ns))
    ss_zeta   = 1d0
    !
    allocate(ss_Sz(Ns), ss_Sz_ineq(Nineq,Nss))
    allocate(ss_Op(Ns), ss_Op_ineq(Nineq,Nss))
    allocate(ss_SzSz(Nlat,4,Norb,Norb), ss_SzSz_ineq(Nineq,4,Norb,Norb))
    ss_Sz     = 0d0; ss_Sz_ineq   = 0d0
    ss_Op     = 0d0; ss_Op_ineq   = 0d0
    ss_SzSz   = 0d0; ss_SzSz_ineq = 0d0
    !
    allocate(ss_Hk(Ns,Ns,Nk))
    ss_Hk = zero
    !
    allocate(ss_Wtk(Ns,Ns,Nk))
    ss_Wtk= 0d0
    !
    allocate(ss_Hloc(Ns,Ns))
    ss_Hloc = zero
  end subroutine ss_setup_structure


  subroutine ss_setup_dimensions()
    !the total number of degrees of freedom. Nlat=#inequivalent sites.
    !-TODO: Norb be the total number of orbitals P+D. Now it is only correlated orbitals.
    Nlso = Nlat*Norb*Nspin
    !ordering of the degrees of freedom in H(k).
    !THE SS ORDERING IS ASSUMED TO BE ALWAYS [Norb,Nlat,Nspin]
    DefOrder  = [character(len=5) :: "Norb","Nlat","Nspin"]
    nDefOrder = [Norb,Nlat,Nspin]
    !
    Nineq = size(ss_ineq2ilat)
    !
    Ns   = 2*Nlat*Norb          !total number of parameters
    Nss  = 2*Norb
    !
  end subroutine ss_setup_dimensions



  subroutine ss_init_params()
    logical                          :: IOfile
    real(8),dimension(:),allocatable :: params
    integer                          :: Len
    inquire(file=trim(Pfile)//trim(ss_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       len = file_length(trim(Pfile)//trim(ss_file_suffix)//".restart")
       allocate(params(len))
       call read_array(trim(Pfile)//trim(ss_file_suffix)//".restart",params)
       ss_lambda = params(1:Ns)
       ss_zeta   = params(Ns+1:2*Ns)
       if(size(params)==2*Ns+1)xmu=params(2*Ns+1)
    else
       ss_lambda = -ss_lambda0
       ss_zeta   = 1d0
    endif
  end subroutine ss_init_params




  function ss_reorder_hk(Huser,UserOrder) result(Hord)
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb) :: Huser
    character(len=*),dimension(3)                         :: UserOrder
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb) :: Hord
    integer,dimension(3)                                  :: Ivec,Jvec
    integer,dimension(3)                                  :: UserIndex
    integer,dimension(3)                                  :: nUserOrder
    integer                                               :: iord,jord,iuser,juser,i
    !
    !Construct an index array corresponding to the User ordering.
    !This is a permutation of the default ordering [1,2,3].
    !For each entry in Default Order we look for the position of the
    !corresponding entry in User Order using Fortran findloc.
    !If 0 entries exist, corresponding components are not found. stop. 
    do i=1,3     
       UserIndex(i:i)=findloc(UserOrder,DefOrder(i))
    enddo
    if(any(UserIndex==0))then
       print*,"SS_REORDER_HK ERROR: wrong entry in UserIndex at: ",findloc(UserIndex,0)
       stop
    endif
    !
    !From UserIndex we can re-order the dimensions array to get the User dimensions array 
    nUserOrder=indx_reorder(nDefOrder,UserIndex)
    !
    if(any(UserIndex/=[1,2,3]))then
       !Re-order
       do iord=1,Nlso
          Ivec  = i2indices(iord,nDefOrder)     !get components in default ordering
          Jvec  = indx_reorder(Ivec,UserIndex)      !reorder according to User ordering using UserIndex
          iuser = indices2i(Jvec,nUserOrder) !get corresponding total index in user ordering
          do jord=1,Nlso
             Ivec  = i2indices(jord,nDefOrder)
             Jvec  = indx_reorder(Ivec,UserIndex)
             juser = indices2i(Jvec,nUserOrder)
             !
             Hord(iord,jord) = Huser(iuser,juser)
             !
          enddo
       enddo
       !
    else
       !copy that
       Hord = Huser
    endif
    return
  end function ss_reorder_hk


  function ss_reorder_bands(Huser,UserOrder) result(Hord)
    real(8),dimension(Nspin*Nlat*Norb) :: Huser
    character(len=*),dimension(3)      :: UserOrder
    real(8),dimension(Nspin*Nlat*Norb) :: Hord
    integer,dimension(3)               :: Ivec,Jvec
    integer,dimension(3)               :: UserIndex
    integer,dimension(3)               :: nUserOrder
    integer                            :: iord,iuser,i
    !
    !Construct an index array corresponding to the User ordering.
    !This is a permutation of the default ordering [1,2,3].
    !For each entry in Default Order we look for the position of the
    !corresponding entry in User Order using Fortran findloc.
    !If 0 entries exist, corresponding components are not found. stop. 
    do i=1,3     
       UserIndex(i:i)=findloc(UserOrder,DefOrder(i))
    enddo
    if(any(UserIndex==0))then
       print*,"SS_REORDER_HK ERROR: wrong entry in UserIndex at: ",findloc(UserIndex,0)
       stop
    endif
    !
    !From UserIndex we can re-order the dimensions array to get the User dimensions array 
    nUserOrder=indx_reorder(nDefOrder,UserIndex)
    !
    if(any(UserIndex/=[1,2,3]))then
       !Re-order
       do iord=1,Nlso
          Ivec  = i2indices(iord,nDefOrder)     !get components in default ordering
          Jvec  = indx_reorder(Ivec,UserIndex)      !reorder according to User ordering using UserIndex
          iuser = indices2i(Jvec,nUserOrder) !get corresponding total index in user ordering
          !
          Hord(iord) = Huser(iuser)
          !
       enddo
    else
       !Copy
       Hord = Huser
    endif
    return
  end function ss_reorder_bands





  subroutine ss_spin_symmetry_nlso(array)
    real(8),dimension(Ns)    :: array
    array(Nlat*Norb+1:) = array(1:Nlat*Norb)
  end subroutine ss_spin_symmetry_nlso

  subroutine ss_spin_symmetry_nn(array)
    real(8),dimension(Nlat,Nss) :: array
    integer                     :: ilat
    do ilat=1,Nlat
       array(ilat,Norb+1:) = array(ilat,1:Norb)
    enddo
  end subroutine ss_spin_symmetry_nn



  subroutine ss_spread_array(Ain,Aout)
    real(8),dimension(2*Nlat*Norb) :: Ain
    real(8),dimension(Nlat,2*Norb) :: Aout
    integer                        :: ispin,ilat,iorb,io,ii
    do ispin=1,2
       do ilat=1,Nlat
          do iorb=1,Norb
             ii = Indices2i([iorb,ilat,ispin],nDefOrder) !io = iorb + (ispin-1)*Norb
             io = Indices2i([iorb,ispin],[Norb,2])
             !
             Aout(ilat,io) = Ain(ii)
             !
          enddo
       enddo
    enddo
  end subroutine ss_spread_array



  function indx_reorder(Ain,Index)  result(Aout)
    integer,dimension(:)         :: Ain
    integer,dimension(size(Ain)) :: Index
    integer,dimension(size(Ain)) :: Aout
    integer                        :: i
    forall(i=1:size(Ain))Aout(Index(i)) = Ain(i)
  end function indx_reorder


  function indices2i(ivec,Nvec) result(istate)
    integer,dimension(:)          :: ivec
    integer,dimension(size(ivec)) :: Nvec
    integer                       :: istate,i
    istate=ivec(1)
    do i=2,size(ivec)
       istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
    enddo
  end function indices2i

  function i2indices(istate,Nvec) result(ivec)
    integer                       :: istate
    integer,dimension(:)          :: Nvec
    integer,dimension(size(Nvec)) :: Ivec
    integer                       :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end function i2indices




END MODULE SS_SETUP
