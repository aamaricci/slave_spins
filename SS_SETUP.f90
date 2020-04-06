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

  interface ss_user2ss
     module procedure :: ss_user2ss_vec
     module procedure :: ss_user2ss_mat
  end interface ss_user2ss

  interface ss_ss2user
     module procedure :: ss_ss2user_vec
     module procedure :: ss_ss2user_mat
  end interface ss_ss2user


  public :: ss_setup_structure
  public :: ss_init_params
  !
  public :: ss_user2ss
  public :: ss_ss2user
  !
  public :: ss_indices2i
  public :: ss_i2indices
  public :: ss_indx_reorder
  !
  public :: ss_pack_array
  public :: ss_unpack_array
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
    allocate(ss_Weiss(Nlat,Nss), ss_Weiss_ineq(Nineq,Nss))
    allocate(ss_C(Nlat,Nss),ss_C_ineq(Nineq,Nss))
    allocate(ss_Dens(Nlat,Nss), ss_Dens_ineq(Nineq,Nss))
    allocate(ss_Lambda0(Nlat,Nss), ss_Lambda0_ineq(Nineq,Nss))
    allocate( ss_lambda(Nlat,Nss), ss_lambda_ineq(Nineq,Nss) )
    allocate( ss_zeta(Nlat,Nss),   ss_zeta_ineq(Nineq,Nss) )
    allocate( ss_Sz(Nlat,Nss), ss_Sz_ineq(Nineq,Nss))
    allocate( ss_Op(Nlat,Nss), ss_Op_ineq(Nineq,Nss))
    ss_weiss  = 0d0; ss_weiss_ineq  = 0d0
    ss_c      = 0d0; ss_c_ineq      = 0d0
    ss_dens   = 0d0; ss_dens_ineq   = 0d0
    ss_lambda0= 0d0; ss_lambda0_ineq= 0d0
    ss_lambda = 0d0;ss_lambda_ineq = 0d0
    ss_zeta   = 1d0;ss_zeta_ineq   = 1d0
    ss_Sz     = 0d0; ss_Sz_ineq    = 0d0
    ss_Op     = 0d0; ss_Op_ineq    = 0d0
    !
    allocate( ss_SzSz(Nlat,4,Norb,Norb), ss_SzSz_ineq(Nineq,4,Norb,Norb))
    ss_SzSz   = 0d0; ss_SzSz_ineq  = 0d0
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
    !ordering of the degrees of freedom in H(k).
    !THE SS ORDERING IS ASSUMED TO BE ALWAYS [Norb,Nlat,Nspin]
    DefOrder  = [character(len=5) :: "Norb","Nlat","Nspin"]
    nDefOrder = [Norb,Nlat,Nspin]
    !
    Nineq = size(ss_ineq2ilat)
    !
    !the total number of degrees of freedom. Nlat=#inequivalent sites.
    !-TODO: Norb be the total number of orbitals P+D. Now it is only correlated orbitals.
    Nlso = Nlat*Norb*Nspin
    Niso = Nineq*Norb*Nspin
    !
    Ns   = 2*Nlat*Norb          !total number of parameters
    Nss  = 2*Norb
    !
  end subroutine ss_setup_dimensions



  subroutine ss_init_params()
    logical                          :: IOfile
    real(8),dimension(:),allocatable :: params
    integer                          :: Len
    real(8),dimension(Ns)            :: lambda,zeta
    inquire(file=trim(Pfile)//trim(ss_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       len = file_length(trim(Pfile)//trim(ss_file_suffix)//".restart")
       print*,len,2*Ns
       if( (len/=2*Ns) .AND. (len/=2*Ns+1) )stop "SS_INIT_PARAMS ERROR: len!=2*Ns OR 2*Ns+1"
       allocate(params(len))
       call read_array(trim(Pfile)//trim(ss_file_suffix)//".restart",params)
       ss_lambda = ss_unpack_array(params(1:Ns),Nlat)
       ss_zeta   = ss_unpack_array(params(Ns+1:2*Ns),Nlat)
       if(size(params)==2*Ns+1)xmu=params(2*Ns+1)
    else
       ss_lambda = -ss_lambda0
       ss_zeta   = 1d0
    endif
  end subroutine ss_init_params




  function ss_user2ss_vec(Huser,UserOrder) result(Hss)
    real(8),dimension(Nspin*Nlat*Norb) :: Huser
    character(len=*),dimension(3)      :: UserOrder
    real(8),dimension(Nspin*Nlat*Norb) :: Hss
    integer,dimension(3)               :: Ivec,Jvec
    integer,dimension(3)               :: UserIndex
    integer,dimension(3)               :: nUserOrder
    integer                            :: iss,iuser,i
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
       do iss=1,Nlso
          Ivec  = i2indices(iss,nDefOrder)     !get components in default ordering
          Jvec  = indx_reorder(Ivec,UserIndex)  !reorder according to User ordering using UserIndex
          iuser = indices2i(Jvec,nUserOrder)    !get corresponding total index in user ordering
          !
          Hss(iss) = Huser(iuser)
       enddo
    else
       Hss = Huser
    endif
    return
  end function ss_user2ss_vec

  function ss_user2ss_mat(Huser,UserOrder) result(Hss)
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb) :: Huser
    character(len=*),dimension(3)                         :: UserOrder
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb) :: Hss
    integer,dimension(3)                                  :: Ivec,Jvec
    integer,dimension(3)                                  :: UserIndex
    integer,dimension(3)                                  :: nUserOrder
    integer                                               :: iss,jss,iuser,juser,i
    !
    do i=1,3     
       UserIndex(i:i)=findloc(UserOrder,DefOrder(i))
    enddo
    if(any(UserIndex==0))then
       print*,"SS_REORDER_HK ERROR: wrong entry in UserIndex at: ",findloc(UserIndex,0)
       stop
    endif
    !
    nUserOrder=indx_reorder(nDefOrder,UserIndex)
    !
    if(any(UserIndex/=[1,2,3]))then
       do iss=1,Nlso
          Ivec  = i2indices(iss,nDefOrder)
          Jvec  = indx_reorder(Ivec,UserIndex)
          iuser = indices2i(Jvec,nUserOrder)
          do jss=1,Nlso
             Ivec  = i2indices(jss,nDefOrder)
             Jvec  = indx_reorder(Ivec,UserIndex)
             juser = indices2i(Jvec,nUserOrder)
             !
             Hss(iss,jss) = Huser(iuser,juser)
          enddo
       enddo
    else
       Hss = Huser
    endif
    return
  end function ss_user2ss_mat


  function ss_ss2user_vec(Hss,UserOrder) result(Huser)
    real(8),dimension(Nspin*Nlat*Norb) :: Hss
    character(len=*),dimension(3)      :: UserOrder
    real(8),dimension(Nspin*Nlat*Norb) :: Huser
    integer,dimension(3)               :: Ivec,Jvec
    integer,dimension(3)               :: UserIndex
    integer,dimension(3)               :: nUserOrder
    integer                            :: iss,iuser,i
    !
    do i=1,3     
       UserIndex(i:i)=findloc(UserOrder,DefOrder(i))
    enddo
    if(any(UserIndex==0))then
       print*,"SS_REORDER_HK ERROR: wrong entry in UserIndex at: ",findloc(UserIndex,0)
       stop
    endif
    !
    nUserOrder=indx_reorder(nDefOrder,UserIndex)
    !
    if(any(UserIndex/=[1,2,3]))then
       do iss=1,Nlso
          Ivec  = i2indices(iss,nDefOrder)
          Jvec  = indx_reorder(Ivec,UserIndex)
          iuser = indices2i(Jvec,nUserOrder)
          Huser(iuser) = Hss(iss)
       enddo
    else
       Huser = Hss
    endif
    return
  end function ss_ss2user_vec

  function ss_ss2user_mat(Hss,UserOrder) result(Huser)
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb) :: Hss
    character(len=*),dimension(3)                         :: UserOrder
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb) :: Huser
    integer,dimension(3)                                  :: Ivec,Jvec
    integer,dimension(3)                                  :: UserIndex
    integer,dimension(3)                                  :: nUserOrder
    integer                                               :: iss,jss,iuser,juser,i
    !
    do i=1,3     
       UserIndex(i:i)=findloc(UserOrder,DefOrder(i))
    enddo
    if(any(UserIndex==0))then
       print*,"SS_REORDER_HK ERROR: wrong entry in UserIndex at: ",findloc(UserIndex,0)
       stop
    endif
    !
    nUserOrder=indx_reorder(nDefOrder,UserIndex)
    !
    if(any(UserIndex/=[1,2,3]))then
       do iss=1,Nlso
          Ivec  = i2indices(iss,nDefOrder)
          Jvec  = indx_reorder(Ivec,UserIndex)
          iuser = indices2i(Jvec,nUserOrder)
          do jss=1,Nlso
             Ivec  = i2indices(jss,nDefOrder)
             Jvec  = indx_reorder(Ivec,UserIndex)
             juser = indices2i(Jvec,nUserOrder)
             !
             Huser(iuser,juser) = Hss(iss,jss)
          enddo
       enddo
    else
       Huser = Hss
    endif
    return
  end function ss_ss2user_mat





  subroutine ss_spin_symmetry_nlso(array,Nsite)
    integer                      :: Nsite
    real(8),dimension(Nsite*Nss) :: array
    array(Nsite*Norb+1:) = array(1:Nsite*Norb)
  end subroutine ss_spin_symmetry_nlso

  subroutine ss_spin_symmetry_nn(array,Nsite)
    integer                      :: Nsite
    real(8),dimension(Nsite,Nss) :: array
    integer                      :: isite
    do isite=1,Nsite
       array(isite,Norb+1:) = array(isite,1:Norb)
    enddo
  end subroutine ss_spin_symmetry_nn




  function ss_pack_array(Ain,Nsite) result(Aout)
    integer                         :: Nsite
    real(8),dimension(Nsite,2*Norb) :: Ain
    real(8),dimension(2*Nsite*Norb) :: Aout
    integer                         :: ispin,isite,iorb,io,ii
    do ispin=1,2
       do isite=1,Nsite
          do iorb=1,Norb
             ii = Indices2i([iorb,isite,ispin],[Norb,Nsite,Nspin]) !io = iorb + (ispin-1)*Norb
             io = Indices2i([iorb,ispin],[Norb,2])
             !
             Aout(ii) = Ain(isite,io)
             !
          enddo
       enddo
    enddo
  end function ss_pack_array

  function ss_unpack_array(Ain,Nsite) result(Aout)
    integer                         :: Nsite
    real(8),dimension(2*Nsite*Norb) :: Ain
    real(8),dimension(Nsite,2*Norb) :: Aout
    integer                         :: ispin,isite,iorb,io,ii
    do ispin=1,2
       do isite=1,Nsite
          do iorb=1,Norb
             ii = Indices2i([iorb,isite,ispin],[Norb,Nsite,Nspin]) !io = iorb + (ispin-1)*Norb
             io = Indices2i([iorb,ispin],[Norb,2])
             !
             Aout(isite,io) = Ain(ii)
             !
          enddo
       enddo
    enddo
  end function ss_unpack_array



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
