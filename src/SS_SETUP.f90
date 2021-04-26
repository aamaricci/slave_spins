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
     module procedure :: ss_spin_symmetry_nlso_d
     module procedure :: ss_spin_symmetry_nlso_c
     module procedure :: ss_spin_symmetry_nn_d
     module procedure :: ss_spin_symmetry_nn_c
  end interface ss_spin_symmetry

  interface ss_user2ss
     module procedure :: ss_user2ss_vec
     module procedure :: ss_user2ss_mat
  end interface ss_user2ss

  interface ss_ss2user
     module procedure :: ss_ss2user_vec
     module procedure :: ss_ss2user_mat
  end interface ss_ss2user

  interface ss_pack_array
     module procedure :: ss_pack_array_d
     module procedure :: ss_pack_array_c
     module procedure :: ss_pack_array2_d
     module procedure :: ss_pack_array2_c
  end interface ss_pack_array

  interface ss_unpack_array
     module procedure :: ss_unpack_array_d
     module procedure :: ss_unpack_array_c
     module procedure :: ss_unpack_array2_d
     module procedure :: ss_unpack_array2_c
  end interface ss_unpack_array


  public :: ss_setup_structure
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


  integer :: ispin,isite,jsite,iorb,jorb,io,jo,ii,jj

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
    !
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
    Nso  = Norb*Nspin
    !
    Ns   = 2*Nlat*Norb          !total number of parameters
    Nss  = 2*Norb
    !        
    !
    allocate(ss_C(Nlat,Nso),ss_C_ineq(Nineq,Nso))
    allocate(ss_Dens(Nlat,Nso), ss_Dens_ineq(Nineq,Nso))
    allocate(ss_Lambda0(Nlat,Nso), ss_Lambda0_ineq(Nineq,Nso))
    allocate(ss_lambda(Nlat,Nso), ss_lambda_ineq(Nineq,Nso) )
    allocate(ss_Sz(Nlat,Nso), ss_Sz_ineq(Nineq,Nso))
    !
    allocate(ss_Op(Nlat,Nso), ss_Op_ineq(Nineq,Nso))
    allocate(ss_OdgOp(Nlat,Nso,Nso), ss_OdgOp_ineq(Nineq,Nso,Nso))
    !
    allocate(ss_Heff(Nlat,Nso), ss_Heff_ineq(Nineq,Nso))
    allocate(ss_Jhybr(Nlat,Nso,Nso), ss_Jhybr_ineq(Nineq,Nso,Nso))
    !
    allocate( ss_SzSz(Nlat,4,Norb,Norb), ss_SzSz_ineq(Nineq,4,Norb,Norb))
    ss_c      = 0d0; ss_c_ineq      = 0d0
    ss_dens   = 0d0; ss_dens_ineq   = 0d0
    ss_lambda0= 0d0; ss_lambda0_ineq= 0d0
    ss_lambda = 0d0; ss_lambda_ineq = 0d0
    ss_Sz     = 0d0; ss_Sz_ineq     = 0d0
    !
    ss_Op     = 1d0; ss_Op_ineq     = 1d0
    ss_OdgOp  = 0d0; ss_OdgOp_ineq  = 0d0
    !
    ss_Heff   = 0d0; ss_Heff_ineq   = 0d0
    ss_Jhybr  = 0d0; ss_Jhybr_ineq  = 0d0
    !
    ss_SzSz   = 0d0; ss_SzSz_ineq   = 0d0
    !
    !
    allocate(ss_Hk(Nlso,Nlso,Nk))
    allocate(ss_Hloc(Nlso,Nlso))
    allocate(ss_Hhyb(Nlso,Nlso))
    allocate(ss_Hdiag(Nlso))
    ss_Hk    = zero
    ss_Hloc  = zero
    ss_Hhyb  = zero
    ss_Hdiag = 0d0
    !
    !
    if(is_dos)then
       allocate(ss_Wtk(Nlso,Nk))
    else
       allocate(ss_Wtk(1,Nk))
    end if
    ss_Wtk= 0d0
    !
  end subroutine ss_setup_structure






  function ss_pack_array_d(Ain,Nsite) result(Aout)
    integer                             :: Nsite
    real(8),dimension(Nsite,Nso)        :: Ain
    real(8),dimension(Nspin*Nsite*Norb) :: Aout
    do ispin=1,Nspin
       do isite=1,Nsite
          do iorb=1,Norb
             ii = Indices2i([iorb,isite,ispin],[Norb,Nsite,Nspin]) !io = iorb + (ispin-1)*Norb
             io = Indices2i([iorb,ispin],[Norb,Nspin])
             !
             Aout(ii) = Ain(isite,io)
             !
          enddo
       enddo
    enddo
  end function ss_pack_array_d
  function ss_pack_array_c(Ain,Nsite) result(Aout)
    integer                                :: Nsite
    complex(8),dimension(Nsite,Nso)        :: Ain
    complex(8),dimension(Nspin*Nsite*Norb) :: Aout
    do ispin=1,Nspin
       do isite=1,Nsite
          do iorb=1,Norb
             ii = Indices2i([iorb,isite,ispin],[Norb,Nsite,Nspin])
             io = Indices2i([iorb,ispin],[Norb,Nspin])
             !
             Aout(ii) = Ain(isite,io)
             !
          enddo
       enddo
    enddo
  end function ss_pack_array_c
  function ss_pack_array2_d(Ain,Nsite) result(Aout)
    integer                                              :: Nsite
    real(8),dimension(Nsite,Nso,Nso)                     :: Ain
    real(8),dimension(Nspin*Nsite*Norb,Nspin*Nsite*Norb) :: Aout
    do ispin=1,Nspin
       do isite=1,Nsite
          do iorb=1,Norb
             ii = Indices2i([iorb,isite,ispin],[Norb,Nsite,Nspin])
             io = Indices2i([iorb,ispin],[Norb,Nspin])
             do jorb=1,Norb
                jj = Indices2i([jorb,isite,ispin],[Norb,Nsite,Nspin])
                jo = Indices2i([jorb,ispin],[Norb,Nspin])
                !
                Aout(ii,jj) = Ain(isite,io,jo)
                !
             enddo
          enddo
       enddo
    enddo
  end function ss_pack_array2_d
  function ss_pack_array2_c(Ain,Nsite) result(Aout)
    integer                                                 :: Nsite
    complex(8),dimension(Nsite,Nso,Nso)                     :: Ain
    complex(8),dimension(Nspin*Nsite*Norb,Nspin*Nsite*Norb) :: Aout
    do ispin=1,Nspin
       do isite=1,Nsite
          do iorb=1,Norb
             ii = Indices2i([iorb,isite,ispin],[Norb,Nsite,Nspin])
             io = Indices2i([iorb,ispin],[Norb,Nspin])
             do jorb=1,Norb
                jj = Indices2i([jorb,isite,ispin],[Norb,Nsite,Nspin])
                jo = Indices2i([jorb,ispin],[Norb,Nspin])
                !
                Aout(ii,jj) = Ain(isite,io,jo)
                !
             enddo
          enddo
       enddo
    enddo
  end function ss_pack_array2_c



  function ss_unpack_array_d(Ain,Nsite) result(Aout)
    integer                                :: Nsite
    real(8),dimension(Nspin*Nsite*Norb)    :: Ain
    real(8),dimension(Nsite,Nspin*Norb)    :: Aout
    do ispin=1,Nspin
       do isite=1,Nsite
          do iorb=1,Norb
             ii = Indices2i([iorb,isite,ispin],[Norb,Nsite,Nspin])
             io = Indices2i([iorb,ispin],[Norb,Nspin])
             !
             Aout(isite,io) = Ain(ii)
             !
          enddo
       enddo
    enddo
  end function ss_unpack_array_d
  function ss_unpack_array_c(Ain,Nsite) result(Aout)
    integer                                :: Nsite
    complex(8),dimension(Nspin*Nsite*Norb) :: Ain
    complex(8),dimension(Nsite,Nspin*Norb) :: Aout
    do ispin=1,Nspin
       do isite=1,Nsite
          do iorb=1,Norb
             ii = Indices2i([iorb,isite,ispin],[Norb,Nsite,Nspin])
             io = Indices2i([iorb,ispin],[Norb,Nspin])
             !
             Aout(isite,io) = Ain(ii)
             !
          enddo
       enddo
    enddo
  end function ss_unpack_array_c
  function ss_unpack_array2_d(Ain,Nsite) result(Aout)
    integer                                              :: Nsite
    real(8),dimension(Nspin*Nsite*Norb,Nspin*Nsite*Norb) :: Ain
    real(8),dimension(Nsite,Nspin*Norb,Nspin*Norb)       :: Aout
    do ispin=1,Nspin            !always diagonal in the spin:
       do isite=1,Nsite         !loop over block diagonal: 1-site index only
          do iorb=1,Norb
             ii = Indices2i([iorb,isite,ispin],[Norb,Nsite,Nspin])
             io = Indices2i([iorb,ispin],[Norb,Nspin])
             do jorb=1,Norb
                jj = Indices2i([jorb,isite,ispin],[Norb,Nsite,Nspin])
                jo = Indices2i([jorb,ispin],[Norb,Nspin])
                !
                Aout(isite,io,jo) = Ain(ii,jj)
                !
             enddo
          enddo
       enddo
    enddo
  end function ss_unpack_array2_d
  function ss_unpack_array2_c(Ain,Nsite) result(Aout)
    integer                                                 :: Nsite
    complex(8),dimension(Nspin*Nsite*Norb,Nspin*Nsite*Norb) :: Ain
    complex(8),dimension(Nsite,Nspin*Norb,Nspin*Norb)       :: Aout
    do ispin=1,Nspin            !always diagonal in the spin:
       do isite=1,Nsite         !loop over block diagonal: 1-site index only
          do iorb=1,Norb
             ii = Indices2i([iorb,isite,ispin],[Norb,Nsite,Nspin])
             io = Indices2i([iorb,ispin],[Norb,Nspin])
             do jorb=1,Norb
                jj = Indices2i([jorb,isite,ispin],[Norb,Nsite,Nspin])
                jo = Indices2i([jorb,ispin],[Norb,Nspin])
                !
                Aout(isite,io,jo) = Ain(ii,jj)
                !
             enddo
          enddo
       enddo
    enddo
  end function ss_unpack_array2_c










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








  subroutine ss_spin_symmetry_nlso_d(array,Nsite)
    integer                      :: Nsite
    real(8),dimension(Nsite*Nso) :: array
    if(Nspin>1)array(Nsite*Norb+1:) = array(1:Nsite*Norb)
  end subroutine ss_spin_symmetry_nlso_d

  subroutine ss_spin_symmetry_nlso_c(array,Nsite)
    integer                         :: Nsite
    complex(8),dimension(Nsite*Nso) :: array
    if(Nspin>1)array(Nsite*Norb+1:) = array(1:Nsite*Norb)
  end subroutine ss_spin_symmetry_nlso_c

  subroutine ss_spin_symmetry_nn_d(array,Nsite)
    integer                      :: Nsite
    real(8),dimension(Nsite,Nso) :: array
    integer                      :: isite
    do isite=1,Nsite
       if(Nspin>1)array(isite,Norb+1:) = array(isite,1:Norb)
    enddo
  end subroutine ss_spin_symmetry_nn_d

  subroutine ss_spin_symmetry_nn_c(array,Nsite)
    integer                         :: Nsite
    complex(8),dimension(Nsite,Nso) :: array
    integer                         :: isite
    do isite=1,Nsite
       if(Nspin>1)array(isite,Norb+1:) = array(isite,1:Norb)
    enddo
  end subroutine ss_spin_symmetry_nn_c



















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
