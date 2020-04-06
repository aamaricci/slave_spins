MODULE SS_IO
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  !  
  implicit none
  private


  interface ss_get_dens
     module procedure :: ss_get_dens_NN
     module procedure :: ss_get_dens_Ns
  end interface ss_get_dens

  interface ss_get_zeta
     module procedure :: ss_get_zeta_NN
     module procedure :: ss_get_zeta_Ns
  end interface ss_get_zeta

  interface ss_get_sz
     module procedure :: ss_get_Sz_NN
     module procedure :: ss_get_Sz_Ns
  end interface ss_get_sz

  interface ss_get_lambda
     module procedure :: ss_get_lambda_NN
     module procedure :: ss_get_lambda_Ns
  end interface ss_get_lambda

  interface ss_get_lambda0
     module procedure :: ss_get_lambda0_NN
     module procedure :: ss_get_lambda0_Ns
  end interface ss_get_lambda0

  interface ss_get_self
     module procedure :: ss_get_self_NN
     module procedure :: ss_get_self_Ns
  end interface ss_get_self

  public :: ss_get_dens
  public :: ss_get_zeta
  public :: ss_get_lambda
  public :: ss_get_lambda0  
  public :: ss_get_sz
  public :: ss_get_self
  public :: ss_get_ssHk


  integer                       :: ilat,iorb,ispin,io,Ivec(3)
  character(len=5),dimension(3) :: UserOrder_ = [character(len=5) :: "Norb","Nspin","Nlat"]
  integer,dimension(3)          :: nUserOrder

contains




  subroutine ss_get_dens_Ns(dens,UserOrder)
    real(8),dimension(Ns)                  :: Tmp
    real(8),dimension(Nspin*Nlat*Norb)     :: dens
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp  = ss_pack_array(ss_dens,Nlat)
    dens = ss_ss2user(Tmp(:Nlso),UserOrder_)
  end subroutine ss_get_dens_Ns
  !
  subroutine ss_get_dens_NN(dens,UserOrder)
    real(8),dimension(Ns)                  :: Tmp
    real(8),dimension(:,:,:)               :: dens
    real(8),dimension(Nspin*Nlat*Norb)     :: dens_
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp        = ss_pack_array(ss_dens,Nlat)
    dens_      = ss_ss2user(Tmp(:Nlso),UserOrder_)
    nUserOrder = ss_nOrder(UserOrder_)
    call assert_shape(dens,nUserOrder(3:1:-1),"ss_get_dens","dens")
    do io=1,Nlso
       Ivec = ss_i2indices(io,nUserOrder)       
       dens(Ivec(3),Ivec(2),Ivec(1)) = dens_(io)
    enddo
  end subroutine ss_get_dens_NN





  subroutine ss_get_zeta_Ns(zeta,UserOrder)
    real(8),dimension(Ns)                  :: Tmp
    real(8),dimension(Nspin*Nlat*Norb)     :: zeta
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp  = ss_pack_array(ss_zeta,Nlat)
    zeta = ss_ss2user(Tmp(:Nlso),UserOrder_)
  end subroutine ss_get_zeta_Ns
  !
  subroutine ss_get_zeta_NN(zeta,UserOrder)
    real(8),dimension(Ns)                  :: Tmp
    real(8),dimension(:,:,:)               :: zeta
    real(8),dimension(Nspin*Nlat*Norb)     :: zeta_
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp        = ss_pack_array(ss_zeta,Nlat)
    zeta_      = ss_ss2user(Tmp(:Nlso),UserOrder_)
    nUserOrder = ss_nOrder(UserOrder_)
    call assert_shape(zeta,nUserOrder(3:1:-1),"ss_get_zeta","zeta")
    do io=1,Nlso
       Ivec = ss_i2indices(io,nUserOrder)       
       zeta(Ivec(3),Ivec(2),Ivec(1)) = zeta_(io)
    enddo
  end subroutine ss_get_zeta_NN





  subroutine ss_get_Sz_Ns(Sz,UserOrder)
    real(8),dimension(Ns)                  :: Tmp
    real(8),dimension(Nspin*Nlat*Norb)     :: Sz
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp  = ss_pack_array(ss_Sz,Nlat)
    Sz   = ss_ss2user(Tmp(:Nlso),UserOrder_)
  end subroutine ss_get_Sz_Ns
  !
  subroutine ss_get_Sz_NN(Sz,UserOrder)
    real(8),dimension(Ns)                  :: Tmp
    real(8),dimension(:,:,:)               :: Sz
    real(8),dimension(Nspin*Nlat*Norb)     :: Sz_
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp      = ss_pack_array(ss_Sz,Nlat)
    Sz_      = ss_ss2user(Tmp(:Nlso),UserOrder_)
    nUserOrder = ss_nOrder(UserOrder_)
    call assert_shape(Sz,nUserOrder(3:1:-1),"ss_get_Sz","Sz")
    do io=1,Nlso
       Ivec = ss_i2indices(io,nUserOrder)       
       Sz(Ivec(3),Ivec(2),Ivec(1)) = Sz_(io)
    enddo
  end subroutine ss_get_Sz_NN






  subroutine ss_get_Lambda_Ns(Lambda,UserOrder)
    real(8),dimension(Ns)                  :: Tmp
    real(8),dimension(Nspin*Nlat*Norb)     :: Lambda
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp    = ss_pack_array(ss_Lambda,Nlat)
    Lambda = ss_ss2user(Tmp(:Nlso),UserOrder_)
  end subroutine ss_get_Lambda_Ns
  !
  subroutine ss_get_Lambda_NN(Lambda,UserOrder)
    real(8),dimension(Ns)                  :: Tmp
    real(8),dimension(:,:,:)               :: Lambda
    real(8),dimension(Nspin*Nlat*Norb)     :: Lambda_
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp          = ss_pack_array(ss_Lambda,Nlat)
    Lambda_      = ss_ss2user(Tmp(:Nlso),UserOrder_)
    nUserOrder = ss_nOrder(UserOrder_)
    call assert_shape(Lambda,nUserOrder(3:1:-1),"ss_get_Lambda","Lambda")
    do io=1,Nlso
       Ivec = ss_i2indices(io,nUserOrder)       
       Lambda(Ivec(3),Ivec(2),Ivec(1)) = Lambda_(io)
    enddo
  end subroutine ss_get_Lambda_NN




  subroutine ss_get_Lambda0_Ns(Lambda0,UserOrder)
    real(8),dimension(Ns)                  :: Tmp
    real(8),dimension(Nspin*Nlat*Norb)     :: Lambda0
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp     = ss_pack_array(ss_Lambda0,Nlat)
    Lambda0 = ss_ss2user(Tmp(:Nlso),UserOrder_)
  end subroutine ss_get_Lambda0_Ns
  !
  subroutine ss_get_Lambda0_NN(Lambda0,UserOrder)
    real(8),dimension(Ns)                  :: Tmp
    real(8),dimension(:,:,:)               :: Lambda0
    real(8),dimension(Nspin*Nlat*Norb)     :: Lambda0_
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp        = ss_pack_array(ss_Lambda0,Nlat)
    Lambda0_   = ss_ss2user(Tmp(:Nlso),UserOrder_)
    nUserOrder = ss_nOrder(UserOrder_)
    call assert_shape(Lambda0,nUserOrder(3:1:-1),"ss_get_Lambda0","Lambda0")
    do io=1,Nlso
       Ivec = ss_i2indices(io,nUserOrder)       
       Lambda0(Ivec(3),Ivec(2),Ivec(1)) = Lambda0_(io)
    enddo
  end subroutine ss_get_Lambda0_NN





  subroutine ss_get_self_Ns(self,UserOrder)
    real(8),dimension(Ns)                  :: Tmp,Tmp0
    real(8),dimension(Nlat*Nspin*Norb)     :: self
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp  = ss_pack_array(ss_Lambda,Nlat)
    Tmp0 = ss_pack_array(ss_Lambda0,Nlat)
    Self =  -ss_ss2user(Tmp(:Nlso),UserOrder_) + ss_ss2user(Tmp0(:Nlso),UserOrder_)
  end subroutine ss_get_self_Ns

  subroutine ss_get_self_NN(Self,UserOrder)
    real(8),dimension(Ns)                  :: Tmp,Tmp0
    real(8),dimension(:,:,:)               :: Self
    real(8),dimension(Nspin*Nlat*Norb)     :: Self_
    character(len=*),dimension(3),optional :: UserOrder
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp   = ss_pack_array(ss_Lambda,Nlat)
    Tmp0  = ss_pack_array(ss_Lambda0,Nlat)
    Self_ = -ss_ss2user(Tmp(:Nlso),UserOrder_) + ss_ss2user(Tmp0(:Nlso),UserOrder_)
    nUserOrder = ss_nOrder(UserOrder_)
    call assert_shape(Self,nUserOrder(3:1:-1),"ss_get_Self","Self")
    do io=1,Nlso
       Ivec = ss_i2indices(io,nUserOrder)       
       Self(Ivec(3),Ivec(2),Ivec(1)) = Self_(io)
    enddo
  end subroutine ss_get_self_NN




  subroutine ss_get_ssHk(ssHk,UserOrder)
    real(8),dimension(Ns)                                    :: Tmp,Lambda,Lambda0
    complex(8),dimension(Nspin*Nlat*Norb,Nspin*Nlat*Norb,Nk) :: ssHk
    character(len=*),dimension(3),optional                   :: UserOrder
    complex(8),dimension(Ns,Ns)                              :: Hk_f,diagZ
    integer,dimension(3)                                     :: Ivec,Jvec
    integer,dimension(3)                                     :: UserIndex
    integer,dimension(3)                                     :: nUserOrder
    integer                                                  :: ik,i
    integer                                                  :: iord,jord,iuser,juser
    character(len=5),dimension(3)                            :: UserOrder_
    !
    if(present(UserOrder))UserOrder_ = UserOrder
    !
    Tmp   = ss_pack_array(ss_Zeta,Nlat)
    diagZ = diag( sqrt(Tmp) )
    !
    Lambda  = ss_pack_array(ss_Lambda,Nlat)
    Lambda0 = ss_pack_array(ss_Lambda0,Nlat)
    do ik = 1,Nk 
       Hk_f   = (diagZ .x. ss_Hk(:,:,ik)) .x. diagZ
       Hk_f   = Hk_f + ss_Hloc - diag(lambda) + diag(lambda0)
       !
       ss_Hk(:,:,ik) = ss_ss2user(Hk_f(:Nlso,:Nlso),UserOrder_)
    enddo
  end subroutine ss_get_ssHk






  function ss_Norder(Order) result(nOrder)
    character(len=*),dimension(3)                         :: Order
    integer,dimension(3)                                  :: Index
    integer,dimension(3)                                  :: nOrder
    integer                                               :: i
    !
    do i=1,3     
       Index(i:i)=findloc(Order,DefOrder(i))
    enddo
    if(any(Index==0))then
       print*,"SS_Norder ERROR: wrong entry in Index at: ",findloc(Index,0)
       stop
    endif
    !
    nOrder = ss_indx_reorder(nDefOrder,Index)
  end function ss_Norder


END MODULE SS_IO
