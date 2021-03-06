MODULE SS_SOLVE_SPIN
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  USE SS_SPARSE_MATRIX
  !
  USE SF_TIMER, only: start_timer,stop_timer
  implicit none

  private

  public :: ss_solve_spins


  real(8),allocatable,dimension(:,:),target :: ss_Evecs
  real(8),allocatable,dimension(:)          :: ss_Evals
  integer                                   :: ss_Ndegen
  !
  !< site-resolved quantities, local to this module
  real(8),dimension(:),allocatable          :: ii_lambda
  real(8),dimension(:),allocatable          :: ii_c
  real(8),dimension(:),allocatable          :: ii_Sz
  real(8),dimension(:),allocatable          :: ii_Op
  real(8),dimension(:),allocatable          :: ii_Heff
  real(8),dimension(:,:),allocatable        :: ii_Jhybr
  real(8),dimension(:,:),allocatable        :: ii_OdgOp
  real(8),dimension(:,:,:),allocatable      :: ii_SzSz  
  !
  integer                                   :: Ndim,Nblock
  integer                                   :: istate,jstate,kstate
  integer                                   :: iorb,jorb,ilat,ispin
  integer                                   :: io,il,jo,jl

  type(sparse_matrix_csr)                   :: spHs



contains



  subroutine ss_solve_spins(ineq)
    integer                    :: ineq,ilat
    integer                    :: Nstate
    real(8),dimension(Nso)     :: Lambda,Heff,Const
    real(8),dimension(Nso,Nso) :: Jhybr
    !
    Ndim = 2**Nss
    !
    if(lanc_solve)then
       allocate(Ss_Evals(lanc_Neigen))
       allocate(Ss_Evecs(Ndim,lanc_Neigen))
    else
       allocate(Ss_Evals(Ndim))
       allocate(Ss_Evecs(Ndim,Ndim))
    endif
    !
    ss_Evecs=0d0
    ss_Evals=0d0
    !
    !< get the ilat-th index corresponding to this ineq-site
    ilat = ss_ineq2ilat(ineq)
    !
    !< extract the ilat-th component of ss_lambda, ss_Heff, ss_c, ...
    allocate(ii_lambda(Nss))
    allocate(ii_c(Nss))
    allocate(ii_Sz(Nss))
    allocate(ii_Op(Nss))
    allocate(ii_Heff(Nss))
    allocate(ii_Jhybr(Nss,Nss))
    allocate(ii_OdgOp(Nss,Nss))
    allocate(ii_SzSz(4,Norb,Norb))
    !
    Lambda = ss_Lambda(ilat,:)
    Heff   = ss_Heff(ilat,:)
    Const  = ss_C(ilat,:)
    Jhybr  = ss_Jhybr(ilat,:,:)
    if(Nspin==1)then
       ii_Lambda = [Lambda,Lambda]
       ii_Heff   = [Heff,Heff]
       ii_C      = [Const,Const]
       ii_Jhybr  = kron(pauli_0,one*Jhybr)
    else
       ii_Lambda = Lambda
       ii_Heff   = Heff
       ii_C      = Const
       ii_Jhybr  = Jhybr
    endif
    !
    !< Init the sparse matrix spHs (SParse Hamiltonian Spins)
    call sp_init_matrix(spHs,Ndim)
    !
    !
    !< Build the spin Hamiltonian
    if(master.AND.verbose>3)call start_timer
    call ss_build_Hs()
    if(master.AND.verbose>3)call stop_timer("build_Hs")
    !
    if(lanc_solve)then
       !< Solve with ARPACK Lanczos
       if(master.AND.verbose>3)call start_timer
       Nblock   = min(Ndim,10*lanc_Neigen)
       call sp_eigh(spMatVec_p,ss_Evals,ss_Evecs,Nblock,Nitermax=1000,iverbose=(verbose>5))
       if(master.AND.verbose>3)call stop_timer("sp_eigh")
       Nstate = lanc_Neigen
    else
       !< Solve with Lapack exact diagonalization
       if(master.AND.verbose>3)call start_timer
       call sp_dump_matrix(spHs,ss_Evecs) !dimensions check is done internally
       call eigh(Ss_Evecs,Ss_Evals)
       if(master.AND.verbose>3)call stop_timer("eigh")
       Nstate=Ndim
    endif
    !
    !< Get degeneracy of the groundstate (we are doing T=0)
    ss_Ndegen=1
    do istate=2,Nstate
       if(abs(Ss_Evals(istate)-Ss_Evals(1))<= 1d-10)ss_Ndegen=ss_Ndegen+1
    end do
    zeta_function = dble(ss_Ndegen)
    !
    !< Delete matrix, no more used
    call sp_delete_matrix(spHs)
    !
    !
    !< Get <Sz> and <O>, <O^+>, <O^+O>, ...
    call ss_Spin_observables()
    !
    !< Get Sz-Sz correlations:
    call ss_SpinSpin_correlations()
    !
    !< Copy back into main arrays ss_XYZ at proper position
    ss_Sz_ineq(ineq,:) = ii_Sz(:Nso)
    ss_Op_ineq(ineq,:) = ii_Op(:Nso)
    ss_OdgOp_ineq(ineq,:,:)  = ii_OdgOp(:Nso,:Nso)
    ss_SzSz_ineq(ineq,:,:,:) = ii_SzSz
    !
    if(master.AND.verbose>4)then
       write(*,"(A,12G18.9)")"Sz   =",ss_Sz_ineq(ineq,:)
       write(*,"(A,12G18.9)")"Op   =",ss_Op_ineq(ineq,:)
       if(verbose>5)then
          write(*,"(A)")"OdgOp="
          do io=1,Nso
             write(*,"(100G18.9)")(ss_OdgOp_ineq(ineq,io,jo),jo=1,Nso)
          enddo
       endif
    endif
    !
    deallocate(ii_lambda,ii_Heff,ii_c,ii_Jhybr)
    deallocate(ii_Sz,ii_Op,ii_OdgOp,ii_SzSz)
    deallocate(ss_Evecs,ss_Evals)
    !
  end subroutine ss_solve_spins




  subroutine ss_build_Hs()
    real(8),dimension(Nss)    :: Sz
    real(8),dimension(2,Norb) :: tSz
    real(8)                   :: htmp
    !
    do istate=1,Ndim
       Sz = Bdecomp(istate,Nss) - 0.5d0
       do ispin=1,2
          do iorb=1,Norb
             tSz(ispin,iorb) = Sz(iorb+(ispin-1)*Norb)
          enddo
       enddo
       !
       !Diagonal elements
       htmp = 0d0
       !
       !< sum_{m,s}lambda_{m,s}*(Sz_{m,s}+1/2)
       htmp = htmp + sum(ii_lambda*(Sz+0.5d0))
       !
       !< (U-Ust)/2 * sum_{m} (sum_s Sz_{m,s})**2
       do iorb=1,Norb
          htmp = htmp  + (Uloc(iorb)-Ust)/2d0*sum(tSz(:,iorb))**2
       enddo
       !
       !< Ust/2 * ( sum_{m,s}Sz_{m,s} )**2
       htmp = htmp + Ust/2*sum(tSz)**2
       !
       !< -Jh/2 * sum_s (sum_m Sz_{m,s})**2
       do ispin=1,2
          htmp = htmp - Jh/2d0*sum(tSz(ispin,:))**2
       enddo
       !
       call sp_insert_element(spHs,htmp,Istate,Istate)
       !
       !Non-diagonal elements:
       !< sum_{m,s} h_{m,s} * [c S^+_{m,s} + S^-_{m,s}] + H.c.
       do io=1,Nss
          if(.not.test_Splus(io,Istate))cycle
          call Splus(io,Istate,Jstate)
          htmp = ii_c(io)*ii_Heff(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
       do io=1,Nss
          if(.not.test_Sminus(io,Istate))cycle
          call Sminus(io,Istate,Jstate)
          htmp = ii_Heff(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
       !
       do io=1,Nss
          if(.not.test_Sminus(io,Istate))cycle
          call Sminus(io,Istate,Jstate)
          htmp = ii_c(io)*ii_Heff(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
       do io=1,Nss
          if(.not.test_Splus(io,Istate))cycle
          call Splus(io,Istate,Jstate)
          htmp = ii_Heff(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
    enddo
    !
  end subroutine ss_build_Hs



  subroutine ss_Spin_observables()
    integer                      :: Idegen
    real(8),dimension(Nss)       :: Sz
    real(8)                      :: htmp
    real(8),dimension(:),pointer :: gs_vec
    !
    ii_Sz=0d0
    ii_Op=0d0
    ii_OdgOp=0d0
    do Idegen=1,ss_Ndegen
       gs_vec => Ss_Evecs(:,Idegen)
       !
       do istate=1,Ndim
          Sz = Bdecomp(Istate,Nss) - 0.5d0
          !
          !<Sz>
          ii_Sz = ii_Sz + Sz*gs_vec(istate)**2/zeta_function
          !
          !<O_m> =  c_m <S^+_m> + <S^-_m>
          !c_m <S^+_m>
          do io=1,Nss
             if(.not.test_Sminus(io,Istate))cycle
             call Sminus(io,Istate,Jstate)
             htmp = 1d0
             ii_Op(io) = ii_Op(io) + htmp*gs_vec(Istate)*gs_vec(Jstate)/zeta_function
          enddo
          !<S^-_m>
          do io=1,Nss
             if(.not.test_Splus(io,Istate))cycle
             call Splus(io,Istate,Jstate)
             htmp = ii_c(io)
             ii_Op(io) = ii_Op(io) + htmp*gs_vec(Istate)*gs_vec(Jstate)/zeta_function
          enddo
          !
          !
          if(product(ii_Jhybr)==0d0)then
             !<O^+ O> = J_ij <[c*_i S^-_i + S^+_i] [c_j S^+_j + S^-_j]>
             !        = J_ij [c*_i.c_j <S^-_i.S^+_j>  + &
             !                 c*_i    <S^-_i.S^-_j>  + &
             !                 c_j     <S^+_i.S^+_j>  + &
             !                 1       <S^+_i.S^-_j>  ]
             !
             !"c*_i.c_j <S^-_i S^+_j>"
             do io=1,Nss
                do jo=1,Nss
                   if(.not.test_Sminus(io,Istate))cycle
                   call Sminus(io,Istate,Kstate)
                   if(.not.test_Splus(jo,Kstate))cycle
                   call Splus(jo,Kstate,Jstate)
                   htmp = ii_c(io)*ii_c(jo)*ii_Jhybr(io,jo)
                   ii_OdgOp(io,jo) = ii_OdgOp(io,jo) + &
                        htmp*gs_vec(Istate)*gs_vec(Jstate)/zeta_function
                enddo
             enddo
             !
             !"c*_i <S^-_i S^-_j>"             
             do io=1,Nss
                do jo=1,Nss
                   if(.not.test_Sminus(io,Istate))cycle
                   call Sminus(io,Istate,Kstate)
                   if(.not.test_Sminus(jo,Kstate))cycle
                   call Sminus(jo,Kstate,Jstate)
                   htmp = ii_c(io)*ii_Jhybr(io,jo)
                   ii_OdgOp(io,jo) = ii_OdgOp(io,jo) + &
                        htmp*gs_vec(Istate)*gs_vec(Jstate)/zeta_function
                enddo
             enddo
             !
             !"c_j <S^+_i S^+_j>"
             do io=1,Nss
                do jo=1,Nss
                   if(.not.test_Splus(io,Istate))cycle
                   call Splus(io,Istate,Kstate)
                   if(.not.test_Splus(jo,Kstate))cycle
                   call Splus(jo,Kstate,Jstate)
                   htmp = ii_c(jo)*ii_Jhybr(io,jo)
                   ii_OdgOp(io,jo) = ii_OdgOp(io,jo) + &
                        htmp*gs_vec(Istate)*gs_vec(Jstate)/zeta_function
                enddo
             enddo
             !
             !"<S^+_i S^-_j>"
             do io=1,Nss
                do jo=1,Nss
                   if(.not.test_Splus(io,Istate))cycle
                   call Splus(io,Istate,Kstate)
                   if(.not.test_Sminus(jo,Kstate))cycle
                   call Sminus(jo,Kstate,Jstate)
                   htmp = ii_Jhybr(io,jo)
                   ii_OdgOp(io,jo) = ii_OdgOp(io,jo) + &
                        htmp*gs_vec(Istate)*gs_vec(Jstate)/zeta_function
                enddo
             enddo
          endif
          !
          !
       enddo
    enddo
  end subroutine ss_Spin_observables



  subroutine ss_SpinSpin_Correlations()
    integer                      :: Idegen
    real(8),dimension(Nss)       :: Sz
    real(8),dimension(:),pointer :: gs_vec
    real(8),dimension(Nss,Nss)   :: avSzSz
    !
    if(.not.allocated(ss_Evecs))stop "SS_SPIN_CORR Error: ss_Evecs not allocated"
    !
    avSzSz=0d0
    do Idegen=1,ss_Ndegen
       gs_vec => Ss_Evecs(:,Idegen)
       do istate=1,Ndim
          Sz = Bdecomp(Istate,Nss) - 0.5d0
          avSzSz = avSzSz + outerprod(Sz,Sz)*gs_vec(istate)**2/zeta_function
       enddo
    enddo
    !
    do iorb=1,Norb
       do jorb=1,Norb
          ii_SzSz(1,iorb,jorb) = avSzSz(is_indx(iorb,1),is_indx(jorb,1)) !up-up
          ii_SzSz(2,iorb,jorb) = avSzSz(is_indx(iorb,2),is_indx(jorb,2)) !dw-dw
          ii_SzSz(3,iorb,jorb) = avSzSz(is_indx(iorb,1),is_indx(jorb,2)) !up-dw
          ii_SzSz(4,iorb,jorb) = avSzSz(is_indx(iorb,2),is_indx(jorb,1)) !dw-up
       enddo
    enddo
    !
  contains
    !
    function is_indx(iorb,ispin)
      integer :: iorb,ispin
      integer :: is_indx
      is_indx = iorb + (ispin-1)*Norb
    end function is_indx
    !
  end subroutine ss_SpinSpin_Correlations





  !##################################################################
  !##################################################################
  !MATRIX-VECTOR OPERATION USING SPARSE STRUCTURE
  !##################################################################
  !##################################################################
  subroutine spMatVec_p(Nloc,v,Hv)
    integer                         :: Nloc
    real(8),dimension(Nloc)      :: v
    real(8),dimension(Nloc)      :: Hv
    integer                         :: i,j
    Hv=zero
    do i=1,Nloc
       matmul: do j=1,spHs%row(i)%Size
          Hv(i) = Hv(i) + spHs%row(i)%vals(j)*v(spHs%row(i)%cols(j))
       end do matmul
    end do
  end subroutine spMatVec_p






  !##################################################################
  !##################################################################
  !CREATION / DESTRUCTION OPERATORS
  !##################################################################
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=S^-_pos|in>  OR  |out>=S^+_pos|in> ; 
  !   pos labels the sites
  !+-------------------------------------------------------------------+
  function test_Sminus(pos,in) result(bool)
    integer,intent(in) :: pos
    integer,intent(in) :: in
    logical            :: bool
    bool = btest(in-1,pos-1)
  end function test_Sminus

  function test_Splus(pos,in) result(bool)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    logical            :: bool
    bool = .not.btest(in-1,pos-1)
  end function test_Splus


  subroutine Sminus(pos,in,out)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    integer               :: in_
    in_ = in-1
    if(.not.btest(in_,pos-1))stop "S^- error: S^-_i|...dw_i...>"
    out = ibclr(in_,pos-1)+1
  end subroutine Sminus


  subroutine Splus(pos,in,out)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    integer               :: in_
    in_ = in-1
    if(btest(in_,pos-1))stop "S^+ error: S^+_i|...up_i...>"
    out = ibset(in_,pos-1)+1
  end subroutine Splus






  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !I = 0,...,2**Ns-1
    do l=0,Ntot-1
       busy=btest(i-1,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp




END MODULE SS_SOLVE_SPIN
