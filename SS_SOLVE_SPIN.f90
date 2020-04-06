MODULE SS_SOLVE_SPIN
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  USE SS_SPARSE_MATRIX
  implicit none

  private

  public :: ss_solve_spins


  real(8),allocatable,dimension(:,:),target :: ss_Evecs
  real(8),allocatable,dimension(:)          :: ss_Evals
  integer                                   :: ss_Ndegen
  !
  !< site-resolved quantities, local to this module
  real(8),dimension(:),allocatable          :: ii_lambda
  real(8),dimension(:),allocatable          :: ii_Weiss
  real(8),dimension(:),allocatable          :: ii_c
  real(8),dimension(:),allocatable          :: ii_Sz
  real(8),dimension(:),allocatable          :: ii_Op
  real(8),dimension(:),allocatable          :: ii_Zeta
  real(8),dimension(:,:,:),allocatable      :: ii_SzSz
  !
  integer                                   :: Ndim
  integer                                   :: istate,jstate
  integer                                   :: iorb,jorb,ilat,ispin
  integer                                   :: io,il

contains



  subroutine ss_solve_spins(ineq)
    integer                :: ineq,ilat
    real(8),dimension(Nss) :: avSz,avOO,avZeta
    !
    Ndim = 2**Nss
    !
    allocate(Ss_Evecs(Ndim,Ndim)) ; Ss_Evecs=0d0
    allocate(Ss_Evals(Ndim)) ;Ss_Evals=0d0
    !
    !< get the ilat-th index corresponding to this ineq-site
    ilat = ss_ineq2ilat(ineq)
    !
    !< extract the ilat-th component of ss_lambda, ss_Weiss, ss_c
    !...[1...Norb]_ilat-th_up ...[1...Norb]_ilat-th_dw
    allocate(ii_lambda(Nss))
    allocate(ii_Weiss(Nss))
    allocate(ii_c(Nss))
    allocate(ii_Sz(Nss))
    allocate(ii_Op(Nss))
    allocate(ii_Zeta(Nss))
    allocate(ii_SzSz(4,Norb,Norb))
    !
    ii_Lambda = ss_Lambda(ilat,:)
    ii_Weiss  = ss_Weiss(ilat,:)
    ii_C      = ss_C(ilat,:)
    !
    !Build spin Hamiltonian and write it onto ss_Evecs
    call ss_build_Hs()
    call eigh(Ss_Evecs,Ss_Evals)
    !
    ss_Ndegen=1
    do istate=2,Ns
       if(abs(Ss_Evals(istate)-Ss_Evals(1))<= 1d-10)ss_Ndegen=ss_Ndegen+1
    end do
    zeta_function = dble(ss_Ndegen)
    !
    !< Get <Sz> and <O>--> Z=<O>**2
    call ss_Spin_observables()
    !< Get Sz-Sz correlations:
    call ss_SpinSpin_correlations()
    !< Copy back into main arrays ss_XYZ at proper position
    ss_Sz_ineq(ineq,:) = ii_Sz
    ss_Op_ineq(ineq,:) = ii_Op
    ss_Zeta_ineq(ineq,:) = ii_Zeta
    ss_SzSz_ineq(ineq,:,:,:) = ii_SzSz
    !
    deallocate(ii_lambda,ii_Weiss,ii_c)
    deallocate(ii_Sz,ii_Op,ii_Zeta,ii_SzSz)
    deallocate(ss_Evecs,ss_Evals)
    !
  end subroutine ss_solve_spins




  subroutine ss_build_Hs()
    integer,dimension(Nss)    :: Ivec
    real(8),dimension(Nss)    :: Sz
    real(8),dimension(2,Norb) :: tSz
    real(8)                   :: htmp
    type(sparse_matrix_csr)   :: spHs
    !
    !   
    call sp_init_matrix(spHs,Ndim)
    !
    do istate=1,Ndim
       Ivec = Bdecomp(istate,Nss)
       !
       Sz = 0.5d0 ; where(Ivec==0)Sz=-0.5d0
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
       !< Ust/2 * ( sum_{m,s}Sz_{m,s} )**2
       htmp = htmp + Ust/2*sum(Sz)**2
       !
       !< (U-Ust)/2 * sum_{m} (sum_s Sz_{m,s})**2
       do iorb=1,Norb
          htmp = htmp  + (Uloc(iorb)-Ust)/2d0*sum(tSz(:,iorb))**2
       enddo
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
          if(Sz(io)/=0.5d0)cycle
          call Sminus(io,Istate,Jstate)
          htmp = ii_Weiss(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
       do io=1,Nss
          if(Sz(io)/=-0.5d0)cycle
          call Splus(io,Istate,Jstate)
          htmp = ii_c(io)*ii_Weiss(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
       !
       do io=1,Nss
          if(Sz(io)/=-0.5d0)cycle
          call Splus(io,Istate,Jstate)
          htmp = ii_Weiss(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
       do io=1,Nss
          if(Sz(io)/=0.5d0)cycle
          call Sminus(io,Istate,Jstate)
          htmp = ii_c(io)*ii_Weiss(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
    enddo
    !
    call sp_dump_matrix(spHs,ss_Evecs) !dimensions check is done internally
    call sp_delete_matrix(spHs)
    !
  end subroutine ss_build_Hs



  subroutine ss_Spin_observables()

    integer                      :: Idegen
    integer,dimension(Nss)       :: Ivec
    real(8),dimension(Nss)       :: Sz
    real(8)                      :: htmp
    real(8),dimension(:),pointer :: gs_vec
    !
    ii_Sz=0d0
    ii_Op=0d0
    do Idegen=1,ss_Ndegen
       gs_vec => Ss_Evecs(:,Idegen)
       !
       do istate=1,Ndim
          Ivec = Bdecomp(Istate,Nss)
          Sz = 0.5d0 ; where(Ivec==0)Sz=-0.5d0
          !
          ii_Sz = ii_Sz + Sz*gs_vec(istate)**2/zeta_function
          !
          do io=1,Nss
             if(Sz(io)/=0.5d0)cycle
             call Sminus(io,Istate,Jstate)
             htmp = 1d0
             ii_Op(io) = ii_Op(io) + gs_vec(Jstate)*htmp*gs_vec(Istate)/zeta_function
          enddo
          do io=1,Nss
             if(Sz(io)/=-0.5d0)cycle
             call Splus(io,Istate,Jstate)
             htmp = ii_c(io)
             ii_Op(io) = ii_Op(io) + gs_vec(Jstate)*htmp*gs_vec(Istate)/zeta_function
          enddo
          !
       enddo
    enddo
    ii_zeta = ii_Op**2
  end subroutine ss_Spin_observables



  subroutine ss_SpinSpin_Correlations()
    integer                      :: Idegen
    integer,dimension(Nss)       :: Ivec
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
          Ivec = Bdecomp(Istate,Nss)
          Sz = 0.5d0 ; where(Ivec==0)Sz=-0.5d0
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
  !CREATION / DESTRUCTION OPERATORS
  !##################################################################
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=S^-_pos|in>  OR  |out>=S^+_pos|in> ; 
  !   pos labels the sites
  !+-------------------------------------------------------------------+
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
