MODULE SS_SOLVE_SPIN
  USE SS_VARS_GLOBAl
  USE SS_SETUP
  USE SS_SPARSE_MATRIX
  !
  USE SF_LINALG, only:eigh
  implicit none

  private

  type(sparse_matrix_csr) :: spHs

  public :: ss_solve_spins

contains


  subroutine ss_solve_spins
    integer                                   :: istate,jstate,io,iorb,ispin,Idegen
    integer                                   :: Ndim,Ndegen
    integer,dimension(Ns)                     :: Ivec
    real(8),dimension(Ns)                     :: Sz
    real(8),dimension(2,Norb)                 :: tSz
    real(8)                                   :: htmp
    !
    real(8),allocatable,dimension(:,:),target :: Evecs
    real(8),allocatable,dimension(:)          :: Evals
    real(8),dimension(:),pointer              :: gs_vec
    real(8),dimension(Ns)                     :: avSz,avOO
    !   
    Ndim = 2**Ns    
    call sp_init_matrix(spHs,Ndim)
    !
    do istate=1,Ndim
       Ivec = Bdecomp(istate,Ns)
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
       htmp = htmp + sum(ss_lambda*(Sz+0.5d0))
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
       !
       !Non-diagonal elements:
       !< sum_{m,s} h_{m,s} * [c S^+_{m,s} + S^-_{m,s}] + H.c.
       do io=1,Ns
          if(Sz(io)/=0.5d0)cycle
          call Sminus(io,Istate,Jstate)
          htmp = ss_Weiss(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
       do io=1,Ns
          if(Sz(io)/=-0.5d0)cycle
          call Splus(io,Istate,Jstate)
          htmp = ss_c(io)*ss_Weiss(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
       !
       do io=1,Ns
          if(Sz(io)/=-0.5d0)cycle
          call Splus(io,Istate,Jstate)
          htmp = ss_Weiss(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
       do io=1,Ns
          if(Sz(io)/=0.5d0)cycle
          call Sminus(io,Istate,Jstate)
          htmp = ss_c(io)*ss_Weiss(io)
          call sp_insert_element(spHs,htmp,Istate,Jstate)
       enddo
    enddo



    allocate(Evecs(Ndim,Ndim));Evecs=0d0
    allocate(Evals(Ndim));Evals=0d0
    !
    call sp_dump_matrix(spHs,Evecs)
    call sp_delete_matrix(spHs)
    !
    call eigh(Evecs,Evals)
    !
    Ndegen=1
    do istate=2,Ns
       if(abs(Evals(istate)-Evals(1))<= 1d-10)Ndegen=Ndegen+1
    end do
    zeta_function = dble(Ndegen)
    !
    avSz=0d0
    avOO=0d0
    do Idegen=1,Ndegen
       gs_vec => Evecs(:,Idegen)
       !
       do istate=1,Ndim
          !
          Ivec = Bdecomp(Istate,Ns)
          !
          Sz = 0.5d0 ; where(Ivec==0)Sz=-0.5d0
          !
          avSz = avSz + Sz*gs_vec(istate)**2/zeta_function          
          !
          do io=1,Ns
             if(Sz(io)/=0.5d0)cycle
             call Sminus(io,Istate,Jstate)
             htmp = 1d0
             avOO(io) = avOO(io) + gs_vec(Jstate)*htmp*gs_vec(Istate)/zeta_function
          enddo
          do io=1,Ns
             if(Sz(io)/=-0.5d0)cycle
             call Splus(io,Istate,Jstate)
             htmp = ss_c(io)
             avOO(io) = avOO(io) + gs_vec(Jstate)*htmp*gs_vec(Istate)/zeta_function
          enddo
       enddo
    enddo
    !
    ss_Sz = avSz
    print*,"avOO",avOO
    ss_zeta = avOO**2
    print*,"ps Z",ss_zeta
    !
    if(Nspin==1)call ss_spin_symmetry(ss_zeta)
    !
  end subroutine ss_solve_spins










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
