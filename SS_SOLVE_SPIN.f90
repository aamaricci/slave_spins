MODULE SS_SOLVE_SPIN
  USE SS_VARS_GLOBAl
  USE SS_SETUP
  implicit none

  private

  public :: ss_solve_spins
  public :: ss_spin_corr


contains


  subroutine ss_solve_spins
    integer                                   :: istate,jstate,io,iorb,ispin,Idegen
    integer,dimension(Ns)                     :: Ivec
    real(8),dimension(Ns)                     :: Sz
    real(8),dimension(2,Norb)                 :: tSz
    real(8)                                   :: htmp
    !
    real(8),dimension(:),pointer              :: gs_vec
    real(8),dimension(Ns)                     :: avSz,avOO
    !   
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
    !
    Ss_Evecs=0d0
    Ss_Evals=0d0
    !
    call sp_dump_matrix(spHs,Ss_Evecs)
    call sp_delete_matrix(spHs)
    !
    call eigh(Ss_Evecs,Ss_Evals)
    !
    ss_Ndegen=1
    do istate=2,Ns
       if(abs(Ss_Evals(istate)-Ss_Evals(1))<= 1d-10)ss_Ndegen=ss_Ndegen+1
    end do
    zeta_function = dble(ss_Ndegen)
    !
    ss_Sz=0d0
    ss_Op=0d0
    do Idegen=1,ss_Ndegen
       gs_vec => Ss_Evecs(:,Idegen)
       !
       do istate=1,Ndim
          !
          Ivec = Bdecomp(Istate,Ns)
          !
          Sz = 0.5d0 ; where(Ivec==0)Sz=-0.5d0
          !
          ss_Sz = ss_Sz + Sz*gs_vec(istate)**2/zeta_function          
          !
          do io=1,Ns
             if(Sz(io)/=0.5d0)cycle
             call Sminus(io,Istate,Jstate)
             htmp = 1d0
             ss_Op(io) = ss_Op(io) + gs_vec(Jstate)*htmp*gs_vec(Istate)/zeta_function
          enddo
          do io=1,Ns
             if(Sz(io)/=-0.5d0)cycle
             call Splus(io,Istate,Jstate)
             htmp = ss_c(io)
             ss_Op(io) = ss_Op(io) + gs_vec(Jstate)*htmp*gs_vec(Istate)/zeta_function
          enddo
       enddo
    enddo
    !
    ss_zeta = ss_Op**2
    !
    if(Nspin==1)call ss_spin_symmetry(ss_zeta)
    !
  end subroutine ss_solve_spins



  subroutine ss_spin_corr
    integer                      :: istate,iorb,jorb,Idegen
    integer,dimension(Ns)        :: Ivec
    real(8),dimension(Ns)        :: Sz
    real(8),dimension(:),pointer :: gs_vec
    real(8),dimension(Ns,Ns)     :: avSzSz
    !
    if(.not.allocated(ss_Evecs))stop "SS_SPIN_CORR Error: ss_Evecs not allocated"
    !
    avSzSz=0d0
    !
    do Idegen=1,ss_Ndegen
       gs_vec => Ss_Evecs(:,Idegen)
       !
       do istate=1,Ndim
          Ivec = Bdecomp(Istate,Ns)
          Sz = 0.5d0 ; where(Ivec==0)Sz=-0.5d0
          avSzSz = avSzSz + outerprod(Sz,Sz)*gs_vec(istate)**2/zeta_function
       enddo
    enddo
    !
    do iorb=1,Norb
       do jorb=1,Norb
          ss_SzSz(1,iorb,jorb) = avSzSz(is_indx(iorb,1),is_indx(jorb,1)) !up-up
          ss_SzSz(2,iorb,jorb) = avSzSz(is_indx(iorb,2),is_indx(jorb,2)) !dw-dw
          ss_SzSz(3,iorb,jorb) = avSzSz(is_indx(iorb,1),is_indx(jorb,2)) !up-dw
          ss_SzSz(4,iorb,jorb) = avSzSz(is_indx(iorb,2),is_indx(jorb,1)) !dw-up
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
  end subroutine ss_spin_corr














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
