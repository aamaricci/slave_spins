MODULE SS_SOLVE_FERMION
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  !
  USE SF_OPTIMIZE,only: brentq,fzero
  USE SF_MISC,    only: sort_array
  implicit none

  private


  public :: ss_solve_fermions
  public :: ss_solve_fermions_Ef
  public :: ss_solve_lambda0

  
  real(8),parameter           :: mch=1d-6


contains


  subroutine ss_solve_fermions()
    complex(8),dimension(Ns,Ns) :: Hk_f,Uk_f,Eweiss,diagZ,Wtk,diagR
    real(8),dimension(Ns)       :: sq_zeta,lambda,lambda0,rhoDiag,Ek_f
    real(8),dimension(Ns,Ns)    :: rhoK
    integer                     :: ik,iorb,jorb,ispin,io,jo,indx
    logical                     :: bool
    !
    bool = (Ns==Nspin*Norb)
    !
    if(Nspin==1)then
       call ss_spin_symmetry(ss_zeta)
       call ss_spin_symmetry(ss_lambda)
    endif
    !
    lambda  = ss_lambda
    lambda0 = ss_lambda0
    if(any(ss_zeta<0d0))then
       print*,ss_zeta
       stop "ERROR in ss_solve_fermions: any(ss_zeta)<0"
    endif
    sq_zeta = sqrt(ss_zeta)
    diagZ   = diag(sq_zeta)
    !
    Eweiss  = 0d0
    ss_dens = 0d0
    !
    do ik = 1,Nk 
       Hk_f   = (diagZ .x. ss_Hk(:,:,ik)) .x. diagZ
       Uk_f   = Hk_f + ss_Hloc - xmu*eye(Ns)  - diag(lambda) + diag(lambda0)
       call eigh(Uk_f,Ek_f)
       diagR  = diag(fermi(Ek_f, beta))
       RhoK   = (Uk_f .x. diagR) .x. (conjg(transpose(Uk_f)))
       Eweiss = Eweiss + ss_Hk(:,:,ik)*RhoK*ss_Wtk(:,:,ik) !element wise product
       ss_dens= ss_dens + diagonal(RhoK*ss_Wtk(:,:,ik)) !element wise product
    enddo
    if(Nspin==1)call ss_spin_symmetry(ss_dens)
    !
    ! Get H_{a,s} = \sum_{b} sqrt(Z_{b,s})* sum_k H_{a,s, b,s}*\rho_{a,s, b,s}
    !             = \sum_{b} sqrt(Z_{b,s})* Eweiss_{a,s,b,s}
    ss_weiss= 0d0
    do ispin=1,2
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          do jorb=1,Norb
             jo = jorb + (ispin-1)*Norb
             ss_weiss(io) = ss_weiss(io) + sq_zeta(jo)*Eweiss(io,jo)
          enddo
       enddo
    enddo
    ! Get C = ( n_{l,s}*(1-n_{l,s}))**{-1/2} - 1, at half-filling C=1
    ss_c  = 1d0/(sqrt(ss_dens*(1d0-ss_dens))+mch) - 1d0
  end subroutine ss_solve_fermions









  subroutine ss_solve_fermions_Ef()
    complex(8),dimension(Ns,Ns) :: Hk_f,Uk_f,Eweiss,diagZ,Wtk,diagR
    real(8),dimension(Ns)       :: sq_zeta,lambda,lambda0
    integer                     :: ik,iorb,jorb,ispin,io,jo,indx
    logical                     :: bool
    real(8),dimension(Ns,Nk)    :: eK
    real(8),dimension(Ns)       :: rhoDiag,Ek_f
    real(8),dimension(Ns,Ns,Nk) :: rhoK
    real(8),dimension(Nk*Ns)    :: Ek_all
    integer                     :: stride,N_electrons,N_index
    integer,dimension(Nk*Ns)    :: Ek_indx
    real(8)                     :: Efilling
    integer                     :: Nall
    !
    if(is_dos)stop "This method GG_xyx works only for H(k) on Bravais lattices."
    !
    bool = (Ns==Nspin*Norb)
    !
    if(Nspin==1)call ss_spin_symmetry(ss_zeta)
    !
    lambda  = ss_lambda
    lambda0 = ss_lambda0
    sq_zeta = sqrt(ss_zeta)
    diagZ   = diag(sq_zeta)
    !
    stride  = 0
    Eweiss  = 0d0
    ss_dens = 0d0
    do ik = 1,Nk 
       Hk_f = (diagZ.x.ss_Hk(:,:,ik)) .x. diagZ
       !
       Uk_f = Hk_f + ss_Hloc - diag(lambda) + diag(lambda0)
       !
       call eigh(Uk_f,Ek_f)
       !
       eK(:,ik)     = Ek_f
       rhoK(:,:,ik) = Uk_f
       !
       Ek_all(stride+1:stride+Ns) = Ek_f
       stride = stride+Ns
    enddo
    !
    call sort_array(Ek_all,Ek_indx)
    !
    indx  = ceiling(filling*Nk)
    xmu   = Ek_all(indx)
    !
    do ik = 1,Nk 
       rhoDiag = fermi(eK(:,ik)-xmu, beta)
       diagR   = diag(rhoDiag)
       Uk_f    = rhoK(:,:,ik)
       rhoK(:,:,ik) = (Uk_f .x. diagR) .x. (conjg(transpose(Uk_f)))
       !
       Eweiss = Eweiss + ss_Hk(:,:,ik)*rhoK(:,:,ik)*ss_Wtk(:,:,ik) !element wise product
       ss_dens= ss_dens + diagonal(rhoK(:,:,ik)*ss_Wtk(:,:,ik)) !element wise product
       !
    enddo
    if(Nspin==1)call ss_spin_symmetry(ss_dens)
    !
    ! Get H_{a,s} = \sum_{b} sqrt(Z_{b,s})* sum_k H_{a,s, b,s}*\rho_{a,s, b,s}
    !             = \sum_{b} sqrt(Z_{b,s})* Eweiss_{a,s,b,s}
    ss_weiss= 0d0
    do ispin=1,2
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          do jorb=1,Norb
             jo = jorb + (ispin-1)*Norb
             ss_weiss(io) = ss_weiss(io) + sq_zeta(jo)*Eweiss(io,jo)
          enddo
       enddo
    enddo
    ! Get C = ( n_{l,s}*(1-n_{l,s}))**{-1/2} - 1, at half-filling C=1
    ss_c  = 1d0/(sqrt(ss_dens*(1d0-ss_dens))+mch) - 1d0
  end subroutine ss_solve_fermions_Ef



  subroutine ss_solve_lambda0()
    complex(8),dimension(Ns,Ns) :: Hk_f,Uk_f,Eweiss,diagZ,Wtk,diagR
    real(8),dimension(Ns)       :: sq_zeta,lambda
    integer                     :: ik,iorb,jorb,ispin,io,jo,indx,unit
    logical                     :: bool
    real(8),dimension(Ns,Nk)    :: eK
    real(8),dimension(Ns)       :: rhoDiag,Ek_f
    real(8),dimension(Ns,Ns,Nk) :: rhoK
    real(8),parameter           :: mch=1d-5
    real(8),dimension(Nk*Ns)    :: Ek_all
    integer                     :: stride,N_electrons,N_index
    integer,dimension(Nk*Ns)    :: Ek_indx
    real(8)                     :: mu0,Dmin,Dmax
    integer                     :: Nall,info
    !
    bool = (Ns==Nspin*Norb)
    !
    stride = 0
    do ik = 1,Nk 
       Uk_f = ss_Hk(:,:,ik) + ss_Hloc
       call eigh(Uk_f,Ek_f)
       eK(:,ik)     = Ek_f
       rhoK(:,:,ik) = Uk_f
       Ek_all(stride+1:stride+Ns) = Ek_f
       !
       stride = stride+Ns
    enddo

    Dmin = minval(Ek)
    Dmax = maxval(Ek)
    mu0 = Dmin
    call fzero(get_dens,mu0,Dmax,info,rguess=Dmin+0.5d0*(Dmax-Dmin))
    if(info/=1)then
       write(*,*)"ERROR ss_get_lambda0: fzero returned info>1 ",info
       stop
    endif
    if(verbose>3)write(*,"(A6,12G18.9)")"mu0  =",mu0

    Eweiss  = 0d0
    ss_dens = 0d0
    do ik = 1,Nk 
       diagR        = diag(fermi(eK(:,ik)-mu0, beta))
       Uk_f         = rhoK(:,:,ik)
       rhoK(:,:,ik) = (Uk_f .x. diagR) .x. (conjg(transpose(Uk_f)))
       Eweiss       = Eweiss  + ss_Hk(:,:,ik)*rhoK(:,:,ik)*ss_Wtk(:,:,ik) !element wise product
       ss_dens      = ss_dens + diagonal(rhoK(:,:,ik)*ss_Wtk(:,:,ik))     !element wise product
    enddo
    if(Nspin==1)call ss_spin_symmetry(ss_dens)
    if(verbose>2)write(*,"(A6,12G18.9)")"N0   =",ss_dens,sum(ss_dens),filling
    !
    !< Get H_{a,s} = \sum_{b}sum_k H_{a,s, b,s}*\rho_{a,s, b,s}
    ss_weiss= 0d0
    do ispin=1,2
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          do jorb=1,Norb
             jo = jorb + (ispin-1)*Norb
             ss_weiss(io) = ss_weiss(io) + Eweiss(io,jo)
          enddo
       enddo
    enddo
    !
    !< Get Lambda0 = -2* h0_{m,s}*[n0_{m,s}-0.5]/[n0_{m,s}*(1-n0_{m,s})]
    ss_lambda0 = -2d0*ss_Weiss*(ss_dens-0.5d0)/(ss_dens*(1d0-ss_dens)+mch)
    xmu        = 2*ss_lambda0(1)+mu0

    open(free_unit(unit),file="lambda0.ss")
    write(unit,*)ss_lambda0
    close(unit)

  contains

    function get_dens(mu) result(dens)
      real(8)                     :: mu
      real(8)                     :: dens
      real(8),dimension(Ns)       :: ndens(Ns)
      real(8),dimension(Ns,Ns)    :: Rho
      complex(8),dimension(Ns,Ns) :: Uk_f,diagRho
      ndens = 0d0
      do ik = 1,Nk 
         diagRho = diag(fermi(eK(:,ik)-mu, beta))
         Uk_f    = rhoK(:,:,ik)
         Rho     = (Uk_f .x. diagRho) .x. (conjg(transpose(Uk_f)))
         ndens   = ndens + diagonal(Rho*ss_Wtk(:,:,ik))     !element wise product
      enddo
      if(verbose>4)write(*,"(A9,3G18.9)")"Ef,N0   =",mu,sum(ndens),filling
      dens = sum(ndens)-filling
    end function get_dens

  end subroutine ss_solve_lambda0


END MODULE SS_SOLVE_FERMION



