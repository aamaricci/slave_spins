MODULE SS_SOLVE_FERMION
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  !
  USE SF_OPTIMIZE,only: brentq,fzero
  USE SF_MISC,    only: sort_array
  implicit none

  private

  public :: ss_solve_lambda0
  public :: ss_solve_fermions
  ! public :: ss_solve_fermions_Ef


  real(8),parameter :: mch=1d-6
  integer           :: iorb,jorb,ispin,io,jo,ilat,jlat,ineq

contains




  subroutine ss_solve_fermions()
    complex(8),dimension(Ns,Ns) :: Hk_f
    complex(8),dimension(Ns,Ns) :: Uk_f
    real(8),dimension(Ns)       :: Ek_f
    !
    complex(8),dimension(Ns,Ns) :: Eweiss
    complex(8),dimension(Ns,Ns) :: diagZ,diagR
    complex(8),dimension(Ns,Ns) :: rhoK
    real(8),dimension(Ns)       :: sq_zeta
    real(8),dimension(Ns)       :: lambda,lambda0,dens,weiss,const,zeta,Op
    real(8),dimension(Ns)       :: rhoDiag
    integer                     :: ik,i,j
    real(8),dimension(Nlat,Nss) :: TmpDens
    !
    if(Nspin==1)then
       call ss_spin_symmetry(ss_zeta,Nlat)
       call ss_spin_symmetry(ss_lambda,Nlat)
    endif
    !
    lambda  = ss_pack_array(ss_Lambda,Nlat)
    lambda0 = ss_pack_array(ss_Lambda0,Nlat)
    zeta    = ss_pack_array(ss_Zeta,Nlat)
    Op      = ss_pack_array(ss_Op,Nlat)
    if(any(zeta<0d0))then
       where(abs(zeta)< 10*solve_tolerance)zeta=0d0
       if(sum(abs(zeta)) > 1d-1)then
          print*,"WARNING:",zeta
          zeta=abs(zeta)
          ! stop "ERROR in ss_solve_fermions: any(ss_zeta)<0"
          ! else
          !    zeta=abs(zeta)
       endif
    endif
    sq_zeta = sqrt(zeta)
    diagZ   = diag(sq_zeta)
    ! diagZ = diag(Op)
    !
    Eweiss  = 0d0
    dens    = 0d0
    do ik = 1,Nk 
       Hk_f   = (diagZ .x. ss_Hk(:,:,ik)) .x. diagZ
       !forall(i=1:Ns,j=1:Ns)Hk_f(i,j) = Op(i)*ss_Hk(i,j,ik)*Op(j) !(diagZ .x. ss_Hk(:,:,ik)) .x. diagZ
       Uk_f   = Hk_f + ss_Hloc - xmu*eye(Ns)  - diag(lambda) + diag(lambda0)
       call eigh(Uk_f,Ek_f)
       diagR  = diag(step_fermi(Ek_f))
       RhoK   = (Uk_f .x. diagR) .x. (conjg(transpose(Uk_f)))
       Eweiss = Eweiss + ss_Hk(:,:,ik)*RhoK*ss_Wtk(:,:,ik) !element wise product
       dens   = dens + diagonal(RhoK*ss_Wtk(:,:,ik)) !element wise product
    enddo
    if(Nspin==1)call ss_spin_symmetry(dens,Nlat)
    !
    ! Get H_{a,s} = \sum_{b} sqrt(Z_{b,s})* sum_k H_{a,s, b,s}*\rho_{a,s, b,s}
    !             = \sum_{b} sqrt(Z_{b,s})* Eweiss_{a,s,b,s}
    weiss = 0d0
    do ispin=1,2
       do ilat=1,Nlat
          do iorb=1,Norb
             io = ss_Indices2i([iorb,ilat,ispin],nDefOrder) !io = iorb + (ispin-1)*Norb
             do jlat=1,Nlat
                do jorb=1,Norb
                   jo = ss_Indices2i([jorb,jlat,ispin],nDefOrder) !jo = jorb + (ispin-1)*Norb
                   weiss(io) = weiss(io) + sq_zeta(jo)*Eweiss(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    if(Nspin==1)call ss_spin_symmetry(weiss,Nlat)
    !
    ! Get C = ( n_{l,s}*(1-n_{l,s}))**{-1/2} - 1, at half-filling C=1
    const  = 1d0/(sqrt(dens*(1d0-dens))+mch) - 1d0
    !
    ss_Dens = ss_unpack_array(dens,Nlat)
    ss_C    = ss_unpack_array(const,Nlat)
    ss_Weiss= ss_unpack_array(weiss,Nlat)
    !
    do ineq=1,Nineq
       ilat = ss_ineq2ilat(ineq)
       ss_Dens_ineq(ineq,:)  = ss_Dens(ilat,:)
       ss_C_ineq(ineq,:)     = ss_C(ilat,:)
       ss_Weiss_ineq(ineq,:) = ss_Weiss(ilat,:)
    enddo
    !
  end subroutine ss_solve_fermions






  subroutine ss_solve_lambda0()
    complex(8),dimension(Ns,Ns)    :: Uk_f,Eweiss,diagRho,Rho
    complex(8),dimension(Ns,Ns,Nk) :: rhoK
    real(8),dimension(Ns,Nk)       :: eK
    real(8),dimension(Ns)          :: rhoDiag,Ek_f,lambda0,dens,weiss
    integer                        :: ik,unit
    integer                        :: stride
    real(8)                        :: mu0,Dmin,Dmax
    integer                        :: info
    real(8),dimension(Nlat,Nss)    :: TmpArray
    !
    do ik = 1,Nk
       Uk_f = ss_Hk(:,:,ik) + ss_Hloc
       call eigh(Uk_f,Ek_f)
       eK(:,ik)     = Ek_f
       rhoK(:,:,ik) = Uk_f
    enddo
    !
    Dmin = minval(Ek)
    Dmax = maxval(Ek)
    mu0 = Dmin
    call fzero(get_dens,mu0,Dmax,info)!,rguess=Dmin+0.5d0*(Dmax-Dmin))
    if(info/=1)then
       write(*,*)"ERROR ss_get_lambda0: fzero returned info>1 ",info
       stop
    endif
    if(verbose>3)write(*,"(A6,12G18.9)")"mu0  =",mu0
    !
    !
    Eweiss = 0d0    
    Dens   = 0d0
    do ik = 1,Nk 
       diagRho      = diag(step_fermi(eK(:,ik)-mu0))
       Uk_f         = rhoK(:,:,ik)
       rhoK(:,:,ik) = (Uk_f .x. diagRho) .x. (conjg(transpose(Uk_f)))
       Eweiss       = Eweiss  + ss_Hk(:,:,ik)*rhoK(:,:,ik)*ss_Wtk(:,:,ik) !element wise product
       Dens         = Dens + diagonal(rhoK(:,:,ik)*ss_Wtk(:,:,ik))        !element wise product
    enddo
    if(Nspin==1)call ss_spin_symmetry(Dens,Nlat)
    ss_Dens = ss_unpack_array(Dens,Nlat)
    !
    if(verbose>2)then
       do ilat=1,Nlat
          write(*,"(A6,12G18.9)")"N0   =",ss_Dens(ilat,:Nspin*Norb),&
               sum(ss_dens(Ilat,:))*(3-Nspin),filling
       enddo
    endif
    !
    !< Get H_{a,s} = \sum_{b}sum_k H_{a,s, b,s}*\rho_{a,s, b,s}
    Weiss = 0d0
    do ispin=1,2
       do ilat=1,Nlat
          do iorb=1,Norb
             io = ss_Indices2i([iorb,ilat,ispin],nDefOrder) !io = iorb + (ispin-1)*Norb
             do jlat=1,Nlat
                do jorb=1,Norb
                   jo = ss_Indices2i([jorb,jlat,ispin],nDefOrder) !jo = jorb + (ispin-1)*Norb
                   Weiss(io) = Weiss(io) + Eweiss(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !< Get Lambda0 = -2* h0_{m,s}*[n0_{m,s}-0.5]/[n0_{m,s}*(1-n0_{m,s})]
    lambda0 = -2d0*Weiss*(Dens-0.5d0)/(Dens*(1d0-Dens)+mch)
    xmu     =  2*sum(lambda0)/Ns + mu0!2*ss_lambda0(1)+mu0
    !
    ss_Lambda0 = ss_unpack_array(lambda0,Nlat)
    !
    do ilat=1,Nlat
       open(free_unit(unit),file="lambda0_site"//str(ilat,4)//".ss")
       write(unit,*)ss_lambda0(ilat,:)
       close(unit)
    enddo
    !
  contains

    function get_dens(mu) result(dens)
      real(8),intent(in)          :: mu
      real(8)                     :: dens
      real(8),dimension(Ns)       :: ndens
      ndens = 0d0
      do ik = 1,Nk 
         diagRho = diag(step_fermi(eK(:,ik)-mu))
         Uk_f    = rhoK(:,:,ik)
         Rho     = (Uk_f .x. diagRho) .x. (conjg(transpose(Uk_f)))
         ndens   = ndens + diagonal(Rho*ss_Wtk(:,:,ik))     !element wise product
      enddo
      if(verbose>4)write(*,"(A9,3G18.9)")"Ef,N0   =",mu,sum(ndens),filling
      dens = sum(ndens)-filling
    end function get_dens
  end subroutine ss_solve_lambda0





  !+------------------------------------------------------------------+
  !PURPOSE  : calculate step function
  !+------------------------------------------------------------------+
  elemental function step_fermi(x) result(step)
    real(8),intent(in)          :: x
    real(8)                     :: step
    ! step=1d0
    ! ! if(x>=0.d0)step=0d0
    ! if(x>0.d0)step=0d0
    step = fermi(x,beta)
  end function step_fermi

END MODULE SS_SOLVE_FERMION





! subroutine ss_solve_fermions_Ef()
!   complex(8),dimension(Ns,Ns) :: Hk_f,Uk_f,Eweiss,diagZ,Wtk,diagR
!   real(8),dimension(Ns)       :: sq_zeta,lambda,lambda0
!   integer                     :: ik,iorb,jorb,ispin,io,jo,indx
!   logical                     :: bool
!   real(8),dimension(Ns,Nk)    :: eK
!   real(8),dimension(Ns)       :: rhoDiag,Ek_f
!   real(8),dimension(Ns,Ns,Nk) :: rhoK
!   real(8),dimension(Nk*Ns)    :: Ek_all
!   integer                     :: stride,N_electrons,N_index
!   integer,dimension(Nk*Ns)    :: Ek_indx
!   real(8)                     :: Efilling
!   integer                     :: Nall
!   !
!   if(is_dos)stop "This method GG_xyx works only for H(k) on Bravais lattices."
!   !
!   bool = (Ns==Nspin*Norb)
!   !
!   if(Nspin==1)call ss_spin_symmetry(ss_zeta)
!   !
!   lambda  = ss_lambda
!   lambda0 = ss_lambda0
!   sq_zeta = sqrt(ss_zeta)
!   diagZ   = diag(sq_zeta)
!   !
!   stride  = 0
!   Eweiss  = 0d0
!   ss_dens = 0d0
!   do ik = 1,Nk 
!      Hk_f = (diagZ.x.ss_Hk(:,:,ik)) .x. diagZ
!      !
!      Uk_f = Hk_f + ss_Hloc - diag(lambda) + diag(lambda0)
!      !
!      call eigh(Uk_f,Ek_f)
!      !
!      eK(:,ik)     = Ek_f
!      rhoK(:,:,ik) = Uk_f
!      !
!      Ek_all(stride+1:stride+Ns) = Ek_f
!      stride = stride+Ns
!   enddo
!   !
!   call sort_array(Ek_all,Ek_indx)
!   !
!   indx  = ceiling(filling*Nk)
!   xmu   = Ek_all(indx)
!   !
!   do ik = 1,Nk 
!      rhoDiag = fermi(eK(:,ik)-xmu, beta)
!      diagR   = diag(rhoDiag)
!      Uk_f    = rhoK(:,:,ik)
!      rhoK(:,:,ik) = (Uk_f .x. diagR) .x. (conjg(transpose(Uk_f)))
!      !
!      Eweiss = Eweiss + ss_Hk(:,:,ik)*rhoK(:,:,ik)*ss_Wtk(:,:,ik) !element wise product
!      ss_dens= ss_dens + diagonal(rhoK(:,:,ik)*ss_Wtk(:,:,ik)) !element wise product
!      !
!   enddo
!   if(Nspin==1)call ss_spin_symmetry(ss_dens)
!   !
!   ! Get H_{a,s} = \sum_{b} sqrt(Z_{b,s})* sum_k H_{a,s, b,s}*\rho_{a,s, b,s}
!   !             = \sum_{b} sqrt(Z_{b,s})* Eweiss_{a,s,b,s}
!   ss_weiss= 0d0
!   do ispin=1,2
!      do ilat=1,Nlat
!         do iorb=1,Norb
!            io = ss_Indices2i([iorb,ilat,ispin],nDefOrder) !io = iorb + (ispin-1)*Norb
!            do jlat=1,Nlat
!               do jorb=1,Norb
!                  jo = ss_Indices2i([jorb,jlat,ispin],nDefOrder) !jo = jorb + (ispin-1)*Norb
!                  ss_weiss(io) = ss_weiss(io) + sq_zeta(jo)*Eweiss(io,jo)
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
!   !
!   ! Get C = ( n_{l,s}*(1-n_{l,s}))**{-1/2} - 1, at half-filling C=1
!   ss_c  = 1d0/(sqrt(ss_dens*(1d0-ss_dens))+mch) - 1d0
! end subroutine ss_solve_fermions_Ef



