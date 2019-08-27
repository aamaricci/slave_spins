MODULE SS_SOLVE_FERMION
  USE SS_INPUT_VARS
  USE SS_SETUP
  !
  USE SF_LINALG,  only: kron,eigh,diag,diagonal,operator(.x.)
  USE SF_SPECIAL, only: fermi,heaviside,step
  USE SF_MISC,    only: sort_array
  USE SF_IOTOOLS, only: save_array
  implicit none

  private


  public :: ss_solve_fermions


contains


  subroutine ss_solve_fermions()
    complex(8),dimension(Ns,Ns) :: Hk_f,Uk_f,Eweiss,diagZ,Wtk,diagR
    real(8),dimension(Ns)       :: sq_zeta,lambda
    integer                     :: ik,iorb,jorb,ispin,io,jo,indx
    logical                     :: bool
    real(8),dimension(Ns,Nk)    :: eK
    real(8),dimension(Ns)       :: rhoDiag,Ek_f
    real(8),dimension(Ns,Ns,Nk) :: rhoK
    real(8),parameter           :: mch=1d-5
    real(8),dimension(Nk*Ns)    :: Ek_all
    integer                     :: stride,N_electrons,N_index
    integer,dimension(Nk*Ns)    :: Ek_indx
    real(8)                     :: Efilling
    integer                     :: Nall
    !
    bool = (Ns==Nspin*Norb)
    !
    if(Nspin==1)call ss_spin_symmetry(ss_zeta)
    !
    lambda  = ss_lambda
    sq_zeta = sqrt(ss_zeta)
    diagZ   = diag(sq_zeta)
    !
    stride  = 0
    Eweiss  = 0d0
    ss_dens = 0d0
    do ik = 1,Nk 
       Hk_f = (diagZ.x.ss_Hk(:,:,ik)) .x. diagZ
       !
       Uk_f = Hk_f + ss_Hloc - diag(lambda)
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
    !< Get chemical potential: two approaches if is_bethe=F,T
    call sort_array(Ek_all,Ek_indx)
    !
    select case(is_bethe)
    case(.true.)                !DOS like Bethe needs special treatment 
       Efilling = 0d0
       ik_loop: do ik=1,Nk
          do io=1,Ns
             indx = io + (ik-1)*Ns
             Efilling = Efilling + ss_Wtk(io,io,ik) !this does not include Hloc so it is onyl for degenerate bands
             if(Efilling >= filling)exit ik_loop
          enddo
       enddo ik_loop
       indx=indx-1
    case (.false.)              !
       indx = ceiling(filling*Nk)
    end select

    print*,"Indx=",indx,Nk*Ns

    xmu = Ek_all(indx)
    print*,"Ef=",xmu
    !
    do ik = 1,Nk 
       ! rhoDiag = step(ek(:,ik)-xmu,.false.)
       rhoDiag = fermi(eK(:,ik)-xmu, beta)
       diagR   = diag(rhoDiag)
       Uk_f    = rhoK(:,:,ik)
       !
       select case(ss_Hdiag)
       case(.true.)
          rhoK(:,:,ik) = diagR
       case(.false.)
          rhoK(:,:,ik) = (Uk_f .x. diagR) .x. (conjg(transpose(Uk_f)))
       end select
       !
       Eweiss = Eweiss + ss_Hk(:,:,ik)*rhoK(:,:,ik)*ss_Wtk(:,:,ik) !element wise product
       ss_dens= ss_dens + diagonal(rhoK(:,:,ik)*ss_Wtk(:,:,ik)) !element wise product
       !
    enddo
    if(Nspin==1)call ss_spin_symmetry(ss_dens)
    print*,"dens=",ss_dens,sum(ss_dens),filling
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
    print*,"C = ",ss_c
  end subroutine ss_solve_fermions



END MODULE SS_SOLVE_FERMION
