MODULE SS_SOLVE_FERMION
  USE SS_INPUT_VARS
  USE SS_SETUP
  !
  USE SF_LINALG,  only: kron,eigh,diag,diagonal
  USE SF_SPECIAL, only: fermi
  USE SF_MISC,    only: sort
  implicit none

  private


  public :: ss_solve_fermions




contains


  subroutine ss_solve_fermions()
    complex(8),dimension(Ns,Ns) :: Hk_f,Uk_f,Hk,Eweiss,Hloc
    real(8),dimension(Ns)       :: sq_zeta,lambda
    integer                     :: ik,iorb,jorb,ispin,io,jo
    logical                     :: bool
    real(8),dimension(Ns,Nk)    :: eK
    real(8),dimension(Ns)       :: rhoDiag,Ek_f
    real(8),dimension(Ns,Ns,Nk) :: rhoK
    real(8),parameter           :: mch=1d-5
    real(8),dimension(Nk*Ns)    :: Ek_all
    integer                     :: stride,N_electrons
    integer,dimension(Nk*Ns)    :: Ek_indx
    real(8)                     :: Nall
    !
    bool = (Ns==Nspin*Norb)

    lambda  = ss_lambda
    print*,"lambda",lambda
    sq_zeta = sqrt(ss_zeta)
    print*,"ss_zeta",ss_zeta
    Eweiss  = 0d0
    !
    ss_dens = 0d0

    stride  = 0
    do ik = 1,Nk 
       select case(bool)
       case(.true.)
          Hk   = ss_Hk(:,:,ik)
          Hloc = ss_Hloc
       case(.false.)
          Hk   = kron(pauli_0,ss_Hk(:,:,ik))
          Hloc = kron(pauli_0,one*ss_Hloc)
       end select
       !
       !Hk_f = (diag(sq_zeta).mat_product.Hk).mat_product.diag(sq_zeta)
       Hk_f = matmul( matmul(diag(sq_zeta), Hk), diag(sq_zeta) ) !diag(sqrt(z)) . Htmp . diag(sqrt(z))
       !
       Uk_f = Hk_f + Hloc - diag(lambda)
       !
       call eigh(Uk_f,Ek_f)
       !
       eK(:,ik)     = Ek_f
       rhoK(:,:,ik) = Uk_f
       !
       Ek_all(stride+1:stride+Ns) = Ek_f
       stride = stride+Ns
    enddo

    call indexx(Nk*Ns,Ek_all,Ek_indx) !sort Ek_all

    ! print*,"Nk =",Nk,Ns*Nk
    ! print*,"filling*Nk",filling*Nk
    N_electrons = ceiling(filling*Nk)
    print*,"N =",N_electrons

    Nall = 0d0
    get_dens: do ik=1,Nk*Ns
       Nall = Nall + 1d0!fermi(Ek_all(Ek_indx(ik)),beta)
       if(Nall > N_electrons)exit get_dens
    enddo get_dens
    ! print*,"ciao",Nall,ik-1

    xmu = Ek_all(Ek_indx(ik-1)) !Ek_all(ik-1)

    print*,"mu=",xmu

    do ik = 1,Nk 
       !
       rhoDiag = fermi(eK(:,ik)-xmu, beta)
       Uk_f    = rhoK(:,:,ik)
       select case(ss_Hdiag)
       case(.true.)
          rhoK(:,:,ik) = diag(rhoDiag)
       case(.false.)
          rhoK(:,:,ik) = matmul( matmul(Uk_f, diag(rhoDiag)), conjg(transpose(Uk_f)) )
       end select
       !
       select case(bool)
       case(.true.)
          Hk = ss_Hk(:,:,ik)
       case(.false.)
          Hk = kron(pauli_0,ss_Hk(:,:,ik))
       end select
       Eweiss = Eweiss + Hk(:,:)*rhoK(:,:,ik)*ss_Wtk(ik)
       ss_dens= ss_dens+ diagonal(rhoK(:,:,ik))*ss_Wtk(ik)
       !
    enddo
    print*,ss_dens,sum(ss_dens)


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

    print*,ss_c
  end subroutine ss_solve_fermions



END MODULE SS_SOLVE_FERMION
