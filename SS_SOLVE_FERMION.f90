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


  real(8),parameter                  :: mch=1d-6
  integer                            :: iorb,jorb,ispin,io,jo,ilat,jlat,ineq


contains

  subroutine ss_solve_fermions()
    complex(8),dimension(Nlso,Nlso) :: Hk_f
    complex(8),dimension(Nlso,Nlso) :: Uk_f
    real(8),dimension(Nlso)         :: Ek_f
    real(8),dimension(Nlso,Nlso)    :: Wtk
    !
    complex(8),dimension(Nlso,Nlso) :: Eweiss,Eweiss_tmp
    real(8),dimension(Nlso)         :: dens,dens_tmp
    complex(8),dimension(Nlso,Nlso) :: diagR
    complex(8),dimension(Nlso,Nlso) :: rhoK
    real(8),dimension(Ns)           :: lambda,lambda0
    real(8),dimension(Ns)           :: const
    real(8),dimension(Ns)           :: Op
    real(8),dimension(Ns)           :: rhoDiag
    !
    complex(8),dimension(Ns)        :: weiss
    integer                         :: ik,i,j,N

    !
    if(Nspin==1)then
       call ss_spin_symmetry(ss_lambda0,Nlat)
       call ss_spin_symmetry(ss_lambda,Nlat)
       call ss_spin_symmetry(ss_Op,Nlat)
    endif
    !
    lambda  = ss_pack_array(ss_Lambda,Nlat)
    lambda0 = ss_pack_array(ss_Lambda0,Nlat)
    Op      = ss_pack_array(ss_Op,Nlat)
    !
    !
#ifdef _MPI
    if(check_MPI())then
       mpi_rank=get_rank_MPI()
       mpi_size=get_size_MPI()
    else
       mpi_rank=0
       mpi_size=1
    endif
#endif
    !
    Eweiss_tmp  = zero;Eweiss=zero
    dens_tmp    = 0d0 ;dens=0d0
    do ik=1+mpi_rank,Nk,mpi_size
       forall(i=1:Nlso,j=1:Nlso)Hk_f(i,j) = Op(i)*ss_Hk(i,j,ik)*Op(j)
       Uk_f   = Hk_f + diag(ss_Hdiag)-xmu*eye(Nlso)-diag(lambda(:Nlso))+diag(lambda0(:Nlso))
       call eigh(Uk_f,Ek_f)
       diagR  = diag(step_fermi(Ek_f))
       RhoK   = (Uk_f .x. diagR) .x. (conjg(transpose(Uk_f)))
       Wtk    = ss_Wtk(1,ik) ; if(is_dos)Wtk = diag(ss_Wtk(:,ik))
       Eweiss_tmp = Eweiss_tmp + ss_Hk(:,:,ik)*RhoK*Wtk
       dens_tmp   = dens_tmp + diagonal(RhoK*Wtk)
    enddo
#ifdef _MPI
    if(check_MPI())then
       call AllReduce_MPI(MPI_COMM_WORLD,Eweiss_tmp,Eweiss)
       call AllReduce_MPI(MPI_COMM_WORLD,dens_tmp,dens)
    else
       Eweiss = Eweiss_tmp
       dens   = dens_tmp
    endif
#else
    Eweiss = Eweiss_tmp
    dens   = dens_tmp
#endif
    !
    ! Get H_{a,s} = \sum_{b} sqrt(Z_{b,s})* sum_k H_{a,s, b,s}*\rho_{a,s, b,s}
    !             = \sum_{b} sqrt(Z_{b,s})* Eweiss_{a,s,b,s}
    weiss = 0d0
    do ispin=1,Nspin
       do ilat=1,Nlat
          do iorb=1,Norb
             io = ss_Indices2i([iorb,ilat,ispin],[Norb,Nlat,Nspin])
             do jlat=1,Nlat
                do jorb=1,Norb
                   jo = ss_Indices2i([jorb,jlat,ispin],[Norb,Nlat,Nspin])
                   weiss(io) = weiss(io) + Op(jo)*Eweiss(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    ! Get C = (n_{l,s}*(1-n_{l,s}))**{-1/2} - 1, at half-filling C=1
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
    if(Nspin==1)then
       call ss_spin_symmetry(ss_Dens_ineq,Nineq)
       call ss_spin_symmetry(ss_C_ineq,Nineq)
       call ss_spin_symmetry(ss_Weiss_ineq,Nineq)
    endif
    !
  end subroutine ss_solve_fermions





  subroutine ss_solve_lambda0()
    complex(8),dimension(Ns,Ns)      :: Uk_f,Eweiss,diagRho,Rho,Eweiss_tmp
    complex(8),dimension(Ns,Ns,Nk)   :: rhoK,rhoK_tmp
    real(8),dimension(Ns,Nk)         :: eK,eK_tmp
    real(8),dimension(Ns,Ns)         :: Wtk
    real(8),dimension(Ns)            :: rhoDiag,Ek_f
    real(8),dimension(Ns)            :: lambda0,dens,dens_tmp
    complex(8),dimension(Ns)         :: weiss
    integer                          :: ik,unit,N
    integer                          :: stride
    real(8)                          :: mu0,Dmin,Dmax
    integer                          :: info
    logical                          :: IOfile
    integer                          :: Len
    real(8),dimension(:),allocatable :: params
    !
#ifdef _MPI
    if(check_MPI())then
       mpi_rank=get_rank_MPI()
       mpi_size=get_size_MPI()
    endif
#endif
    !
    inquire(file=trim(Pfile)//"0"//trim(ss_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       len = file_length(trim(Pfile)//"0"//trim(ss_file_suffix)//".restart")
       if( (len/=2*Ns+1) )stop "SS_SOLVE_LAMBDA0 ERROR: len!=2*Ns OR 2*Ns+1"
       allocate(params(len))
       call read_array(trim(Pfile)//"0"//trim(ss_file_suffix)//".restart",params)
       lambda0 = params(1:Ns)
       dens    = params(Ns+1:2*Ns)
       mu0     = params(2*Ns+1)
       !
    else
       !
       ek_tmp=0d0   ;ek=0d0
       rhoK_tmp=zero;rhoK=zero
       do ik=1+mpi_rank,Nk,mpi_size
          Uk_f = ss_Hk(:,:,ik) + diag(ss_Hdiag)
          call eigh(Uk_f,Ek_f)
          eK_tmp(:,ik)     = Ek_f
          rhoK_tmp(:,:,ik) = Uk_f
       enddo
#ifdef _MPI
       if(check_MPI())then
          call AllReduce_MPI(MPI_COMM_WORLD,eK_tmp,eK)
          call AllReduce_MPI(MPI_COMM_WORLD,rhoK_tmp,rhoK)
       else
          eK   = eK_tmp
          rhoK = rhoK_tmp
       endif
#else
       eK   = eK_tmp
       rhoK = rhoK_tmp
#endif
       !
       if(filling/=0d0)then
          Dmin = minval(Ek); Dmax = maxval(Ek); mu0 = Dmin
          call fzero(get_dens,mu0,Dmax,info,rguess=Dmin+0.5d0*(Dmax-Dmin))
          if(info/=1)then
             write(*,*)"ERROR ss_get_lambda0: fzero returned info>1 ",info
             stop
          endif
          if(master.AND.verbose>3)write(*,"(A6,12G18.9)")"mu0  =",mu0
       else
          mu0 = xmu
       endif
       !
       Eweiss_tmp  = zero;Eweiss=zero
       dens_tmp    = 0d0 ;dens=0d0
       do ik=1+mpi_rank,Nk,mpi_size
          diagRho      = diag(step_fermi(eK(:,ik)-mu0))
          Uk_f         = rhoK(:,:,ik)
          rhoK(:,:,ik) = (Uk_f .x. diagRho) .x. (conjg(transpose(Uk_f)))
          Wtk = ss_Wtk(1,ik) ; if(is_dos)Wtk = diag(ss_Wtk(:,ik))
          Eweiss_tmp = Eweiss_tmp + ss_Hk(:,:,ik)*rhoK(:,:,ik)*Wtk
          dens_tmp   = dens_tmp + diagonal(rhoK(:,:,ik)*Wtk)
       enddo
#ifdef _MPI
       if(check_MPI())then
          call AllReduce_MPI(MPI_COMM_WORLD,Eweiss_tmp,Eweiss)
          call AllReduce_MPI(MPI_COMM_WORLD,dens_tmp,dens)
       else
          Eweiss = Eweiss_tmp
          dens   = dens_tmp
       endif
#else
       Eweiss = Eweiss_tmp
       dens   = dens_tmp
#endif       
       !
       !< Get H_{a,s} = \sum_{b}sum_k H_{a,s, b,s}*\rho_{a,s, b,s}
       Weiss = 0d0
       do ispin=1,2
          do ilat=1,Nlat
             do iorb=1,Norb
                io = ss_Indices2i([iorb,ilat,ispin],nDefOrder)
                do jlat=1,Nlat
                   do jorb=1,Norb
                      jo = ss_Indices2i([jorb,jlat,ispin],nDefOrder)
                      Weiss(io) = Weiss(io) + Eweiss(io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
       !< Get Lambda0 = -2* h0_{m,s}*[n0_{m,s}-0.5]/[n0_{m,s}*(1-n0_{m,s})]
       lambda0 = -2d0*abs(Weiss)*(Dens-0.5d0)/(Dens*(1d0-Dens)+mch)
       !where(dens==1d0.OR.dens==0d0)lambda0=0d0
    endif

    ss_Dens = ss_unpack_array(Dens,Nlat)
    ss_Lambda0 = ss_unpack_array(lambda0,Nlat)
    if(Nspin==1)then
       call ss_spin_symmetry(ss_Dens,Nlat)
       call ss_spin_symmetry(ss_Lambda0,Nlat)
    endif
    !
    xmu = mu0 !+ 2d0*sum(lambda0)/Ns
    !if(filling/=0d0)xmu =  mu0 + 2d0*sum(lambda0)/Ns
    !
    if(master)then
       if(verbose>2)then
          do ineq=1,Nineq
             ilat = ss_ineq2ilat(ineq)
             write(*,"(A7,12G18.9)")"Lam0  =",ss_lambda0(ilat,:Nspin*Norb)
             write(*,"(A6,12G18.9)")"N0    =",ss_Dens(ilat,:Nspin*Norb),sum(ss_dens(Ilat,:Nspin*Norb))*(3-Nspin),filling
          enddo
       endif
       do ilat=1,Nlat
          open(free_unit(unit),file="lambda0_site"//str(ilat,4)//".ss")
          write(unit,*)ss_lambda0(ilat,:)
          close(unit)
          !
          open(free_unit(unit),file="N0_site"//str(ilat,4)//".ss")
          write(unit,*)ss_dens(ilat,:)
          close(unit)
       enddo
       call save_array(trim(Pfile)//"0"//trim(ss_file_suffix)//".ss",[lambda0,dens,mu0])
    endif
    !

  contains
    !
    function get_dens(mu) result(dens)
      real(8),intent(in)          :: mu
      real(8)                     :: dens
      real(8),dimension(Ns)       :: ndens,ndens_tmp
      ndens_tmp = 0d0
      ndens     = 0d0
      do ik=1+mpi_rank,Nk,mpi_size
         diagRho = diag(step_fermi(eK(:,ik)-mu))
         Uk_f    = rhoK(:,:,ik)
         Rho     = (Uk_f .x. diagRho) .x. (conjg(transpose(Uk_f)))
         Wtk     = ss_Wtk(1,ik) ; if(is_dos)Wtk = diag(ss_Wtk(:,ik))
         ndens_tmp = ndens_tmp + diagonal(Rho*Wtk)
      enddo
#ifdef _MPI
      if(check_MPI())then
         call AllReduce_MPI(MPI_COMM_WORLD,ndens_tmp,ndens)
      else
         ndens   = ndens_tmp
      endif
#else
      ndens   = ndens_tmp
#endif
      if(master.AND.verbose>3)write(*,"(A9,3G18.9)")"Ef,N0   =",mu,sum(ndens),filling
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







