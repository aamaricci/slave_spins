MODULE SS_SOLVE_FERMION
  USE SS_VARS_GLOBAL
  USE SS_SETUP
  !
  USE SF_LINALG, only: outerprod
  USE SF_OPTIMIZE,only: brentq,fzero
  USE SF_MISC,    only: sort_array
  implicit none

  private

  public :: ss_solve_lambda0
  public :: ss_update_lambda0
  public :: ss_solve_fermions


  real(8),parameter                  :: mch=1d-6
  integer                            :: iorb,jorb,ispin,jspin,io,jo,ilat,jlat,ineq,ii,jj


contains




  subroutine ss_solve_lambda0()
    complex(8),dimension(Nlso,Nlso)    :: Uk_f
    real(8),dimension(Nlso)            :: Ek_f
    complex(8),dimension(Nlso,Nlso)    :: Rho
    real(8),dimension(Nlso,Nlso)       :: Wtk
    !
    real(8),dimension(Nlso,Nk)         :: eK,eK_tmp
    complex(8),dimension(Nlso,Nlso,Nk) :: rhoK,rhoK_tmp
    !
    complex(8),dimension(Nlso,Nlso)    :: FdgF,FdgF_tmp
    complex(8),dimension(Nlso,Nlso)    :: Hfk,Hfk_tmp
    complex(8),dimension(Nlso,Nlso)    :: Hfloc
    !
    real(8),dimension(Nlso)            :: dens
    real(8),dimension(Nlso)            :: lambda0
    complex(8),dimension(Nlso)         :: heff
    complex(8),dimension(Nlso,Nlso)    :: jhybr    
    !
    real(8)                            :: mu0,Dmin,Dmax
    integer                            :: ik,unit,N
    integer                            :: info
    logical                            :: IOfile
    integer                            :: Len
    real(8),dimension(:),allocatable   :: params
    !
#ifdef _MPI
    if(check_MPI())then
       mpi_rank=get_rank_MPI()
       mpi_size=get_size_MPI()
    endif
#endif
    !
    inquire(file=trim(Pfile)//"0"//trim(ss_file_suffix)//".restart",exist=IOfile)
    FILE_READ: if(IOfile)then
       len = file_length(trim(Pfile)//"0"//trim(ss_file_suffix)//".restart")
       if( (len/=2*Nlso+1) )stop "SS_SOLVE_LAMBDA0 ERROR: len!=2*Ns OR 2*Ns+1"
       allocate(params(len))
       call read_array(trim(Pfile)//"0"//trim(ss_file_suffix)//".restart",params)
       lambda0 = params(1:Nlso)
       dens    = params(Nlso+1:2*Nlso)
       xmu     = params(2*Nlso+1)
       deallocate(params)
       !
       ss_Dens    = ss_unpack_array(Dens,Nlat)
       ss_Lambda0 = ss_unpack_array(lambda0,Nlat)
       !
    else
       !
       !
       !< Solve non-interacting problem and store solution
       ek_tmp=0d0   ;ek=0d0
       rhoK_tmp=zero;rhoK=zero
       do ik=1+mpi_rank,Nk,mpi_size
          Uk_f = ss_Hk(:,:,ik) + ss_Hloc + ss_Hhyb + diag(ss_Hdiag)
          call eigh(Uk_f,Ek_f)
          eK_tmp(:,ik)     = Ek_f
          rhoK_tmp(:,:,ik) = Uk_f
       enddo
       !< MPI dispatch 
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
       !< Solve non-linear problem to find correct chemical potential for
       !< given filling, if any.  
       if(filling/=0d0)then
          Dmin = minval(Ek); Dmax = maxval(Ek); mu0 = Dmin
          call fzero(get_dens,mu0,Dmax,info,rguess=Dmin+0.5d0*(Dmax-Dmin))
          if(info/=1)then
             write(*,*)"ERROR ss_get_lambda0: fzero returned info>1 ",info
             stop
          endif
          if(master.AND.verbose>3)write(*,"(A6,12G18.9)")"mu0  =",mu0
          xmu = mu0
       endif
       !
       !
       !< Get Rho^{ab} = sum_k \rho^{ab}_k = sum_k <fdg^a_k f^b_k>
       !< Get h^{ab} = sum_k H^{ab}_k*\rho^{ab}_k
       Hfk_tmp  = 0d0 ; Hfk = 0d0
       FdgF_tmp = 0d0 ; FdgF= 0d0       
       do ik=1+mpi_rank,Nk,mpi_size
          Wtk  = ss_Wtk(1,ik) ; if(is_dos)Wtk = diag(ss_Wtk(:,ik))
          Uk_f = rhoK(:,:,ik)             !Rotation matrix
          Ek_f = step_fermi(eK(:,ik)-xmu) !bands> Fermi function
          forall(io=1:Nlso,jo=1:Nlso)&    !build up <f^+f>
               Rho(io,jo) = sum(Uk_f(jo,:)*conjg(Uk_f(io,:))*EK_f)
          !< Get Rho
          FdgF_tmp = FdgF_tmp + Rho*Wtk
          !< Get inter-cell contribution to effective/Weiss field
          Hfk_tmp  = Hfk_tmp + (ss_Hk(:,:,ik)*Rho)*Wtk !dreal(ss_Hk(:,:,ik)*Rho)*Wtk
       enddo
       !< MPI reduction where applicable
#ifdef _MPI
       if(check_MPI())then
          call AllReduce_MPI(MPI_COMM_WORLD,FdgF_tmp,FdgF)            
          call AllReduce_MPI(MPI_COMM_WORLD,Hfk_tmp,Hfk)
       else
          FdgF = FdgF_tmp
          Hfk  = Hfk_tmp
       endif
#else
       FdgF = FdgF_tmp
       Hfk  = Hfk_tmp
#endif
       !
       !< Get occupation n^a = Rho^{aa}
       dens = dreal(diagonal(FdgF))
       !
       !< Get intra-cell, non-local contribution to Weiss field h'
       !< h'^{ab} = sum_k H^{ab}_{loc,\mu_i!=\nu_i} \rho^{ab}_k
       !          = H^{ab}_{loc,\mu!=\nu} Rho^{ab}
       Hfloc= ss_Hloc*FdgF
       !
       !< sum up Weiss field contributions. 
       !< H^a = \sum_{b} [h^{ab} + h'^{ab}]
       Heff = zero
       do ispin=1,Nspin
          do ilat=1,Nlat
             do iorb=1,Norb
                io = ss_Indices2i([iorb,ilat,ispin],[Norb,Nlat,Nspin])
                do jspin=1,Nspin
                   do jlat=1,Nlat
                      do jorb=1,Norb
                         jo = ss_Indices2i([jorb,jlat,jspin],[Norb,Nlat,Nspin])
                         Heff(io) = Heff(io) + Hfk(io,jo) + Hfloc(io,jo)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
       !< Get intra-cell, local off-diagonal contribution to Weiss field
       !< J^ab = sum_k H^{a!=b}_{loc} \rho^{ab}_k
       !       = H^{a!=b}_{loc} Rho^{ab}
       Jhybr = ss_Hhyb*FdgF
       !
       !< Get Lambda0 = -2*|h0|*[n0-0.5]/[n0*(1-n0)]
       lambda0 = -2d0*abs(heff)*(Dens-0.5d0)/(Dens*(1d0-Dens)+mch)
       !
       !
       !< Dump to SS global parameters
       ss_Dens    = ss_unpack_array(Dens,Nlat)
       ss_Lambda0 = ss_unpack_array(lambda0,Nlat)
       ss_Heff    = ss_unpack_array(heff,Nlat)
       ss_Jhybr   = ss_unpack_array(jhybr,Nlat)
    endif FILE_READ
    !

    if(master.AND.verbose>2)then
       do ineq=1,Nineq
          ilat = ss_ineq2ilat(ineq)
          write(*,"(A7,12G18.9)")"Lam0  =",ss_lambda0(ilat,:)
          write(*,"(A6,12G18.9)")"N0    =",ss_Dens(ilat,:)
          if(verbose>3)write(*,"(A6,12G18.9)")"Weiss =",ss_Heff(ilat,:)
       enddo
    endif
    !
    if(master)then
       do ilat=1,Nlat
          open(free_unit(unit),file="lambda0_site"//str(ilat,4)//".ss")
          write(unit,*)ss_lambda0(ilat,:)
          close(unit)
          !
          open(free_unit(unit),file="N0_site"//str(ilat,4)//".ss")
          write(unit,*)ss_dens(ilat,:)
          close(unit)
       enddo
       call save_array(trim(Pfile)//"0"//trim(ss_file_suffix)//".ss",[lambda0,Dens,xmu])
    endif
    !
  contains
    !
    !
    function get_dens(mu) result(dens)
      real(8),intent(in)      :: mu
      real(8)                 :: dens
      real(8),dimension(Nlso) :: ndens,ndens_tmp
      real(8),dimension(Nlso) :: Rtmp,wt,Efvec
      ndens_tmp = 0d0
      ndens     = 0d0
      do ik=1+mpi_rank,Nk,mpi_size
         Wt    = ss_Wtk(1,ik) ; if(is_dos)Wt = ss_Wtk(:,ik)
         Efvec = step_fermi(eK(:,ik)-mu)
         forall(io=1:Nlso)&
              Rtmp(io) = sum( abs(rhoK(io,:,ik))**2*Efvec*wt )
         ndens_tmp = ndens_tmp + Rtmp*(3-Nspin)
      enddo
      !
#ifdef _MPI
      if(check_MPI())then
         call AllReduce_MPI(MPI_COMM_WORLD,ndens_tmp,ndens)
      else
         ndens   = ndens_tmp
      endif
#else
      ndens   = ndens_tmp
#endif
      !
      if(master.AND.verbose>3)write(*,"(A9,3G18.9)")"Ef,N0   =",mu,sum(ndens),filling
      !
      dens = sum(ndens)-filling
      return
    end function get_dens
  end subroutine ss_solve_lambda0



  subroutine ss_update_lambda0()
    real(8),dimension(Nlso)         :: dens
    real(8),dimension(Nlso)         :: lambda0
    real(8),dimension(Nlso)         :: Op
    real(8),dimension(Nlso)         :: Heff
    !< Get Lambda0 = -2*|h0|*[n0-0.5]/[n0*(1-n0)]
    Dens    = ss_pack_array(ss_Dens,Nlat)
    Op      = ss_pack_array(ss_Op,Nlat)
    Heff    = ss_pack_array(ss_Heff,Nlat)
    !
    lambda0 = -2d0*Op*abs(heff)*(Dens-0.5d0)/(Dens*(1d0-Dens)+mch)
    !
    ss_lambda0 = ss_unpack_array(lambda0,Nlat)
  end subroutine ss_update_lambda0



  subroutine ss_solve_fermions()
    complex(8),dimension(Nlso,Nlso) :: Hk_f
    complex(8),dimension(Nlso,Nlso) :: Uk_f
    real(8),dimension(Nlso)         :: Ek_f
    complex(8),dimension(Nlso,Nlso) :: Rho
    real(8),dimension(Nlso,Nlso)    :: Wtk
    !
    complex(8),dimension(Nlso,Nlso) :: FdgF,FdgF_tmp
    complex(8),dimension(Nlso,Nlso) :: Hfk,Hfk_tmp
    complex(8),dimension(Nlso,Nlso) :: Hfloc
    !
    real(8),dimension(Nlso)         :: dens
    real(8),dimension(Nlso)         :: lambda
    real(8),dimension(Nlso)         :: lambda0
    real(8),dimension(Nlso)         :: const      
    real(8),dimension(Nlso)         :: Op
    complex(8),dimension(Nlso)      :: heff
    complex(8),dimension(Nlso,Nlso) :: jhybr
    complex(8),dimension(Nlso,Nlso) :: OdgOp
    !
    integer                         :: ik,i,j,N,io,jo
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
    lambda  = ss_pack_array(ss_Lambda,Nlat)
    lambda0 = ss_pack_array(ss_Lambda0,Nlat)
    Op      = ss_pack_array(ss_Op,Nlat)
    OdgOp   = ss_pack_array(ss_OdgOp,Nlat)
    !
    if(master.AND.verbose>5)then
       write(*,"(A,100G18.9)")"Op   =",Op
       write(*,"(A)")"OdgOp="
       do io=1,Nlso
          write(*,"(100G11.4)")( OdgOp(io,jo),jo=1,Nlso )
       enddo
    endif
    !
    !< Get Rho^{ab} = sum_k \rho^{ab}_k = sum_k <fdg^a_k f^b_k>
    !< Get h^{ab} = sum_k H^{ab}_k*\rho^{ab}_k
    Hfk_tmp  = 0d0 ;Hfk =0d0
    FdgF_tmp = 0d0 ;FdgF=0d0
    do ik=1+mpi_rank,Nk,mpi_size
       Wtk= ss_Wtk(1,ik) ; if(is_dos)Wtk = diag(ss_Wtk(:,ik))
       !
       !< Construct the renormalized Hamiltonian:
       !< H(k):non-local,inter-cell and *Hloc:non-local intra-cell as
       !< <O_a>.[H_{ab}(k)+Hloc_{ab}].<O_b>
       !< Hhyb:local,intra-cell, off-diagonal as
       !< <O^+_a O_b>.Hhyb_{ab}
       !< Hdiag:local,intra-cell, diagonal remains identical (couples to density).
       Hk_f  = ( ss_Hk(:,:,ik) + ss_Hloc )*outerprod(Op,Op) + ss_Hhyb*OdgOp 
       Uk_f  = Hk_f + diag(ss_Hdiag - lambda + lambda0)
       call eigh(Uk_f,Ek_f)
       forall(io=1:Nlso,jo=1:Nlso)& !build up <f^+f>
            Rho(io,jo) = sum(Uk_f(jo,:)*conjg(Uk_f(io,:))*step_fermi(Ek_f-xmu))
       !< Get Rho
       FdgF_tmp = FdgF_tmp  + Rho*Wtk
       !< Get inter-cell contribution to effective field
       Hfk_tmp  = Hfk_tmp   + (ss_Hk(:,:,ik)*Rho*Wtk) !dreal(ss_Hk(:,:,ik)*Rho*Wtk)
    enddo
    !< MPI reduction where applicable
#ifdef _MPI
    if(check_MPI())then
       call AllReduce_MPI(MPI_COMM_WORLD,Hfk_tmp,Hfk)
       call AllReduce_MPI(MPI_COMM_WORLD,FdgF_tmp,FdgF)
    else
       Hfk  = Hfk_tmp
       FdgF = FdgF_tmp
    endif
#else
    Hfk  = Hfk_tmp
    FdgF = FdgF_tmp
#endif
    !
    !< Get occupation n^a = Rho^{aa}
    dens = diagonal(FdgF)
    !
    !< BUILD UP EFFECTIVE FIELDS FOR THE SPIN PROBLEM:
    !< Get intra-cell, non-local contribution to Weiss field h'
    !< h'^{ab} = sum_k H^{ab}_{loc,\mu_i!=\nu_i} \rho^{ab}_k
    !          = H^{ab}_{loc,\mu!=\nu} Rho^{ab}
    Hfloc= ss_Hloc*FdgF
    !
    !< sum up Weiss field contributions. 
    !< H^a = \sum_{b} [h^{ab} + h'^{ab}]
    Heff = 0d0
    do ispin=1,Nspin
       do ilat=1,Nlat
          do iorb=1,Norb
             io = ss_Indices2i([iorb,ilat,ispin],[Norb,Nlat,Nspin])
             do jspin=1,Nspin
                do jlat=1,Nlat
                   do jorb=1,Norb
                      jo = ss_Indices2i([jorb,jlat,jspin],[Norb,Nlat,Nspin])
                      Heff(io) = Heff(io) + Op(jo)*Hfk(io,jo) + Op(jo)*Hfloc(io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !< Get intra-cell, local off-diagonal contribution to Weiss field
    !< J^ab = sum_k H^{a!=b}_{loc} \rho^{ab}_k
    !       = H^{a!=b}_{loc} Rho^{ab}
    Jhybr = ss_Hhyb*FdgF
    ! !>>DEBUG
    ! print*,"Hhyb:"
    ! do io=1,Nlso
    !    write(*,"(100G11.2)")( ss_Hhyb(io,jo),jo=1,Nlso )
    ! enddo
    ! print*,"FdgF:"
    ! do io=1,Nlso
    !    write(*,"(100G11.2)")( FdgF(io,jo),jo=1,Nlso )
    ! enddo
    ! !<<DEBUG
    !
    !
    !< Get C = (n_{l,s}*(1-n_{l,s}))**{-1/2} - 1, at half-filling C=1
    Const  = 1d0/(sqrt(Dens*(1d0-Dens))+mch) - 1d0
    !
    !< Dump to SS global parameters
    ss_Dens = ss_unpack_array(dens,Nlat)
    ss_C    = ss_unpack_array(const,Nlat)
    ss_Heff = ss_unpack_array(Heff,Nlat)
    ss_Jhybr= ss_unpack_array(Jhybr,Nlat)
    do ineq=1,Nineq
       ilat = ss_ineq2ilat(ineq)
       ss_Dens_ineq(ineq,:)    = ss_Dens(ilat,:)
       ss_C_ineq(ineq,:)       = ss_C(ilat,:)
       ss_Heff_ineq(ineq,:)    = ss_Heff(ilat,:)
       ss_Jhybr_ineq(ineq,:,:) = ss_Jhybr(ilat,:,:)
    enddo
    !
    if(master.AND.verbose>3)then
       do ineq=1,Nineq
          write(*,"(A,12G18.9)")"C    =",ss_C_ineq(ineq,:)
          write(*,"(A,12G18.9)")"N    =",ss_Dens_ineq(ineq,:)
          write(*,"(A,12G18.9)")"Heff =",ss_Heff_ineq(ineq,:)
          if(verbose>4)then
             write(*,"(A)")"Jhyb ="
             do io=1,Nso
                write(*,"(100G11.2)")( ss_Jhybr_ineq(ineq,io,jo),jo=1,Nso )
             enddo
          endif
       enddo
    endif
    !
  end subroutine ss_solve_fermions








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







