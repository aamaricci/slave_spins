c Code DFT+SS. 
      integer i,j,l,ii,jj,ll 
      integer indx
      integer iorb

      integer whatdoyourun
      integer whattypeofJh
      integer whattypeofdiago
      integer ifdoublecounting
      real*8 Threshold

      integer ciclo
      logical check

c      real*8 EtotS
c      real*8 EtotF
c      real*8 EtotFS
c      real*8 Edc

c variables related to Fermionic hamiltonian:

      integer num_wann,num_wann_ss,nspin,ndim
      parameter (num_wann=12,num_wann_ss=12,nspin=2,ndim=3)

      integer Nsite,isite,indx_orb_site

      integer loop_kpt
      integer num_kpts
      parameter (num_kpts=1728)
      real*8 kpt_latt(ndim,num_kpts)

      integer loop_rpt
      integer nrpts
      parameter (nrpts=617)
      real*8 irvec(ndim,nrpts)
      real*8 rdotk
      real*8 a,b,factor_hr
      complex*16 ham_r(num_wann,num_wann,nrpts)
      complex*16 ham_r2(num_wann,num_wann,nrpts)

      integer num_wann_ss_ind
      parameter (num_wann_ss_ind=3)
      integer if_sym
      integer if_sym_updn
      integer index_sym(num_wann_ss)

      complex*16 Jijsigma(num_wann,num_wann,nrpts,nspin)
      complex*16 fdagf(num_wann,num_wann,nrpts)

      real*8 fnorm
      complex*16 factexp

      real lambda(num_wann_ss)
      real lambdawrk(num_wann_ss)
      real e_on_site(num_wann)
      real*8 shifteps(num_wann_ss)

      integer Npara
c      real x(num_wann_ss),f(num_wann_ss)
      real x(num_wann_ss_ind),f(num_wann_ss_ind)

      real*8 n_state(num_wann)

      complex*16 sumaux,sumauxb

      complex*16 ham_k(num_wann,num_wann)
      real*8 eneW 

c ZHEEVX 
      integer neigval
      integer lwork
      integer iwork(5*(num_wann))
      integer ifail(num_wann),info
      parameter (lwork=2*(num_wann))
      double precision  rwork(7*(num_wann))
      double precision  vl,vu
      complex*16 work(lwork)
      double precision eigval(num_wann)
      double precision eigenvalues_k(num_kpts,num_wann)
      double precision eigenvalues(num_wann*num_kpts)
      complex*16 Z(num_wann,num_wann)
      complex*16 Z_k(num_kpts
     &                        ,num_wann,num_wann)

      integer e_ind(num_wann*num_kpts)
      integer Ntot
      integer N_e_k
      double precision EF

      real filling 
      real*8 navhartree

      integer ie,mne,if_Ldos
      parameter (mne=1001)
      real*8 LdosE(num_wann,mne),emin,emax,sumd
      real*8 dosE(mne)
      real*8 egrid(mne)
      real*8 gamma
      real*8 pi

c variables related to Spin hamiltonian:

      integer N_up,N_dn
      integer Ns_up,Ns_dn
      integer istate,Nstate
      integer N_site_tot
      integer jorb
      integer J_dn
      integer J_up
      integer count_dn_ini,count_dn_fin,count_up_ini,count_up_fin
      integer Ip_up,Ip_dn
      integer k,iii,jjj,pos
      integer indx_iorb

      real*8 Uprime,Uvalue,Jh
      logical bool,bool1,bool2,bool3,bool4

      real*8 sumup,sumdn,sumup1,sumup2,sumdn1,sumdn2,sum1,sum2,sum3
      real*8 sumtot,sumaux1,sumaux2
      complex*16 sum1aux,sum2aux,sum3aux,sum4aux

      parameter (N_site_tot=3)                      ! # orbitals in the single site
                                                    ! it must be: N_site_tot=num_wann_ss
      complex*16 Oi(N_site_tot,nspin)               ! <O>
      complex*16 Opi(N_site_tot,nspin)              ! <O^+> 
      complex*16 mOi(N_site_tot*nspin)
      complex*16 mOpi(N_site_tot*nspin)
      complex*16 mOi_orb(num_wann_ss*nspin)
      complex*16 mOpi_orb(num_wann_ss*nspin)
      complex*16 mOi_orb_ini(num_wann_ss*nspin)     ! initial <O> 
      complex*16 mOpi_orb_ini(num_wann_ss*nspin)    ! initial <O^+>

      real*8 cv
      complex*16 cm(N_site_tot,nspin)               ! gauge variables

      complex*16 factorHM

      integer NmaxNs_up,NmaxNs_dn
      parameter (NmaxNs_up=2**N_site_tot,NmaxNs_dn=2**N_site_tot) ! # of up/dn states in Hilbert space
      integer I_up(NmaxNs_up),I_dn(NmaxNs_dn)       ! possible states in Hilbert space
 
      integer NmaxNstatesHS
      parameter (NmaxNstatesHS = 2**(2*N_site_tot)) ! dimension of Hilbert space

      complex*16 Sz(nspin,N_site_tot)               ! <Sz>
      real mSz(nspin*N_site_tot)
      real mSz_orb(num_wann_ss*nspin)

      complex*16 SzSz(2*nspin,N_site_tot,N_site_tot)
      complex*16 mSzSz_iorb_jorb(2*nspin,num_wann_ss,num_wann_ss)
      character*8 x1,x2
      character*30 filename

      complex*16 sumsum

c ZHEEVX LAPACK:
      integer neig
      integer llwork
      integer iiwork(5*(NmaxNstatesHS))
      integer iifail(NmaxNstatesHS),iinfo
      parameter (llwork=2*(NmaxNstatesHS))
      double precision  rrwork(7*(NmaxNstatesHS)),ivl,ivu
      complex*16 wwork(llwork)
      double precision eig(NmaxNstatesHS)
      complex*16 ZZ(NmaxNstatesHS,NmaxNstatesHS)
      complex*16 ZZ_gs(NmaxNstatesHS)
      complex*16 haccaS(NmaxNstatesHS,NmaxNstatesHS)

      integer Ndegen

c LANCZOS:
      integer iterL,Niter,Nerr,NiterL
      integer Ncheck
      parameter (Ncheck = 50)
      integer NmaxiterL
      parameter (NmaxiterL = 200)
      double precision size,norm
      double precision eThresholdLanczos,diff,esave(NmaxiterL)
      double precision TA(NmaxiterL),TB(NmaxiterL+1),TA1,TB1
      double precision Zaux(NmaxiterL,NmaxiterL)
      double precision diag(NmaxiterL),subdiag(NmaxiterL)
      complex*16 Qlin(NmaxNstatesHS),Qin(NmaxNstatesHS)
      complex*16 tmp(NmaxNstatesHS)
      complex*16 Qsave(NmaxNstatesHS)
      complex*16 ZA(NmaxNstatesHS)
c Random generator variables:
      real rand
      double precision sirand1,sirand2
      external rand

c time variables:
      integer*4 today(3),now(3)

c Common variables:

      common / / ciclo,whatdoyourun
     &           ,if_Ldos,if_sym,if_sym_updn
     &           ,index_sym
     &           ,ifdoublecounting,navhartree
     &           ,Uvalue,Jh,filling
     &           ,Threshold,eThresholdLanczos
     &           ,whattypeofJh,whattypeofdiago
     &           ,factor_hr,shifteps
     &           ,n_state,Jijsigma
     &           ,mSz,mOi,mOpi
     &           ,mSz_orb,mOi_orb,mOpi_orb
     &           ,mOi_orb_ini,mOpi_orb_ini
     &           ,SzSz
     &           ,mSzSz_iorb_jorb
     &           ,kpt_latt,ham_r2,irvec
c     &           ,EtotF,EtotS,EtotFS,Edc
