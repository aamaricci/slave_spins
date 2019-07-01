c     Code DFT+SS. 
      integer i,j,l,ii,jj,ll 
      integer indx
      integer iorb

      real*8 Threshold

      integer ciclo
      logical check

      
      real*8 thop

c     variables related to Fermionic hamiltonian:

      integer num_wann,nspin,ndim
      parameter (num_wann=1,nspin=2,ndim=2)

      integer loop_kpt
      integer num_kpts
      integer Nkx,Nky,Nkz
      parameter(Nkx=1500,Nky=1500,Nkz=1)
      parameter (num_kpts=Nkx*Nky*Nkz)
      real*8 kpt_latt(ndim,num_kpts)


      complex*16 Jijsigma(num_wann,num_wann,nspin)
      complex*16 fdagf(num_wann,num_wann)

      real*8 fnorm
      complex*16 factexp

      real lambda(num_wann)
      real lambdawrk(num_wann)

      integer Npara
      real x(num_wann),f(num_wann)

      real*8 n_state(num_wann)

      complex*16 sumaux,sumauxb

      complex*16 hdiag(num_wann,num_wann)
      complex*16 ham_k(num_wann,num_wann,num_kpts)
      complex*16 ham_k_in(num_wann,num_wann,num_kpts)
      real*8 eneW,enew2


c     ZHEEVX 
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
     &     ,num_wann,num_wann)

      integer e_ind(num_wann*num_kpts)
      integer Ntot
      integer N_e_k
      double precision EF

      real filling 

      integer ie,mne,if_Ldos
      parameter (mne=1024)
      real*8 LdosE(num_wann,mne),emin,emax,sumd
      real*8 dosE(mne)
      real*8 egrid(mne)
      real*8 gamma
      real*8 pi

c     variables related to Spin hamiltonian:
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

      real*8 Uprime,Uvalue
      logical bool,bool1,bool2,bool3,bool4

      real*8 sumup,sumdn,sumup1,sumup2,sumdn1,sumdn2,sum1,sum2,sum3
      real*8 sumtot,sumaux1,sumaux2
      complex*16 sum1aux,sum2aux,sum3aux,sum4aux

      parameter (N_site_tot=num_wann) ! # orbitals in the single site
!     it must be: N_site_tot=num_wann
      complex*16 Oi(N_site_tot,nspin) ! <O>
      complex*16 Opi(N_site_tot,nspin) ! <O^+> 
      complex*16 mOi(N_site_tot*nspin)
      complex*16 mOpi(N_site_tot*nspin)
      complex*16 mOi_orb(num_wann*nspin)
      complex*16 mOpi_orb(num_wann*nspin)
      complex*16 mOi_orb_ini(num_wann*nspin) ! initial <O> 
      complex*16 mOpi_orb_ini(num_wann*nspin) ! initial <O^+>

      real*8 cv
      complex*16 cm(N_site_tot,nspin) ! gauge variables

      complex*16 factorHM

      integer NmaxNs_up,NmaxNs_dn
      parameter (NmaxNs_up=2**N_site_tot,NmaxNs_dn=2**N_site_tot) ! # of up/dn states in Hilbert space
      integer I_up(NmaxNs_up),I_dn(NmaxNs_dn) ! possible states in Hilbert space
      
      integer NmaxNstatesHS
      parameter (NmaxNstatesHS = 2**(2*N_site_tot)) ! dimension of Hilbert space

      complex*16 Sz(nspin,N_site_tot) ! <Sz>
      real mSz(nspin*N_site_tot)
      real mSz_orb(num_wann*nspin)

      character*8 x1,x2
      character*30 filename

      complex*16 sumsum

c     ZHEEVX LAPACK:
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

c     time variables:
      integer*4 today(3),now(3)

c     Common variables:

      common/variables/ciclo,if_ldos,Uvalue,filling,thop,
     &     Threshold,ham_k_in,
     &     n_state,Jijsigma,
     &     mSz,mOi,mOpi,
     &     mSz_orb,mOi_orb,mOpi_orb,
     &     mOi_orb_ini,mOpi_orb_ini

