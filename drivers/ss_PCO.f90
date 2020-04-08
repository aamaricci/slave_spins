program ss_PCO
  USE SLAVE_SPINS
  !
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                     :: Nk,Nkpath,Nkx,Npts
  integer                                     :: i,j,k,ik,iorb,jorb,io,ispin
  real(8),dimension(:,:),allocatable          :: kpath,Self
  complex(8),dimension(:,:,:),allocatable     :: Hk
  real(8),dimension(:),allocatable            :: Wtk,Zeta
  complex(8),dimension(:,:),allocatable       :: Hloc
  complex(8),dimension(:,:,:,:,:),allocatable :: Greal,Gmats
  character(len=32)                           :: w90_file
  real(8)                                     :: ef

  call parse_input_variable(Nkx,"NKX","inputPCO.conf",default=10)
  call parse_input_variable(nkpath,"NKPATH","inputPCO.conf",default=500)
  call parse_input_variable(w90_file,"W90_FILE","inputPCO.conf",default="W90_hr_bulk.w90")
  call ss_read_input('inputPCO.conf')

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")


  !>Direct lattice basis vector
  call TB_set_ei([1d0,0d0,0d0],[0d0,1d0,0d0],[0d0,0d0,1d0])
  call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])

  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nk=Nkx*Nkx*Nkx ;   write(*,*) "Using Nk_total="//txtfy(Nk)

  call TB_w90_setup(reg(w90_file),nlat=1,nspin=Nspin,norb=Norb,verbose=.true.)
  call TB_w90_FermiLevel([Nkx,Nkx,Nkx],dble(Nspin*Norb),Ef)


  allocate(Hk(Nspin*Norb,Nspin*Norb,Nk))
  allocate(Wtk(Nk))
  allocate(Hloc(Nspin*Norb,Nspin*Norb))
  !
  call TB_build_model(Hk,Nspin*Norb,[Nkx,Nkx,Nkx])
  Wtk = 1d0/Nk
  Hloc= sum(Hk(:,:,:),dim=3)/Nk
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc,"Hloc_PCO.dat")

  !solve along the standard path in the 3D BZ.
  Npts = 8
  Nk=(Npts-1)*Nkpath
  allocate(kpath(Npts,3))
  kpath(1,:)=[0,0,0]*pi
  kpath(2,:)=[1,0,0]*pi
  kpath(3,:)=[1,1,0]*pi
  kpath(4,:)=[0,0,0]*pi
  kpath(5,:)=[1,1,1]*pi
  kpath(6,:)=[1,0,1]*pi
  kpath(7,:)=[0,0,1]*pi
  kpath(8,:)=[0,0,0]*pi
  call TB_Solve_model(TB_w90_model,Nspin*Norb,kpath,Nkpath,&
       colors_name=[red1,green1,blue1],&
       points_name=[character(len=20) ::'G', 'X', 'M', 'G', 'R', 'A', 'Z','G'],&
       file="Eigenbands_PCO")


  !> SS INIT:
  call ss_solve(Hk)


  !  TODO: get the GF of the problem by convoluting the spin-spin correlation function
  ! with the fermionic (free) green's function. 

  call ss_get_Hf(Hk)
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  call get_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,Lreal))
  call print_gf_realaxis(Greal,"w90Gloc",iprint=1)

  allocate(Zeta(Nspin*Norb))
  allocate(Self(Nspin*Norb,Nspin*Norb))
  call ss_get_zeta(zeta)
  call ss_get_self(self)


  call TB_w90_Zeta(zeta)
  call TB_w90_Self(self)
  call TB_Solve_model(TB_w90_model,Nspin*Norb,kpath,Nkpath,&
       colors_name=[red1,green1,blue1],&
       points_name=[character(len=20) ::'G', 'X', 'M', 'G', 'R', 'A', 'Z','G'],&
       file="zEigenbands_PCO")

  call TB_w90_delete()


end program SS_PCO


