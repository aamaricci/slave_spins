program ss_PCO
  USE SLAVE_SPINS
  !
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                 :: Nk,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb,io,ispin
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk
  complex(8),dimension(:,:),allocatable   :: Hloc
  character(len=32)                       :: w90_file

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

  !METHOD 1 (setup W90 --> use internal W90 model)
  call TB_w90_setup(reg(w90_file),nlat=1,nspin=Nspin,norb=Norb)

  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nk=Nkx*Nkx*Nkx ;   write(*,*) "Using Nk_total="//txtfy(Nk)
  !
  allocate(Hk(Nspin*Norb,Nspin*Norb,Nk))
  allocate(Wtk(Nk))
  allocate(Hloc(Nspin*Norb,Nspin*Norb))
  !
  call TB_build_model(Hk,TB_w90_model,Nspin*Norb,[Nkx,Nkx,Nkx])
  Wtk = 1d0/Nk
  Hloc= sum(Hk(:,:,:),dim=3)/Nk
  where(abs(Hloc)<1d-6)Hloc=zero
  forall(ik=1:Nk)Hk(:,:,ik) = Hk(:,:,ik) - Hloc !set the Fermi level to 0d0
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
  call ss_init(Hk,Wtk)

  !> SS SOLVE:
  call ss_solve()



  call TB_w90_delete()


end program SS_PCO


