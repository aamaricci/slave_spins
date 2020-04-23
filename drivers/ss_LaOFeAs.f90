program ss_LaOFeAs
  USE SLAVE_SPINS
  !
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none


  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts,Nlso
  integer                                 :: i,j,k,ik,ilat,iorb,jorb,io,ispin
  real(8),dimension(3)                    :: e1,e2,e3
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk
  complex(8),dimension(:,:),allocatable   :: Hloc
  real(8),allocatable                     :: Dens(:),Zeta(:), Self(:)
  character(len=20)                       :: w90file,InputFile
  real(8)                                 :: ef=1d0
  logical                                 :: master=.true.,bool

#ifdef _MPI
  call init_MPI
  master = get_master_MPI()
#endif


  call parse_cmd_variable(InputFile,"INPUTFILE",default="inputLaOFeAs.conf")
  call parse_input_variable(w90file,"w90file",InputFile,default="pwscf_LaOFeAs.w90")
  call parse_input_variable(Nkx,"NKX",InputFile,default=10)
  call parse_input_variable(nkpath,"NKPATH",InputFile,default=500)
  call ss_read_input(reg(InputFile))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  Nlso = Nlat*Nspin*Norb

  !>Direct lattice basis vector
  e1 = [7.625044946193571d0,0d0,0d0]
  e2 = [0d0,7.625044946193571d0,0d0]
  e3 = [0d0,0d0,16.51809612755341d0]  
  call TB_set_ei(e1,e2,e3)
  call TB_build_bk(verbose=.true.)

  !METHOD 1 (setup W90 --> use internal W90 model)
  call TB_w90_setup(reg(w90file),nlat=Nlat,nspin=Nspin,norb=Norb,verbose=.true.)
  call TB_w90_FermiLevel([Nkx,Nkx,Nkx],filling,Ef)

  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx*Nkx ;   write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nlso,Nlso,Nktot))
  ! call set_store_size(200)
  ! inquire(file="LaOFeAs_hk.in",exist=bool)
  ! if(bool)then
  !    call read_array("LaOFeAs_hk.in",Hk)
  ! else

  call start_timer
  call TB_build_model(Hk,TB_w90_model,Nlso,[Nkx,Nkx,Nkx])
  call stop_timer("TB_build_model")
  ! call save_array("LaOFeAs_hk.in",Hk)
  ! endif


  allocate(Hloc(Nlso,Nlso))
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  if(master)call TB_write_Hloc(Hloc,"w90Hloc.dat")

  !SOLVE SS:
  call ss_solve(Hk,ineq_sites=[1,1])


  !Retrieve Zeta and ReSigma(0)=lambda0-lambda
  allocate(Zeta(Nlat*Nspin*Norb))
  allocate(Self(Nlat*Nspin*Norb))
  call ss_get_zeta(zeta)
  call ss_get_Self(self)

  !Push em to Wannier 90 setup
  call TB_w90_Zeta(zeta)
  call TB_w90_Self(diag(self))

  ! !solve along a path in the 3D BZ.
  Npts = 9
  Nk=(Npts-1)*Nkpath
  allocate(kpath(Npts,3))
  kpath(1,:)=[0.5d0,0.5d0,0d0]
  kpath(2,:)=[0.5d0,0.5d0,0.5d0]
  kpath(3,:)=[0d0,0d0,0d0]
  kpath(4,:)=[0.5d0,0d0,0d0]
  kpath(5,:)=[0.5d0,0.5d0,0d0]
  kpath(6,:)=[0d0,0d0,0d0]
  kpath(7,:)=[0d0,0d0,0.5d0]
  kpath(8,:)=[0.5d0,0d0,0.5d0]
  kpath(9,:)=[0.5d0,0.5d0,0.5d0]

  !Solve for the renormalized bands:
  if(master)call TB_Solve_model(TB_w90_model,Nlso,kpath,Nkpath,&
       colors_name=[black,red,red,green,blue,black,red,red,green,blue],&
       points_name=[character(len=40) ::'M', 'R', 'G', 'X', 'M', 'G', 'Z','A', 'R'],&
       file="zBands_LaOFeAs",iproject=.true.)


  call TB_w90_delete()


#ifdef _MPI
  call finalize_MPI()
#endif

end program ss_LaOFeAs
