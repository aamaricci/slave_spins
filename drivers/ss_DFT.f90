program ss_DFT
  USE SLAVE_SPINS
  !
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                 :: Nktot,Nkpath,Nkvec(3),Npts,Nlso
  integer                                 :: i,j,k,ik,ilat,iorb,jorb,io,ispin
  real(8),dimension(3)                    :: e1,e2,e3
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk
  complex(8),dimension(:,:),allocatable   :: Hloc
  real(8),allocatable                     :: Dens(:),Zeta(:), Self(:),Tmp(:)
  character(len=60)                       :: w90file,InputFile,latfile,kpathfile,ineqfile,hkfile
  character(len=40),allocatable           :: points_name(:)
  real(8)                                 :: ef
  logical                                 :: FSflag,EFflag
  logical                                 :: master=.true.,bool
  integer                                 :: unit
  integer,allocatable,dimension(:)        :: ineq_sites

#ifdef _MPI
  call init_MPI
  master = get_master_MPI()
#endif


  call parse_cmd_variable(InputFile,"INPUTFILE",default="inputDFT.conf")
  call parse_input_variable(w90file,"w90file",InputFile,default="hij.conf")
  call parse_input_variable(hkfile,"hkfile",InputFile,default="hk.conf")
  call parse_input_variable(latfile,"latfile",InputFile,default="lat.conf")
  call parse_input_variable(kpathfile,"kpathfile",InputFile,default="kpath.conf")
  call parse_input_variable(ineqfile,"ineqfile",InputFile,default="ineq.conf")
  call parse_input_variable(FSflag,"FSflag",InputFile,default=.false.)
  call parse_input_variable(EFflag,"EFflag",InputFile,default=.false.)
  call parse_input_variable(Nkvec,"NKVEC",InputFile,default=[10,10,10])
  call parse_input_variable(nkpath,"NKPATH",InputFile,default=500)
  call ss_read_input(reg(InputFile))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  Nlso = Nlat*Nspin*Norb


  !METHOD 1 (setup W90 --> use internal W90 model)
  inquire(file=reg(latfile),exist=bool)
  if(bool)then
     open(free_unit(unit),file=reg(latfile))
     read(unit,*)e1
     read(unit,*)e2
     read(unit,*)e3
     close(unit)
  else
     e1 = [1d0,0d0,0d0]
     e2 = [0d0,1d0,0d0]
     e3 = [0d0,0d0,1d0]
  endif
  call TB_set_ei(e1,e2,e3)
  call TB_build_bk(verbose=.true.)
  call start_timer
  call TB_w90_setup(reg(w90file),nlat=Nlat,nspin=Nspin,norb=Norb,verbose=.true.)
  call stop_timer("TB_w90_setup")
  !
  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  inquire(file=reg(hkfile),exist=bool)
  if(bool)then
     call TB_read_hk(Hk,reg(hkfile),Nkvec)
     if(size(Hk,1)/=Nlat*Nspin*Norb)stop "ss_DFT error: wrong size in Hk as read from file"
  else
     write(*,"(A)")"get Fermi Level"
     call start_timer  
     call TB_w90_FermiLevel(Nkvec,filling,Ef)
     call stop_timer("TB_w90_FermiLevel")
     !
     Nktot=product(Nkvec)
     allocate(Hk(Nlso,Nlso,Nktot))
     call TB_set_dos_lreal(256)
     call start_timer
     call TB_build_model(Hk,Nlso,Nkvec,wdos=.false.)
     call TB_write_hk(reg(hkfile),Nkvec)
     call stop_timer("TB_build_model")
  endif


  write(*,*)"Using Nk_total="//str(size(Hk,3))

  allocate(Hloc(Nlso,Nlso))
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  if(master)call TB_write_Hloc(Hloc,"w90Hloc.dat")


  allocate(ineq_sites(Nlat));ineq_sites=1
  inquire(file=reg(ineqfile),exist=bool)
  if(bool)then
     open(free_unit(unit),file=reg(ineqfile))
     do i=1,Nlat
        read(unit,*)ineq_sites(i)
        write(*,"(A,I5,A,I5)")"Site",i,"corresponds to ",ineq_sites(i)
     enddo
     close(unit)
  endif

  !SOLVE SS:
  call start_timer
  call ss_solve(Hk,ineq_sites=ineq_sites,UserOrder=[character(len=5)::"Norb","Nlat","Nspin"])
  call stop_timer("SS SOLUTION")



  !Retrieve Zeta and ReSigma(0)=lambda0-lambda
  allocate(Zeta(Nlat*Nspin*Norb))
  allocate(Self(Nlat*Nspin*Norb))
  call ss_get_zeta(zeta)
  call ss_get_Self(self)
  call save_array("renorm.save",[zeta,self])

  !Push em to Wannier 90 setup
  call TB_w90_Zeta(zeta)
  call TB_w90_Self(diag(self))
  deallocate(zeta,self)

  ! !solve along a path in the 3D BZ.
  inquire(file=reg(kpathfile),exist=bool)
  if(bool)then
     Npts = file_length(reg(kpathfile))
     allocate(kpath(Npts,3))
     allocate(points_name(Npts))
     open(free_unit(unit),file=reg(kpathfile))
     do i=1,Npts
        read(unit,*)points_name(i),kpath(i,:)
     enddo
     close(unit)
  else
     Npts = 9
     allocate(kpath(Npts,3),points_name(Npts))
     kpath(1,:)=[0.5d0,0.5d0,0d0]
     kpath(2,:)=[0.5d0,0.5d0,0.5d0]
     kpath(3,:)=[0d0,0d0,0d0]
     kpath(4,:)=[0.5d0,0d0,0d0]
     kpath(5,:)=[0.5d0,0.5d0,0d0]
     kpath(6,:)=[0d0,0d0,0d0]
     kpath(7,:)=[0d0,0d0,0.5d0]
     kpath(8,:)=[0.5d0,0d0,0.5d0]
     kpath(9,:)=[0.5d0,0.5d0,0.5d0]
     points_name=[character(len=40) ::'M', 'R', 'G', 'X', 'M', 'G', 'Z','A', 'R']
  endif
  !Solve for the renormalized bands:
  if(master)call TB_Solve_model(TB_w90_model,Nlso,kpath,Nkpath,&
       colors_name=[black,red,green,blue,magenta,black,red,green,blue,magenta],&
       points_name=points_name,& 
       file="zBands_ssDFT",iproject=.true.)



  if(FSflag)then
     inquire(file='renorm.save',exist=bool)
     if(bool)then
        allocate(tmp(2*Nlat*Nspin*Norb))
        call read_array("zeta_self.restart",tmp)
        allocate(Zeta(Nlat*Nspin*Norb))
        allocate(Self(Nlat*Nspin*Norb))
        zeta = tmp(:Nlat*Nspin*Norb)
        self = tmp(Nlat*Nspin*Norb+1:)
        call TB_w90_Zeta(zeta)
        call TB_w90_Self(diag(self))
        deallocate(zeta,self)
     endif
     call TB_FSurface(Nlso,0d0,Nkvec(1:2),&
          colors_name=[black,red,red,green,blue],&
          file='FS_ssDFT',cutoff=1d-1,Niter=3,Nsize=2)
  endif



  call TB_w90_delete()


#ifdef _MPI
  call finalize_MPI()
#endif

end program ss_DFT
