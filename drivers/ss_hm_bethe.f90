program ss_bethe
  USE SLAVE_SPINS
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  implicit none
  !
  integer                                     :: Le,Nlso,iorb,ilat,io
  real(8),dimension(5)                        :: Wbethe,Dbethe
  !  
  real(8),dimension(:,:),allocatable          :: Dbands
  real(8),dimension(:,:),allocatable          :: Ebands
  real(8),dimension(:),allocatable            :: H0
  real(8),dimension(:),allocatable            :: de,dens
  !
  real(8),dimension(:),allocatable            :: Wband
  !
  character(len=32)                           :: finput



  call parse_cmd_variable(finput,"FINPUT",default='inputSS_BETHE.conf')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wbethe,"WBETHE",finput,default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(Dbethe,"DBETHE",finput,default=[0d0,0d0,0d0,0d0,0d0])
  !
  call ss_read_input(trim(finput))

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Nlat,"nlat")
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  if(Nspin/=1.OR.Norb>5.OR.Nlat>2)stop "Wrong setup from input file: Nspin!=1 OR Norb>5 OR Nlat>2"
  Nlso=Nlat*Nspin*Norb


  allocate(Ebands(Nlso,Le))
  allocate(Dbands(Nlso,Le))
  allocate(H0(Nlso))
  allocate(Wband(Nlso))
  allocate(de(Nlso))
  !
  do ilat=1,Nlat
     do iorb=1,Norb
        io = iorb + (ilat-1)*Norb
        Wband(io) = Wbethe(iorb)
        H0(io)    = Dbethe(iorb)
     enddo
  enddo
  !
  do io=1,Nlso
     Ebands(io,:) = linspace(-Wband(io),Wband(io),Le,mesh=de(io))
     Dbands(io,:) = dens_bethe(Ebands(io,:),Wband(io))*de(io)
     call splot("dos_l"//str(io)//".ss",Ebands(io,:)+H0(io),Dbands(io,:)/de(io))
  enddo
  call TB_write_Hloc(one*diag(H0))
  !

  call ss_solve( Ebands,Dbands,Hloc=H0,ineq_sites=(/(1,ilat=1,Nlat)/) )


end program ss_bethe
