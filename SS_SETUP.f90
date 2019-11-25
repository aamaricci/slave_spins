MODULE SS_SETUP
  USE SS_VARS_GLOBAL
  !
  USE SF_IOTOOLS
  implicit none
  private


  public :: ss_setup_structure
  public :: ss_init_params
  public :: ss_spin_symmetry

contains


  subroutine ss_setup_dimensions()
    Ns = 2*Norb
  end subroutine ss_setup_dimensions



  subroutine ss_setup_structure()
    !
    if(Nspin>2)stop "ss_setup error: Nspin>2"
    if(Norb>5)stop "ss_setup error: Norb>5. Fix me at SS_SETUP.f90/line.25"
    !
    call ss_setup_dimensions()
    !
    allocate(ss_lambda(Ns))
    allocate(ss_lambda0(Ns))
    allocate(ss_zeta(Ns))
    allocate(ss_weiss(Ns))
    allocate(ss_c(Ns))
    ss_lambda = 0d0
    ss_lambda0= 0d0
    ss_zeta   = 1d0
    ss_weiss  = 0d0
    ss_c      = 0d0
    !
    allocate(ss_dens(Ns))
    ss_dens  = 0d0
    !
    allocate(ss_Sz(Ns))
    ss_Sz = 0d0
    !
    allocate(ss_Hk(Ns,Ns,Nk))
    ss_Hk = zero
    !
    allocate(ss_Wtk(Ns,Ns,Nk))
    ss_Wtk= 0d0
    !
    allocate(ss_Hloc(Ns,Ns))
    ss_Hloc = zero
    !
    Nso = Nspin*Norb
    !    
  end subroutine ss_setup_structure



  subroutine ss_init_params()
    logical                          :: IOfile
    real(8),dimension(:),allocatable :: params
    inquire(file=trim(Pfile)//trim(ss_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       allocate(params(2*Ns))
       call read_array(trim(Pfile)//trim(ss_file_suffix)//".restart",params)
       ss_lambda = params(1:Ns)
       ss_zeta   = params(Ns+1:2*Ns)
    else
       ss_lambda = -ss_lambda0
       ss_zeta   = 1d0
    endif
  end subroutine ss_init_params



  subroutine ss_spin_symmetry(array)
    real(8),dimension(Ns)    :: array
    array(Norb+1:) = array(1:Norb)
  end subroutine ss_spin_symmetry


END MODULE SS_SETUP
