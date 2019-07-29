MODULE SS_SETUP
  USE SS_INPUT_VARS
  !
  USE SF_IOTOOLS
  implicit none
  private


  public :: ss_setup_structure
  public :: ss_init_params

contains


  subroutine ss_setup_dimensions()
    Ns = 2*Norb
  end subroutine ss_setup_dimensions



  subroutine ss_setup_structure()
    !
    call ss_setup_dimensions()
    !
    allocate(ss_lambda(Ns))
    allocate(ss_zeta(Ns))
    allocate(ss_weiss(Ns))
    allocate(ss_c(Ns))
    ss_lambda= 0d0
    ss_zeta  = 1d0
    ss_weiss = 0d0
    ss_c     = 0d0
    !
    allocate(ss_dens(Ns))
    ss_dens  = 0d0
    !
    allocate(ss_Sz(Ns))
    ss_Sz = 0d0
    !
    allocate(ss_Hk(Nspin*Norb,Nspin*Norb,Nk))
    ss_Hk = zero
    !
    allocate(ss_Wtk(Nk))
    ss_Wtk= 0d0
    !
    allocate(ss_Hloc(Nspin*Norb,Nspin*Norb))
    ss_Hloc = zero
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
    endif

  end subroutine ss_init_params


END MODULE SS_SETUP
