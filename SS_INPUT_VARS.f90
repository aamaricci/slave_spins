MODULE SS_INPUT_VARS
  USE SS_VARS_GLOBAL
  USE SF_VERSION
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str
  implicit none

  !GIT VERSION
  include "revision.inc"  !this file is generated at compilation time in the Makefile


  !input variables
  !=========================================================
  integer              :: Norb                !Norb =# of impurity orbitals
  integer              :: Nspin               !Nspin=# spin degeneracy (max 2)
  real(8),dimension(5) :: Uloc                !local interactions
  real(8)              :: Ust                 !intra-orbitals interactions
  real(8)              :: Jh                  !J_Hund: Hunds' coupling constant 
  real(8)              :: Jx                  !J_X: coupling constant for the spin-eXchange interaction term
  real(8)              :: Jp                  !J_P: coupling constant for the Pair-hopping interaction term 
  real(8)              :: xmu                 !chemical potential
  real(8)              :: beta                !inverse temperature
  real(8)              :: eps                 !broadening
  real(8)              :: wini,wfin           !
  integer              :: nloop               !max convergence loop variables
  integer              :: Nsuccess            !
  integer              :: verbose          !


  !Some parameters for function dimension:
  !=========================================================
  integer              :: Lmats
  integer              :: Lreal

  !LOG AND Hamiltonian UNITS
  !=========================================================
  character(len=100)   :: Pfile
  integer,save         :: LOGfile



contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine ss_read_input(INPUTunit,comm)
#ifdef _MPI
    USE MPI
    USE SF_MPI
#endif
    character(len=*) :: INPUTunit
    integer,optional :: comm
    logical          :: master=.true.
    integer          :: i,rank=0
#ifdef _MPI
    if(present(comm))then
       master=get_Master_MPI(comm)
       rank  =get_Rank_MPI(comm)
    endif
#endif
    !
    !Store the name of the input file:
    ss_input_file=str(INPUTunit)
    !
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of impurity orbitals (max 5).")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy (max 2)")
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=[2d0,0d0,0d0,0d0,0d0],comment="Values of the local interaction per orbital (max 5)")
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="S-E coupling")
    call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="P-H coupling")
    call parse_input_variable(beta,"BETA",INPUTunit,default=1000.d0,comment="Inverse temperature, at T=0 is used as a IR cut-off.")
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")
    call parse_input_variable(nloop,"NLOOP",INPUTunit,default=100,comment="Max number of DMFT iterations.")
    call parse_input_variable(nsuccess,"NSUCCESS",INPUTunit,default=1,comment="Number of successive iterations below threshold for convergence")
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=5000,comment="Number of Matsubara frequencies.")
    call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of real-axis frequencies.")
    call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
    call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")
    call parse_input_variable(Pfile,"PFILE",INPUTunit,default="parameters",comment="File where to retrieve/store the parameters.")
    call parse_input_variable(verbose,"VERBOSE",INPUTunit,default=3,comment="Verbosity level: 0=almost nothing --> 5:all. Really: all")

#ifdef _MPI
    if(present(comm))then
       if(.not.master)then
          LOGfile=1000-rank
          open(LOGfile,file="stdOUT.rank"//str(rank)//".ed")
          do i=1,get_Size_MPI(comm)
             if(i==rank)write(*,"(A,I0,A,I0)")"Rank ",rank," writing to unit: ",LOGfile
          enddo
       endif
    endif
#endif
    !
    !
    if(master)then
       call print_input()
       call save_input(INPUTunit)
       call scifor_version()
       call code_version(revision)
    endif
    !Act on the input variable only after printing.
    !In the new parser variables are hard-linked into the list:
    !any change to the variable is immediately copied into the list... (if you delete .ed it won't be printed out)
    call substring_delete(Pfile,".restart")
    call substring_delete(Pfile,".ss")
  end subroutine ss_read_input




  subroutine substring_delete (s,sub)
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop


END MODULE SS_INPUT_VARS
