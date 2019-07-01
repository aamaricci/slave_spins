MODULE SS_SETUP
  USE SS_INPUT_VARS
  implicit none
  private

  public :: ss_get_dimension


contains


  function ss_get_dimension() result(Nss)
    integer :: Nss
    Nss = 2*Norb
  end function ss_get_dimension
  

END MODULE SS_SETUP
