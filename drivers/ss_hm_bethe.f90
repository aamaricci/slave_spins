program ss_bethe
  USE SLAVE_SPINS
  implicit none

  integer :: Nlambda
  real(8),dimension(:),allocatable :: lambda
  

  call ss_read_input('inputSS.conf')


  Nlambda = ss_get_dimension()
  



end program ss_bethe
