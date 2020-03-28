program ss_bethe
  USE SLAVE_SPINS
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: Hk
  real(8)                                 :: ts(5)
  real(8)                                 :: Mh(5)
  integer                                 :: Nkx,Nktot,Nso


  call parse_input_variable(ts,"ts","inputSS.conf",default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(Mh,"Mh","inputSS.conf",default=[0d0,0d0,0d0,0d0,0d0])
  call parse_input_variable(Nkx,"Nkx","inputSS.conf",default=20)
  call ss_read_input('inputSS.conf')

  Nktot = Nkx**2

  Nso = Nspin*Norb
  if(Nso>5)stop "Nso > 5"

  allocate(Hk(Nso,Nso,Nktot))

  call TB_set_bk([pi2,0d0],[0d0,pi2])
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx])

  !init SS
  call ss_init(Hk)

  call ss_solve()


contains

  function hk_model(kvec,N) result(hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N,i
    complex(8),dimension(N,N) :: Hk
    real(8)                   :: kx,ky,cx,cy
    kx = kvec(1)
    ky = kvec(2)
    cx = cos(kx)
    cy = cos(ky)
    hk = zero
    do i=1,N
       hk(i,i) = Mh(i)-2d0*ts(i)*(cx+cy)
    enddo
  end function hk_model
end program ss_bethe
