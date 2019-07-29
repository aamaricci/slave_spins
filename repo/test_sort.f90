program test_sort
  USE SCIFOR
  implicit none
  integer,parameter    :: L=5*10**7
  real(8),dimension(L) :: array,acopy
  integer,dimension(L) :: aindx
  real(8)              :: t_start,t_stop
  integer              :: seed,i

  seed = 5449489
  call mt_init(seed)
  call mt_random(array)
  acopy = array

  call cpu_time(t_start)
  call sort_array(acopy,aindx)
  call cpu_time(t_stop)
  print*,"Time:",t_stop-t_start
  do i=1,10
     print*,acopy(i),aindx(i)
  enddo

  call cpu_time(t_start)
  call sort(acopy,L)
  call cpu_time(t_stop)
  print*,"Time:",t_stop-t_start
  do i=1,10
     print*,acopy(i)
  enddo


  acopy = array
  call cpu_time(t_start)
  call quicksort_d(acopy,aindx,L)
  call cpu_time(t_stop)
  print*,"Time:",t_stop-t_start
  do i=1,10
     print*,acopy(i),aindx(i)
  enddo

  acopy = array
  call cpu_time(t_start)
  call indexx(L,acopy,aindx)
  call cpu_time(t_stop)
  print*,"Time:",t_stop-t_start
  do i=1,10
     print*,acopy(i),aindx(i),acopy(aindx(i))

  enddo
end program test_sort
