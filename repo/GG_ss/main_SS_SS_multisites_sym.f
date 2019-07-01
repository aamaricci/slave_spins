      PROGRAM SlaveSpins_multisites
      implicit none
      include 'main_SS_SS_multisites_sym.h'  
c input:       
      open(2,file='INPUT.in')    ! INPUT file
       read(2,*)Uvalue           ! U value
       read(2,*)Jh               ! Jh value
       read(2,*)whattypeofJh     ! S-F and P-H included (1) or not (0)
       read(2,*)filling          ! filling 
       read(2,*)Threshold        ! Threshold for the self-consistent loops
       read(2,*)whatdoyourun     ! Type of Run
       read(2,*)ifdoublecounting ! if doublecounting is taken into account
       read(2,*)whattypeofdiago  ! Lapack (0), Lanczos (1)
       if(whattypeofdiago.eq.1)then
       read(2,*)eThresholdLanczos! Threshold for the Lanczos
       endif
       read(2,*)if_sym
       read(2,*)if_sym_updn
      close(2)

      if_Ldos=0

      write(*,*)' '
      write(*,*)' U,Jh value = ',Uvalue,Jh
      write(*,*)' '
      write(*,*)' S-F and P-H in (1) or not (0) '
     &          ,whattypeofJh
      write(*,*)' ' 
      write(*,*)' filling = ',filling
      write(*,*)' '
      write(*,*)' Threshold = ',Threshold
      write(*,*)' '
      write(*,*)' Type of run : ',whatdoyourun
      write(*,*)' '
      write(*,*)' Double Counting = ',ifdoublecounting
      write(*,*)' '
      write(*,*)' Diagonalization Lapack (0), Lanczos (1) = '
     &          ,whattypeofdiago
      write(*,*)' '
      write(*,*)' Threshold Lanczos = ',eThresholdLanczos
      write(*,*)' '
      write(*,*)' Symmetry = ',if_sym
      write(*,*)' '
      write(*,*)' Symmetry UP/DN = ',if_sym_updn

      if(if_sym.eq.0)then
       if(num_wann_ss_ind.ne.num_wann_ss)then
        write(*,*)' check for consistency ind symmetry '
        stop
       endif
      endif

c DFT input:
c read num_kpts k-points:

      write(*,*)' '
      write(*,*)' k - mesh in *'
      open(2,file='kpoints.scf')
      do loop_kpt=1,num_kpts
       read(2,*)(kpt_latt(j,loop_kpt),j=1,ndim) 
      enddo
      close(2)
      write(*,*)' k - mesh in **'

c read DFT hoppings between unit cells and orbitals 
      write(*,*)' '
      write(*,*)' h lda n,m,R *'

      do loop_rpt=1,nrpts
       do i=1,num_wann
        do j=1,num_wann
        ham_r2(i,j,loop_rpt)=dcmplx(0.d0,0.d0)
        enddo
       enddo
      enddo

      open(2,file='hr.dat')
      read(2,*)factor_hr  
      write(*,*)' Reduced or Un-Reduced unit cell ',factor_hr
      do loop_rpt=1,nrpts
       do i=1,num_wann
        do j=1,num_wann
         read(2,*)irvec(1,loop_rpt)         
     &           ,irvec(2,loop_rpt)
     &           ,irvec(3,loop_rpt)
     &           ,jj,ii,a,b
         ham_r2(i,j,loop_rpt)=dcmplx(a,b)   
        enddo
       enddo
      enddo
      close(2)
      write(*,*)' h lda n,m,R **'
      write(*,*)

c symmetries
      if(if_sym.eq.0)then

       Npara=num_wann_ss_ind

       do iorb=1,num_wann_ss
        index_sym(iorb)=iorb
       enddo

       write(*,*)' symmetries OFF '
       write(*,*)' # independent lambda = ',Npara
       write(*,*)
       write(*,*)' description symmetry : '
       write(*,*)' orb index sym orb '
       do iorb=1,num_wann_ss
        write(*,*)iorb,index_sym(iorb)
       enddo

      endif

      if(if_sym.eq.1)then

       Npara=num_wann_ss_ind

       open(2,file='sym.in')
       do iorb=1,num_wann_ss
        read(2,*)index_sym(iorb)
       enddo
       close(2)

       write(*,*)' symmetries ON '
       write(*,*)' # independent lambda = ',Npara
       write(*,*)
       write(*,*)' description symmetry : '
       write(*,*)' orb index sym orb '
       do iorb=1,num_wann_ss
        write(*,*)iorb,index_sym(iorb)
       enddo

      endif

c initial guess:
c read initial Lagrange multipliers lambda
c read initial <O>

      do iorb=1,num_wann_ss
       lambda(iorb)=0.0d0 
       mOi_orb(iorb)=cmplx(0.d0,0.d0)
      enddo

      write(*,*)''
      write(*,*)' l and <O> parameters '
      open(2,file='parameters')
      do iorb=1,num_wann_ss
       read(2,*)lambda(iorb)
      enddo
      do iorb=1,num_wann_ss
       read(2,*)mOi_orb(iorb)
      enddo
      do iorb=1,num_wann_ss
       read(2,*)shifteps(iorb)
      enddo
      if(ifdoublecounting.eq.1)then
      read(2,*)navhartree
      endif
      close(2)
c write all the initial parameters:
      write(*,*)''
      write(*,*)' initial l parameters '
      do iorb=1,num_wann_ss
       write(*,*)lambda(iorb)
      enddo
      write(*,*)''
      write(*,*)' initial <O> parameters '
      do iorb=1,num_wann_ss
       write(*,*)mOi_orb(iorb)
      enddo
      write(*,*)''
      write(*,*)' shift on-site energies '
      do iorb=1,num_wann_ss
       write(*,*)shifteps(iorb)
      enddo
c 
      do j=1,num_wann_ss
       mOi_orb(j+num_wann_ss)=mOi_orb(j)
      enddo
      do j=1,nspin*num_wann_ss
       mOpi_orb(j)=mOi_orb(j)
      enddo
c 
      do j=1,num_wann_ss
       mOi_orb_ini(j)=mOi_orb(j)
      enddo
      do j=1,num_wann_ss
       mOi_orb_ini(j+num_wann_ss)=mOi_orb(j)
      enddo

      do j=1,num_wann_ss
       mOpi_orb_ini(j)=mOpi_orb(j)
      enddo
      do j=1,num_wann_ss
       mOpi_orb_ini(j+num_wann_ss)=mOpi_orb(j)
      enddo

c broyden subroutine (based on Broyden`s method) changes the lambda_i to get the constrain n=S+1/2 satisfied.
       write(*,*)' parameters passed to the Broyden '
      do iorb=1,num_wann_ss_ind
       x(iorb)=lambda(iorb)
       write(*,*)iorb,x(iorb)
      enddo
       write(*,*)

      if(if_sym.eq.0)then
       write(*,*)' '
       write(*,*)' broyden will work with # lambda_i ',Npara
       write(*,*)' *** no symmetries applied **** '
       write(*,*)' '
      endif
      if(if_sym.eq.1)then
       write(*,*)' '
       write(*,*)' broyden will work with # lambda_i ',Npara
       write(*,*)' *** symmetries applied **** '
       write(*,*)' '
      endif
c 
c      Npara=num_wann_ss

      ciclo=0

      write(*,*)''
      write(*,*)' call broydn '
      write(*,*)''
c 
c this subroutine is taken from Numerical Recipes.
c 
      call broydn(x,Npara,check)
      write(*,*)''
      write(*,*)' out  broydn '
      write(*,*)''
      write(*,*)' spin-slaves equations solved'
      write(*,*)''

      if (check) then
        write(*,*)' Convergence problems !!!!'
      endif

c output 
c write lambda and <O>
c      do iorb=1,num_wann_ss
c       lambda(iorb)=x(iorb)
c      enddo

       write(*,*)' parameters passed from the Broyden '
      do iorb=1,num_wann_ss
       lambda(iorb)=x(index_sym(iorb))
      enddo
       write(*,*)

      write(*,*)''
      write(*,*)' final l parameters '
      do iorb=1,num_wann_ss
       write(*,*)lambda(iorb)
      enddo
      write(*,*)''

      open(2,file='parameters')
      do iorb=1,num_wann_ss
       write(2,*)lambda(iorb)
      enddo
      do iorb=1,num_wann_ss
       write(2,*)mOi_orb(iorb)
      enddo
      do iorb=1,num_wann_ss
       write(2,*)shifteps(iorb)
      enddo
      if(ifdoublecounting.eq.1)then
      write(2,*)navhartree
      endif
      close(2)
c output code:
      open(2,file='Z_n_U_Jh')
       if(num_wann_ss.eq.1)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &        ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.eq.2)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.eq.3)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.eq.4)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.eq.5)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.eq.6)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.eq.7)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.eq.8)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.eq.9)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.eq.10)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.ge.11)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.ge.12)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.ge.13)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.ge.14)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.ge.15)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.ge.16)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.ge.17)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.ge.18)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif  
       if(num_wann_ss.ge.19)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.ge.20)then
       write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5
     &         ,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &         Uvalue,Jh,(real(mOi_orb(iorb))**2,iorb=1,num_wann_ss)
     &         ,(n_state(iorb),iorb=1,num_wann_ss)
       endif
       if(num_wann_ss.ge.21)then
        write(2,*)' NOT IMPLEMENTED !!! '
        STOP
       endif
      close(2)
c 
      if(if_sym.eq.0)then
       Nsite=int(real(num_wann_ss/N_site_tot))
      endif
      if(if_sym.eq.1)then
       Nsite=int(real(num_wann_ss_ind/N_site_tot))
      endif

      indx_orb_site=0
      do isite=1,Nsite
       do iorb=1,N_site_tot
         write(x1,'(I5.5)')iorb
         write(x2,'(I5.5)')isite
         filename='SzSz'//x1//x2//''
         open(700,file=filename,form='formatted',status='unknown')
         write(700,"(100F18.12)")Uvalue,Jh
     &        ,(real(mSz_orb(jorb+indx_orb_site)),jorb=1,N_site_tot)
     &        ,((real(mSzSz_iorb_jorb(1,iorb+indx_orb_site
     &               ,jorb+indx_orb_site))),jorb=1,N_site_tot)
     &        ,((real(mSzSz_iorb_jorb(3,iorb+indx_orb_site
     &               ,jorb+indx_orb_site))),jorb=1,N_site_tot)
         close(700) 
       enddo
      indx_orb_site=indx_orb_site+N_site_tot
      enddo
c  
      if_Ldos=1
      call fermionic_hamiltonian(x)
c


      stop
      END

      SUBROUTINE funcv(Npara,x,f)
      implicit none
      include 'main_SS_SS_multisites_sym.h'
      integer iter,jter
      integer infoconv1
      integer Nmaxiter
c     integer Nsite,isite,indx_orb_site
      parameter (Nmaxiter=25)
      complex*16 f_aux1(Nmaxiter,nspin*num_wann_ss)  
      complex*16 mOi_orb_aux(nspin*num_wann_ss)
      complex*16 mOpi_orb_aux(nspin*num_wann_ss)
      complex*16 mSz_orb_aux(nspin*num_wann_ss)

      ciclo=ciclo+1
      write(*,*)'***************'
      write(*,*)'**** Cycle ****',ciclo 
      write(*,*)'***************'

      do i=1,nspin*num_wann_ss
         do iter=1,Nmaxiter
            f_aux1(iter,i)=0.d0  
         enddo
      enddo
c     same starting guess at all the cicles
c     up
      do iorb=1,num_wann_ss
         mOi_orb(iorb)=mOi_orb_ini(iorb)
         mOpi_orb(iorb)=mOpi_orb_ini(iorb)
      enddo
c     dn
      do iorb=1,num_wann_ss
         mOi_orb(iorb+num_wann_ss)=mOi_orb_ini(iorb+num_wann_ss)
         mOpi_orb(iorb+num_wann_ss)=mOpi_orb_ini(iorb+num_wann_ss)
      enddo
c     
      do i=1,num_wann
         n_state(i)=0.d0      
      enddo
c     
      do i=1,nspin*num_wann_ss
         mSz_orb(i)=0.d0          
      enddo
c     
      write(*,*)''
      write(*,*)' Self-Consistency between Hs and Hf starts '
      write(*,*)''

      do iter=1,Nmaxiter 

         write(*,*)'============== '
         write(*,*)' Iteration # = ',iter
         write(*,*)'============== '
         write(*,*)' '
         write(*,*)' l in funcv (i) '
c     
         do i=1,num_wann_ss
            write(*,*)x(index_sym(i))
         enddo
c     
         write(*,*)' '
         write(*,*)' call H F'
         write(*,*)' '
         call fermionic_hamiltonian(x) 
c     
         write(*,*)' '
         write(*,*)' call H S'
         write(*,*)' '

         if(if_sym.eq.0)then

            Nsite=int(real(num_wann_ss/N_site_tot))

            write(*,*)' symmetry OFF'
            write(*,*)' Nsites = ',Nsite
            write(*,*)

         endif

         if(if_sym.eq.1)then

            write(*,*)' symmetry ON'
            Nsite=int(real(num_wann_ss_ind/N_site_tot))

            write(*,*)' Nsites = ',Nsite
            write(*,*)

         endif

         do i=1,2*nspin
            do iorb=1,num_wann
               do jorb=1,num_wann
                  mSzSz_iorb_jorb(i,iorb,jorb)=0.d0
               enddo
            enddo
         enddo
c     up             
c     do iorb=1,num_wann_ss
c     mOi_orb_aux(iorb)=mOi_orb(iorb)
c     mOpi_orb_aux(iorb)=mOpi_orb(iorb)
c     mSz_orb_aux(iorb)=mSz_orb(iorb)
c     enddo
c     c dn
c     do iorb=1,num_wann_ss
c     mOi_orb_aux(iorb+num_wann_ss)=mOi_orb(iorb+num_wann_ss)
c     mOpi_orb_aux(iorb+num_wann_ss)=mOpi_orb(iorb+num_wann_ss)
c     mSz_orb_aux(iorb+num_wann_ss)=mSz_orb(iorb+num_wann_ss)
c     enddo
c     

         indx_orb_site=0

         do isite=1,Nsite

            write(*,*)' ============='
            write(*,*)' site = ',isite
            write(*,*)' ============='

            call spin_hamiltonian(x,isite,indx_orb_site)

c     up
            do iorb=1,N_site_tot
               mOi_orb_aux(iorb+indx_orb_site)=mOi(iorb)
               mOpi_orb_aux(iorb+indx_orb_site)=mOpi(iorb)
               mSz_orb_aux(iorb+indx_orb_site)=mSz(iorb)
            enddo
            write(*,*)
            write(*,*)' <O^+> up (update) ' 
            do iorb=1,num_wann_ss
               write(*,*)mOi_orb_aux(iorb)
            enddo
            write(*,*)
            write(*,*)' <Sz> up (update) '
            do iorb=1,num_wann_ss
               write(*,*)mSz_orb_aux(iorb)
            enddo
            write(*,*)
c     dn
            do iorb=1,N_site_tot
               mOi_orb_aux(iorb+indx_orb_site+num_wann_ss)=mOi(iorb+N_site_tot)
               mOpi_orb_aux(iorb+indx_orb_site+num_wann_ss)=
     &              mOpi(iorb+N_site_tot)
               mSz_orb_aux(iorb+indx_orb_site+num_wann_ss)=mSz(iorb+N_site_tot)
            enddo    
            write(*,*)
            write(*,*)' <O^+> dn (update) '
            do iorb=1,num_wann_ss
               write(*,*)mOi_orb_aux(iorb+num_wann_ss)
            enddo
            write(*,*)
            write(*,*)' <Sz> dn (update) '
            do iorb=1,num_wann_ss
               write(*,*)mSz_orb_aux(iorb+num_wann_ss)
            enddo
            write(*,*)                               
c     
            do i=1,2*nspin
               do iorb=1,N_site_tot
                  do jorb=1,N_site_tot
                     mSzSz_iorb_jorb(i,iorb+indx_orb_site,jorb+indx_orb_site)=
     &                    SzSz(i,iorb,jorb)
                  enddo
               enddo
            enddo
c     
            indx_orb_site=indx_orb_site+N_site_tot

         enddo              
c     
c     
c     do iorb=1,num_wann_ss
c     mOi_orb(iorb)=mOi_orb_aux(iorb)
c     mOpi_orb(iorb)=mOpi_orb_aux(iorb)
c     mSz_orb(iorb)=mSz_orb_aux(iorb)
c     enddo
c     
c     do iorb=1,num_wann_ss
c     mOi_orb(iorb+num_wann_ss)=mOi_orb_aux(iorb+num_wann_ss)
c     mOpi_orb(iorb+num_wann_ss)=mOpi_orb_aux(iorb+num_wann_ss)
c     mSz_orb(iorb+num_wann_ss)=mSz_orb_aux(iorb+num_wann_ss)
c     enddo
c     
c     
         do iorb=1,num_wann_ss
            mOi_orb(iorb)=mOi_orb_aux(index_sym(iorb))
            mOpi_orb(iorb)=mOpi_orb_aux(index_sym(iorb))
            mSz_orb(iorb)=mSz_orb_aux(index_sym(iorb))
         enddo
c     
         do iorb=1,num_wann_ss
            mOi_orb(iorb+num_wann_ss)=mOi_orb_aux(index_sym(iorb)
     &           +num_wann_ss)
            mOpi_orb(iorb+num_wann_ss)=mOpi_orb_aux(index_sym(iorb)
     &           +num_wann_ss)
            mSz_orb(iorb+num_wann_ss)=mSz_orb_aux(index_sym(iorb)
     &           +num_wann_ss)
         enddo
c     
         write(*,*)''
         write(*,*)' solution of delta-functions = 0 (NM, SP cases) '
         write(*,*)' <n>=Sz+1/2 '
         write(*,*)''
c     up
         write(*,*)' spin = 1 '
         do i=1,num_wann_ss
            write(*,*)n_state(i),(mSz_orb(i)+0.5)
            write(*,*)n_state(i)-(mSz_orb(i)+0.5)
            write(*,*)
         enddo
         write(*,*)''
c     dn
         write(*,*)' spin = 2 '
         do i=1,num_wann_ss
            write(*,*)n_state(i),(mSz_orb(i+num_wann_ss)+0.5)
            write(*,*)n_state(i)-(mSz_orb(i+num_wann_ss)+0.5)
            write(*,*)
         enddo
         write(*,*)''

         indx=0
         do j=1,nspin*num_wann_ss
            indx=indx+1 
            f_aux1(iter,indx)=mOi_orb(j)             
         enddo

         write(*,*)
         write(*,*)' function for the convergence : '
         do j=1,nspin*num_wann_ss
            write(*,*)f_aux1(iter,j)
         enddo
         write(*,*)

         write(*,*)'<O^+_i><O_i>'
         do j=1,nspin*num_wann_ss
            write(*,*)mOi_orb(j)*mOpi_orb(j)
         enddo
         write(*,*)

         write(*,*)' <O_i> '
         write(*,*)''
         do i=1,nspin*num_wann_ss
            write(*,*)mOi_orb(i)
         enddo
         write(*,*)''

         write(*,*)' <O^+_i> '
         write(*,*)''
         do i=1,nspin*num_wann_ss
            write(*,*)mOpi_orb(i)
         enddo
         write(*,*)''

         if(iter.gt.1)then
            infoconv1=0
            
            do  i=1,nspin*num_wann_ss               
               write(*,*)' abs ',i,abs(f_aux1(iter,i)-f_aux1(iter-1,i))
     &              ,Threshold
               if(abs(f_aux1(iter,i)-f_aux1(iter-1,i))
     &              .lt.Threshold)then
                  
                  infoconv1=infoconv1+1
                  
                  write(*,*)' conv check ',i
     &                 ,abs(f_aux1(iter,i)-f_aux1(iter-1,i))
                  
               endif
               
            enddo
            write(*,*)
            if(infoconv1.eq.nspin*num_wann_ss)then
               goto 1
            endif
         endif
         
      enddo 

 1    write(*,*)''
      write(*,*)' Self-Consistency between Hs and Hf ends '
      write(*,*)''

      write(*,*)''
      write(*,*)' Self-Consistent Equations '
      write(*,*)''
c     
      do i=1,num_wann_ss_ind
         f(i)=n_state(i)-(mSz_orb(i)+0.5)     
         write(*,*)'functions for the broyden ',f(i)
      enddo
      write(*,*)

      do i=1,num_wann_ss
         write(*,*)' occupancies in funcv ',n_state(i)
         write(*,*)' <Sz>+0.5 ',mSz_orb(i)+0.5
         write(*,*)' deviation ',n_state(i)-(mSz_orb(i)+0.5)
      enddo
      write(*,*)
c     
      return
      END

      SUBROUTINE fermionic_hamiltonian(x)
      implicit none
      include 'main_SS_SS_multisites_sym.h'
c      integer Nsite,isite
c 
c      do j=1,num_wann_ss
c       lambdawrk(j)=x(j)
c      enddo

       write(*,*)' '
      if(if_sym.eq.0)then
       write(*,*)' *** symmetry OFF *** '
       write(*,*)' orb ind orb '
       do iorb=1,num_wann_ss
        write(*,*)iorb,index_sym(iorb)
       enddo
       write(*,*)
      endif
      if(if_sym.eq.1)then
       write(*,*)' *** symmetry ON *** '
       write(*,*)' orb ind orb '
       do iorb=1,num_wann_ss
        write(*,*)iorb,index_sym(iorb)
       enddo
       write(*,*)
      endif

      do j=1,num_wann_ss
       lambdawrk(j)=x(index_sym(j))
      enddo
c 
      do loop_rpt=1,nrpts
       do i=1,num_wann
        do j=1,num_wann
        ham_r(i,j,loop_rpt)=dcmplx(0.d0,0.d0)
        enddo
       enddo
      enddo
c 
      do loop_rpt=1,nrpts
       do i=1,num_wann
        do j=1,num_wann
         ham_r(i,j,loop_rpt)=ham_r2(i,j,loop_rpt)
        enddo
       enddo
      enddo
c 
      do loop_rpt=1,nrpts
       if(irvec(1,loop_rpt).eq.0.and.
     &    irvec(2,loop_rpt).eq.0.and.
     &    irvec(3,loop_rpt).eq.0)then
       do i=1,num_wann

          e_on_site(i)=real(ham_r(i,i,loop_rpt))

          ham_r(i,i,loop_rpt)=dcmplx(0.d0,0.d0) 

       enddo
       endif
      enddo
c 
      do loop_rpt=1,nrpts

       do i=1,num_wann
        do j=1,num_wann

         if(i.le.num_wann_ss.and.j.le.num_wann_ss)then

         ham_r(i,j,loop_rpt)=
     &   ham_r(i,j,loop_rpt)*
     &   mOi_orb(i)*mOpi_orb(j)

          if(irvec(1,loop_rpt).eq.0.and.
     &       irvec(2,loop_rpt).eq.0.and.
     &       irvec(3,loop_rpt).eq.0)then
          write(*,*)' Z iorb, Z jorb ',i,j
          write(*,*)mOi_orb(i)*mOpi_orb(j)
          endif

         elseif(i.le.num_wann_ss.and.j.gt.num_wann_ss)then

         ham_r(i,j,loop_rpt)=
     &   ham_r(i,j,loop_rpt)*
     &   mOi_orb(i)

          if(irvec(1,loop_rpt).eq.0.and.
     &       irvec(2,loop_rpt).eq.0.and.
     &       irvec(3,loop_rpt).eq.0)then
          write(*,*)' Z iorb ',i,j
          write(*,*)mOi_orb(i)
          endif

         elseif(i.gt.num_wann_ss.and.j.le.num_wann_ss)then

         ham_r(i,j,loop_rpt)=
     &   ham_r(i,j,loop_rpt)*
     &   mOi_orb(j)

          if(irvec(1,loop_rpt).eq.0.and.
     &       irvec(2,loop_rpt).eq.0.and.
     &       irvec(3,loop_rpt).eq.0)then
          write(*,*)' Z jorb ',i,j
          write(*,*)mOi_orb(j)
          endif

         elseif(i.gt.num_wann_ss.and.j.gt.num_wann_ss)then

         ham_r(i,j,loop_rpt)=
     &   ham_r(i,j,loop_rpt)

          if(irvec(1,loop_rpt).eq.0.and.
     &       irvec(2,loop_rpt).eq.0.and.
     &       irvec(3,loop_rpt).eq.0)then
          write(*,*)' Z=0 ',i,j
          write(*,*)
          endif

         endif

        enddo
       enddo

      enddo
c 
      write(*,*)' '
      write(*,*)' Type of run : ',whatdoyourun
      write(*,*)' '

      do loop_rpt=1,nrpts
       if(irvec(1,loop_rpt).eq.0.and.
     &    irvec(2,loop_rpt).eq.0.and.
     &    irvec(3,loop_rpt).eq.0)then
        do i=1,num_wann_ss

          if(whatdoyourun.eq.0)then
          ham_r(i,i,loop_rpt)=
     &    ham_r(i,i,loop_rpt)
     &    +dcmplx(e_on_site(i),0.d0)
          endif
          if(whatdoyourun.eq.1)then
          ham_r(i,i,loop_rpt)=
     &    ham_r(i,i,loop_rpt)
     &    -dcmplx(lambdawrk(i),0.d0)
     &    +dcmplx(e_on_site(i),0.d0)
     &    +dcmplx(shifteps(i),0.d0)
          endif

        enddo

        do i=num_wann_ss+1,num_wann

          ham_r(i,i,loop_rpt)=
     &    ham_r(i,i,loop_rpt)
     &    +dcmplx(e_on_site(i),0.d0)

        enddo

       endif
      enddo
c
c      Edc=0.d0
      if(ifdoublecounting.eq.1)then

      write(*,*)
      write(*,*)' double counting '
      write(*,*)

      Uprime=Uvalue-2*Jh

      navhartree=navhartree-0.5d0

      do loop_rpt=1,nrpts
       if(irvec(1,loop_rpt).eq.0.and.
     &    irvec(2,loop_rpt).eq.0.and.
     &    irvec(3,loop_rpt).eq.0)then
        do i=1,num_wann_ss
          ham_r(i,i,loop_rpt)=
     &    ham_r(i,i,loop_rpt)
     &    -dcmplx(Uvalue*navhartree,0.d0)
     &    -dcmplx(Uprime
     &            *navhartree*(N_site_tot-1),0.d0)
     &    -dcmplx((Uprime-Jh)
     &            *navhartree*(N_site_tot-1),0.d0)

          write(*,*)
          write(*,*)i,
     &    -dcmplx(Uvalue*navhartree,0.d0)
     &    -dcmplx(Uprime
     &            *navhartree*(N_site_tot-1),0.d0)
     &    -dcmplx((Uprime-Jh)
     &            *navhartree*(N_site_tot-1),0.d0)
          write(*,*)

          write(*,*)'U'
          write(*,*)i,
     &    -dcmplx(Uvalue*navhartree,0.d0)
          write(*,*)

          write(*,*)'U`'
          write(*,*)i,
     &    -dcmplx(Uprime
     &            *navhartree*(N_site_tot-1),0.d0)
          write(*,*)

          write(*,*)'U`-Jh'
          write(*,*)i,
     &    -dcmplx((Uprime-Jh)
     &            *navhartree*(N_site_tot-1),0.d0)
          write(*,*)

        enddo
       endif
      enddo

c      Nsite=int(real(num_wann_ss/N_site_tot))
c
c      do isite=1,Nsite
c
c       Edc=Edc+N_site_tot*Uvalue*navhartree**2
c       k=0
c       do i=1,N_site_tot
c        do j=i+1,N_site_tot
c        k=k+1
c        enddo
c       enddo
c       write(*,*)
c       write(*,*)k
c       write(*,*)
c       Edc=Edc+nspin*Uprime*k*navhartree**2
c       Edc=Edc+nspin*(Uprime-Jh)*k*navhartree**2
c
c       Edc=Edc+N_site_tot*Uvalue*navhartree
c       Edc=Edc+nspin*Uprime*k*navhartree
c       Edc=Edc+nspin*(Uprime-Jh)*k*navhartree
c
c      enddo


      endif

      write(*,*)
      write(*,*)' Diagonalization Fermionic problem '
      call idate(today)
      call itime(now)
      write ( *, 1007 )  today(2), today(1), today(3), now
 1007 format ( '  Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
         write(*,*)

c Diagonalization Fermionic problem
      do loop_kpt=1,num_kpts 

        do i=1,num_wann
         do j=1,num_wann
          ham_k(i,j)=(0.d0,0.d0)
         enddo
        enddo
c 
        do i=1,num_wann
         do j=1,num_wann
          do loop_rpt=1,nrpts
           rdotk=0.d0
            do ll=1,ndim
             rdotk=rdotk+factor_hr*(kpt_latt(ll,loop_kpt)
     &                        *(irvec(ll,loop_rpt)))
            enddo
              ham_k(i,j)=ham_k(i,j)+ham_r(i,j,loop_rpt)
     &               *dcmplx(cos(rdotk),-sin(rdotk))
          enddo
         enddo
        enddo
c
        call zheevx('V'
     &             ,'A'
     &             ,'U'
     &             ,num_wann
     &             ,ham_k
     &             ,num_wann
     &             ,vl,vu
     &             ,1
     &             ,num_wann
     &             ,0.0
     &             ,neigval,eigval
     &             ,Z,num_wann
     &             ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

        if(info.ne.0)then
         write(*,*)' Failure in ZHEEVX. INFO =',info 
         stop
        endif
c 
        do l=1,num_wann
         eigenvalues(num_wann*(loop_kpt-1)+l)=eigval(l)
         eigenvalues_k(loop_kpt,l)=eigval(l)
        enddo
c
        do ll=1,num_wann
        do ii=1,num_wann
         Z_k(loop_kpt,ii,ll)=Z(ii,ll)
        enddo
        enddo
c
      enddo 

      write(*,*)
      write(*,*)' End Diagonalization Fermionic problem '
      call idate(today)
      call itime(now)
      write ( *, 1008 )  today(2), today(1), today(3), now
 1008 format ( '  Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
      write(*,*)
c 
      write(*,*)' '
      Ntot=num_wann*num_kpts
      write(*,*)' # states ',Ntot
      N_e_k=ifix(filling*real(Ntot)+0.5)
      write(*,*)' # states filled ',N_e_k
c
      call indexx(Ntot,eigenvalues,e_ind)
c 
      EF=eigenvalues(e_ind(N_e_k))
c
      write(*,*)' '
      write(*,*)' EF = ',EF 
c 
      eneW=0.d0
      do i=1,N_e_k 
       eneW=eneW+eigenvalues(e_ind(i))
      enddo
c
      write(*,*)''
      write(*,*)' ene tot ',eneW/real(num_kpts)
c 
      do ll=1,num_wann
       n_state(ll)=0.d0
      enddo
      do ii=1,num_wann
       do loop_kpt=1,num_kpts
        do ll=1,num_wann
        if(eigenvalues_k(loop_kpt,ll).lt.EF)then
         n_state(ii)=n_state(ii)+abs(Z_k(loop_kpt,ii,ll))**2.
        endif
        enddo
       enddo
      enddo

      write(*,*)''
      do ll=1,num_wann
       write(*,*)' Occupancies : ',ll,n_state(ll)/real(num_kpts)
       n_state(ll)=n_state(ll)/real(num_kpts)
      enddo

c 
      write(*,*)''
      sumaux=0.d0
      do ll=1,num_wann
      sumaux=sumaux+n_state(ll)
      enddo
      write(*,*)' Sum Occupancies : ',sumaux
c      
      write(*,*)' '
      sumaux=0.d0
      do ll=1,num_wann_ss
      sumaux=sumaux+n_state(ll)
      enddo
      navhartree=(real(sumaux)/real(num_wann_ss))
      write(*,*)' n0 averaged Hartree = ',navhartree
c
      write(*,*)''
      write(*,*)' W ',abs(eigenvalues(e_ind(1))
     &                   -eigenvalues(e_ind(Ntot)))
     &                   ,eigenvalues(e_ind(1))
     &                   ,eigenvalues(e_ind(Ntot))
      write(*,*)

      if(if_Ldos.eq.1)then

      write(*,*)''
      write(*,*)' LDOS '
      write(*,*)''

      gamma=0.01
      write(*,*)' smearing ',gamma
      write(*,*)''
      pi=acos(-1.d0)

      do ii=1,num_wann
      do ie=1,mne
       LdosE(ii,ie)=0.d0
      enddo
      enddo

      emin=eigenvalues(e_ind(1))-2.d0
      emax=eigenvalues(e_ind(Ntot))+2.d0

      do i=1,mne
       egrid(i)=emin+real(i-1)*(emax-emin)/real(mne-1)
      enddo

      do ii=1,num_wann
       do ie=1,mne
        LdosE(ii,ie)=0.d0
        do loop_kpt=1,num_kpts
         do ll=1,num_wann
           LdosE(ii,ie)=LdosE(ii,ie)+
     &     (abs(Z_k(loop_kpt,ii,ll))**2)
     &     *gamma/pi/
     &     (gamma**2
     &     +((egrid(ie))
     &     -(eigenvalues_k(loop_kpt,ll)))**2)
         enddo
        enddo
       enddo
      enddo

      write(*,*)
      write(*,*)' DOS '
      write(*,*)''
      sumd=0.d0
      do ie=1,mne
       do ii=1,num_wann
        sumd=sumd+LdosE(ii,ie)/real(num_kpts)
       enddo
      enddo
      write(*,*)''
      write(*,*)' int dos(E) de = ',sumd*abs((emax-emin)/real(mne))

      do ii=1,num_wann
       write(x1,'(I5.5)')ii
       filename='Ldos'//x1//''
       open(700,file=filename,form='formatted',status='unknown')
       do ie=1,mne
         write(700,*)egrid(ie),LdosE(ii,ie)/real(num_kpts)
       enddo
       close(700)
      enddo

      goto 3000

      endif

c calculations of < f^+ f >
      write(*,*)
      write(*,*)' < f^+ f > '
      call idate(today)
      call itime(now)
      write ( *, 1009 )  today(2), today(1), today(3), now
 1009 format ( '  Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
      write(*,*)

      fnorm=1.d0/dfloat(num_kpts)

      do loop_rpt=1,nrpts
       do i=1,num_wann
        do l=i,num_wann 
         
         fdagf(i,l,loop_rpt)=0.d0

         sumaux=cmplx(0.d0,0.d0)
         do loop_kpt=1,num_kpts

          rdotk=0.d0
          do ll=1,ndim
           rdotk=rdotk+factor_hr*(kpt_latt(ll,loop_kpt)
     &                      *(irvec(ll,loop_rpt)))
          enddo

           factexp=dcmplx(cos(rdotk),-sin(rdotk))

          do jj=1,num_wann
           if(eigenvalues_k(loop_kpt,jj).le.EF)then
            sumaux=sumaux
     &     +dconjg(Z_k(loop_kpt,i,jj))*Z_k(loop_kpt,l,jj)
     &     *factexp
           endif
       
          enddo

         enddo

         fdagf(i,l,loop_rpt)=sumaux*fnorm

        enddo
       enddo
      enddo

      write(*,*)
      write(*,*)
      write(*,*)' < f^+ f > calculated '
      call idate(today)
      call itime(now)
      write ( *, 1010 )  today(2), today(1), today(3), now
 1010 format ( '  Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
      write(*,*)
c
      sumaux=0.d0
      do loop_rpt=1,nrpts
       do i=1,num_wann
        do j=i+1,num_wann
        sumaux=sumaux+fdagf(i,j,loop_rpt)
     &                           *ham_r(i,j,loop_rpt)
        enddo
       enddo
      enddo
c 
      sumaux=2.d0*sumaux
c 
      do loop_rpt=1,nrpts
       do i=1,num_wann
        sumaux=sumaux+fdagf(i,i,loop_rpt)
     &                           *ham_r(i,i,loop_rpt)
       enddo
      enddo
      write(*,*)' check ene tot \sum < f^+_i f_j > = ',sumaux
      write(*,*)
c 
      write(*,*)' t < f^+ f > '
      write(*,*)

      do loop_rpt=1,nrpts
       do i=1,num_wann
        do j=1,num_wann
        ham_r(i,j,loop_rpt)=dcmplx(0.d0,0.d0)
        enddo
       enddo
      enddo
      do loop_rpt=1,nrpts
       do i=1,num_wann
        do j=1,num_wann
         ham_r(i,j,loop_rpt)=ham_r2(i,j,loop_rpt)
        enddo
       enddo
      enddo
c 
      write(*,*)' '
      do loop_rpt=1,nrpts
       if(irvec(1,loop_rpt).eq.0.and.
     &    irvec(2,loop_rpt).eq.0.and.
     &    irvec(3,loop_rpt).eq.0)then
       do i=1,num_wann

          write(*,*)' i,e ',i,ham_r(i,i,loop_rpt)
          ham_r(i,i,loop_rpt)=dcmplx(0.d0,0.d0)

       enddo
       endif
      enddo
      write(*,*)' '

      do loop_rpt=1,nrpts
       do i=1,num_wann
        do j=i,num_wann
        Jijsigma(i,j,loop_rpt,1)=fdagf(i,j,loop_rpt)
     &                         *ham_r(i,j,loop_rpt)
        Jijsigma(i,j,loop_rpt,2)=fdagf(i,j,loop_rpt)
     &                         *ham_r(i,j,loop_rpt)
        enddo
       enddo
      enddo

 3000 continue

      return
      END


      SUBROUTINE spin_hamiltonian(x,isite,indx_orb_site)
      implicit none
      include 'main_SS_SS_multisites_sym.h'
c      integer indx_orb_site,isite
 
c      do iorb=1,N_site_tot
c       lambda(iorb)=x(iorb+indx_orb_site)
c      enddo

       write(*,*)' '
      if(if_sym.eq.0)then
       write(*,*)' *** symmetry OFF *** '
       write(*,*)' orb ind orb '
       do iorb=1,N_site_tot
        write(*,*)iorb,iorb+indx_orb_site,
     &            index_sym(iorb+indx_orb_site)
       enddo
       write(*,*)
      endif
      if(if_sym.eq.1)then
       write(*,*)' *** symmetry ON *** '
       write(*,*)' orb ind orb '
       do iorb=1,N_site_tot
        write(*,*)iorb,iorb+indx_orb_site,
     &            index_sym(iorb+indx_orb_site)
       enddo
       write(*,*)
      endif

      do iorb=1,N_site_tot
       lambda(iorb)=x(index_sym(iorb+indx_orb_site))
      enddo

      Uprime=Uvalue-2*Jh
c 
      write(*,*)''
      write(*,*)' generating states HS'
      write(*,*)' dn'
      J_dn = 0
      do i=1,2**N_site_tot 
        J_dn = J_dn + 1  
        I_dn(J_dn) = i   
      enddo

      write(*,*)' up'
      J_up = 0
      do j=1,2**N_site_tot
        J_up=J_up+1     
        I_up(J_up) = j  
      enddo

      write(*,*)''
      write(*,*)'',J_up,J_dn
      Ns_up=J_up
      Ns_dn=J_dn
      write(*,*)''
      write(*,*)' Hilbert Space dimensions '
     &         ,Ns_up,Ns_up,2**(2*N_site_tot) ! 2**(2*N_site_tot) = # states in Hilbert space 

      Nstate=2**(2*N_site_tot)
c 
      write(*,*)''
      write(*,*)' <i|H|j>'

      do i=1,Nstate
       do j=1,Nstate
       haccaS(i,j)=cmplx(0.d0,0.d0)
       enddo
      enddo

      write(*,*)''
      write(*,*)' diagonal terms '

c lambda terms:
c for such term see Eq. 14 http://arxiv.org/pdf/cond-mat/0503764

      do i=1,Ns_dn     
       do j=1,Ns_up    

       sumtot=0.d0

        sumaux1=0.d0
        do jorb=1,N_site_tot 
         sumup=0.d0
         bool1 = btest(I_up(j)-1,jorb-1)
         if(bool1.eqv..TRUE.)then          
          sumup=sumup+0.5d0             
         endif
         if(bool1.eqv..FALSE.)then         
          sumup=sumup-0.5d0
         endif
          sumaux1=sumaux1+lambda(jorb)*(sumup+0.5d0) 
        enddo

        sumaux2=0.d0
        do jorb=1,N_site_tot
         sumdn=0.d0
         bool2 = btest(I_dn(i)-1,jorb-1)
         if(bool2.eqv..TRUE.)then
          sumdn=sumdn+0.5d0
         endif
         if(bool2.eqv..FALSE.)then
          sumdn=sumdn-0.5d0
         endif
          sumaux2=sumaux2+lambda(jorb)*(sumdn+0.5d0) 
        enddo

        sumtot=sumtot+(sumaux1+sumaux2)

       ii = Ns_up*(i-1)+j              

       haccaS(ii,ii) = haccaS(ii,ii)
     &                +cmplx(sumtot,0.d0)

       enddo
      enddo
c first term Eq. 18 of the paper "Orbital selective Mott Transition in multi-band systems..." 
c http://arxiv.org/pdf/cond-mat/0503764
      do i=1,Ns_dn
       do j=1,Ns_up

       sum1=0.d0

        sumup=0.d0
        do jorb=1,N_site_tot
        bool1 = btest(I_up(j)-1,jorb-1)
        if(bool1.eqv..TRUE.)then 
         sumup=sumup+0.5d0
        endif
        if(bool1.eqv..FALSE.)then 
         sumup=sumup-0.5d0
        endif
        enddo

        sumdn=0.d0
        do jorb=1,N_site_tot
        bool2 = btest(I_dn(i)-1,jorb-1)
        if(bool2.eqv..TRUE.)then 
         sumdn=sumdn+0.5d0
        endif
        if(bool2.eqv..FALSE.)then 
         sumdn=sumdn-0.5d0
        endif
        enddo

        sum1=sum1+(sumup+sumdn)**2

       ii = Ns_up*(i-1)+j

       haccaS(ii,ii) = haccaS(ii,ii)
     &                +cmplx(0.5d0*Uprime*sum1,0.d0)

       enddo
      enddo

      if(Jh.gt.0)then
c second term Eq. 18 of the paper "Orbital selective Mott Transition in multi-band systems..."
c http://arxiv.org/pdf/cond-mat/0503764

      do i=1,Ns_dn
       do j=1,Ns_up

        sum2=0.d0

        do jorb=1,N_site_tot

         sumup=0.d0
         bool1 = btest(I_up(j)-1,jorb-1)
         if(bool1.eqv..TRUE.)then
          sumup=sumup+0.5d0
         endif
         if(bool1.eqv..FALSE.)then
          sumup=sumup-0.5d0
         endif

         sumdn=0.d0
         bool2 = btest(I_dn(i)-1,jorb-1)
         if(bool2.eqv..TRUE.)then
          sumdn=sumdn+0.5d0
         endif
         if(bool2.eqv..FALSE.)then
          sumdn=sumdn-0.5d0
         endif

         sum2=sum2+(sumup+sumdn)**2

        enddo

        ii = Ns_up*(i-1)+j
        haccaS(ii,ii) = haccaS(ii,ii)
     &                 +cmplx(Jh*sum2,0.d0)

       enddo
      enddo

c third term Eq. 18 of the paper "Orbital selective Mott Transition in multi-band systems..."
c http://arxiv.org/pdf/cond-mat/0503764
      do i=1,Ns_dn
       do j=1,Ns_up

        sum3=0.d0

         sumup=0.d0
         do jorb=1,N_site_tot
          bool1 = btest(I_up(j)-1,jorb-1)
          if(bool1.eqv..TRUE.)then
           sumup=sumup+0.5d0
          endif
          if(bool1.eqv..FALSE.)then
           sumup=sumup-0.5d0
          endif
         enddo
         sumup=sumup**2

         sumdn=0.d0
         do jorb=1,N_site_tot
          bool2 = btest(I_dn(i)-1,jorb-1)
          if(bool2.eqv..TRUE.)then
           sumdn=sumdn+0.5d0
          endif
          if(bool2.eqv..FALSE.)then
           sumdn=sumdn-0.5d0
          endif
         enddo
         sumdn=sumdn**2

         sum3=sum3+sumup+sumdn

        ii = Ns_up*(i-1)+j
        haccaS(ii,ii) = haccaS(ii,ii)
     &                 +cmplx(-0.5d0*Jh*sum3,0.d0)
       enddo
      enddo

      write(*,*)''
      write(*,*)' <i|H|j>'
      write(*,*)' off-diagonal terms '

      write(*,*)''
      write(*,*)' S-F(Spin-Flip) and P-H(Pair-Hopping) terms '

      write(*,*)whattypeofJh 

      if(whattypeofJh.eq.1)then

      write(*,*)''
      write(*,*)' S-F(Spin-Flip) '

c S-F:
      do i=1,Ns_dn                        
       do j=1,Ns_up

        count_dn_ini=0                    
        do jorb=1,N_site_tot
        bool1 = btest(I_dn(i)-1,jorb-1)
        if(bool1.eqv..TRUE.)then
         count_dn_ini=count_dn_ini+1
        endif
        enddo
        count_up_ini=0
        do jorb=1,N_site_tot
        bool1 = btest(I_up(j)-1,jorb-1)   
        if(bool1.eqv..TRUE.)then
         count_up_ini=count_up_ini+1
        endif
        enddo

         do ii=1,Ns_dn
          do jj=1,Ns_up

        count_dn_fin=0                     
        do jorb=1,N_site_tot
        bool1 = btest(I_dn(ii)-1,jorb-1)
        if(bool1.eqv..TRUE.)then
         count_dn_fin=count_dn_fin+1
        endif
        enddo
        count_up_fin=0                     
        do jorb=1,N_site_tot
        bool1 = btest(I_up(jj)-1,jorb-1)
        if(bool1.eqv..TRUE.)then
         count_up_fin=count_up_fin+1
        endif
        enddo

        if(count_dn_ini.eq.count_dn_fin)then 
        if(count_up_ini.eq.count_up_fin)then 

             do iorb=1,N_site_tot
              do jorb=iorb+1,N_site_tot

              bool1 = btest(I_up(j)-1,iorb-1)
              bool2 = btest(I_up(j)-1,jorb-1)
              bool3 = btest(I_dn(i)-1,iorb-1)
              bool4 = btest(I_dn(i)-1,jorb-1)

              if((bool1.EQV..FALSE.).and.     
     &           (bool2.EQV..TRUE.).and.
     &           (bool3.EQV..TRUE.).and.
     &           (bool4.EQV..FALSE.))then

                Ip_up=IBSET(IBCLR(I_up(j)-1,jorb-1)      
     &                                     ,iorb-1)+1
                Ip_dn=IBSET(IBCLR(I_dn(i)-1,iorb-1)
     &                                     ,jorb-1)+1

                if(I_up(jj).eq.Ip_up)then                
                if(I_dn(ii).eq.Ip_dn)then

                 iii=Ns_up*(i-1)+j
                 jjj=Ns_up*(ii-1)+jj

                 haccaS(min(iii,jjj),max(iii,jjj))=     
     &           haccaS(min(iii,jjj),max(iii,jjj))
     &           +cmplx(-Jh,0.d0)

                endif
                endif
              endif

              enddo
             enddo

           endif
           endif

          enddo
         enddo
       enddo
      enddo

      write(*,*)''
      write(*,*)' P-H(Pair-Hopping) '
c P-H:
      do i=1,Ns_dn      
       do j=1,Ns_up    

        count_dn_ini=0
        do jorb=1,N_site_tot
        bool1 = btest(I_dn(i)-1,jorb-1)
        if(bool1.eqv..TRUE.)then
         count_dn_ini=count_dn_ini+1
        endif
        enddo
        count_up_ini=0
        do jorb=1,N_site_tot
        bool1 = btest(I_up(j)-1,jorb-1)
        if(bool1.eqv..TRUE.)then
         count_up_ini=count_up_ini+1
        endif
        enddo

         do ii=1,Ns_dn
          do jj=1,Ns_up

        count_dn_fin=0
        do jorb=1,N_site_tot
        bool1 = btest(I_dn(ii)-1,jorb-1)
        if(bool1.eqv..TRUE.)then
         count_dn_fin=count_dn_fin+1
        endif
        enddo
        count_up_fin=0
        do jorb=1,N_site_tot
        bool1 = btest(I_up(jj)-1,jorb-1)
        if(bool1.eqv..TRUE.)then
         count_up_fin=count_up_fin+1
        endif
        enddo

        if(count_dn_ini.eq.count_dn_fin)then
        if(count_up_ini.eq.count_up_fin)then

             do iorb=1,N_site_tot
              do jorb=iorb+1,N_site_tot

              bool1 = btest(I_up(j)-1,iorb-1)
              bool2 = btest(I_up(j)-1,jorb-1)
              bool3 = btest(I_dn(i)-1,iorb-1)
              bool4 = btest(I_dn(i)-1,jorb-1)

              if((bool1.EQV..FALSE.).and.
     &           (bool2.EQV..TRUE.).and.
     &           (bool3.EQV..FALSE.).and.
     &           (bool4.EQV..TRUE.))then

                Ip_up=IBSET(IBCLR(I_up(j)-1,jorb-1)
     &                                     ,iorb-1)+1
                Ip_dn=IBSET(IBCLR(I_dn(i)-1,jorb-1)
     &                                     ,iorb-1)+1

               if(I_up(jj).eq.Ip_up)then
               if(I_dn(ii).eq.Ip_dn)then

                 iii=Ns_up*(i-1)+j
                 jjj=Ns_up*(ii-1)+jj

                 haccaS(min(iii,jjj),max(iii,jjj))=
     &           haccaS(min(iii,jjj),max(iii,jjj))
     &           +cmplx(-Jh,0.d0)

                endif
                endif
              endif

              enddo
             enddo

           endif
           endif
 
          enddo
         enddo
       enddo
      enddo

      endif ! if whattypeofJh

      endif 

      write(*,*)' '
      write(*,*)' c_m,sigma in Hs (spin polarized version) '

      do i=1,nspin
       do j=1,N_site_tot
       cv=1.d0/sqrt(n_state(j+indx_orb_site)
     &            *(1-n_state(j+indx_orb_site)))-1.d0
       cm(j,i)=cmplx(cv,0.d0)
       write(*,*)cm(j,i)
       enddo
      enddo
      write(*,*)' '

      write(*,*)' '
      write(*,*)' O^+_i <O_i> terms UP '
      write(*,*)' '

      do i=1,Ns_dn  
       do j=1,Ns_up 
 
        iii=Ns_up*(i-1)+j  

        do iorb=1,N_site_tot  
                            
         indx_iorb=iorb+indx_orb_site

         bool1 = btest(I_up(j)-1,iorb-1) 
c S^+_i                                  
         if(bool1.EQV..FALSE.)then       

         Ip_up=IBSET(I_up(j)-1,iorb-1)+1 

         do ii=1,Ns_dn                   
          do jj=1,Ns_up                  

          if(Ip_up.eq.I_up(jj))then      
          if(I_dn(i).eq.I_dn(ii))then    

          jjj=Ns_up*(ii-1)+jj            

                  do loop_rpt=1,nrpts     

                  do jorb=1,num_wann_ss

               if(indx_iorb.le.jorb)then     

               factorHM=Jijsigma(indx_iorb,jorb,loop_rpt,1)
     &         *mOi_orb(jorb)

                if(jjj.gt.iii)then
 
                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM

                else

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

               else                      

               factorHM=Jijsigma(jorb,indx_iorb,loop_rpt,1)
     &         *mOi_orb(jorb) 

                if(jjj.gt.iii)then
 
                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM

                else

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

               endif

                  enddo

                  enddo
 
          endif
          endif

          enddo
         enddo

         endif
c S^-_i c^*_i,up                         
         if(bool1.EQV..TRUE.)then        

         Ip_up=IBCLR(I_up(j)-1,iorb-1)+1

         do ii=1,Ns_dn
          do jj=1,Ns_up

          if(Ip_up.eq.I_up(jj))then
          if(I_dn(i).eq.I_dn(ii))then

          jjj=Ns_up*(ii-1)+jj

                  do loop_rpt=1,nrpts

                  do jorb=1,num_wann_ss

               if(indx_iorb.le.jorb)then

               factorHM=Jijsigma(indx_iorb,jorb,loop_rpt,1)
     &         *CONJG(cm(iorb,1))
     &         *mOi_orb(jorb)

                if(jjj.gt.iii)then

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM

                else

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

               else

               factorHM=Jijsigma(jorb,indx_iorb,loop_rpt,1)
     &         *CONJG(cm(iorb,1))
     &         *mOi_orb(jorb)

                if(jjj.gt.iii)then

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM
               
                else
               
                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)
     
                endif

               endif

                  enddo

                  enddo

          endif
          endif

          enddo
         enddo

         endif
c
        enddo

       enddo
      enddo


      write(*,*)' '
      write(*,*)' O^+_i <O_j> terms DN '
      write(*,*)' '

      do i=1,Ns_dn
       do j=1,Ns_up

        iii=Ns_up*(i-1)+j

        do iorb=1,N_site_tot

         indx_iorb=iorb+indx_orb_site

         bool1 = btest(I_dn(i)-1,iorb-1)
c S^+_i
         if(bool1.EQV..FALSE.)then

         Ip_dn=IBSET(I_dn(i)-1,iorb-1)+1

         do ii=1,Ns_dn
          do jj=1,Ns_up

          if(Ip_dn.eq.I_dn(ii))then
          if(I_up(j).eq.I_up(jj))then

          jjj=Ns_up*(ii-1)+jj

                  do loop_rpt=1,nrpts

                  do jorb=1,num_wann_ss

               if(indx_iorb.le.jorb)then

               factorHM=Jijsigma(indx_iorb,jorb,loop_rpt,2)
     &         *mOi_orb(jorb+num_wann_ss)

                if(jjj.gt.iii)then

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM

                else

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

               else

               factorHM=Jijsigma(jorb,indx_iorb,loop_rpt,2)
     &         *mOi_orb(jorb+num_wann_ss)

                if(jjj.gt.iii)then

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM
     
                else
               
                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

               endif

                  enddo

                  enddo

          endif
          endif

          enddo
         enddo

         endif
c S^-_i c^*_i,dn
         if(bool1.EQV..TRUE.)then

         Ip_dn=IBCLR(I_dn(i)-1,iorb-1)+1

         do ii=1,Ns_dn
          do jj=1,Ns_up

          if(Ip_dn.eq.I_dn(ii))then
          if(I_up(j).eq.I_up(jj))then

          jjj=Ns_up*(ii-1)+jj

                  do loop_rpt=1,nrpts

                  do jorb=1,num_wann_ss

               if(indx_iorb.le.jorb)then

               factorHM=Jijsigma(indx_iorb,jorb,loop_rpt,2)
     &         *CONJG(cm(iorb,2))
     &         *mOi_orb(jorb+num_wann_ss)

                if(jjj.gt.iii)then

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM

                else

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

               else

               factorHM=Jijsigma(jorb,indx_iorb,loop_rpt,2)
     &         *CONJG(cm(iorb,2))
     &         *mOi_orb(jorb+num_wann_ss)

                if(jjj.gt.iii)then

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM

                else

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

               endif

                  enddo

                  enddo

          endif
          endif

          enddo
         enddo

         endif
c
        enddo

       enddo
      enddo

c

c dp model 
      if(num_wann_ss.ne.num_wann)then

      write(*,*)
      write(*,*)' hybridization ss UP'
      write(*,*)

c UP sector

      do i=1,Ns_dn  
       do j=1,Ns_up 

        iii=Ns_up*(i-1)+j  

        do iorb=1,N_site_tot

         indx_iorb=iorb+indx_orb_site

         bool1 = btest(I_up(j)-1,iorb-1) 
c S^+_i                                  
         if(bool1.EQV..FALSE.)then       

         Ip_up=IBSET(I_up(j)-1,iorb-1)+1 

         do ii=1,Ns_dn                   
          do jj=1,Ns_up                  

          if(Ip_up.eq.I_up(jj))then      
          if(I_dn(i).eq.I_dn(ii))then    

          jjj=Ns_up*(ii-1)+jj            

                  do loop_rpt=1,nrpts     

                  do jorb=num_wann_ss+1,num_wann      

               factorHM=Jijsigma(indx_iorb,jorb,loop_rpt,1)

                if(jjj.gt.iii)then

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM

                else

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

                 enddo

                 enddo

          endif
          endif 

          enddo
         enddo

         endif

c S^-_i c^*_i,up                         
         if(bool1.EQV..TRUE.)then        

         Ip_up=IBCLR(I_up(j)-1,iorb-1)+1

         do ii=1,Ns_dn
          do jj=1,Ns_up

          if(Ip_up.eq.I_up(jj))then
          if(I_dn(i).eq.I_dn(ii))then

          jjj=Ns_up*(ii-1)+jj

                 do loop_rpt=1,nrpts

                  do jorb=num_wann_ss+1,num_wann

               factorHM=Jijsigma(indx_iorb,jorb,loop_rpt,1)
     &         *CONJG(cm(iorb,1))

                if(jjj.gt.iii)then

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM

                else

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

                 enddo

                 enddo

          endif 
          endif 

          enddo
         enddo

         endif

        enddo 

       enddo 
      enddo 

      write(*,*)
      write(*,*)' hybridization ss DN'
      write(*,*)

c DN sector

      do i=1,Ns_dn  
       do j=1,Ns_up 

        iii=Ns_up*(i-1)+j  

        do iorb=1,N_site_tot

         indx_iorb=iorb+indx_orb_site

         bool1 = btest(I_dn(i)-1,iorb-1) 
c S^+_i                                  
         if(bool1.EQV..FALSE.)then       

         Ip_dn=IBSET(I_dn(i)-1,iorb-1)+1 

         do ii=1,Ns_dn                   
          do jj=1,Ns_up                  

          if(Ip_dn.eq.I_dn(ii))then      
          if(I_up(j).eq.I_up(jj))then    

          jjj=Ns_up*(ii-1)+jj            

                  do loop_rpt=1,nrpts     

                  do jorb=num_wann_ss+1,num_wann      

               factorHM=Jijsigma(indx_iorb,jorb,loop_rpt,2)

                if(jjj.gt.iii)then

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM

                else

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

                 enddo

                 enddo

          endif
          endif 

          enddo
         enddo

         endif

c S^-_i c^*_i,up                         
         if(bool1.EQV..TRUE.)then        

         Ip_dn=IBCLR(I_dn(i)-1,iorb-1)+1

         do ii=1,Ns_dn
          do jj=1,Ns_up

          if(Ip_dn.eq.I_dn(ii))then
          if(I_up(j).eq.I_up(jj))then

          jjj=Ns_up*(ii-1)+jj

                  do loop_rpt=1,nrpts

                  do jorb=num_wann_ss+1,num_wann

               factorHM=Jijsigma(indx_iorb,jorb,loop_rpt,2)
     &         *CONJG(cm(iorb,2))

                if(jjj.gt.iii)then

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +factorHM

                else

                haccaS(min(iii,jjj),max(iii,jjj))=
     &          haccaS(min(iii,jjj),max(iii,jjj))
     &          +CONJG(factorHM)

                endif

                 enddo

                 enddo

          endif 
          endif 

          enddo
         enddo

         endif

        enddo 

       enddo 
      enddo 

      endif

c diagonalization:

      write(*,*)''
      write(*,*)' Diagonalization SS problem'
      call idate(today)
      call itime(now)
      write ( *, 1030 )  today(2), today(1), today(3), now
 1030 format ( '  Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
      write(*,*)

       if(whattypeofdiago.eq.0)then
        write(*,*)
        write(*,*)' LAPACK ZHEEVX'
        write(*,*)

       call zheevx('V'
     &            ,'A'
     &            ,'U'
     &            ,Nstate
     &            ,haccaS
     &            ,Nstate
     &            ,ivl,ivu
     &            ,1
     &            ,Nstate
     &            ,0.d0
     &            ,neig,eig
     &            ,ZZ,Nstate
     &            ,WWORK,LLWORK,RRWORK,IIWORK,IIFAIL,IINFO)

      if(iinfo.ne.0)then
       write(*,*)' Failure in ZHEEVX. INFO =',iinfo
       stop
      endif

      write(*,*)''
      write(*,*)' End Diagonalization SS problem'
      call idate(today)
      call itime(now)
      write ( *, 1031 )  today(2), today(1), today(3), now
 1031 format ( '  Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
      write(*,*)

      Ndegen=1    
      do i=2,Nstate
       if(eig(i).eq.eig(1))then
        Ndegen=Ndegen+1
       endif
      enddo
      write(*,*)' '
      write(*,*)' degenerecy GS Hs ',Ndegen
      write(*,*)' '

      write(*,*)' '
      write(*,*)' eigenvalues Hs '
      write(*,*)' '

      do i=1,Nstate
       write(*,*)i,eig(i)
      enddo
      write(*,*)' '
c
       do i=1,Ns_dn
        do j=1,Ns_up

        ii = Ns_up*(i-1)+j

        ZZ_gs(ii)=ZZ(ii,1)

        enddo
       enddo

      endif 

      if(whattypeofdiago.eq.1)then
       write(*,*)
       write(*,*)' Lanczos HLZPACK '
       write(*,*)

       do i=1,Nstate
       do j=1,Nstate

        if(i.gt.j)then
        haccaS(i,j)=CONJG(haccaS(j,i))
        endif

       enddo
       enddo

       call Lanczos(haccaS,diag,ZA) 

      write(*,*)''
      write(*,*)' End Diagonalization SS problem'
      call idate(today)
      call itime(now)
      write ( *, 1032 )  today(2), today(1), today(3), now
 1032 format ( '  Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
      write(*,*)


      Ndegen=1
      do i=2,Nstate
       if(diag(i).eq.diag(1))then
        Ndegen=Ndegen+1
       endif
      enddo
      write(*,*)' '
      write(*,*)' degenerecy GS Hs ',Ndegen
      write(*,*)' '

      do i=1,NmaxiterL
       write(*,*)i,diag(i)
      enddo
      write(*,*)' '

       do i=1,Ns_dn
        do j=1,Ns_up

        ii = Ns_up*(i-1)+j

        ZZ_gs(ii)=ZA(ii)

        enddo
       enddo

      endif

c <Sz>

      write(*,*)''
      write(*,*)' calculation <Sz>'
      write(*,*)''

      do iorb=1,N_site_tot 

       Sz(1,iorb)=cmplx(0.d0,0.d0)
       Sz(2,iorb)=cmplx(0.d0,0.d0)

       sum1aux=cmplx(0.d0,0.d0)
       sum2aux=cmplx(0.d0,0.d0)

       do i=1,Ns_dn       
        do j=1,Ns_up     

        ii = Ns_up*(i-1)+j 

        sumup=0.d0
        bool1 = btest(I_up(j)-1,iorb-1)
        if(bool1.eqv..TRUE.)then        
         sumup=sumup+0.5d0
        endif
        if(bool1.eqv..FALSE.)then       
         sumup=sumup-0.5d0
        endif

        sum1aux=sum1aux+sumup*CONJG(ZZ_gs(ii))*ZZ_gs(ii)

        sumdn=0.d0
        bool2 = btest(I_dn(i)-1,iorb-1)
        if(bool2.eqv..TRUE.)then 
         sumdn=sumdn+0.5d0
        endif
        if(bool2.eqv..FALSE.)then 
         sumdn=sumdn-0.5d0
        endif

        sum2aux=sum2aux+sumdn*CONJG(ZZ_gs(ii))*ZZ_gs(ii)

        enddo
       enddo

       Sz(1,iorb)=sum1aux
       Sz(2,iorb)=sum2aux
      enddo

      write(*,*)''
      write(*,*)' mean values <Sz> (spin comp 1, 2)'
      do iorb=1,N_site_tot
       write(*,*)iorb
       write(*,*)Sz(1,iorb),Sz(2,iorb)
      enddo
      write(*,*)''

      do i=1,nspin*N_site_tot
       mSz(i)=0.d0
      enddo
c up
      indx=0
      do iorb=1,N_site_tot
      indx=indx+1
      mSz(indx)=real(Sz(1,iorb)) 
      enddo
c dn
      do iorb=1,N_site_tot
      indx=indx+1
      mSz(indx)=real(Sz(2,iorb))
      enddo

      write(*,*)''

      write(*,*)''
      write(*,*)' calculation <O> '
      write(*,*)''

c <O_i>
      do iorb=1,N_site_tot 

       Oi(iorb,1)=cmplx(0.d0,0.d0)
       Oi(iorb,2)=cmplx(0.d0,0.d0)

       sum1aux=cmplx(0.d0,0.d0) 

       do j=1,Ns_up        
        do jj=1,Ns_up      

              bool1 = btest(I_up(j)-1,iorb-1)
c S^-_i
              if(bool1.EQV..TRUE.)then         

               Ip_up=IBCLR(I_up(j)-1,iorb-1)+1 

               if(Ip_up.eq.I_up(jj))then

                 do k=1,Ns_dn          

                  iii=Ns_up*(k-1)+j    
                  jjj=Ns_up*(k-1)+jj   

                  sum1aux=sum1aux+
     &            +CONJG(ZZ_gs(jjj))*ZZ_gs(iii)

                 enddo

               endif
              endif
c S^+_i c_i,up
              if(bool1.EQV..FALSE.)then       

               Ip_up=IBSET(I_up(j)-1,iorb-1)+1 

               if(Ip_up.eq.I_up(jj))then

                 do k=1,Ns_dn                 

                  iii=Ns_up*(k-1)+j          
                  jjj=Ns_up*(k-1)+jj        

                  sum1aux=sum1aux+
     &            +cm(iorb,1)
     &            *CONJG(ZZ_gs(jjj))*ZZ_gs(iii)

                 enddo

               endif
              endif

        enddo
       enddo

       Oi(iorb,1)=sum1aux
c 
       sum2aux=cmplx(0.d0,0.d0)

       do i=1,Ns_dn      
        do ii=1,Ns_dn    
c S^-_i
             bool1 = btest(I_dn(i)-1,iorb-1)

              if(bool1.EQV..TRUE.)then   

               Ip_dn=IBCLR(I_dn(i)-1,iorb-1)+1 

               if(Ip_dn.eq.I_dn(ii))then

                 do k=1,Ns_up            

                  iii=Ns_up*(i-1)+k      
                  jjj=Ns_up*(ii-1)+k     

                  sum2aux=sum2aux+
     &            +CONJG(ZZ_gs(jjj))*ZZ_gs(iii)

                 enddo

               endif
              endif
c S^+_i c_i,dn
              if(bool1.EQV..FALSE.)then   

               Ip_dn=IBSET(I_dn(i)-1,iorb-1)+1  

               if(Ip_dn.eq.I_dn(ii))then

                 do k=1,Ns_up             

                  iii=Ns_up*(i-1)+k       
                  jjj=Ns_up*(ii-1)+k      

                  sum2aux=sum2aux+
     &            +cm(iorb,2)
     &            *CONJG(ZZ_gs(jjj))*ZZ_gs(iii)

                 enddo

               endif
              endif

        enddo
       enddo

       Oi(iorb,2)=sum2aux

       enddo
    
       write(*,*)' '
       write(*,*)' mean values <O_i> '
       do iorb=1,N_site_tot
        write(*,*)Oi(iorb,1),Oi(iorb,2)
       enddo
       write(*,*)' '
c
       indx=0
       do i=1,nspin
       do j=1,N_site_tot
       indx=indx+1
       mOi(indx)=Oi(j,i) 
       enddo
       enddo
       write(*,*)' '

c <O^+_i>
      do iorb=1,N_site_tot  

       Opi(iorb,1)=cmplx(0.d0,0.d0)
       Opi(iorb,2)=cmplx(0.d0,0.d0)

       sum1aux=cmplx(0.d0,0.d0)

       do j=1,Ns_up         
        do jj=1,Ns_up       

              bool1 = btest(I_up(j)-1,iorb-1)
c S^+_i
              if(bool1.EQV..FALSE.)then         

               Ip_up=IBSET(I_up(j)-1,iorb-1)+1  

               if(Ip_up.eq.I_up(jj))then

                 do k=1,Ns_dn        

                  iii=Ns_up*(k-1)+j   
                  jjj=Ns_up*(k-1)+jj  

                  sum1aux=sum1aux+
     &            +CONJG(ZZ_gs(jjj))*ZZ_gs(iii)

                 enddo

               endif
              endif
c S^-_i c^*_i,up
              if(bool1.EQV..TRUE.)then          

               Ip_up=IBCLR(I_up(j)-1,iorb-1)+1  

               if(Ip_up.eq.I_up(jj))then

                 do k=1,Ns_dn          

                  iii=Ns_up*(k-1)+j    
                  jjj=Ns_up*(k-1)+jj   

                  sum1aux=sum1aux+
     &            +CONJG(cm(iorb,1))
     &            *CONJG(ZZ_gs(jjj))*ZZ_gs(iii)

                 enddo

               endif
              endif

        enddo
       enddo

       Opi(iorb,1)=sum1aux
c
       sum2aux=cmplx(0.d0,0.d0)

       do i=1,Ns_dn             
        do ii=1,Ns_dn           
c S^+_i
             bool1 = btest(I_dn(i)-1,iorb-1) 

              if(bool1.EQV..FALSE.)then        

               Ip_dn=IBSET(I_dn(i)-1,iorb-1)+1 

               if(Ip_dn.eq.I_dn(ii))then     

                 do k=1,Ns_up           

                  iii=Ns_up*(i-1)+k     
                  jjj=Ns_up*(ii-1)+k    

                  sum2aux=sum2aux+
     &            +CONJG(ZZ_gs(jjj))*ZZ_gs(iii)

                 enddo

               endif
              endif
c S^-_i c^*_i,dn
              if(bool1.EQV..TRUE.)then         

               Ip_dn=IBCLR(I_dn(i)-1,iorb-1)+1 

               if(Ip_dn.eq.I_dn(ii))then       

                 do k=1,Ns_up                  
 
                  iii=Ns_up*(i-1)+k            
                  jjj=Ns_up*(ii-1)+k           

                  sum2aux=sum2aux+
     &            +CONJG(cm(iorb,2))
     &            *CONJG(ZZ_gs(jjj))*ZZ_gs(iii)

                 enddo

               endif
              endif

        enddo
       enddo

       Opi(iorb,2)=sum2aux

       enddo
    
       write(*,*)' '
       write(*,*)' mean values <O^+_i> '
       do iorb=1,N_site_tot
        write(*,*)Opi(iorb,1),Opi(iorb,2)
       enddo
c
       indx=0
       do i=1,nspin
       do j=1,N_site_tot
       indx=indx+1
       mOpi(indx)=Opi(j,i)     
       enddo
       enddo
       write(*,*)' '


       if(if_sym_updn.eq.1)then
        write(*,*)
        write(*,*)' SYMMETRIES ARE APPLIED BETWEEN UP and DN '
        indx=0
        do i=1,nspin
        do j=1,N_site_tot
        indx=indx+1
        mOi(indx)=Oi(j,1)
        enddo
        enddo
        write(*,*)' '
 
        indx=0
        do i=1,nspin
        do j=1,N_site_tot
        indx=indx+1
        mOpi(indx)=Opi(j,1)
        enddo
        enddo
        write(*,*)' '
       endif

c spin-spin correlation functions:
      do iorb=1,N_site_tot
      do jorb=1,N_site_tot

       SzSz(1,iorb,jorb)=cmplx(0.d0,0.d0)
       SzSz(2,iorb,jorb)=cmplx(0.d0,0.d0)
       SzSz(3,iorb,jorb)=cmplx(0.d0,0.d0)
       SzSz(4,iorb,jorb)=cmplx(0.d0,0.d0)

       sum1aux=cmplx(0.d0,0.d0)  !up-up
       sum2aux=cmplx(0.d0,0.d0)  !dn-dn
       sum3aux=cmplx(0.d0,0.d0)  !up-dn
       sum4aux=cmplx(0.d0,0.d0)  !dn-up

       do i=1,Ns_dn
        do j=1,Ns_up

        ii = Ns_up*(i-1)+j

        sumup1=0.d0
        bool1 = btest(I_up(j)-1,iorb-1)
        if(bool1.eqv..TRUE.)then
         sumup1=sumup1+0.5d0
        endif
        if(bool1.eqv..FALSE.)then
         sumup1=sumup1-0.5d0
        endif

        sumup2=0.d0
        bool1 = btest(I_up(j)-1,jorb-1)
        if(bool1.eqv..TRUE.)then
         sumup2=sumup2+0.5d0
        endif
        if(bool1.eqv..FALSE.)then
         sumup2=sumup2-0.5d0
        endif

        sum1aux=sum1aux+sumup1*sumup2*CONJG(ZZ_gs(ii))*ZZ_gs(ii)

        sumdn1=0.d0
        bool2 = btest(I_dn(i)-1,iorb-1)
        if(bool2.eqv..TRUE.)then
         sumdn1=sumdn1+0.5d0
        endif
        if(bool2.eqv..FALSE.)then
         sumdn1=sumdn1-0.5d0
        endif

        sumdn2=0.d0
        bool2 = btest(I_dn(i)-1,jorb-1)
        if(bool2.eqv..TRUE.)then
         sumdn2=sumdn2+0.5d0
        endif
        if(bool2.eqv..FALSE.)then
         sumdn2=sumdn2-0.5d0
        endif

        sum2aux=sum2aux+sumdn1*sumdn2*CONJG(ZZ_gs(ii))*ZZ_gs(ii)

        sumup1=0.d0
        bool2 = btest(I_up(j)-1,iorb-1)
        if(bool2.eqv..TRUE.)then
         sumup1=sumup1+0.5d0
        endif
        if(bool2.eqv..FALSE.)then
         sumup1=sumup1-0.5d0
        endif

        sumdn2=0.d0
        bool2 = btest(I_dn(i)-1,jorb-1)
        if(bool2.eqv..TRUE.)then
         sumdn2=sumdn2+0.5d0
        endif
        if(bool2.eqv..FALSE.)then
         sumdn2=sumdn2-0.5d0
        endif

        sum3aux=sum3aux+sumup1*sumdn2*CONJG(ZZ_gs(ii))*ZZ_gs(ii)

        sumdn1=0.d0
        bool2 = btest(I_dn(i)-1,iorb-1)
        if(bool2.eqv..TRUE.)then
         sumdn1=sumdn1+0.5d0
        endif
        if(bool2.eqv..FALSE.)then
         sumdn1=sumdn1-0.5d0
        endif

        sumup2=0.d0
        bool2 = btest(I_up(j)-1,jorb-1)
        if(bool2.eqv..TRUE.)then
         sumup2=sumup2+0.5d0
        endif
        if(bool2.eqv..FALSE.)then
         sumup2=sumup2-0.5d0
        endif

        sum4aux=sum4aux+sumdn1*sumup2*CONJG(ZZ_gs(ii))*ZZ_gs(ii)

        enddo
       enddo

       SzSz(1,iorb,jorb)=sum1aux
       SzSz(2,iorb,jorb)=sum2aux
       SzSz(3,iorb,jorb)=sum3aux
       SzSz(4,iorb,jorb)=sum4aux

      enddo
      enddo

      write(*,*)''
      write(*,*)' SzSz '
      write(*,*)''
      do iorb=1,N_site_tot
      do jorb=1,N_site_tot
       write(*,*)iorb,jorb,real(SzSz(1,iorb,jorb))
       write(*,*)iorb,jorb,real(SzSz(2,iorb,jorb))
       write(*,*)iorb,jorb,real(SzSz(3,iorb,jorb))
       write(*,*)iorb,jorb,real(SzSz(4,iorb,jorb))
      enddo
      enddo
      write(*,*)''

      return
      END

      SUBROUTINE Lanczos(haccaS,diag,ZA)
      implicit none
      include 'main_SS_SS_multisites_sym.h'

      do i=1,NmaxNstatesHS
       sirand1=DBLE(rand())
       sirand2=DBLE(rand())
       Qin(i)=cmplx(sirand1,sirand2)
       Qlin(i)=cmplx(0.d0,0.d0)
       Qsave(i)=Qin(i)
      enddo

      TB(1)=0

      NiterL=0
      do iterL = 1,NmaxiterL

       call LANCZOSstep(iterL
     &                 ,haccaS
     &                 ,NmaxNstatesHS
     &                 ,TA1,TB1
     &                 ,Qin,Qlin)

c       write(*,*)'iteration Lanczos (E_gs)',iterL

       NiterL=NiterL+1
c       write(*,*)'TA1*',TA1
c       write(*,*)'TB1*',TB1
       TA(iterL)=TA1
       TB(iterL+1)=TB1

c        do i=1,NiterL
c          do j=1,NiterL
c             Zaux(j,i)=0.d0
c          enddo
c          Zaux(i,i)=1.d0
c        enddo
c        do i=1,NiterL
c          diag(i)=TA(i)
c        end do
c        do i=2,NiterL
c          subdiag(i)=TB(i)
c        end do
c        call tql2(NmaxiterL,NiterL,diag,subdiag,Zaux,Nerr)
c        write(*,*)'normal return = 0'
c        write(*,*)Nerr
c        print *,'Output from diagonalization:'
c        print *,'------> eigenvalues <-------'
c        do i=1,NiterL
c         write(*,*)diag(i)
c        enddo
       if(NiterL.ge.Ncheck)then
        esave(NiterL-(Ncheck-1))=diag(1)
         if(NiterL.ge.(Ncheck+1))then
          diff=esave(NiterL-(Ncheck-1))-esave(NiterL-(Ncheck-1)-1)
c          write(*,*)'test deltaE = ',diff
          if(abs(diff).le.eThresholdLanczos)then
           goto 400
          endif
         endif
       endif

      enddo
 400  write(*,*)'----------------------------- '
      write(*,*)'      Lanczos converged       '
      write(*,*)'----------------------------- '

      if(NiterL.eq.NmaxiterL)then
       write(*,*)'too many iteration '
       stop
      endif

      do i=1,NiterL
         do j=1,NiterL
            Zaux(j,i)=0.d0
         enddo
         Zaux(i,i)=1.d0
      enddo

      do i=1,NiterL
         diag(i)=TA(i)
      end do
      do i=2,NiterL
         subdiag(i)=TB(i)
      end do

      call tql2(NmaxiterL,NiterL,diag,subdiag,Zaux,Nerr)
c      write(*,*)'normal return = 0'
c      write(*,*)Nerr
      if(Nerr.ne.0)then
       write(*,*)' TQL2 has problem !'
       stop
      endif
      write(*,*)    
c      print *,'Output from diagonalization:'
      print *,'------> eigenvalues <-------'
      do i=1,NiterL
       write(*,*)diag(i)
      enddo

      write(*,*)
      write(*,*)' recalling Lanczos '
      write(*,*)
      do i=1,NmaxNstatesHS
       Qin(i)=Qsave(i)
       Qlin(i)=cmplx(0.d0,0.d0)
      enddo

      TB(1)=0

      do i=1,NmaxNstatesHS
       ZA(i)=0.d0
      enddo

      do iterL = 1,NiterL  !<---------------

       call LANCZOSstep(iterL
     &                 ,haccaS
     &                 ,NmaxNstatesHS
     &                 ,TA1,TB1
     &                 ,Qin,Qlin)

c       write(*,*)'(re-)iteration Lanczos (Z) ',iterL

c       write(*,*)'TA1*',TA1
c       write(*,*)'TB1*',TB1
       TA(iterL)=TA1
       TB(iterL+1)=TB1

        do l=1,NmaxNstatesHS
         ZA(l)=ZA(l)+Qlin(l)*Zaux(iterL,1)   !<---------------------
        enddo

      enddo
c      write(*,*)'Lanczos converged ** '
c      write(*,*)'----------------------------- '
c      write(*,*)'----------------------------- '

c      do i=1,NiterL
c         do j=1,NiterL
c            Zaux(j,i)=0.d0
c         enddo
c         Zaux(i,i)=1.d0
c      enddo
c
c      do i=1,NiterL
c         diag(i)=TA(i)
c      end do
c      do i=2,NiterL
c         subdiag(i)=TB(i)
c      end do
c
c      call tql2(NmaxiterL,NiterL,diag,subdiag,Zaux,Nerr)
c      write(*,*)'normal return = 0'
c      write(*,*)Nerr
c      if(Nerr.ne.0)then
c       write(*,*)' TQL2 has problem !'
c       stop
c      endif
c      print *,'Output from diagonalization:'
c      print *,'------> (re-)eigenvalues <-------'
c      do i=1,NiterL
c       write(*,*)diag(i)
c      enddo

      size=0.d0
      do i=1,NmaxNstatesHS
       size=size+CONJG(ZA(i))*ZA(i)
      enddo
      do i=1,NmaxNstatesHS
       ZA(i)=ZA(i)/sqrt(size)   !eigenvector ground state with norm 1
      enddo

      return
      END


       SUBROUTINE LANCZOSstep(iter
     &                 ,haccaS,Nmax
     &                 ,TA1,TB1
     &                 ,Qin,Qlin)
       implicit none
       integer i,j
       integer iter,N
       integer Nmax
       complex*16 haccaS(Nmax,Nmax)
       complex*16 Qin(Nmax),Qlin(Nmax)
       complex*16 tmp(Nmax),Qh(Nmax)
       complex*16 sum
       double precision size,TA1,TB1

       N=Nmax

       if(iter.eq.1)then

        size=0.d0
        do i=1,N
         size=size+CONJG(Qin(i))*Qin(i)
        enddo
 
        do i=1,N
         Qin(i)=Qin(i)/sqrt(size)
        enddo

        size=0.d0
        do i=1,N
         size=size+CONJG(Qin(i))*Qin(i)
        enddo
c        write(*,*)'initial vector norm =',size

        TB1=0

       endif

       do i=1,N
        tmp(i)=cmplx(0.d0,0.d0)
        Qh(i)=cmplx(0.d0,0.d0)
       enddo

c
          DO I=1,N
          SUM=0
          DO J=1,N
            SUM=SUM+haccaS(I,J)*Qin(J)
          ENDDO
          tmp(I)=SUM
          ENDDO
c
         DO j=1,N
          tmp(j)=tmp(j)-TB1*Qlin(j)
         ENDDO

        size=0.d0
        do i=1,N
         size=size+CONJG(Qin(i))*tmp(i)
        enddo
         TA1=size
c         write(*,*)'TA1',TA1

         DO j=1,N
          Qh(j)=tmp(j)-TA1*Qin(j)
         ENDDO

        size=0.d0
        do i=1,N
         size=size+CONJG(Qh(i))*Qh(i)
        enddo
         TB1=size
         TB1=sqrt(TB1)
c         write(*,*)'TB1',TB1

         DO j=1,n
            Qlin(j)=Qin(j)
            Qin(j)=Qh(j)/TB1
         ENDDO


       return
       END

C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: tql2
C TYPE   : subroutine
C PURPOSE: Compute the eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix
C I/O    : 
C VERSION:
C COMMENT: Minimally modified SLATEC routine from netlib.
C          To learn about netlib, send an otherwise empty
C          message to netlib@research.att.com
C          containing 'send index' in the subject header)
C          The WWW address of netlib is
C          http://netlib.att.com/netlib/search.html
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
*DECK TQL2
      SUBROUTINE TQL2 (NM, N, D, E, Z, IERR)
C***BEGIN PROLOGUE  TQL2
C***PURPOSE  Compute the eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A5, D4C2A
C***TYPE      SINGLE PRECISION (TQL2-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TQL2,
C     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
C     Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
C     The eigenvectors of a FULL SYMMETRIC matrix can also
C     be found if  TRED2  has been used to reduce this
C     full matrix to tridiagonal form.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, Z, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        D contains the diagonal elements of the symmetric tridiagonal
C          matrix.  D is a one-dimensional DOUBLE PRECISION array, 
C          dimensioned D(N).
C
C        E contains the subdiagonal elements of the symmetric
C          tridiagonal matrix in its last N-1 positions.  E(1) is
C          arbitrary.  E is a one-dimensional DOUBLE PRECISION array, 
C          dimensioned E(N).
C
C        Z contains the transformation matrix produced in the
C          reduction by  TRED2, if performed.  If the eigenvectors
C          of the tridiagonal matrix are desired, Z must contain
C          the identity matrix.  Z is a two-dimensional DOUBLE PRECISION array,
C          dimensioned Z(NM,N).
C
C      On Output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct but
C          unordered for indices 1, 2, ..., IERR-1.
C
C        E has been destroyed.
C
C        Z contains orthonormal eigenvectors of the symmetric
C          tridiagonal (or full) matrix.  If an error exit is made,
C          Z contains the eigenvectors associated with the stored
C          eigenvalues.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  TQL2
C
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      DOUBLE PRECISION D(*),E(*),Z(NM,*)
      DOUBLE PRECISION B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      DOUBLE PRECISION PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = PYTHAG(P,1.0D0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0D0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0D0)
            E(I+1) = S * E(I) * R
            S = 1.0D0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: pythag 
C TYPE   : function
C PURPOSE: compute dsqrt(a**2+b**2) without overflow or
C          destructive underflow
C I/O    : 
C VERSION:
C COMMENT: 
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      double precision function pythag(a,b)
      double precision a,b
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end

      REAL FUNCTION RAND()
C
C  This function returns a pseudo-random number for each invocation.
C  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
C  standard number generator whose Pascal code appears in the article:
C
C     Park, Steven K. and Miller, Keith W., "Random Number Generators:
C     Good Ones are Hard to Find", Communications of the ACM,
C     October, 1988.
C
      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773,
     +           MOMDMP=2836)
C
      COMMON  /SEED/JSEED,IFRST
      INTEGER HVLUE, LVLUE, TESTV, NEXTN
      SAVE    NEXTN
C
      IF (IFRST .EQ. 0) THEN
        NEXTN = JSEED
        IFRST = 1
      ENDIF
C
      HVLUE = NEXTN / MOBYMP
      LVLUE = MOD(NEXTN, MOBYMP)
      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
      IF (TESTV .GT. 0) THEN
        NEXTN = TESTV
      ELSE
        NEXTN = TESTV + MODLUS
      ENDIF
      RAND = REAL(NEXTN)/REAL(MODLUS)
C
      RETURN
      END
      BLOCKDATA RANDBD
      COMMON /SEED/JSEED,IFRST
C
      DATA JSEED,IFRST/123456789,0/
C
      END
