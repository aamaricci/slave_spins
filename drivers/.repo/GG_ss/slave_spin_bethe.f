      PROGRAM SlaveSpins_multisites
      implicit none      
      include 'slave_spin_hm.h'  

      integer ix,iy
      real*8 kx,ky
      real*8 a,b
      

c     input:       
      open(2,file='INPUT.in')   ! INPUT file
      read(2,*)Uvalue           ! U value
      read(2,*)thop             ! hopping
      read(2,*)filling          ! filling 
      read(2,*)Threshold        ! Threshold for the self-consistent loops
      close(2)
      write(*,*)' Uvalue = ',Uvalue
      write(*,*)' '
      write(*,*)' t_hop = ',thop
      write(*,*)' '
      write(*,*)' filling = ',filling
      write(*,*)' '
      write(*,*)' Threshold = ',Threshold
      write(*,*)' '

      if_ldos=0

c!!!!!TO BE MODIFIED:
c     DFT input:
c     read num_kpts k-points:
c$$$  write(*,*)' k - mesh in *'
c$$$  open(2,file='kpoints.scf')
c$$$  do loop_kpt=1,num_kpts
c$$$  read(2,*)(kpt_latt(j,loop_kpt),j=1,ndim)
c$$$  enddo
c$$$  close(2)
c$$$  write(*,*)' k - mesh in **'
c$$$  
c$$$  do loop_kpt=1,num_kpts
c$$$  do i=1,num_wann
c$$$  do j=1,num_wann
c$$$  ham_k_in(i,j,loop_kpt)=dcmplx(0.d0,0.d0)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  
c$$$  open(2,file='hk.dat')
c$$$  do loop_kpt=1,num_kpts
c$$$  do i=1,num_wann
c$$$  do j=1,num_wann
c$$$  read(2,*)a,b
c$$$  ham_k_in(i,j,loop_kpt)=dcmplx(a,b)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  close(2)
c$$$  c!!!!!!!


      pi=acos(-1.d0)
      loop_kpt=0
      do ix=1,Nkx
         kx = -pi + 2.d0*pi*dble(ix-1)/dble(Nkx)
         do iy=1,Nky
            ky = -pi + 2.d0*pi*dble(iy-1)/dble(Nky)

            loop_kpt=loop_kpt+1

            do i=1,num_wann
               ham_k_in(i,i,loop_kpt)=               
     &              -2.d0*thop*(cos(kx)+cos(ky))
            enddo
         enddo
      enddo

      Npara=num_wann

c     initial guess:
c     read initial Lagrange multipliers lambda
c     read initial <O>
      do iorb=1,num_wann
         lambda(iorb)=0.0d0 
         mOi_orb(iorb)=cmplx(0.d0,0.d0)
      enddo

      write(*,*)''
      write(*,*)' l and <O> parameters '
      open(2,file='parameters')
      do iorb=1,num_wann
         read(2,*)lambda(iorb)
      enddo
      do iorb=1,num_wann
         read(2,*)mOi_orb(iorb)
      enddo
      close(2)

c     write all the initial parameters:
      write(*,*)''
      write(*,*)' initial l parameters '
      do iorb=1,num_wann
         write(*,*)lambda(iorb)
      enddo
      write(*,*)''
      write(*,*)' initial <O> parameters '
      do iorb=1,num_wann
         write(*,*)mOi_orb(iorb)
      enddo
      write(*,*)''


c     Initialize spin down from spin up
      do j=1,num_wann
         mOi_orb(j+num_wann)=mOi_orb(j)
      enddo
      do j=1,nspin*num_wann
         mOpi_orb(j)=mOi_orb(j)
      enddo
!     <O>
      do j=1,num_wann
         mOi_orb_ini(j)=mOi_orb(j)
      enddo
      do j=1,num_wann
         mOi_orb_ini(j+num_wann)=mOi_orb(j)
      enddo
!     <O+>
      do j=1,num_wann
         mOpi_orb_ini(j)=mOpi_orb(j)
      enddo
      do j=1,num_wann
         mOpi_orb_ini(j+num_wann)=mOpi_orb(j)
      enddo


c     broyden subroutine (based on Broyden`s method) changes the lambda_i to get the constrain n=S+1/2 satisfied.
      write(*,*)' parameters passed to the Broyden '
      do iorb=1,num_wann
         x(iorb)=lambda(iorb)
         write(*,*)iorb,x(iorb)
      enddo
      write(*,*)

      ciclo=0
      write(*,*)''
      write(*,*)' call broydn '
      write(*,*)''
      call broydn(x,Npara,check)
      write(*,*)''
      write(*,*)' out  broydn '
      write(*,*)''
      write(*,*)' spin-slaves equations solved'
      write(*,*)''

      if (check) then
         write(*,*)' Convergence problems !!!!'
      endif

c     output 
      write(*,*)' parameters passed from the Broyden '
      do iorb=1,num_wann
         lambda(iorb)=x(iorb)
      enddo
      write(*,*)

      write(*,*)''
      write(*,*)' final l parameters '
      do iorb=1,num_wann
         write(*,*)lambda(iorb)
      enddo
      write(*,*)''

      open(2,file='parameters')
      do iorb=1,num_wann
         write(2,*)lambda(iorb)
      enddo
      do iorb=1,num_wann
         write(2,*)mOi_orb(iorb)
      enddo
      close(2)

c     output code:
      open(2,file='Z_n_U_Jh')
      if(num_wann.eq.1)then
         write(2,'(1f10.5,1f10.5,1f10.5,1f10.5)')
     &        Uvalue,(real(mOi_orb(iorb))**2,iorb=1,num_wann)
     &        ,(n_state(iorb),iorb=1,num_wann)
      endif
      if(num_wann.eq.2)then
         write(2,'(1f10.5,1f10.5,1f10.5,1f10.5,1f10.5,1f10.5)')
     &        Uvalue,(real(mOi_orb(iorb))**2,iorb=1,num_wann)
     &        ,(n_state(iorb),iorb=1,num_wann)
      endif
      close(2)
c     

c$$$      if_Ldos=1
c$$$      call fermionic_hamiltonian(x)
      

      stop
      end




c##################################################################
c##################################################################
c##################################################################
c##################################################################
c##################################################################










      SUBROUTINE funcv(Npara,x,f)
      implicit none
      include 'slave_spin_hm.h'
      integer iter,jter
      integer infoconv1
      integer Nmaxiter
c     integer Nsite,isite,indx_orb_site
      parameter (Nmaxiter=50)
      complex*16 f_aux1(Nmaxiter,nspin*num_wann)  
      complex*16 mOi_orb_aux(nspin*num_wann)
      complex*16 mOpi_orb_aux(nspin*num_wann)
      complex*16 mSz_orb_aux(nspin*num_wann)

      ciclo=ciclo+1
      write(*,*)'***************'
      write(*,*)'**** Cycle ****',ciclo 
      write(*,*)'***************'

      do i=1,nspin*num_wann
         do iter=1,Nmaxiter
            f_aux1(iter,i)=0.d0  
         enddo
      enddo

c     same starting guess at all the cicles
c     up
      do iorb=1,num_wann
         mOi_orb(iorb)=mOi_orb_ini(iorb)
         mOpi_orb(iorb)=mOpi_orb_ini(iorb)
      enddo
c     dn
      do iorb=1,num_wann
         mOi_orb(iorb+num_wann)=mOi_orb_ini(iorb+num_wann)
         mOpi_orb(iorb+num_wann)=mOpi_orb_ini(iorb+num_wann)
      enddo



c     set density to zero
      do i=1,num_wann
         n_state(i)=0.d0      
      enddo
c     set Sz to zero
      do i=1,nspin*num_wann
         mSz_orb(i)=0.d0          
      enddo

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
         do i=1,num_wann
            write(*,*)x(i)
         enddo


c     solve the fermionic part: out -> n_state
         write(*,*)' '
         write(*,*)' call H F'
         write(*,*)' '
         call fermionic_hamiltonian(x) 

c     solve the spin part: out -> mSz_orb
         write(*,*)' '
         write(*,*)' call H S'
         write(*,*)' '
         call spin_hamiltonian(x)
         write(*,*)''
         write(*,*)' solution of delta-functions = 0 (NM, SP cases) '
         write(*,*)' <n>=Sz+1/2 '
         write(*,*)''

c     up
         write(*,*)' spin = 1 '
         do i=1,num_wann
            write(*,*)n_state(i),(mSz_orb(i)+0.5d0)
            write(*,*)n_state(i)-(mSz_orb(i)+0.5d0)
            write(*,*)
         enddo
         write(*,*)''
c     dn
         write(*,*)' spin = 2 '
         do i=1,num_wann
            write(*,*)n_state(i),(mSz_orb(i+num_wann)+0.5d0)
            write(*,*)n_state(i)-(mSz_orb(i+num_wann)+0.5d0)
            write(*,*)
         enddo
         write(*,*)''


         indx=0
         do j=1,nspin*num_wann
            indx=indx+1 
            f_aux1(iter,indx)=mOi_orb(j)             
         enddo

         write(*,*)' function for the convergence : '
         do j=1,nspin*num_wann
            write(*,*)f_aux1(iter,j)
         enddo
         write(*,*)

         write(*,*)'<O^+_i><O_i>'
         do j=1,nspin*num_wann
            write(*,*)mOi_orb(j)*mOpi_orb(j)
         enddo
         write(*,*)

         write(*,*)' <O_i> '
         write(*,*)''
         do i=1,nspin*num_wann
            write(*,*)mOi_orb(i)
         enddo
         write(*,*)''

         write(*,*)' <O^+_i> '
         write(*,*)''
         do i=1,nspin*num_wann
            write(*,*)mOpi_orb(i)
         enddo
         write(*,*)''


         if(iter.gt.1)then
            infoconv1=0
            do  i=1,nspin*num_wann
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

            if(infoconv1.eq.nspin*num_wann)then
               goto 1
            endif
         endif
      enddo                     !end the self-consistency loop between fermionic and spin.

 1    write(*,*)''
      write(*,*)' Self-Consistency between Hs and Hf ends '
      write(*,*)''

      write(*,*)''
      write(*,*)' Self-Consistent Equations '
      write(*,*)''
      do i=1,num_wann
         f(i)=n_state(i)-(mSz_orb(i)+0.5d0)     
         write(*,*)'functions for the broyden ',f(i)
      enddo
      write(*,*)

      do i=1,num_wann
         write(*,*)' occupancies in funcv ',n_state(i)
         write(*,*)' <Sz>+0.5 ',mSz_orb(i)+0.5
         write(*,*)' deviation ',n_state(i)-(mSz_orb(i)+0.5)
      enddo
      write(*,*)
c     
      return
      end





c##################################################################
c##################################################################
c##################################################################
c##################################################################
c##################################################################




      SUBROUTINE fermionic_hamiltonian(x)
      implicit none
      include 'slave_spin_hm.h'
      
      do j=1,num_wann
         lambdawrk(j)=x(j)
      enddo

c     
      do loop_kpt=1,num_kpts
         do i=1,num_wann
            do j=1,num_wann
               ham_k(i,j,loop_kpt)=
     &              ham_k_in(i,j,loop_kpt)*mOi_orb(i)*mOpi_orb(j)
            enddo
            ham_k(i,i,loop_kpt)=ham_k(i,i,loop_kpt)-
     &           dcmplx(lambdawrk(i),0.d0)
         enddo
      enddo
      

      write(*,*)
      write(*,*)' Diagonalization Fermionic problem '
      write(*,*)
c     Diagonalization Fermionic problem
      do loop_kpt=1,num_kpts 
         do i=1,num_wann
            do j=1,num_wann
               hdiag(i,j)=ham_k(i,j,loop_kpt)
            enddo
         enddo     

         call zheevx('V'
     &        ,'A'
     &        ,'U'
     &        ,num_wann
     &        ,hdiag
     &        ,num_wann
     &        ,vl,vu
     &        ,1
     &        ,num_wann
     &        ,0.0
     &        ,neigval,eigval
     &        ,Z,num_wann
     &        ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

         if(info.ne.0)then
            write(*,*)' Failure in ZHEEVX. INFO =',info 
            stop
         endif
         
         do l=1,num_wann
            eigenvalues(num_wann*(loop_kpt-1)+l)=eigval(l)
            eigenvalues_k(loop_kpt,l)=eigval(l)
         enddo
         
         do ll=1,num_wann
            do ii=1,num_wann
               Z_k(loop_kpt,ii,ll)=Z(ii,ll)
            enddo
         enddo
      enddo 
      write(*,*)
      write(*,*)' End Diagonalization Fermionic problem '
      write(*,*)

      

      write(*,*)' '
      Ntot=num_wann*num_kpts
      write(*,*)' # states ',Ntot
      N_e_k=ifix(filling*real(Ntot)+0.5)
      write(*,*)' # states filled ',N_e_k
      
      call indexx(Ntot,eigenvalues,e_ind)
      
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
                  n_state(ii)=n_state(ii)+abs(Z_k(loop_kpt,ii,ll))**2.d0
               endif
            enddo
         enddo
      enddo

      write(*,*)''
      do ll=1,num_wann
         write(*,*)' Occupancies : ',ll,n_state(ll)/real(num_kpts)
         n_state(ll)=n_state(ll)/real(num_kpts)
      enddo

      
      write(*,*)''
      sumaux=0.d0
      do ll=1,num_wann
         sumaux=sumaux+n_state(ll)
      enddo
      write(*,*)' Sum Occupancies : ',sumaux
c     

      write(*,*)''
      write(*,*)' W ',abs(eigenvalues(e_ind(1))
     &     -eigenvalues(e_ind(Ntot)))
     &     ,eigenvalues(e_ind(1))
     &     ,eigenvalues(e_ind(Ntot))
      write(*,*)


      if(if_ldos.eq.1)then
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
     &                    (abs(Z_k(loop_kpt,ii,ll))**2)
     &                    *gamma/pi/
     &                    (gamma**2
     &                    +((egrid(ie))
     &                    -(eigenvalues_k(loop_kpt,ll)))**2)
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

         stop

      endif





c     calculations of < f^+ f >
      write(*,*)
      write(*,*)' < f^+ f > '
      write(*,*)
      fnorm=1.d0/dfloat(num_kpts)
      do i=1,num_wann
         do j=1,num_wann 
            fdagf(i,j)=0.d0
            do loop_kpt=1,num_kpts
               sumaux=cmplx(0.d0,0.d0)
               do l=1,num_wann
                  if(eigenvalues_k(loop_kpt,l).le.EF)then
                     sumaux = sumaux + 
     1                    dconjg(Z_k(loop_kpt,i,l))*
     2                    Z_k(loop_kpt,j,l)
                  endif
               enddo
               fdagf(i,j)=fdagf(i,j)+ham_k_in(i,j,loop_kpt)*sumaux*fnorm
            enddo

         enddo
      enddo

      enew2=0.d0
      do i=1,Num_wann
         do j=1,num_wann
            enew2=enew2+fdagf(i,j)
         enddo
      enddo
      write(*,*)""
      write(*,*)"EneK 2 types:",eneW/real(num_kpts),enew2
     &     ,abs(eneW/real(num_kpts)-enew2)
      write(*,*)""
      write(*,*)' < f^+ f > calculated '
      write(*,*)

      
      do i=1,num_wann
         do j=i,num_wann
            Jijsigma(i,j,1)=fdagf(i,j)*thop
            Jijsigma(i,j,2)=fdagf(i,j)*thop
         enddo
      enddo
      return
      end






c##################################################################
c##################################################################
c##################################################################
c##################################################################
c##################################################################





      SUBROUTINE spin_hamiltonian(x)
      implicit none
      include 'slave_spin_hm.h'
      
      do iorb=1,num_wann        !N_site_tot
         lambda(iorb)=x(iorb)
      enddo


      Uprime=Uvalue             !-2*Jh
      print*,"Uprime=",Uprime

!     Setup the Hilbert Space for the Spin problem
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
      write(*,*)'',J_up,J_dn
      Ns_up=J_up
      Ns_dn=J_dn
      write(*,*)' Hilbert Space dimensions '
     &     ,Ns_up,Ns_dn,2**(2*N_site_tot) ! 2**(2*N_site_tot) = # states in Hilbert space 

      Nstate=2**(2*N_site_tot)



      write(*,*)''
      write(*,*)' <i|H|j>'
      do i=1,Nstate
         do j=1,Nstate
            haccaS(i,j)=cmplx(0.d0,0.d0)
         enddo
      enddo


!     DIAGONAL TERMS IN H:
      write(*,*)''
      write(*,*)' diagonal terms '
c     lambda terms:
c     for such term see Eq. 14 http://arxiv.org/pdf/cond-mat/0503764

      do i=1,Ns_dn     
         do j=1,Ns_up    

            sumtot=0.d0

            sumaux1=0.d0
            do jorb=1,num_wann  !N_site_tot 
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
     &           +cmplx(sumtot,0.d0)
         enddo
      enddo



c     first term Eq. 18 of the paper "Orbital selective Mott Transition in multi-band systems..." 
c     http://arxiv.org/pdf/cond-mat/0503764
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
     &           +cmplx(0.5d0*Uprime*sum1,0.d0)

         enddo
      enddo




!     Phases/gauges
      write(*,*)' '
      write(*,*)' c_m,sigma in Hs (spin polarized version) '
      do i=1,nspin
         do j=1,num_wann        !N_site_tot
            cv=1.d0/sqrt(n_state(j)*(1-n_state(j)))-1.d0
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
            iii=Ns_up*(i-1)+j   !full HS state index IN/BRA
            
            do iorb=1,num_wann
               
               bool1 = btest(I_up(j)-1,iorb-1) 
               
               
!     c     S^+_i                                  
               if(bool1.EQV..FALSE.)then       
                  Ip_up=IBSET(I_up(j)-1,iorb-1)+1 
                  
                  do ii=1,Ns_dn                   
                     do jj=1,Ns_up                  
                        
                        if(Ip_up.eq.I_up(jj))then      
                           if(I_dn(i).eq.I_dn(ii))then    
                              
                              jjj=Ns_up*(ii-1)+jj !full HS state index OUT/KET
                              
                              
                              do jorb=1,num_wann
                                 
                                 if(iorb.le.jorb)then !both therms are calculated here <O+>O + <O>O+
                                    factorHM=Jijsigma(iorb,jorb,1)
     &                                   *mOi_orb(jorb)
                                    if(jjj.gt.iii)then                                          
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +factorHM
                                    else
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +CONJG(factorHM)
                                    endif
                                    
                                 else                      
                                    
                                    factorHM=Jijsigma(jorb,iorb,1)
     &                                   *mOi_orb(jorb) 
                                    if(jjj.gt.iii)then                                       
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +factorHM
                                    else
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +CONJG(factorHM)
                                    endif
                                 endif
                                 
                              enddo
                              
                           endif
                        endif
                        
                     enddo
                  enddo
                  
               endif
               
               
               
               
               
!     c     S^-_i c^*_i,up                         
               if(bool1.EQV..TRUE.)then        
                  
                  Ip_up=IBCLR(I_up(j)-1,iorb-1)+1
                  
                  do ii=1,Ns_dn
                     do jj=1,Ns_up
                        
                        if(Ip_up.eq.I_up(jj))then
                           if(I_dn(i).eq.I_dn(ii))then
                              
                              jjj=Ns_up*(ii-1)+jj
                              
                              
                              do jorb=1,num_wann
                                 if(iorb.le.jorb)then
                                    factorHM=Jijsigma(iorb,jorb,1)
     &                                   *CONJG(cm(iorb,1))
     &                                   *mOi_orb(jorb)
                                    if(jjj.gt.iii)then
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +factorHM
                                    else
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +CONJG(factorHM)
                                    endif
                                 else
                                    factorHM=Jijsigma(jorb,iorb,1)
     &                                   *CONJG(cm(iorb,1))
     &                                   *mOi_orb(jorb)
                                    if(jjj.gt.iii)then
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +factorHM
                                    else
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +CONJG(factorHM)
                                       
                                    endif
                                 endif
                              enddo
                              
                           endif
                        endif
                        
                     enddo
                  enddo
                  
               endif

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
               
               bool1 = btest(I_dn(i)-1,iorb-1)
!     c     S^+_i
               if(bool1.EQV..FALSE.)then
                  
                  Ip_dn=IBSET(I_dn(i)-1,iorb-1)+1
                  
                  do ii=1,Ns_dn
                     do jj=1,Ns_up
                        
                        if(Ip_dn.eq.I_dn(ii))then
                           if(I_up(j).eq.I_up(jj))then
                              
                              jjj=Ns_up*(ii-1)+jj
                              
                              do jorb=1,num_wann
                                 
                                 if(iorb.le.jorb)then
                                    factorHM=Jijsigma(iorb,jorb,2)
     &                                   *mOi_orb(jorb+num_wann)
                                    if(jjj.gt.iii)then
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +factorHM
                                    else
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +CONJG(factorHM)
                                    endif
                                 else
                                    factorHM=Jijsigma(jorb,iorb,2)
     &                                   *mOi_orb(jorb+num_wann)
                                    if(jjj.gt.iii)then
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +factorHM                                       
                                    else                                      
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +CONJG(factorHM)
                                    endif                                    
                                 endif
                                 
                              enddo
                              
                           endif
                        endif
                        
                     enddo
                  enddo
                  
               endif
               
               
!     c     S^-_i c^*_i,dn
               if(bool1.EQV..TRUE.)then
                  
                  Ip_dn=IBCLR(I_dn(i)-1,iorb-1)+1
                  
                  do ii=1,Ns_dn
                     do jj=1,Ns_up
                        
                        if(Ip_dn.eq.I_dn(ii))then
                           if(I_up(j).eq.I_up(jj))then
                              
                              jjj=Ns_up*(ii-1)+jj
                              
                              do jorb=1,num_wann
                                 if(iorb.le.jorb)then
                                    factorHM=Jijsigma(iorb,jorb,2)
     &                                   *CONJG(cm(iorb,2))
     &                                   *mOi_orb(jorb+num_wann)
                                    if(jjj.gt.iii)then
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +factorHM
                                    else
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +CONJG(factorHM)
                                    endif
                                 else
                                    factorHM=Jijsigma(jorb,iorb,2)
     &                                   *CONJG(cm(iorb,2))
     &                                   *mOi_orb(jorb+num_wann)
                                    if(jjj.gt.iii)then
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +factorHM
                                    else
                                       haccaS(min(iii,jjj),
     &                                      max(iii,jjj))=
     &                                      haccaS(min(iii,jjj),
     &                                      max(iii,jjj))
     &                                      +CONJG(factorHM)
                                    endif                                    
                                 endif
                                 
                              enddo
                              
                           endif
                        endif
                        
                     enddo
                  enddo
                  
               endif
               
            enddo
            
         enddo
      enddo
      
      
      do i=1,Nstate
         write(*,'(100F12.6)')(haccaS(i,j),j=1,Nstate)
      enddo


c     diagonalization:
      write(*,*)''
      write(*,*)' Diagonalization SS problem'
      write(*,*)

      write(*,*)
      write(*,*)' LAPACK ZHEEVX'
      write(*,*)
      call zheevx('V'
     &     ,'A'
     &     ,'U'
     &     ,Nstate
     &     ,haccaS
     &     ,Nstate
     &     ,ivl,ivu
     &     ,1
     &     ,Nstate
     &     ,0.d0
     &     ,neig,eig
     &     ,ZZ,Nstate
     &     ,WWORK,LLWORK,RRWORK,IIWORK,IIFAIL,IINFO)


      if(iinfo.ne.0)then
         write(*,*)' Failure in ZHEEVX. INFO =',iinfo
         stop
      endif

      write(*,*)''
      write(*,*)' End Diagonalization SS problem'
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

      
      do i=1,Ns_dn
         do j=1,Ns_up
            ii = Ns_up*(i-1)+j
            ZZ_gs(ii)=ZZ(ii,1)
         enddo
      enddo








c     <Sz>
      write(*,*)''
      write(*,*)' calculation <Sz>'
      write(*,*)''
      do iorb=1,num_wann        !N_site_tot 
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
c     up
      indx=0
      do iorb=1,N_site_tot
         indx=indx+1
         mSz_orb(indx)=real(Sz(1,iorb)) 
      enddo
c     dn
      do iorb=1,N_site_tot
         indx=indx+1
         mSz_orb(indx)=real(Sz(2,iorb))
      enddo







c     <O_i>
      write(*,*)''
      write(*,*)' calculation <O> '
      write(*,*)''
      do iorb=1,num_wann        !N_site_tot 
         Oi(iorb,1)=cmplx(0.d0,0.d0)
         Oi(iorb,2)=cmplx(0.d0,0.d0)
         sum1aux=cmplx(0.d0,0.d0) 
         do j=1,Ns_up        
            do jj=1,Ns_up      
               bool1 = btest(I_up(j)-1,iorb-1)
c     S^-_i
               if(bool1.EQV..TRUE.)then         
                  Ip_up=IBCLR(I_up(j)-1,iorb-1)+1 
                  if(Ip_up.eq.I_up(jj))then
                     do k=1,Ns_dn          
                        iii=Ns_up*(k-1)+j    
                        jjj=Ns_up*(k-1)+jj   
                        sum1aux=sum1aux+
     &                       +CONJG(ZZ_gs(jjj))*ZZ_gs(iii)
                     enddo
                  endif
               endif
c     S^+_i c_i,up
               if(bool1.EQV..FALSE.)then       
                  Ip_up=IBSET(I_up(j)-1,iorb-1)+1 
                  if(Ip_up.eq.I_up(jj))then
                     do k=1,Ns_dn                 
                        iii=Ns_up*(k-1)+j          
                        jjj=Ns_up*(k-1)+jj        
                        sum1aux=sum1aux+
     &                       +cm(iorb,1)
     &                       *CONJG(ZZ_gs(jjj))*ZZ_gs(iii)
                     enddo
                  endif
               endif
            enddo
         enddo

         Oi(iorb,1)=sum1aux
         sum2aux=cmplx(0.d0,0.d0)
         do i=1,Ns_dn      
            do ii=1,Ns_dn    
c     S^-_i
               bool1 = btest(I_dn(i)-1,iorb-1)
               if(bool1.EQV..TRUE.)then   
                  Ip_dn=IBCLR(I_dn(i)-1,iorb-1)+1 
                  if(Ip_dn.eq.I_dn(ii))then
                     do k=1,Ns_up            
                        iii=Ns_up*(i-1)+k      
                        jjj=Ns_up*(ii-1)+k     
                        sum2aux=sum2aux+
     &                       +CONJG(ZZ_gs(jjj))*ZZ_gs(iii)
                     enddo
                  endif
               endif
c     S^+_i c_i,dn
               if(bool1.EQV..FALSE.)then   
                  Ip_dn=IBSET(I_dn(i)-1,iorb-1)+1  
                  if(Ip_dn.eq.I_dn(ii))then
                     do k=1,Ns_up             
                        iii=Ns_up*(i-1)+k       
                        jjj=Ns_up*(ii-1)+k      
                        sum2aux=sum2aux+
     &                       +cm(iorb,2)
     &                       *CONJG(ZZ_gs(jjj))*ZZ_gs(iii)
                     enddo
                  endif
               endif
            enddo
         enddo
         Oi(iorb,2)=sum2aux
      enddo
      


!     Store the variables to be passed to the calling routine
      write(*,*)' '
      write(*,*)' mean values <O_i> '
      do iorb=1,num_wann        !N_site_tot
         write(*,*)Oi(iorb,1),Oi(iorb,2)
      enddo
      write(*,*)' '
      indx=0
      do i=1,nspin
         do j=1,num_wann        !N_site_tot
            indx=indx+1
            mOi_orb(indx)=Oi(j,i) 
         enddo
      enddo
      write(*,*)' '





c     <O^+_i>
      do iorb=1,num_wann        !N_site_tot  
         Opi(iorb,1)=cmplx(0.d0,0.d0)
         Opi(iorb,2)=cmplx(0.d0,0.d0)
         sum1aux=cmplx(0.d0,0.d0)
         do j=1,Ns_up         
            do jj=1,Ns_up       
               bool1 = btest(I_up(j)-1,iorb-1)
c     S^+_i
               if(bool1.EQV..FALSE.)then         
                  Ip_up=IBSET(I_up(j)-1,iorb-1)+1  
                  if(Ip_up.eq.I_up(jj))then
                     do k=1,Ns_dn        
                        iii=Ns_up*(k-1)+j   
                        jjj=Ns_up*(k-1)+jj  
                        sum1aux=sum1aux+
     &                       +CONJG(ZZ_gs(jjj))*ZZ_gs(iii)
                     enddo
                  endif
               endif
c     S^-_i c^*_i,up
               if(bool1.EQV..TRUE.)then          
                  Ip_up=IBCLR(I_up(j)-1,iorb-1)+1  
                  if(Ip_up.eq.I_up(jj))then
                     do k=1,Ns_dn          
                        iii=Ns_up*(k-1)+j    
                        jjj=Ns_up*(k-1)+jj   
                        sum1aux=sum1aux+
     &                       +CONJG(cm(iorb,1))
     &                       *CONJG(ZZ_gs(jjj))*ZZ_gs(iii)
                     enddo
                  endif
               endif
            enddo
         enddo

         Opi(iorb,1)=sum1aux
         sum2aux=cmplx(0.d0,0.d0)
         do i=1,Ns_dn             
            do ii=1,Ns_dn           
c     S^+_i
               bool1 = btest(I_dn(i)-1,iorb-1) 
               if(bool1.EQV..FALSE.)then        
                  Ip_dn=IBSET(I_dn(i)-1,iorb-1)+1 
                  if(Ip_dn.eq.I_dn(ii))then     
                     do k=1,Ns_up           
                        iii=Ns_up*(i-1)+k     
                        jjj=Ns_up*(ii-1)+k    
                        sum2aux=sum2aux+
     &                       +CONJG(ZZ_gs(jjj))*ZZ_gs(iii)
                     enddo
                  endif
               endif
c     S^-_i c^*_i,dn
               if(bool1.EQV..TRUE.)then         
                  Ip_dn=IBCLR(I_dn(i)-1,iorb-1)+1 
                  if(Ip_dn.eq.I_dn(ii))then       
                     do k=1,Ns_up                  
                        iii=Ns_up*(i-1)+k            
                        jjj=Ns_up*(ii-1)+k           
                        sum2aux=sum2aux+
     &                       +CONJG(cm(iorb,2))
     &                       *CONJG(ZZ_gs(jjj))*ZZ_gs(iii)
                     enddo
                  endif
               endif
            enddo
         enddo
         Opi(iorb,2)=sum2aux
      enddo
      

!     Store the variables to be passed to the calling routine
      write(*,*)' '
      write(*,*)' mean values <O^+_i> '
      do iorb=1,num_wann        !N_site_tot
         write(*,*)Opi(iorb,1),Opi(iorb,2)
      enddo
      indx=0
      do i=1,nspin
         do j=1,num_wann        !N_site_tot
            indx=indx+1
            mOpi_orb(indx)=Opi(j,i)     
         enddo
      enddo
      write(*,*)' '

      return
      end




