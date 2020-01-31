      program SSpin_grad_magn

      include 'param.dat'
      logical check
      character*80 xyz 
      real*8 fvec(totpar),x(totpar),xsave(totpar)
      complex*16 SS(M,M,Ns),Sz(M,Ns),S(M,Ns)
      real*8 ddiag(totpar)
      real*8 fjac(totpar,totpar),r((totpar*(totpar+1))/2)
      real*8 qtf(totpar)
      real*8 wa1(totpar),wa2(totpar),wa3(totpar),wa4(totpar)

      external funcv

ccccccccccc reads the input parameters 
      open(2,file="input",form="formatted",status="unknown")

      do L=1,M
         do L2=1,M
           tmm(L,L2)=0.d0
         enddo
      enddo
      
      read(2,'(a80)') xyz !! interaction strength 
      read(2,*) Umin, Umax, Ustep, ratioJU
      read(2,'(a80)') xyz !! temperature, total density
      read(2,*) Beta, ntot
      read(2,'(a80)') xyz !! the bare levels
      read(2,*) (epsm(L),L=1,M)
      read(2,'(a80)') xyz
      do L=1,M  !! (only the lower triangle of) the hopping matrix
         read(2,*)(tmm(L,L2),L2=1,L)
      enddo

      do L=1,M
         do L2=1,L-1
            tmm(L2,L)=tmm(L,L2)
         enddo
      enddo

cccccccccccc  reads the seed for self-consistent parameters
      open(3,file="selfparam",form="formatted",status="unknown")

c     read the guess, if present
      do is=1,Ns
         do L=1,M
            read(3,*,err=122,end=122)(hm(L,L2,is),L2=1,M)
         enddo
      enddo
      do L=1,M
         read(3,*,err=122,end=122)(lam(L,is),is=1,Ns)
      enddo
      read(3,*,err=122,end=122) mu
      goto 124

 122  mu=0.d0
      do L=1,M
         do is=1,Ns
            do L2=1,M
               hm(L,L2,is)=dcmplx(-0.4244*tmm(L,L2),0.d0)
               if (L.eq.L2) hm(L,L2,is)=dcmplx(0.4244*tmm(L,L2),0.d0)
            enddo
         lam(L,is)=0.d0
         mu=mu+ntot*epsm(L)/(2*L)
         enddo
      enddo
 124  continue
      do L=1,M
         do is=1,Ns
         cm(L,is)=dcmplx(1/sqrt(ntot/(2*M)-(ntot/(2*M))**2)-1,0.d0)
         enddo
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      print*,"________________________________________________________"
      print*,"___________________PARAMETERS___________________________"
      print*,"M,Ns,"
      print*,M,Ns
      print*,"________________________________________________________"
      print*,"Umin, Umax, Ustep, ratioJU" 
      print*, Umin, Umax, Ustep, ratioJU
      print*,"Beta, ntot"
      print*, Beta, ntot
      print*,"epsm(1...M)"
      print*, (epsm(L),L=1,M)
      print*,"hopping matrix"
      do L=1,M
         print*,(tmm(L,L2),L2=1,M)
      enddo
      print*," "
      print*,"___________________GUESS___________________________"
      print*,"________________________________________________________"
      print*,"hm(1....M,1....M), spin up"
      do L=1,M
         print*,(hm(L,L2,1),L2=1,M)
      enddo
      print*,"hm(1....M,1....M), spin down"
      do L=1,M
         print*,(hm(L,L2,Ns),L2=1,M)
      enddo
      print*,"lam(1...M) spin up"
      print*, (lam(L,1),L=1,M)
      print*,"lam(1...M) spin down"
      print*, (lam(L,Ns),L=1,M)
      print*,"mu"
      print*,mu
      print*,"cm(1,..M) spin up"
      print*,(cm(L,1),L=1,M)
      print*,"cm(1,..M) spin down"
      print*,(cm(L,Ns),L=1,M)
      print*,"________________________________________________________"


C***********************************************************************
C******* we solve first the non-interacting case, to fix both the phase
C******* ph(L,is) of the gauge cm, and the shift of the bare levels
C***********************************************************************
      U=0.d0
      init=.true.
      print*,'(U=',U,')'
      print111,' '
c     here we assign the selfconsistent parameters to be adjusted
      icount=0
      do is=1,Ns
         do L1=1,M
            do L2=1,L1
               icount=icount+1
               x(icount)=dreal(hm(L1,L2,is))
               ddiag(icount)=min(1000.d0,1/tmm(L1,L2))
            enddo
         enddo
      enddo
      do is=1,Ns
         do L=1,M
            icount=icount+1
            x(icount)=lam(L,is)
            ddiag(icount)=1.d0
         enddo
      enddo
      icount=icount+1
      x(icount)=mu
      ddiag(icount)=1.d0

c      do i=1,icount
c         print*,x(i)
c      enddo
c      print*,""

      if (icount.ne.totpar) stop "icount differs from totpar"
      n=0
c      call newt(x,totpar,check)
c      call broydn(x,totpar,check)
cccccccc parameters for hybrd (the rest is in param.dat) ccccccccccccc
      ml=totpar-1
      mmu=totpar-1
      mode=1
      nprint=-1
      ldfjac=totpar
      lr=(totpar*(totpar+1))/2
      epsfcn2=epsfcn
 234  call hybrd(funcv,totpar,x,fvec,xtol,maxfev,ml,mmu,epsfcn2,ddiag,
     *                 mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,
     *                 qtf,wa1,wa2,wa3,wa4)
c     here we reconvert the variational parameters
      icount=0
      do is=1,Ns
         do L1=1,M
            do L2=1,L1
               icount=icount+1
               hm(L1,L2,is)=dcmplx(x(icount),0.d0)
               if (L2.lt.L1) hm(L2,L1,is)=dconjg(hm(L1,L2,is))
            enddo
         enddo
      enddo
      do is=1,Ns
         do L=1,M
            icount=icount+1
            lam(L,is)=x(icount)
         enddo
      enddo
      icount=icount+1
      mu=x(icount)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(88,*) U,nfev,info
      call flush (88)
      if (info.eq.5) then
c         print*, char(7)      
        epsfcn2=epsfcn2*10.d0
         write(88,*) "shift epsfcn to",epsfcn2
         call flush (88)
         goto 234
      endif
      write(88,*)(fvec(i),i=1,icount)

c      write(98,*)U,n,mu
c      call flush (98)
c      if (check) then
c         print*,"check=",check
c         stop
c      endif

      print*
      print*,"physical eps", (epsm(L),L=1,M)
      do L=1,M
c         print*,epsm(L),lam(L,1)
c         epsm(L)=epsm(L)-lam(L,1)
         do is=1,Ns
            ph(L,is)=dcmplx(1.d0,0.d0) !!!!to be changed
         enddo
c        print*,epsm(L),lam(L,1)
      enddo
      print*,"  scaled eps", (epsm(L),L=1,M)

C*********** loop on U ************************************************
      init=.false.
      nU=int((Umax-Umin)/Ustep)
      if (nU.lt.0) then
         nU=abs(nU)
         Ustep=-abs(Ustep)
      endif
      do IU=0,nU
         U=Umin+IU*Ustep
         J=U*ratioJU
         print*,' '
         print*,'---------------------------------------',
     &'----------------------------'
         print*,'***************************************',
     &'****************************'
         print*,'---------------------------------------',
     &'----------------------------'
         print*,' '
         do L=1,M
            rewind(30+L)
         enddo

C***********************************************************************
         print*,'(U=',U,')'
         print111,' '
c     here we assign the selfconsistent parameters to be adjusted
         icount=0
         do is1=1,Ns
            do L1=1,M
               do L2=1,L1
                  icount=icount+1
                  x(icount)=dreal(hm(L1,L2,is1))
               enddo
            enddo
         enddo
         do is=1,Ns
            do L=1,M
               icount=icount+1
               x(icount)=lam(L,is)
            enddo
         enddo
         icount=icount+1
         x(icount)=mu

         do i=1,icount
            xsave(i)=x(i)
         enddo
         epsfcn2=epsfcn
         write(88,*)"prima di hybrd",(fvec(i),i=1,icount)
 235     rewind(97)
c         rewind(99)
         write(97,*)"prima di hybrd",(x(i),i=1,icount)
         write(99,*)"prima di hybrd",(fvec(i),i=1,icount)
         do i=1,icount
            x(i)=xsave(i)
         enddo
         n=0
c         call newt(x,totpar,check)
c         call broydn(x,totpar,check)
         call hybrd(funcv,totpar,x,fvec,xtol,maxfev,ml,mmu,epsfcn2,
     *        ddiag,mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,
     *                 qtf,wa1,wa2,wa3,wa4)
c     here we reconvert the variational parameters
      icount=0
      do is=1,Ns
         do L1=1,M
            do L2=1,L1
               icount=icount+1
               hm(L1,L2,is)=dcmplx(x(icount),0.d0)
               if (L2.lt.L1) hm(L2,L1,is)=dconjg(hm(L1,L2,is))
            enddo
         enddo
      enddo
      do is=1,Ns
         do L=1,M
            icount=icount+1
            lam(L,is)=x(icount)
         enddo
      enddo
      icount=icount+1
      mu=x(icount)
         

      write(88,*) U,nfev,info
      call flush (88)
      if (info.eq.5) then
!         print*, char(7)
         write(88,*)(fvec(i),i=1,icount)
         epsfcn2=epsfcn2*10.d0
         write(88,*) "shift epsfcn to",epsfcn2
         call flush (88)
         goto 235
      endif

c         write(98,*)U,n,mu,lam(1,1)
c         if (check) then
c            print*,"check=",check
c            stop
c         endif
         

c     write the selfconsistent converged parameters (file to be used as guess)
          rewind(3)
          do is=1,Ns
             do L=1,M
                write(3,*)(hm(L,L2,is),L2=1,M)
             enddo
          enddo
          do L=1,M
             write(3,*)(lam(L,is),is=1,Ns)
          enddo
          write(3,*)mu
          write(3,*)" "
          call flush (3)

C*******calculation of physical quantities******************************
          print*,''
          print*,'___________________________________'
          print*,'___________________________________'
          print*,"_____Converged_____________________"
          print*,'___________________________________'
          do L=1,M
             print*,'Z(',L,')=',(dreal(Z(L,L,is)),is=1,Ns)
             print*,'n(',L,')=',(nm(L,is),is=1,Ns)
             print*,'diag hm(',L,')=',(hm(L,L,is),is=1,Ns)
             do L1=1,L-1
                print*,'non-diag hm(',L,L1,')=',(hm(L,L1,is),is=1,Ns)
             enddo
             print*,'lam(',L,')=',(lam(L,is),is=1,Ns)
          enddo
          print*,'ntot=',ntot,'   mu=',mu
          print*,'___________________________________'
          print*,' '
          
          do L=1,M
             write(20+L,*)U,(dreal(Z(L,L,is)),nm(L,is),dreal(cm(L,is))
     &            ,is=1,Ns)
             call flush (20+L)
          enddo

          print*,"##################################################"
 1000  enddo
       
       
 111   format(a1,$)
 333   format(a3,I1,a1,$)
       
       stop
       end

C***********************************************************************
C***********************************************************************

      subroutine funcv(Mdummy,x,fvec,iflag)
c      subroutine funcv(Mdummy,x,fvec)

      include 'param.dat'
      real*8 fvec(totpar),x(totpar)
      complex*16 SS(M,M,Ns),Sz(M,Ns),S(M,Ns)
      complex*16 hmf(M,M,Ns)

      pi=dacos(-1.d0)

      n=n+1
      print111,']'
 111   format(a1,$)
      
c     here we reconvert the variational parameters
      icount=0
      do is=1,Ns
         do L1=1,M
            do L2=1,L1
               icount=icount+1
               hm(L1,L2,is)=dcmplx(x(icount),0.d0)
               if (L2.lt.L1) hm(L2,L1,is)=dconjg(hm(L1,L2,is))
            enddo
         enddo
      enddo
      do is=1,Ns
         do L=1,M
            icount=icount+1
            lam(L,is)=x(icount)
         enddo
      enddo
      icount=icount+1
      mu=x(icount)


      write(97,*)n,(x(i),i=1,icount)
      call flush(97)
c      do i=1,icount
c         print*,x(i)
c      enddo
c      stop

C*****Diagonalization of the Spin hamiltonian****************************
      if(verb) print*,'diag'
c      call Diagonalization(U,J,hm,lam,cm)
      call Diagonalization
         
C***********************************************************************
C**Calculation of the average values of S(L,is), Sz(L,is) and SS(L1,L2,is1)
      if(verb) print*
      if(verb) print*,"Calcolo i valori medi"
c      call AverageValues(Beta,Sz,S,SS,cm)
      call AverageValues(Sz,S,SS)
      
      do is=1,Ns
         do L=1,M
            Z(L,L,is)=dreal(S(L,is))**2+dimag(S(L,is))**2
c     print*,'Z(',L,' ',is,')=',Z(L,L,is)
            do L1=1,L-1
               Z(L,L1,is)=SS(L,L1,is)
               Z(L1,L,is)=dconjg(Z(L,L1,is))
            enddo
         enddo
      enddo
      if(verb) then 
         do L=1,M
            print*,'S(',L,')=',(S(L,is), is =1,Ns)
            print*,'Sz(',L,')=',(Sz(L,is), is =1,Ns)
            do L1=1,M
               print*,'SS(',L,L1,')=',(SS(L,L1,is),is=1,Ns)
            enddo
            print*,'Z(',L,')=',(Z(L,L,is), is =1,Ns)
         enddo
         print*
      endif


C***** Fermionic problem *************************************************
c      dmu=2.d0/dfloat(100)
c      do iiii=1,99
c         mu=-0.999+iiii*dmu
c      call Fermionic(Z,lam,mu,nm,hmf)
      call Fermionic(hmf)
c      enddo


c      print*,nm(1,1)
      
      do is=1,Ns
         do L=1,M
            do L2=1,M
               if (L.eq.L2) then 
                  hmf(L,L,is)=hmf(L,L,is)*S(L,is)
               else
c                  hmf(L,L2,is)=hmf(L,L2,is)*tmm(L,L2)
               endif
            enddo
         enddo
      enddo

      do L=1,M
         do is=1,Ns
c            cm(L,is)=1/sqrt(nm(L,is)-nm(L,is)**2+epsfcn)-1
         enddo
         if(verb) print*,"dopo fermionic"
         if(verb) print*,"L=",L,"  cm="
     &        ,(cm(L,is),is=1,ns),"  nm=",(nm(L,is),is=1,ns)
         if(verb) print*,"hmf= ",(hmf(L,L,is),is=1,ns)
c         write (60+L,*)(lam(L,is),is=1,ns)
      enddo

c      do L=1,M
c        write(30+L,*)n,(dreal(S(L,is)),is=1,ns),
c     &        (dreal(Sz(L,is)),is=1,ns),(lam(L,is),is=1,ns)
c     &        ,(dreal(hm(L,L,is)),is=1,ns),mu
c        call flush(30+L)
c      enddo

      if (verb) print*,' '

c      print*,hmf(1,1,1)
c      stop

C***********************************************************************
C************** self-consistency equations **************************** 
          icount=0
          sum_n=0.d0
          do is=1,Ns
             do L=1,M
                icount=icount+1
                fvec(icount)=Sz(L,is)+0.5d0-nm(L,is)
                sum_n=sum_n+nm(L,is)/dfloat(Ns)*2
                do L2=1,L
                   icount=icount+1
                   fvec(icount)=hm(L,L2,is)-hmf(L,L2,is)
c                   print*,"ecco gli hm ",hm(L,L2,is),hmf(L,L2,is)
                enddo
             enddo
          enddo          
          icount=icount+1
          fvec(icount)=sum_n-ntot

          if (icount.ne.totpar) stop "icount diff from totpar - fvec"

c          do i=1,icount
c             fvec(i)=fvec(i)*100.d0
c          enddo

          if (verb) then
             print*,"*****************************"
             print*,"selfconsistent equations"
             print*,"number of parameters = ",Mdummy,totpar
             print*,"number of equations = ", icount
             print*,"*****************************"
             do i=1,icount
                print*,fvec(i)
             enddo
             write(99,*)(x(i),i=1,icount),(fvec(i),i=1,icount)
             
             print*,"no aspetta,",sum_n,ntot
             
             print*,"*****************************"
          endif

c          print*,"CCCCCCCCC---COMMON---------CCCCCCCC"
c          print*,hm,cm,ph,Z
c          print*,tmm,lam,nm,epsm
c          print*,U,J,ntot,mu,Beta,n


          return
          end
      
C***********************************************************************
C***********************************************************************

      Subroutine Diagonalization()

      include 'param.dat'
      logical want_to_see_the_matrix
      integer ket(2*M),bra(2*M),differ(2*M),Sorb(M)
      real*8 lagr
      complex*16 hmstar(M,M,Ns),cmstar(M,Ns)
      complex*16 H(totst,totst),EVEC(totst,totst)
      real*8 E(totst)
      complex*16 WORK(2*totst-1)
      real*8 WORK2(3*totst-2)
      common /diag/ E,EVEC

c*****************writing matrix elements********************************
c      do is=1,ns
c      do L=1,M
c      do L2=1,M
c      print*,U,((hm(L,L2,is),L2=1,M),is=1,Ns),(lam(L,is),L=1,M),J
c      enddo
c      enddo
c      enddo

      do L=1,M
         do is=1,Ns 
            do L2=1,M
               hmstar(L,L2,is)=dconjg(hm(L,L2,is))
               if(verb) then
                  print*,"hm(",L,",",L2,",",is,")=",hm(L,L2,is)
                  print*,"hmstar(",L,",",L2,",",is,")=",
     &                 hmstar(L,L2,is)
               endif 
            enddo  !L2
            cmstar(L,is)=dconjg(cm(L,is))
         enddo !is
         if(verb) then
            print*,"cm(",L,")=",(cm(L,is),is=1,Ns)
            print*,"cmstar(",L,")=",(cmstar(L,is),is=1,Ns)
            print*,"in diag, lam(",L,")=",(lam(L,is),is=1,Ns)
         endif
      enddo !L

      do II=1,totst
         do JJ=1,totst
            H(II,JJ)=0.d0
            do L=1,2*M  !! L goes from 1 to 2*M (even if Ns=1, obviously)
               ket(L)=mod(int((II-1)/2**(L-1)),2)
               bra(L)=mod(int((JJ-1)/2**(L-1)),2)
            enddo
            if (II.eq.JJ) then ! --------- diagonal part
               sum=0.d0
               do L=1,2*M
                  sum=sum+ket(L)-0.5
               enddo
               sum_up=0.d0
               sum_down=0.d0
               do L=1,M
                  sum_up=sum_up+ket(2*L-1)-0.5
                  sum_down=sum_down+ket(2*L)-0.5
               enddo
               H(II,JJ)=0.5*(U-2*J)*sum*sum
               H(II,JJ)=H(II,JJ)-0.5*J*(sum_up*sum_up+sum_down*sum_down) 
               sum=0.d0
               lagr=0.d0
               do L=1,M
                  Sorb(L)=ket(2*L-1)+ket(2*L)-1
                  sum=sum+Sorb(L)**2
                  if (Ns.eq.1) then
                     lagr=lagr+lam(L,1)*(Sorb(L)+1)
                  else
                     do is=1,Ns
                        lagr=lagr+lam(L,is)*ket(2*L-2+is)
                     enddo
                  endif
               enddo
               H(II,JJ)=H(II,JJ)+lagr+J*sum
c               print*,sum, lagr,Sorb(1),Sorb(2) 
            else  ! --------------- non-diagonal part
               diff=0.d0
cccc diff: We find which bra and ket differ by only one occupied/empty site 
cccc      (for hopping-like terms), or by two sites (for spin-flip type terms)
               do L=1,2*M 
                  diff=diff+abs(ket(L)-bra(L))
                  differ(L)=bra(L)-ket(L) !+1 means creation operator applied
                  if (abs(ket(L)-bra(L)).eq.1) then 
                     id2=id
                     id=L
                     in2=in !if diff=2, this index is the penultimate
                     in=int((L+1)/2) ! keeps track of which site and spin
                     ins2=ins     ! state differ.
                     ins=mod(L+1,2)+1
                  endif
               enddo
               iss=ins !this is used in the arrays and keeps into account
               iss2=ins2  ! the spin symmetry breaking or not
               if (Ns.eq.1) then
                  iss=1
                  iss2=1
               endif
c               print*,II,JJ,diff,diff_up,diff_do
cccccccccc   intra-band hm, terms (h*_mm S_m +h_mm S*_m) 
               if (diff.eq.1) then 
                  if (differ(id).eq.1) then 
c                     print*,hm(1,1,1)
                     H(II,JJ)=hmstar(in,in,iss)*cm(in,iss)
     &                    +hm(in,in,iss)
                  else if (differ(id).eq.-1) then 
                     H(II,JJ)=hm(in,in,iss)*cmstar(in,iss)
     &                    +hmstar(in,in,iss)
                  endif           
               endif           
ccccccc   interband terms (h_mn S*_m S_n +hc) (from the on-site hybridization)
               if ((diff.eq.2).and.(ins.eq.ins2)) then !here we use ins
                  if((differ(id).eq.1).and.(differ(id2).eq.-1))then
                     H(II,JJ)=H(II,JJ)+hm(in,in2,iss)+
     &                    hmstar(in,in2,iss)*cmstar(in2,iss)*cm(in,iss)
                  else if((differ(id).eq.1).and.(differ(id2).eq.1))then
                     H(II,JJ)=H(II,JJ)+hm(in,in2,iss)*cm(in2,iss)+
     &                    hmstar(in,in2,iss)*cm(in,iss)
                  elseif((differ(id).eq.-1).and.(differ(id2).eq.-1))then
                     H(II,JJ)=H(II,JJ)+hm(in,in2,iss)*cmstar(in,iss)+
     &                    hmstar(in,in2,iss)*cmstar(in2,iss)
                  else if((differ(id).eq.-1).and.(differ(id2).eq.1))then
                     H(II,JJ)=H(II,JJ)+hmstar(in,in2,iss)+
     &                    hm(in,in2,iss)*cm(in2,iss)*cmstar(in,iss)
                  endif
               endif
cccccccccccccccccc spin flip term cccccccccccccccccccccccccccccccccccccc
               do L1=1,M-1
                  do L2=L1+1,M
                     if (((ket(2*L1-1).eq.1).and.(ket(2*L1).eq.0).and.
     &                    (ket(2*L2-1).eq.0).and.(ket(2*L2).eq.1).and.
     &                    (bra(2*L1-1).eq.0).and.(bra(2*L1).eq.1).and.
     &                    (bra(2*L2-1).eq.1).and.(bra(2*L2).eq.0)).or.
     &                    ((ket(2*L1-1).eq.0).and.(ket(2*L1).eq.1).and.
     &                    (ket(2*L2-1).eq.1).and.(ket(2*L2).eq.0).and.
     &                    (bra(2*L1-1).eq.1).and.(bra(2*L1).eq.0).and.
     &                  (bra(2*L2-1).eq.0).and.(bra(2*L2).eq.1))) then
                        iflag=1
                        do L3=1,M
                           if ((L3.ne.L1).and.(L3.ne.L2)) then
                              if ((ket(2*L3-1).ne.bra(2*L3-1))
     &                             .or.(ket(2*L3).ne.bra(2*L3))) then
                                 iflag=0
                                 goto 101
                              endif
                           endif
                        enddo
 101                    if (iflag.eq.1)     H(II,JJ)=H(II,JJ)-J
                     endif
                  enddo
               enddo
ccccccccccccccccc pair hopping term cccccccccccccccccccccccccccccccccccc
               do L1=1,M-1
                  do L2=L1+1,M
                     if (((ket(2*L1-1).eq.1).and.(ket(2*L1).eq.1).and.
     &                    (ket(2*L2-1).eq.0).and.(ket(2*L2).eq.0).and.
     &                    (bra(2*L1-1).eq.0).and.(bra(2*L1).eq.0).and.
     &                    (bra(2*L2-1).eq.1).and.(bra(2*L2).eq.1)).or.
     &                    ((ket(2*L1-1).eq.0).and.(ket(2*L1).eq.0).and.
     &                    (ket(2*L2-1).eq.1).and.(ket(2*L2).eq.1).and.
     &                    (bra(2*L1-1).eq.1).and.(bra(2*L1).eq.1).and.
     &                  (bra(2*L2-1).eq.0).and.(bra(2*L2).eq.0))) then
                        iflag=1
                        do L3=1,M
                           if ((L3.ne.L1).and.(L3.ne.L2)) then
                              if ((ket(2*L3-1).ne.bra(2*L3-1))
     &                             .or.(ket(2*L3).ne.bra(2*L3))) then
                                 iflag=0
                                 goto 102
                              endif
                           endif
                        enddo
 102                    if (iflag.eq.1)     H(II,JJ)=H(II,JJ)-J
                     endif
                  enddo
               enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            endif ! diagonal or not            
         enddo  ! state JJ
c         if (verb) print*,(ket(L),L=1,2*M)
       enddo    ! state II

C********************check of the hamiltonian matrix********************
      want_to_see_the_matrix=.false.
      if (want_to_see_the_matrix) then
         print*,''
         do II=1,totst
            do JJ=1,totst
               print333,H(II,JJ)
               if (mod(jj,2).eq.0)  then
                  print334,'|'
               endif
               print334,' '
            enddo
            print*,' '
            if (mod(II,2).eq.0)  then
               print*,'-------------------------
     &              ------------------------------'
            endif
         enddo
         
 333  format ((f5.2,f5.2),$)
 334  format (a2,$)

      endif

C****************diagonalization*****************************************

c      call zheev('V','A','U',totst,H,totst,vl,vu,il,iu,2*dlamch('S'),
c     &     idum,E,ZEVEC,totst,WORK,2*totst-1,RWORK,IWORK,IFAIL,INFO)

      call zheev('V','U',totst,H,totst,E,WORK,2*totst-1,WORK2,INFO)
      
      if (INFO.ne.0) then 
         print*,'INFO=',INFO
         stop 'error in the diagonalization routine'
      endif

      do II=1,totst
         do JJ=1,totst
            EVEC(II,JJ)=H(II,JJ)
            if (want_to_see_the_matrix) then
               print*,'EVEC(',II,',',JJ,')=',evec(II,JJ)
            endif
         enddo
      enddo

      return
      end

C************************************************************************
C************************************************************************

      subroutine AverageValues(Sz,S,SS)

      include 'param.dat'
      complex*16 Sz(M,ns),SS(M,M,Ns)
      complex*16 S(M,ns),Xs,Xz
      complex*16 EVEC(totst,totst)      
      real*8 E(totst)

      common /diag/ E,EVEC

      PartF=0.d0
      do II=1,totst
         PartF=PartF+dexp(-Beta*(E(II)-E(1)))
      enddo
                 ! <Sx(L)>, <Sz(L)>

      do L=1,M
      do is=1,ns
         Sz(L,is)=dcmplx(0.d0,0.d0)
         S(L,is)=dcmplx(0.d0,0.d0)
         do II=1,totst
            Xz=dcmplx(0.d0,0.d0)
            Xs=dcmplx(0.d0,0.d0)
            do JJ=1,totst
c               occ=mod(int((JJ-1)/2**(2*L-1)),2) ! 2*L --> paramagnetic case
               occ=mod(int((JJ-1)/2**(2*L-1-is+1)),2)
               if (occ.eq.1) then
c                  JJflip=JJ-2**(2*L-1) !para
                  JJflip=JJ-2**(2*L-1-is+1)
               else
c                  JJflip=JJ+2**(2*L-1) !para
                  JJflip=JJ+2**(2*L-1-is+1)
               endif
               Xz=Xz+(occ-0.5d0)*Evec(JJ,II)**2
               if (occ.eq.1) then   
               Xs=Xs+Evec(JJ,II)*Evec(JJflip,II)
               else
               Xs=Xs+Evec(JJ,II)*Evec(JJflip,II)*cm(L,is)
               endif
            enddo
            Sz(L,is)=Sz(L,is)+Xz*dexp(-Beta*(E(II)-E(1)))/PartF
            S(L,is)=S(L,is)+Xs*dexp(-Beta*(E(II)-E(1)))/PartF
         enddo   
      enddo
      enddo
                ! mean values <S*(L)S(L')>
      do L1=1,M
      do L2=1,M
      do is=1,Ns
         SS(L1,L2,is)=dcmplx(0.d0,0.d0)
         do II=1,totst
            Xs=0.d0
            do JJ=1,totst
c               occ1=mod(int((JJ-1)/2**(2*L1-1)),2) !2*L --> paramagnetic case
c               occ2=mod(int((JJ-1)/2**(2*L2-1)),2) !2*L --> paramagnetic case
               occ1=mod(int((JJ-1)/2**(2*L1-1-is+1)),2)
               occ2=mod(int((JJ-1)/2**(2*L2-1-is+1)),2)
               if (occ1.eq.1) then   
                  JJflip1=JJ-2**(2*L1-1-is+1)
                  coeff1=1.d0
               else
                  JJflip1=JJ+2**(2*L1-1-is+1)
                  coeff1=cm(L1,is)
               endif
               if (occ2.eq.1) then !we apply S* on the bra   
                  JJflip2=JJ-2**(2*L2-1-is+1)
                  coeff2=1.d0 
               else
                  JJflip2=JJ+2**(2*L2-1-is+1)
                  coeff2=dconjg(cm(L2,is))
               endif
               Xs=Xs+Evec(JJflip1,II)*Evec(JJflip2,II)*coeff1*coeff2
            enddo
            SS(L1,L2,is)=SS(L1,L2,is)
     &           +Xs*dexp(-Beta*(E(II)-E(1)))/PartF
         enddo
      enddo
      enddo
      enddo

ccccccccccccccccc check!

      if (abs(SS(1,2,1)-SS(2,1,1)).gt.1.d-10) then 
         print*,SS(1,2,1),SS(2,1,1)
         stop "error in the average values"
      endif
c      print*,"SxSx",SxSx(1,2)

      return
      end

C************************************************************************
C************************************************************************
      
      subroutine Fermionic(hmf)

      include 'param.dat'
      complex*16 hmf(M,M,Ns),yy
      complex*16 H(M,M)
      real*8 E(M),D(M)
      complex*16 WORK(2*M-1)
      real*8 WORK2(3*M-2)



      pi=dacos(-1.d0)

      do is=1,Ns
         do L1=1,M
               nm(L1,is)=0.d0
            do L2=1,M
               hmf(L1,L2,is)=dcmplx(0.d0,0.d0)
            enddo
         enddo
      enddo

      do L=1,M
         D(L)=2*tmm(L,L)
      enddo
   
      de=2.d0/dfloat(ne)

      xnorm=0.d0
      do ie=1,ne
         en=-1.d0+de*ie        
         if (abs(en).lt.1.d0) then
            dos=2/pi*dsqrt(1.d0-en**2) 
         else
            dos=0.d0
         endif
         xnorm=xnorm+de*dos
      enddo

      do is=1,Ns
         
         do ie=1,ne
            en=-1.d0+de*ie      !! the energy must be rescaled with bwidth,
            if (abs(en).lt.1.d0) then
               dos=2/pi*dsqrt(1.d0-en**2)/xnorm !   dos*de  must not.
            else
               dos=0.d0
            endif
            do L=1,M
               do L2=1,M
                  if (L2.lt.L) then
                     H(L,L2)=Z(L,L2,is)*tmm(L,L2)
                  else if (L2.eq.L) then
                     H(L,L)=epsm(L)-lam(L,is)-mu+D(L)*Z(L,L,is)*en  
                  else
                     H(L,L2)=dconjg(H(L2,L))
                  endif
               enddo
c               print*,(H(L,L2),L2=1,M)
            enddo
            
ccccccc diagonalization of the M-orbital matrix for a given k (a given energy)
            call zheev('V','L',M,H,M,E,WORK,2*M-1,WORK2,INFO)
            if (INFO.ne.0) then 
               print*,'INFO=',INFO
               stop 'error in the fermionic diagonalization routine'
            endif
c            print*,epsm(1),lam(1,1),mu,D(1),Z(1,1,1),en,E(1)
ccc diagonal n.n. mean-field \sum_j Zt_ij<f^+_ia f_ja> and density <f^+_ia f_ia>
            do L1=1,M
               xx=0.d0
               do L2=1,M
                  xx=xx+abs(H(L1,L2))**2/(dexp(Beta*E(L2))+1)
c                  print*,abs(H(L1,L2)) **2,xx,1/(dexp(Beta*E(L2))+1)
               enddo
               nm(L1,is)=nm(L1,is)+de*dos*xx
               hmf(L1,L1,is)=hmf(L1,L1,is)+de*dos*en*D(L1)*xx
            enddo

ccccc off-diagonal local mean-field <f^+_ia f_ib>
            do L1=1,M
               do L2=1,L1-1
                  yy=dcmplx(0.d0,0.d0)
                  do L=1,M
                     yy=yy+dconjg(H(L1,L))*H(L2,L)/(dexp(Beta*E(L))+1)
                  enddo
                  hmf(L1,L2,is)=hmf(L1,L2,is)+de*dos*yy
               enddo
            enddo
            
         enddo  !ie

            do L1=1,M
               do L2=1,L1-1
                  hmf(L1,L2,is)=hmf(L1,L2,is)*tmm(L1,L2)
                  hmf(L2,L1,is)=dconjg(hmf(L1,L2,is))
               enddo
            enddo
      enddo !is
      
c      print*,"$$$$$$$$$$$$$$$$$$$$$$$$"
c      print*,epsm(1),lam(1,1),mu,tmm(1,1)
c      print*,nm(1,1),D(1),Z(1,1,1),hmf(1,1,1)
c      write(91,*)epsm(1),lam(1,1),D(1)
c      call flush (91)
c      print*,"$$$$$$$$$$$$$$$$$$$$$$$$"

      return
      end


C************************************************************************
C************************************************************************
