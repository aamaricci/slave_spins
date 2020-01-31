!	This subroutine solve the problem of the slave fermion f for the MF 
!	analysis of the Slave-Spin problem (single-site MF approximation)
! 	Here we solve the fermion Hamiltonian for the NON-HYBRIDIZED case. 
!	It takes as input the number of the orbitals, the guess value of the 
!	spectral weight Z_{m,sigma} and the one of the Lagrange multiplier 
!	Lambda_{m,sigma} with m= 1...N and sigma=-,+
!	 
!	STRUCTURE:
!	This program is composed by two main parts: 
!	1) Computation of the Fermi level
!	2) Computation of the average of nf, Ef over the ground state for all 
!	the 2N_orb levels.

!	OUTPUT
!	The output of the program are :
!	*) <nf_ms> = 1/N sum_k nf_ms(k)
!	*) <E_ms> = 1/N sum_k E_ms(k) nf_ms(k)
!	for each 2N_orb levels.

	subroutine fermion_equiv(rho,N_orb,N_slave,t_ms,Lambda_ms,Z_ms,nf_ms,E_ms,c_ms)

!	%%%%%%%%%%%%%%%%%%%%%%%%%%    DECLARATIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%
	implicit none
	
!	INPUT DATA
	integer :: N_orb, N_slave
!	N_orb is the number of the orbitals, 
!	N_slave=2*N_orb is the number of slave particle
	integer :: N_el
	real*8 :: rho	
!	rho is the number of the electrons for site 
!	N_el= rho*nsites is the total number of electrons
!	nsites is fixed by the coarse graining of the BZ i.e. nsites=nk**2
	real*8, dimension(N_slave) :: t_ms
!	t_ms are arrays of dimension 2*N_orb, the pedix ms stays for 
!	the orbital and spin labels {m, sigma}. Model Data
	real*8, dimension(N_slave) :: Lambda_ms, Z_ms
!	Lambda_ms, Z_ms are arrays of dimension 2*N_orb, the pedix ms stays for 
!	the orbital and spin labels {m, sigma}. Guesses
	
!	OUTPUT DATA 
	real*8 :: Ef
	real*8, dimension(N_slave) :: nf_ms, E_ms, c_ms
!	Output arrays of dim = N_slave that contains 
!	*) <nf_ms> = 1/N sum_k nf_ms(k)
!	*) <E_ms> = 1/N sum_k E_ms(k) nf_ms(k)
!	*) [sqrt(nf_ms*(1-nf_ms))**(-1)] - 1 
	real*8, dimension(N_slave) :: sum_e, sum_n
	real*8, dimension(N_orb) :: nf_m, E_m	
	real*8 :: nf,  E
!	INDICES AND COUNTERS
	integer, parameter :: nk=301 ! number of k-vector of the BZ
	real*8, parameter :: pi=atan(1.d0)*4.d0
	real*8, parameter :: eps = 2.d0*pi/nk ! k-stepsize
	real*8, parameter :: delta = 10.d-5 ! high-doping renormalization for c_ms
	integer :: j,i,l, m, kl, full, ind,ind_ord, i_ord
	real*8 :: kx, ky !NB: k id in a unit (a lattice parameter)
	integer :: nkmax, nmax


!	ENERGY VARIABLES
	real*8, allocatable, dimension(:) :: Ek,Eigenk,Ek_ord 
	integer, allocatable, dimension(:) :: i_orden
!	Ek is the array of the eigenstates dim = nmax= N_orb*nk**2
!	This array must be ordered in order to establish the state of filling 
!	of the eigenstates.
!	i_orden is the arrays we use for the sorting procedure
!	Ek_ord is the array in which we can optionally save the energies 
!	in increasing order to save the Ef value

!	%%%%%%%%%%%%%%%%%%%%%%    ASSIGNMENT & ALLOCATION    %%%%%%%%%%%%%%%%%%%	
	nkmax=nk**2
	N_el = rho*nkmax 
	nmax=nkmax*N_slave 
!	Eigenvalues and indices array/matrix
	allocate(Ek(nmax))
	allocate(Eigenk(nmax))
	allocate(Ek_ord(nmax))
	allocate(i_orden(nmax))

!	write(*,*)'INPUT'
!	write(*,*)'L',  Lambda_ms(:)
!	write(*,*)'Z', Z_ms(:)
!	%%%%%%%%%%%%%%%%%%%%%    COMPUTING THE FERMI LEVEL   %%%%%%%%%%%%%%%%%%%
!	Filling the array Ek and the matrix Ek_ms
!	Example: 2dimension, all s-symmetry NN-hopping
	kl = 1
	do l=1,N_slave
		sum_e=0.d0
		sum_n=0.d0
		kx=-pi
		ky=-pi
		do i=1,nk
			do j=1,nk 
				Eigenk(kl) = - 2.d0*t_ms(l)*(dcos(kx)+ dcos(ky)) 
				Ek(kl) = Z_ms(l)*Eigenk(kl) - Lambda_ms(l) 
				kl = kl + 1
			   	ky=ky+eps
		   end do
		   kx=kx+eps
		   ky=-pi
		end do
	end do !l-loop

!	Initialization of the Indices Array (nmax = N_slave*nk**2)
 	do ind=1,nmax 
         i_orden(ind)=ind
    end do 
!	Sorting Subroutine    
	call  indexx(nmax, Ek, i_orden)

!	Ordening of the Eigenvalues Arrays	
!	do ind=1,nmax
! 	   ind_ord = i_orden(ind)
! 	   Ek_ord(ind) = Ek(ind_ord)
!    end do     
!	Assignation of the Fermi Energy 
!	(we save it just in case we would like to know it but actually we don't need
!	this information)   
!	i_ord =  i_orden(N_el)
!	Ef =  Ek(i_ord) 

!	%%%%%%%%%%%%%%%%%%%%%      COMPUTING THE OUTPUT     %%%%%%%%%%%%%%%%%%%%	
!	This cycle looks only the filled states and save the sum over k of each l 
!	The empty states, not considered here, do not contribute to the sums, 
!	indeed if nf_k(l) = 0 also Eigenk(l)*nf_k(l)=0 
	sum_e= 0.d0
	sum_n=0.d0
    do full= 1,N_el 
    	i_ord =  i_orden(full)
        l= 1 + (i_ord - 1)/nkmax
        sum_e(l)= sum_e(l) + Eigenk(i_ord)
        sum_n(l)=sum_n(l) + 1 
    end do
    do l=1,N_slave 
		E_ms(l)= sum_e(l)/nkmax
		nf_ms(l)= sum_n(l)/nkmax
		c_ms(l)=(sqrt(nf_ms(l)*(1-nf_ms(l))))**(-1.d0) - 1.d0 
      end do

!	PARAMAGNETIC STATES
!	Let us require the state not to be polarized
!	do m=1,N_orb
!	  E_m(m)= (E_ms(m*2-1)  + E_ms(m*2))/2.d0
!	  E_ms(m*2)=E_m(m)
!	  E_ms(m*2-1)=E_m(m)
!	  nf_m(m)= (nf_ms(m*2-1)  + nf_ms(m*2))/2.d0 
!	  nf_ms(m*2)=nf_m(m)
!	  nf_ms(m*2-1)=nf_m(m)
!	end do

!	ALL EQUIVALENT ORBITALS 
	E=0.d0
	nf=0.d0
	do l=1,N_slave 
	  E= E + E_ms(l)
	  nf= nf + nf_ms(l)
	end do
	do l=1,N_slave 
	  E_ms(l)=E/N_slave
	  nf_ms(l)=nf/N_slave
	end do   
	
    do l=1,N_slave 
	  c_ms(l)=(sqrt( nf_ms(l)*(1-nf_ms(l))) + delta )**(-1.d0) - 1.d0 	
	  !we use delta to avoid unphisical divergencies as n goes to 0
    end do	

	deallocate(Ek, Eigenk, Ek_ord, i_orden)
	end subroutine fermion_equiv

!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	%%%%%%%%%%%%%%%%%%%%%%%%%%    SUBROUTINES    %%%%%%%%%%%%%%%%%%%%%%%%%%%
!	Sorting 	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	SUBROUTINE indexx(n,arr,indx)
	implicit real*8 (a-h,o-z)
	dimension  arr(n),indx(n)
	PARAMETER (M=7,NSTACK=50)
	dimension istack(NSTACK)
      
	do 11 j=1,n
	  indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
	END