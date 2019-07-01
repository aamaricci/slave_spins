!	OUT HF VERSION 
!	This program solves the problem of the pseudospin Hamiltonian in the 
!	Slave-Spin approach within the single-site MF approximation. 
!	It takes as input the number of the orbitals, the guess value of the 
!	Weiss field h_{m,sigma} and the one of the Lagrange multiplier 
!	Lambda_{m,sigma} with m= 1...N and sigma=-,+
 
!	STRUCTURE:
!	This program is composed by three parts: 
!	1) Construction of the pseudospin Hamiltonian matrix in the Hilbert 
!	space of the (2N_orbitals)-pseudospins
!	2) Diagonalization of the pseudospin Hamiltonian and selection of the
!	 ground state (i.e. the eigenvector corresponding to the lowest 
!	eigenvalue). For the diagonalization procedure we use the LAPACK routine 
!	dsyev.f via the Fortran 90 interface diaham. 
!	NB: Diagonalization subroutine for a REAL & SYMMETRIC matrix.
!	3) Computation of the average of Sz, Sx over the ground state for all 
!	the 2N_orb pseudospins.

!	METHOD: The BIT MAPPING
!	The Hilbert space of 2N_orb pseudospin admits 2**{2N_orbitals} states 
!	corresponding to the possible arrangements of all the pseudospins. 
!	Thus every states can be written as a vector of 2N_orbitals components: 
!	|Sz_{m,sigma}> = |Sz_{1-},Sz_{1+},...Sz_{N_orb,-},Sz_{N_orb,+}> 
!	Using the mapping: 
!	Sz up <-> TRUE		Sz down <-> FALSE
!	We can associate every |Sz_{m,sigma}> vector to a 2N_orb-bit string i.e.
!	every state to an integer number I 
!	|I>=|bit_{n}> 
!	with the convention:  
!	n	|	m      sigma
!	0	|	1	-
!	1	|	1	+
!	2	|	2	-
!	3	|	2	+
!	...	|	...	...
!	2N_orb-1|	N_orb	+
!
!	Example1: 
!	Given the state |I>, the value of Sz_{M,S} is given by the value of the 
!	bit in position pos((M-1)*2) if S= -, pos((M-1)*2 + 1) if S= + .
!	Example2: 
!	Given the state |I>, the value of the bit in position pos(i) gives the 
!	value of Sz of the pseudospin in the orbital m=(i/2)+1 and with spin 
!	sigma = -/+ if i is even/odd respectively. 
!
!	Following this procedure the states are already ordered simply counting 
!	from |0> up to |2N_orbital-1> (-1 since we start from 0). Moreover the 
!	integers can be use as labels for the Hamiltonian elements always having 
!	in mind that the first line(column) j=1 is given by the state <0|(|0>)  
!	i.e. j=I+1. 
!	NB: This mapping allows us to treat a problem of N-spin with N <= 31 
!	since a standard integer has 31 bits available (the 32-th is useless 
!	since is the sign of I). In our multiorbila case where N= 2N_orbitals 
!	this means N_orbitals <=15. 
!	NB: In what follows we use these intrinsic bit manipulation functions: 
!	BTEST(I, POS)	.TRUE. if the position number POS of I is 1
!	IEOR(I, J)		gives an INTEGER whose bits are the logical exclusive OR 
!					between I and J (i.e. the bits series of this new integer
!					has .FALSE. in each position but the ones in which I and J 
!					has a different bit value)
!	IBCHNG(I,POS)	reverses the value of the bit POS in the integer I, 
!					gives back a new INTEGER (this function needs the intel 
!					compiler)

!	OUTPUT
!	The output of the program are :
!	*)  <Sz_{m, sigma}>  in order to compare 
!	    <nf_{m, sigma}> vs <Sz_{m, sigma}> + 1/2
!	*)  <Sx_{m, sigma}> in order to compute the new Weiss field 
!		h_ms Z_{m, sigma} ~ <Sx_{m, sigma}> 
!		and the new q.ple weight, i.e.
!		Z_{m, sigma} = 4 <O_{m, sigma}>^2
!	    if there is no orbital hybridization 
!	
	subroutine pspin_corr(U,Up,JH,N_orb,N_pspins,Lambda_ms,h_ms,O_ms,OT_ms,c_ms,Sz0_ms, flagCorr, corr_1d1d, corr_1d1u, corr_1d2d,corr_1d2u)

!	%%%%%%%%%%%%%%%%%%%%%%%%%%    DECLARATIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%
	implicit none

!	INPUT DATA
!	The interaction couplings are given by:
	real*8 :: U ! intraorbital Hubbard term
	real*8 :: Up! interorbital Hubbard term
	real*8 :: JH! Hund coupling 
	integer :: N_orb,N_pspins,N_states
!	Given the number of the orbitals (N_orb) the number of possible states 
!	is given by N_states= 2**(2*N_orb) where every states is defined by the 
!	configuration of N_pspins=2*N_orb pseudospins
	real*8, dimension(N_pspins) :: Lambda_ms, h_ms
!	Lambda_ms, h_ms are arrays of dimension 2*N_orb, the pedix ms stays for 
!	the orbital and spin labels {m, sigma}

!	PSPIN HAMILTONIAN VARIABLES
	real*8, allocatable, dimension(:,:) :: Ham
!	Ham is the matrix N_states x N_states 
	real*8, allocatable, dimension(:) :: gs
!	gs is the array of dimension N_states and allows us to define the ground 
!	states in terms of the N_states base vector as |GS> = sum_j gs(j) |j-1> 
!	with j = 1,... N_states 

!	OUTPUT DATA 
	real*8, dimension(N_pspins) :: Sz0_ms, O_ms, OT_ms, c_ms
    real*8 :: corr_1d1d, corr_1d1u, corr_1d2d, corr_1d2u
!	Again here the pedix ms stays for the orbital and spin labels {m, sigma}
!	Sz0_ms, O_ms are arrays of dimension N_pspins. They contain the 
!	average values on the GS of Sz0_ms, O_ms, transpose(O_ms) respectively.
!	c_ms is the array containing [sqrt(nf_ms*(1-nf_ms))**(-1)] - 1 
	integer :: flagCorr
!	INDICES AND COUNTERS
	integer :: j,k,l

!	%%%%%%%%%%%%%%%%%%%%%%%%%%    OPENING FILES    %%%%%%%%%%%%%%%%%%%%%%%%%
!	open(unit=10,file='input.dat') li prende non gli serve piu leggerli 
!	open(unit=11,file='Ham.dat')

!	%%%%%%%%%%%%%%%%%%%%%%    ASSIGNMENT & ALLOCATION    %%%%%%%%%%%%%%%%%%%
!	Defines N_states and allacates memory for the Hamiltonian	
	N_states= 2**(2*N_orb) 
	allocate(Ham(N_states,N_states))
	allocate(gs(N_states))
	
!	%%%%%%%%%%%%%%%%%    1.HAMITONIAN MATRIX CONSTRUCTION    %%%%%%%%%%%%%%%
!	write(*,*)'imput', N_orb,N_pspins,N_states
!	write(*,*)'imput2', Lambda_ms(1), Lambda_ms(2), h_ms(1), U,Up,JH
	call Hambuilder(N_orb,N_pspins,N_states,Lambda_ms,h_ms,c_ms,Ham,U,Up,JH)

!	Check of 'Ham.dat', output of the hamiltonian builder routine 
!	do j=1,N_states
!	    write(*,41) (Ham(j,k),k=1,N_states)
!	    write(11,41) (Ham(j,k),k=1,N_states) 
!	end do
!	close(11)

!	%%%%%%%%%%%%%%%%%%%%%%%%    2.DIAGONALIZATION    %%%%%%%%%%%%%%%%%%%%%%%
!	The subroutine diaham takes as input 
!	*) N_states -> dim of Ham
!	*) Ham matrix 
!	and assigns gs(j)  
	call diaham(N_states,Ham,gs)
!	write(*,*) gs(:)

!	%%%%%%%%%%%%%%%%%%%    3.AVERAGE Sz, Sx OVER THE GS    %%%%%%%%%%%%%%%%%
	call av_s(N_pspins,N_states,gs,Sz0_ms,O_ms,OT_ms,c_ms)
!	write(*,*) 'Sz=',  Sz0_ms(:)
!	write(*,*) 'O=',  O_ms(:)
!	write(*,*) 'OT=',  OT_ms(:)
!	Constraint
!	nf_ms(:) = Sz0_ms(:) + 0.5d0

!	%%%%%%%%%%%%%%%%%%%%%      4.CORRELATIONS          %%%%%%%%%%%%%%%%%%%%%
!	Only at the end
	if(flagCorr.eq.1)then
	call correlation(N_pspins,N_states, gs, Sz0_ms, corr_1d1d, corr_1d1u, corr_1d2d, corr_1d2u)
	end if

!	%%%%%%%%%%%%%%%%%%%%    DEALLOCATION & END PROGRAM    %%%%%%%%%%%%%%%%%%

	40   format(10(1x,f12.6,2x))
	41   format(500(f12.2))

	deallocate(Ham, gs)
	end subroutine pspin_corr

!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	%%%%%%%%%%%%%%%%%%%%%%%%%%    SUBROUTINES    %%%%%%%%%%%%%%%%%%%%%%%%%%%
!	Hamiltonian Builder	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine Hambuilder(N_orb, N_pspins, N_states, Lambda_ms, h_ms, c_ms, Ham, U, Up, JH)
   !	%% DECLARATIONS %%
   implicit none 
   !	Input
   integer :: N_orb, N_pspins, N_states
   real*8, dimension(N_pspins) :: Lambda_ms, h_ms, c_ms
   real*8, dimension(N_states,N_states) :: Ham
   real*8 :: U ! intraorbital Hubbard term
   real*8 :: Up! interorbital Hubbard term
   real*8 :: JH! Hund coupling 

   !	Mute Indices
   integer :: j,k 
   !	j, k = 1, N_states 
   !	*) lines/columns labels of the Hamiltonian matrix i.e. Ham(j,k) 
   !	*) states |j-1>, <k-1|
   integer :: l, i, imax 
   !	l=1,N_pspins is the index the we use to run over the arrays Lambda_ms, 
   !	h_ms, Sz_ms 
   !	i=0,imax is the index that we use to run over the bits of a given state 
   !	imax=N_pspins-1 (since the first bit is in pos(0))
   integer :: m, s
   !	m = 1, N_orb, orbital index;  s= 1,2 spin index: sigma= +,-
   integer :: flips
   !	flips saves the number of TRUE bits coming out from the IEOR(j,k) of two
   !	integer j,k i.e. the number of spin flips between the  two states |j-1>,
   !	|k-1> (see Fcounter function). We use this varable in order to find the 
   !	states connected by the Sx pseudospin operator

   !	Pseudospin Arrays
   integer, allocatable, dimension(:) :: Sz_ms
   !	Sz_ms is an array of dimension N_pspins. It contains the Sz_ms value of 
   !	the pseudospins in a given state
   integer, allocatable, dimension(:) :: Sz_m
   integer, dimension(2) :: Sz_s
   !	*) Sz_m has N_orb components and contains the sum of the pseudospins of 
   !	the two neighbouring boxes of n's bit string i.e SUM_s Sz_ms for each m
   !	*) Sz_s has 2 components and contains the sum over all the orbitals of 
   !	the pseudospins with the same spin i.e. SUM_m Sz_m-  and  SUM_m Sz_m+
   !
   !	Counters
   real*8 :: sum, sum0, sum1, sum2, sum3
   real*8 :: sum_s, sum_p, sum_m
   integer :: flag
   !	Functions
   integer :: Fcounter, spincode

   !	%% ALLOCATIONS %%
   allocate (Sz_ms(N_pspins))
   allocate (Sz_m(N_orb))

   !	%% HAMILTONIAN MAKING LOOP %%
   imax=N_pspins-1 
   !	Max number of bits we need starting from the bit 0-th
   do j=1,N_states !lines of Ham 
      do k=1,N_states !columns of Ham
         sum=0.d0
         sum0=0.d0
         sum1=0.d0
         sum2=0.d0
         sum3=0.d0 
         sum_p=0.d0
         sum_m=0.d0
         sum_s=0.d0
         flips=Fcounter(j-1,k-1, imax)
         !diagonal elements j=k
         if (flips.eq.0) then 
            !	PseudoSpin Array Sz_ms Assigmment
            do i=0,imax
               Sz_ms(i+1) = spincode(BTEST(j-1,i))
            end do
            !	Save the pseudospins'configuration for the state |j-1> in an array 
            !	Sz_ms(l) with l=1, N_pspins
            do l=1,N_pspins
               !	H0_Sz
               sum0= sum0 + Lambda_ms(l)*0.5d0*(Sz_ms(l) + 1.d0)
               !	Hint_U'/2
               sum1= sum1 + 0.5d0*Sz_ms(l) 
            end do
            sum1= sum1**2.d0
            !	PseudoSpin Array Sz_m, Sz_s Assigmment
            do m=1,N_orb
               Sz_m(m)= Sz_ms(m*2-1) + Sz_ms(m*2)
               sum_m= sum_m + Sz_ms(2*m-1)
               sum_p= sum_p + Sz_ms(2*m)
            end do
            Sz_s(1)=sum_m
            Sz_s(2)=sum_p
            !	The array Sz_m contains the sum over the two sigma values of the 
            !	pseudospin of each orbital (i.e. the sum of the two neighbouring bits)
            !	====>	Sz_m= (sum_s Sz_ms)
            !	The array Sz_s contains the sum over the orbitals of the pseudospin with
            !	with s=-, s=+  (i.e. the sum of the bits in the even (sigma=-), and odd
            !	(sigma=+) boxes)
            !	====>	Sz_s= (sum_m Sz_ms)	
            !	Hint_J
            do m=1,N_orb 
               sum2= sum2 + (0.5d0*Sz_m(m))**2.d0
            end do
            !	Hint_J/2
            do s=1,2 
               sum3= sum3 + (0.5d0*Sz_s(s))**2.d0
            end do
            !	All
            sum=sum0 + 0.5d0*Up*sum1 + JH*sum2 - 0.5d0*JH*sum3
            !off-diagonal elements
         else if(flips.eq.1) then 
            !	This loop locate the position of the only one TRUE bit of the integer 
            !	IEOR()
            do i=0,imax
               if(BTEST(IEOR(j-1,k-1),i)) then
                  flag= i
               end if
            end do
            !	H0_Sx
            sum=(1.d0+c_ms(flag+1))*h_ms(flag+1)
            !	NB: the index of the arrays h_ms, c_ms are taken +1 with respect to the 
            !	flag position since the bits string start from i=0, while the arrays start
            !	from l=1 
         else
            sum=0.d0
         end if
         Ham(j,k)=sum
      end do ! k loop
   end do ! j loop

   !	%% CLOSING %%	
   deallocate (Sz_ms, Sz_m)
   return
 end subroutine Hambuilder

 !	Diagonalization	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 !	DIAHAM is the Fortran90 interface for calling the LAPACK diagonalization
 !	subroutine DSYEV.    
 !	Use as input the dimension n of the matrix Ham 
 !	Ham(n,n) = REAL SYMMETRIC matrix to be diagonalized
 !	output: m1(n,n) = orthonormal eigenvectors of Ham          
 !	        eig(n) = eigenvalues of a in ascending order
 !	the ground state (i.e. m1(n,1)) is saved in gs(n)  
 subroutine  diaham(n, Ham, gs)
   !	%% DECLARATIONS %%
   implicit none 
   !	Input
   integer :: n
   real*8, dimension (n,n) :: Ham
   real*8, DIMENSION(n) :: gs

   !	Indices
   integer :: i,l,inf 
   real*8 :: work(n*(3+n/2))

   !	Output
   real*8, allocatable :: m1(:,:),eig(:)
   real*8, allocatable ::	m2(:,:) 
   !	we need m2 only for the diagonalization check 
   !	NB: When added in the iterative procedure it will be eliminated

   !	%% ALLOCATIONS %%
   allocate (m1(n,n))
   allocate (m2(n,n))
   allocate (eig(n))
   m1(:,:)=Ham(:,:)


   !	%% DSYEV CALL %%
   l=n*(3+n/2)
   call dsyev('V','U',n,m1,n,eig,work,l,inf)

   !	%% GS SELECTION %%
   !	Save the eigenvector with the lowest eigen value
   gs(:) = m1(:,1)

   !	%% CHECK OPERATION %%
   !	Save the output in file dia.dat':
   !	*) eigenvalues                                                          
   !	*) eigenvectors (diagonalizing matrix [D])                              
   !	*) the original matrix [M] transformed by [D]; [1/D][M][D] dys.
   !	NB: When added in the iterative procedure this part will be eliminated.

   !	open(20,file='dia.dat')
   !	write(20,*)'Eigenvalues:'
   !	do i=1,n
   !	write(20,10)i,eig(i)
   !	10 format(I3,'   ',f14.8)
   !	end do
   !	write(20,*)
   !	write(20,*)'Eigenvectors:'
   !	do i=1,n
   !	write(20,20)i,m1(:,i)
   !	20 format(i3,'   ',10f14.8)
   !	end do
   !	write(20,*)

   !	m2=matmul(transpose(m1),Ham)
   !	Ham=matmul(m2,m1)

   !	write(20,*)'Transformed matrix (check):'
   !	do i=1,n
   !	write(20,30)Ham(:,i)
   !	30 format(10f14.8)
   !	enddo
   !	write(20,*)
   !	close(20)

   !	%% CLOSING %%
   deallocate(m1); deallocate(m2);  deallocate(eig)
   return
 end subroutine diaham


 function IBCHNG(int,pos) result(jnt)
   integer :: int,jnt
   integer :: pos
   logical :: bool
   bool = btest(int,pos)
   select case (bool)
   case (.true.)                !1
      jnt = ibclr(int,pos)
   case (.false.)               !0
      jnt = ibset(int,pos)
   end select
 end function IBCHNG

 
!	<Sz(O)_ms>_GS	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
!	Our ground state can be written as |GS> = sum_j gs(j) |j> 
!	Thus we have to compute 
!	<GS|Sz_ms|GS> = 
!	sum_{jk} gs(j)*gs(k) <k-1|Sz_ms|j-1> = sum_j gs(j)**2 Sz_ms|_{j-1}
!	where Sz_ms|_{j-1} = +/- (1/2)
!	and
!	<GS|O_ms|GS> = 
!	sum_{jk} gs(k)*gs(j) <k-1|O_ms|j-1> = sum_j gs(j)gs(k*) (1/2 or c/2) 
!	where the state |k*> is simply |j-1> with the spin Sz_ms flipped
!	the presence of the prefactor c depends if we flip a spin + or a spin - 
!	<GS|O^T_ms|GS> = 
!	sum_{jk} gs(k)*gs(j) <k-1|O^T_ms|j-1> = sum_j gs(j)gs(k*) (c/2 or 1/2) 
!	where the state |k*> is simply |j-1> with the spin Sz_ms flipped
!	the presence of the prefactor c depends if we flip a spin - or a spin + 
 subroutine av_s(N_pspins,N_states, gs,sz,o,oT,c)
   implicit none
   interface
      function ibchng(int,pos)
        integer :: int
        integer :: pos
        integer :: ibchng
      end function ibchng
   end interface
   integer :: N_pspins, N_states
   real*8, DIMENSION(N_states) :: gs
   real*8, DIMENSION(N_pspins) :: sz, o, oT, c
   integer :: spincode!, b
   logical :: b
   integer :: j,k,i,imax
   real*8 :: sumo, sumoT, sumz

   !	For each pspin (i=0, imax) we compute the sum over the states |j-1>
   imax=N_pspins-1
   do i=0,imax
      sumz=0.d0
      sumo=0.d0
      sumoT=0.d0
      do j=1,N_states ! run over the base state |j-1> 
         b=BTEST(j-1,i) 
         !			save the the state of the i-th pspin of the integer j-1			
         k = IBCHNG(j-1,i) + 1 
         !			select the state |k-1> as the state with the Sz_i flipped with 
         !			respect to the state |j-1> 
         sumz= sumz + 0.5*spincode(b)*gs(j)**2.d0
         if(b)then
            sumo= sumo + gs(j)*gs(k)
            sumoT= sumoT + c(i+1)*gs(j)*gs(k)
         else 
            sumo= sumo + c(i+1)*gs(j)*gs(k)
            sumoT= sumoT + gs(j)*gs(k)	
         end if
      end do ! j
      !		The averages Sx, O, O^T computed for each i (i.e. ms) are saved in the 
      !		i-th position of the sz, o, oT arrays 
      sz(i+1) = sumz 
      o(i+1) = sumo
      oT(i+1) = sumoT	 	
   end do  ! i (over the bit of the integer j)

   return
 end subroutine av_s

!	CORRELATIONS	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

 subroutine correlation(N_pspins,N_states, gs, sz, corr_1d1d, corr_1d1u, corr_1d2d, corr_1d2u)
   implicit none
   integer :: N_pspins, N_states
   real*8, DIMENSION(N_states) :: gs
   real*8, DIMENSION(N_pspins) :: sz
   integer :: spincode
   logical :: a,b,c,d
   integer :: j,p,k,l, i,imax
   real*8 :: sum1d, sum1u,sum2d, sum2u, sumconst 
   real*8 :: sum1d1d, sum1d1u, sum1d2d, sum1d2u
   real*8 :: corr_1d1d, corr_1d1u, corr_1d2d, corr_1d2u
   sum1d= 0.d0
   sum1u= 0.d0
   sum2d= 0.d0
   sum2u= 0.d0
   sum1d1d= 0.d0	
   sum1d1u= 0.d0
   sum1d2d= 0.d0
   sum1d2u= 0.d0	
   sumconst=0.d0
   i=0 !1d
   k=1 !1u
   p=2 !2d
   l=3 !2u
   do j=1,N_states ! run summing over the base state |j-1> 
      a=BTEST(j-1,i) 
      b=BTEST(j-1,k) 
      c=BTEST(j-1,p)
      d=BTEST(j-1,l) 		
      !		save the the state of the i-th pspin of the integer j-1			
      sum1d= sum1d + 0.5d0*spincode(a)*gs(j)**2.d0
      sum1u= sum1u + 0.5d0*spincode(b)*gs(j)**2.d0
      sum2d= sum2d + 0.5d0*spincode(c)*gs(j)**2.d0
      sum2u= sum2u + 0.5d0*spincode(d)*gs(j)**2.d0		
      sum1d1d= sum1d1d + 0.5d0*spincode(a)*0.5d0*spincode(a)*gs(j)**2.d0	
      sum1d1u= sum1d1u + 0.5d0*spincode(a)*0.5d0*spincode(b)*gs(j)**2.d0
      sum1d2d= sum1d2d + 0.5d0*spincode(a)*0.5d0*spincode(c)*gs(j)**2.d0
      sum1d2u= sum1d2u + 0.5d0*spincode(a)*0.5d0*spincode(d)*gs(j)**2.d0		
      !		sumconst= sumconst + 0.25d0*gs(j)**2.d0		
   end do ! j 
   corr_1d1d = 0.5d0*sum1d + 0.5d0*sum1d+ sum1d1d+0.25d0! sumconst
   corr_1d1u = 0.5d0*sum1d + 0.5d0*sum1u+ sum1d1u+0.25d0! sumconst	
   corr_1d2d = 0.5d0*sum1d + 0.5d0*sum2d+ sum1d2d+0.25d0! sumconst
   corr_1d2u = 0.5d0*sum1d + 0.5d0*sum2u+ sum1d2u+0.25d0! sumconst	

   return
 end subroutine correlation


!	%%%%%%%%%%%%%%%%%%%%%%%%%%%    FUNCTIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
!	Spin Code	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	One bit function. It converts a logical variable b=(T,F) in a spin state 
!	b(+1,-1)  
 integer function spincode(b)
   implicit none 
   logical :: b
   if(b) then
      !	T <-> 1 Pseudospin UP
      spincode =  1
   else 
      !	F <-> -1 Pseudospin DOWN
      spincode= - 1 
   end if
   return 
 end function spincode

 !	Flip Counter	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !	Compares the bits of the integers n, m and counts the number of the flip 
 !	bits (i.e. pseudospins). In order to do this we execute the exclusive or 
 !	(IEOR) between n, m and then sum the bits' values. 
 !	NB: the IEOR is an integer whose bits string has FAlSE every bit except 
 !	the i-postion where the n, m integers have a different bit value.	
 integer function Fcounter(n,m,imax)
   implicit none 
   integer n,m,imax,i,sum
   sum=0
   do i=0,imax
      sum= sum + TRANSFER(BTEST(IEOR(n,m),i),0)
   end do
   Fcounter=sum
   return
 end function Fcounter

