!	BASIC INFORMATION
!	The program gives us the mapping of the slave solution while changing the 
!	coupling interaction U, JH (Up = U-2JH). For a given value of (U, JH) the 
!	program calls the slave routine that solves self consistently the problem 
!	for the two Hamiltonians (the fermion and the pseudospin one)  
!	The pedix  _r of the slave subroutine means that at each call of the pspin
!	routine we are refreshing the set of Weiss fields h_ms using the Sx_ms 
!	computed ad the previous call (see slave_r for details).  

!	We restrict our the phase of space JH/U < 0.33
!	The run over the coupling starts from small values of interactions so that 
!	one can put as guesses of the Lagrange multipliers and quasiparticle weights
!	values close to the metallic solution. As the coupling increase (U+dU, 
!	JH+dJH) we use as input for the new run the self-consistent solution found 
!	in the previous step in the coupling runs but the cases in which we entered 
!	the Mott phases (i.e. all q.ple weights vanished), in this case the program 
!	restart from the output of the first step computed. 

!	INSTRUCTIONS FOR COMPILATION
!	To run the program for different number of orbitals one have to change
!	only the N_orb parameter here and create in the folder data_input the 
!	hoppings file according the convention we generally use:
!	t_{1-}	
!	t_{1+}	
!	t_{2-}	
!	t_{2+}	
!	.....
!	t_{N_orb+}	

!	INPUT - OUTPUT
!	Input data:  hopping_Norb.dat in data_input
!	Output data: Time of each map point, Info about convergence and Mott phase, 
!	the values of Lambda_ms(l), Z_ms(l), nf_ms(l) with l=1, N_slave are saved in 
!	the folder data_output. For details read List_output.txt

!	TIMING
!	The timing of the job is computed using the intrinsic fortran function
!	date_and_time. This subroutine is for me a black box simpling marking 
!	the time at every call, so that at the end it is possible to know the 
!	duration of a job by difference: finish_time - start_time. The results 
!	is in seconds.

!	RANGE OF VALIDITY 
!	The present version of the subroutines works on a 2 dimensional 
!	N_orb-system with hoppings at first nearest neighbors. 
!	Orbital hybridization is not considered. Rotational invariance is assumed 
!	and used in the definition Up=U-2JH. This version of the program is already 
!	written in a complete general way for running with arbitrary N_orb and 
!	arbitrary rho filling.  

!	GENERALIZATIONS
!	1) General hoppings:
!	1a) Anisotropy between x(y)-hoppings, second NN ...:
!	* change of the input file organization and as a consequence of its reading 
!	in Main_UJ.f90 
!	* in fermion.f90 the computation of E_k = 2*t(cos kx+ cos ky ) has to be 
!	computed explicitly via E_k = sum_d exp(ikd) t_d  with d= i-j 
!	1b) Hybridization: 
!	* in fermion.f90 diagonalization of the fermion problem is needed
!	2) Rotational Invariance Broken: 
!	* in pspin.f90 the Hambuilder subroutine has to be modify since the 
!	interacting part of the pspin Hamiltonian has another analithical expression 
	
	program map

!	%%%%%%%%%%%%%%%%%%%%%%%%%%    DECLARATIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%
	implicit none
!	INPUT DATA
	integer, parameter :: N_orb=3
	real*8, parameter :: rho_el=3.00d0
!	rho_el = #  is the number of the electrons for site - fix the doping level 	
!	NB: #= 0 - 2*N_orb empty-complete filled cases	
	integer :: N_slave
!	N_orb is the number of the orbitals, N_slave=2*N_orb is the number of the 
!	slave variables we use (pseudospins and fermions f) 
	real*8, allocatable, dimension(:) :: t_ms
!	hopping parameters array of dimension N_slave=2*N_orb(only intraorbital
!	hoppings)
	real*8, allocatable, dimension(:) :: Lambda_ms, Z_ms, nf_ms
!	Lambda_ms, Z_ms, h_ms are arrays of dimension N_slave=2*N_orb
	real*8 :: corr_1d1d, corr_1d1u, corr_1d2d, corr_1d2u
!	correlators	
	real*8, allocatable, dimension(:) :: FirstLambda_ms, FirstZ_ms
	real*8, allocatable, dimension(:) :: BackLambda_ms, BackZ_ms
!	Back* save the initial guess of Lambda and Z. Every map point restart from 
!	these values.
!	NB: The the pedix ms of the arrays stays for the orbital and spin labels 
!	{m, sigma}
	
!	U, J MAPPING VARIABLES 
!	The interaction couplings are given by:
	real*8 :: U ! intraorbital Hubbard term
	real*8 :: Up! interorbital Hubbard term
	real*8 :: JH! Hund coupling 
	real*8 :: JH_U! Hund coupling in U unit
	real*8 , parameter ::  Umin=0.0d0
!	real*8 , parameter ::  Umin=10.0d0	
!	We start from the weak limit (metallic side) for which we can 
!	easily guess reasonable starting parameters Lambda_ms and Z_ms
	real*8, parameter :: dU=0.5d0 , dJH_U=0.005d0 
!	real*8, parameter :: dU=1.d0 , dJH_U=0.01d0! dJH_U=0.006d0 
!	thickness of the steps 
!	integer, parameter :: nstepsJH_U= 34, nstepsU=29!nstepsU=55
	integer, parameter :: nstepsJH_U=67,  nstepsU=57 
!	we span the space JH/U = 0 - 0.33
	integer :: nU, nJH
!	mute indices for the cycles

!	Convergency procedure parameters
!	*)Auxiliary for the pseudospins loop
	real*8 :: alpha
!	fix the stepsize for the Lambda changing procedure inside the slave routine
!	*)Auxiliary for the pseudospins-fermions output-input 
	real*8 :: mix
!	fix the degree of mixing parameter for the Lambda vector 
!	NB: Probably these values should be optimized depending on U, J, N_orb
!	(for detailes about the convergency parameters see the comments within the
!	slave routine)
	integer :: flagUc,flagZ
	integer :: step_mott 
!	auxiliary parameters that highlight the proximity of the Mott-phases at JH=0 

!	MUTE INDICES & LABELS
!	Counter
	integer :: l, i
!	Timing Variables
!	Variables for subroutine date_and_time
	integer :: time_array_0(8), time_array_1(8)
	real*8 :: start_time, finish_time
!	Variables for computation
	integer, parameter :: n = 1000000
	double precision :: a(n), b(n), c(n)	
!	Naming Varible
	character (len=200) :: NumOrb
	character (len=200) :: NumEl	

!	%%%%%%%%%%%%%%%%%%%%%%%%%%    OPENING FILES    %%%%%%%%%%%%%%%%%%%%%%%%%
	write(NumOrb,'( I5 )') N_orb
	write(NumEl,'( f8.2 )') rho_el	
!	input	
	open(40,file='data_input/hopping_'//trim(adjustl(NumOrb))//'orb.dat')  
!	output		
	open(45,file='data_output/Zmap_n'//trim(adjustl(NumEl))//'_'//trim(adjustl(NumOrb))//'orb.dat') 
	open(46,file='data_output/Lmap_n'//trim(adjustl(NumEl))//'_'//trim(adjustl(NumOrb))//'orb.dat') 
	open(47,file='data_output/nfmap_n'//trim(adjustl(NumEl))//'_'//trim(adjustl(NumOrb))//'orb.dat')	
	open(60,file='data_output/time_n'//trim(adjustl(NumEl))//'_'//trim(adjustl(NumOrb))//'orb.dat')
	open(90,file='data_output/correlators_n'//trim(adjustl(NumEl))//'_'//trim(adjustl(NumOrb))//'orb.dat') 
		 

!	%%%%%%%%%%%%%%%%%%%%%%    ASSIGNMENT & ALLOCATION    %%%%%%%%%%%%%%%%%%%
!	Input Informations:
!	*) Defines N_slave
!	*) Allocates memory for t_ms, Lambda_ms, Z_ms 
!	*) Acquires the hoppings 
!	*) Defines the starting guesses for Lambda_ms, Z_ms
	N_slave=2*N_orb
	allocate(t_ms(N_slave))
	allocate(Lambda_ms(N_slave))
	allocate(Z_ms(N_slave))
	allocate(nf_ms(N_slave))	
	allocate(BackLambda_ms(N_slave))
	allocate(BackZ_ms(N_slave))
	allocate(FirstLambda_ms(N_slave))
	allocate(FirstZ_ms(N_slave))	
	do l=1,N_slave
	    read(40,*) t_ms(l)
		BackLambda_ms(l)=0.1d-3 
		BackZ_ms(l)=0.999d0
		Lambda_ms(l)=BackLambda_ms(l)
		Z_ms(l)=BackZ_ms(l)
		FirstLambda_ms(l)=BackLambda_ms(l)
	  	FirstZ_ms(l)=BackZ_ms(l)		
	end do
	close(40)

!	CONVERGENCY ISSUE - IMPORTANT PARAMETERS :
!	Mixing Variables
	alpha=0.5d0
	mix=0.5d0

!	MAPPING (U,JH)     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!	Interactions couplings  U, U', J 
!	We run over U up to Uc + DU 
!	In the practice we change the indicator Mott=1 at the criticality and run 
!	for other 5 points of the map
!	NB: Notice that this version of the map runs its peculiar for this case in 
!	which we actually find the Mott phase
	flagUc = 0
	flagZ = 0	
	step_mott=0
	U=Umin
	do nU = 1, nstepsU 
		JH_U=0.d0
		do nJH = 1, nstepsJH_U 
			JH = JH_U*U 
			Up=U-2.d0*JH
!			write(*,*)'NUOVO CICLO %%%%%%%%%%%%%', U, JH_U	
!			Initial guesses
			if(flagZ==N_slave)then
				do l=1,N_slave
					Lambda_ms(l)=FirstLambda_ms(l)
	  				Z_ms(l)=FirstZ_ms(l)
	  			end do
			else
!			Recovery of the output  of the previous nU step 
				if(nJH==1) then
					do l=1,N_slave
	  					Lambda_ms(l)=BackLambda_ms(l)
	  					Z_ms(l)=BackZ_ms(l)
					end do
				end if
			end if

!			Mark the beginning of the self-consisten procedure
			call date_and_time(values=time_array_0)
			start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 + time_array_0 (7) + 0.001 * time_array_0 (8)					

!			CALLING THE SELF CONSISTENT ROUTINE 
!			call slave_r(U,Up,JH,rho_el,N_orb,N_slave,t_ms,Lambda_ms,Z_ms,nf_ms,alpha,mix)	
!			h_ms is refreshed at every pspin loop	
!			call slave_para_r(U,Up,JH,rho_el,N_orb,N_slave,t_ms,Lambda_ms,Z_ms,nf_ms,alpha,mix)
!			forces paramagnetic solutions; h_ms is refreshed at every pspin loop	
			call slave_equiv_r(U,Up,JH,rho_el,N_orb,N_slave,t_ms,Lambda_ms,Z_ms,nf_ms, corr_1d1d, corr_1d1u, corr_1d2d,corr_1d2u, alpha,mix,flagZ)
!			forces equivalent solutions for all the orbitals; h_ms is refreshed at every pspin loop	
!			call slave_f(U,Up,JH,rho_el,N_orb,N_slave,t_ms,Lambda_ms,Z_ms,nf_ms,alpha,mix)	
!			h_ms is fixed in the pspin loop
!			call slave_para_f(U,Up,JH,rho_el,N_orb,N_slave,t_ms,Lambda_ms,Z_ms,nf_ms,alpha,mix)	
!			forces paramagnetic solutions; h_ms is fixed in the pspin loop
!			call slave_equiv_f(U,Up,JH,rho_el,N_orb,N_slave,t_ms,Lambda_ms,Z_ms,nf_ms,alpha,mix)
!			forces equivalent solutions for all the orbitals; h_ms is fixed in the pspin loop						

!			*) U, Up, JH change at every loop's call
!			*) rho_el, N_orb, N_slave are initial data of our problem. 
!			   They are fixed. 
!			*) t_ms is the hopping parameters vector. Is fixed in this loop. 
!			   NB: If we are interested in we could run this program for a 
!			   given set of t_ms in order to find critical values of U/JH and 
!			   repeat changing the t_ms set. 
!			*) Lambda_ms, Z_ms are output arrays of the slave subroutine 
!			   At each map point we use the initial guess values for Lambda_ms, 
!			   Z_ms.
!			*)nf_ms is the array that the fermion routine has to fill
!			*) alpha, mix are the convergency auxilairy parameters.  
!			*) flaZ is the smoking gun of the Mott transition
!			   Z_ms.

!			Mark the end of the self-consistent procedure
			call date_and_time(values=time_array_1)
			finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 + time_array_1 (7) + 0.001 * time_array_1 (8)	

!			Storage of the output for the first step in U, JH 
			if(U==Umin.and.nJH==1 ) then
				do l=1,N_slave
	  				FirstLambda_ms(l)=Lambda_ms(l)
	  				FirstZ_ms(l)=Z_ms(l)
				end do
			end if

!			Check if we are near Uc (i.e. put flagUc=N_slave if all the qp. 
!			spectral weights are smaller then 5.d-3)
			if(nJH==1) then
				do l=1,N_slave
	  				BackLambda_ms(l)= Lambda_ms(l)
	  				BackZ_ms(l)=Z_ms(l)				
				end do ! l = 1,N_slave	
			end if	! (nJH==1)

!			Save the results for each U, JH   
!			write(*,*) 'U=', U,'JH/U=',JH_U, 'JH=', JH				
			write(45,21) U, JH_U, Z_ms(:)
			write(46,21) U, JH_U, Lambda_ms(:)		
			write(47,21) U, JH_U, nf_ms(:)												
			write(60,21) U, JH_U, finish_time - start_time  
			write(90,21) U, JH_U, corr_1d1d, corr_1d1u, corr_1d2d, corr_1d2u
							

!			Closing loops		
			JH_U= JH_U + dJH_U
		end do !JH 
		U = U + dU 
		write(45,*) '   '
		write(46,*) '   '
		write(46,*) '   '		
		write(60,*) '   ' 	
		write(90,*) '   ' 				
	end do !U

	deallocate(t_ms,Lambda_ms,Z_ms,nf_ms,BackLambda_ms,BackZ_ms,FirstLambda_ms,FirstZ_ms)
	close(45)
	close(46)
	close(47)	
	close(60)
	close(90)	
	21   format(20(1x,E12.6,2x))	
	end program 






