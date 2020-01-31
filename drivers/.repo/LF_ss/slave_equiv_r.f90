!	This is the core of the iterative procedure of the Slave-Spin problem 
!	within the single-site MF approximation.
!	In contains two main subroutines one solves the MF fermion Hamiltonian, 
!	the other solve the problem for the pseudospins. 
 
!	STRUCTURE:
!	Convergence Issue
!	*) Exit from the pseudospins loop
!	We want to exit from the pseudospins loop as soon as the constraint
!	conditions has been fulfilled i.e as as soon as for a given solution of the 
!	fermions' problem, nf_ms, the set of Lambda_ms has been arranged in the 
!	pseudospin loop in order to gives us a solution Sz_ms / 
!	nf_ms(i) = Sz_ms(i)  + 1/2  for each component of the arrays. 
!	In order to guaranteed the likeness of the two arrays we need to define a 
!	proper norm (a scalar quantity) of the difference array
!	ndiff(i)=nf_ms(i)-(Sz_ms(i) + 1/2)   and to require that this norm is small 
!	(i.e. smaller than an epsilon)
!	NB: Several choice of this norm are possible (and in principle equivalent). 
!		We'll test:
!		1)The modulus 		
!			norm_ndiff = sum_i abs(ndiff_ms(i))
!		2)Euclidean norm	
!			norm_ndiff = sqrt(sum_i ndiff_ms(i)**2.d0)
!		2)The square of the euclidean norm	
!			norm_ndiff = sum_i ndiff_ms(i)**2.d0)
!	Another possibility is to check this difference component by component 
!	without the use of the norm. NB: this check is requiring always the same 
!	level of precision also increasing the number of orbitals. 
!	*) Exit from the fermions loop
!	We want to exit from the fermions loop as soon as the set of Lambda_ms 
!	obtained from the pesudospins loop interations no longer changes
!
!	OUTPUT:
!	Final data are nf_ms, Lambda_ms(l), Z_ms(l) with l=1, N_slave . 
!	They are saved in files organized as: 
!	*) nf_U#_alpha#_mix#_refresh.dat
!			nf_ms(l)	inters_f
!	*) Lambda_U#_alpha#_mix#_refresh.dat
!			Lambda_ms(l)	inters_f	iters_PS	iters_tot
!	*) Z_U#_alpha#_mix#_refresh.dat
!			Z_ms(l)	inters_f	iters_PS	iters_tot
!	*) Sz_U#_alpha#_mix#_refresh.dat
!			Sz_ms(l)	inters_f	iters_PS	iters_tot

	subroutine slave_equiv_r(U,Up,JH,rho,N_orb,N_slave,t_ms,Lambda_ms,Z_ms,nf_ms, corr_1d1d, corr_1d1u, corr_1d2d, corr_1d2u, alpha,mix,flagZ)	

!	%%%%%%%%%%%%%%%%%%%%%%%%%%    DECLARATIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%
	implicit none
!	MODEL DATA
!	The interaction couplings are given by:
	real*8 :: U ! intraorbital Hubbard term
	real*8 :: Up! interorbital Hubbard term
	real*8 :: JH! Hund coupling 
	real*8 :: rho
!	rho is the number of the electrons for site - fix the doping level  	
	integer :: N_orb, N_slave
!	N_orb is the number of the orbitals, N_slave=2*N_orb is the number of the 
!	slave variables we use (pseudospins and fermions f) 
	real*8, dimension(N_slave) :: t_ms, Lambda_ms, Z_ms
	real*8, dimension(N_slave) :: New_Lambda_ms, Old_Lambda_ms
	real*8, dimension(N_slave) :: New_Z_ms, Old_Z_ms		
!	t_ms, Lambda_ms, Z_ms, h_ms are arrays of dimension 2*N_orb, the pedix ms 
!	stays for the orbital and spin labels {m, sigma}

!	SELF-CONSISTENT ITERATIVE PROCEDURE VARIABLES
!	Pseudospin & Fermion Arrays 
	real*8, dimension(N_slave) :: nf_ms, Sz_ms
!	nf_ms, Sz_ms are arrays containing the solution of the fermion and 
!	pseudospin loop respectively. We use these quantities in order to check the
!	constraint condition <nf_ms> = <Sz_ms> + 1/2
	real*8, dimension(N_slave) :: h_ms, E_ms
	real*8, dimension(N_slave) :: O_ms, OT_ms, c_ms
!	h_ms contains the values of the effective magnetic field for the
!	pseudospin problem. 
!	E_m saves the average value of E*nf on the gs for the fermion problem.
!	O_ms, OT_ms contain the average values on the GS of O_ms, transpose(O_ms) 
!	respectively.
!	c_ms is the array containing [sqrt(nf_ms*(1-nf_ms))**(-1)] - 1 
	real*8, dimension(N_orb) :: nf_m, Ef_m, Lambda_m, Z_m
	real*8 :: corr_1d1d, corr_1d1u, corr_1d2d, corr_1d2u
!	correlators		
!	Auxiliar variables	
	real*8 :: sumeq_Lambda, sumeq_Z,sumeq_n, sumeq_E, sumeq_O, sumeq_OT, sumeq_Sz 

!	Convergency parameters
	real*8, parameter :: eps_PS= 1.d-6, eps_f= 1.d-4
!	the exit condition for the PSpin loop is norm_ndiff < eps_PS
!	the exit condition for the fermion loop is 
!	norm_DZ < eps_f + norm_DLambda < eps_f
!	*)Auxiliary for the pseudospins loop
	real*8, dimension(N_slave) :: ndiff_ms 
!	Array of dimension N_slave=2*N_orb containing the differences between the 
!	output of the fermions loop nf_ms and the one of the pseudospins Sz_ms 
	real*8 :: norm_ndiff, sum_ndiff 
!	norm_ndiff is the norm that we define in order to check the fulfilling of 
!	constraint
	integer :: flag_ndiff
!	is the flag that we can use instead of norm_ndiff in order o check the 
!	fulfilling of constraint for EVERY component of ndiff 
	real*8 :: alpha
!	coefficient for the Lambda changing procedure
!	Lambda_ms(i) = Lambda_ms(i) - alpha*ndiff_ms(i)
!	*)Auxiliary for the pseudospins loop
	real*8, dimension(N_slave) :: DLambda_ms, DZ_ms
	real*8 :: norm_DLambda, sum_DLambda
	real*8 :: norm_DZ, sum_DZ	
	real*8 :: mix 
!	mixing parameter for the Lambda vector 
!	Lambda_ms(i)= mix* New_Lambda_ms(i)  + (1-mix)*Old_Lambda_ms(i) 
!	*)Emergency exit parameters 
	integer, parameter :: itermax_PS = 5000, itermax_f = 5000
!	max_PS, Max_f are the maximum number of steps for the pseudospins' and the 
!	fermions' loop respectively  
	integer :: flagZ, flagCorr
	
!	Mute Indices, Counters & Labels
	integer :: iters_PS, iters_f, true_iters_f, iters_tot,l, i, m
	character (len=200) :: NumOrb	
	character (len=200) :: Hubb
	character (len=200) :: HundU
	character (len=200) :: sizechange
	character (len=200) :: mixing
	character (len=200) :: NumEl		

!	%%%%%%%%%%%%%%%%%%%%%%%%%%    OPENING FILES    %%%%%%%%%%%%%%%%%%%%%%%%%
	write(Hubb,'( f8.4 )') U
	write(HundU,'( f8.4 )') JH/U	
	write(sizechange,'( f8.2 )') alpha
	write(mixing,'( f8.2 )') mix
	write(NumOrb,'( I5 )') N_orb
	write(NumEl,'( f8.2 )') rho	
	
	!open(unit=30,file='data_output/noconv_n'//trim(adjustl(NumEl))//'_'//trim(adjustl(NumOrb))//'orb_equiv_r.dat')
!	We save here the points in which the program didn't find the solution	 
!	open(unit=31,file='data_output/Mott_n'//trim(adjustl(NumEl))//'_'//trim(adjustl(NumOrb))//'orb_equiv_r.dat') 
!	We save here the points in which all the orbitals are below the treshold			 	

!	open(40,file='data_output/iterations/nf_'//trim(adjustl(NumOrb))//'orb_U'//trim(adjustl(Hubb))//'_JHU'//trim(adjustl(HundU))//'_alpha'//trim(adjustl(sizechange))//'_mix'//trim(adjustl(mixing))//'_r.dat')          
!	open(41,file='data_output/iterations/Lambda_'//trim(adjustl(NumOrb))//'orb_U'//trim(adjustl(Hubb))//'_JHU'//trim(adjustl(HundU))//'_alpha'//trim(adjustl(sizechange))//'_mix'//trim(adjustl(mixing))//'_r.dat')
!	open(42,file='data_output/iterations/Z_'//trim(adjustl(NumOrb))//'orb_U'//trim(adjustl(Hubb))//'_JHU'//trim(adjustl(HundU))//'_alpha'//trim(adjustl(sizechange))//'_mix'//trim(adjustl(mixing))//'_r.dat')
!	open(43,file='data_output/iterations/Sz_'//trim(adjustl(NumOrb))//'orb_U'//trim(adjustl(Hubb))//'_JHU'//trim(adjustl(HundU))//'_alpha'//trim(adjustl(sizechange))//'_mix'//trim(adjustl(mixing))//'_r.dat')
!	open(44,file='data_output/iterations/h_'//trim(adjustl(NumOrb))//'orb_U'//trim(adjustl(Hubb))//'_JHU'//trim(adjustl(HundU))//'_alpha'//trim(adjustl(sizechange))//'_mix'//trim(adjustl(mixing))//'_r.dat')
	
!	%%%%%%%%%%%%%%%%%%%%%    SELF-CONSISTENT PROCEDURE   %%%%%%%%%%%%%%%%%%%	
	norm_DLambda= 1.d0
 	norm_DZ= 1.d0
	flagZ=0
	flagCorr=0
	iters_tot=0
	iters_f =0
	do i=1,N_slave
		nf_ms(i)=0.d0
		c_ms(i)=0.d0 
		E_ms(i)=0.d0 
	end do !i-loop 	

!	write(*,*)'here in slave', Z_ms(1), Lambda_ms(1) 
	call fermion_equiv(rho,N_orb,N_slave,t_ms,Lambda_ms,Z_ms,nf_ms,E_ms,c_ms)
 !                      nb: in principle this passage is useless since is already done in the fermion subroutine
                        sumeq_n=0.d0
                        sumeq_Z=0.d0
                        sumeq_E=0.d0
                        do i=1,N_slave
                                sumeq_n=sumeq_n+nf_ms(i)
                                sumeq_Z=sumeq_Z+Z_ms(i)
                                sumeq_E=sumeq_E+E_ms(i)
                        end do
                        do i=1,N_slave
                                nf_ms(i)=sumeq_n/N_slave
                                Z_ms(i)=sumeq_Z/N_slave
                                E_ms(i)=sumeq_E/N_slave
                        end do
	


!	Fermions cycle
	do while (iters_f < itermax_f .and. (norm_DLambda>eps_f .or. norm_DZ>eps_f)) 
	iters_f=iters_f+1
	true_iters_f = iters_f 	
 !		Save the initial Lambda_ms set for the final mixing with the final one 
		do i=1,N_slave
			Old_Lambda_ms(i)=Lambda_ms(i) 
			Old_Z_ms(i)=Z_ms(i) 			
		end do !i-loop 

!		call fermion_equiv(rho,N_orb,N_slave,t_ms,Lambda_ms,Z_ms,nf_ms,E_ms,c_ms)
!		Output of fermion subroutine nf_ms, E_ms
!		Saving data sets
!		write(40,12) iters_f, nf_ms(:) 
!		write(*,12) iters_f, nf_ms(:)
!		write(*,12) iters_f , c_ms(:)		

!		Computation of the input variables for the pspin routine
		do l=1,N_slave
			O_ms(l)= sqrt(Z_ms(l))
	 		h_ms(l)=O_ms(l)*E_ms(l)
		end do 	

		norm_ndiff= 1.d0
		flag_ndiff = 0  
		iters_PS =0	
		
!		Pspins cycle
!		do while (iters_PS < itermax_PS .and. norm_ndiff>eps_PS) 
		do while (iters_PS < itermax_PS .and. flag_ndiff<N_slave) 
			iters_PS=iters_PS+1
			iters_tot= iters_tot + 1

!			At each iteration refresh the Weiss fields value with the Sx_ms 
!			computed in the previous pspin run 
			do l=1,N_slave
	 			h_ms(l)=O_ms(l)*E_ms(l)
			end do
!			Solver of the Pspin Hamiltonian 		
			call pspin_corr(U,Up,JH,N_orb,N_slave,Lambda_ms,h_ms,O_ms,OT_ms,c_ms,Sz_ms,flagCorr, corr_1d1d, corr_1d1u, corr_1d2d,corr_1d2u) 
			
!			ALL EQUIVALENT ORBITALS - Arrange the spins
			sumeq_Sz=0
                        sumeq_O=0
                        sumeq_OT=0
			do i=1,N_slave
				sumeq_Sz=sumeq_Sz+Sz_ms(i)
                                sumeq_O=sumeq_O+O_ms(i)
                                sumeq_OT=sumeq_OT+OT_ms(i)

                       end do
			do i=1,N_slave
				Sz_ms(i)=sumeq_Sz/N_slave
                                O_ms(i)=sumeq_O/N_slave
                                OT_ms(i)=sumeq_OT/N_slave
			end do			
				
!			Constraint Condition
			sum_ndiff=0.d0
			do i=1,N_slave
				ndiff_ms(i) = nf_ms(i) - (Sz_ms(i) + 0.5d0) 
!				Arranging the set of Lambda_ms according to the difference 
!				vector ndiff_ms
!				Lambda_ms(i) = Lambda_ms(i) - alpha*ndiff_ms(i)! this is the simple one	
				Lambda_ms(i) = Lambda_ms(i)*(1-alpha*ndiff_ms(i))! scale-free update
!				Compute the norm
				sum_ndiff = sum_ndiff + abs(ndiff_ms(i))
!				sum_ndiff = sum_ndiff + ndiff_ms(i)**2.d0
!				Check component by component updating flag_ndiff
				if(abs(ndiff_ms(i))<eps_PS)then
					flag_ndiff =flag_ndiff +1
				end if	
			end do !i-loop 
			norm_ndiff=	sum_ndiff	
			
!			Saving data sets
!			write(41,11) iters_f, iters_PS, iters_tot, Lambda_ms(:)
!			write(42,11) iters_f, iters_PS, iters_tot, Z_ms(:)
!			write(43,11) iters_f, iters_PS, iters_tot, Sz_ms(:)
!			write(44,11) iters_f, iters_PS, iters_tot, h_ms(:)
		end do!  pseudospin loop 

!               Let us require the state not to be polarized - PARAMAGNETIC PHASE                                                                 !               do m=1,N_orb                                                                                                                      
!                       Lambda_m(m)= (Lambda_ms(m*2-1)  + Lambda_ms(m*2))/2.d0                                                                    
!                       Lambda_ms(m*2)=Lambda_m(m)                                                                                                
!                       Lambda_ms(m*2-1)=Lambda_m(m)                                                                                              
!                       Z_m(m)= (Z_ms(m*2-1)  + Z_ms(m*2))/2.d0                                                                                   
!                       Z_ms(m*2)=Z_m(m)                                                                                                          
!                       Z_ms(m*2-1)=Z_m(m)                                                                                                        
!               end do                                                                                                                            
!               ALL EQUIVALENT ORBITALS                                                                                                            
                sumeq_Lambda=0.d0
                sumeq_Z=0.d0
                do l=1,N_slave
                        sumeq_Lambda = sumeq_Lambda + Lambda_ms(l)
                        sumeq_Z = sumeq_Z + Z_ms(l)
                     end do
                do l=1,N_slave
                        Lambda_ms(l) =sumeq_Lambda/N_slave
                        Z_ms(l)=sumeq_Z/N_slave
                end do

!		Convergence Test 
		sum_DLambda=0.d0
		sum_DZ=0.d0
		!do i=1,N_slave
			DLambda_ms(:) = Old_Lambda_ms(:) - Lambda_ms(:) 
			DZ_ms(:) = Z_ms(:) - O_ms(:)*OT_ms(:) 
!			Compute the norms of the differences' vectors
!  			sum_DLambda = 	sum_DLambda + abs(DLambda_ms(1))*0.5d0
                        sum_DLambda = abs(DLambda_ms(1))
!  			sum_DZ = 	sum_DZ + abs(DZ_ms(1))
 		        sum_DZ = abs(DZ_ms(1))
		!end do !i-loop 
		norm_DLambda=sum_DLambda
		norm_DZ=sum_DZ	
!		write(*,*)iters_f,iters_PS
!		write(*,*) Z_ms(1) , O_ms(1)*OT_ms(1), norm_DZ
!		write(*,*) Old_Lambda_ms(1),Lambda_ms(1), norm_DLambda
!		write(*,*) nf_ms(1)
		
!		Prepare the input for the fermion call:
!		*) Put limiter if Lambda_ms go beyond the precision 
!		*) Check the Z_ms values and go out if all the components are very small
!		   otherwise problem in the pspin routine could appear due to the 
!		   vanishing of the non interacting part of the pspin hamiltonian
!		*) If we didn't find the PS solution do not refresh Lambda Z, go out
!		   from both the routines (fermion and pseudospin) and save this 
!		   problematic point in the convergence file, otherwise compute the new 
!		  weight Z_ms and the new Lambda using a mixing of the new sets and the
!		  original one Old_Lambda_ms
!		Convergency Issue    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		do i=1,N_slave	
!			New values come out from the pspin routine		
			New_Z_ms(i) = O_ms(i)*OT_ms(i)
			New_Lambda_ms(i)= Lambda_ms(i)
!			Lambda Limiter	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(abs(New_Lambda_ms(i))< 1.d-3)then 
				New_Lambda_ms(i)=0.0d0 
			end if			
!			This avoid the unphysical fluctuation (different for each l) 
!			of the Lambda that could appear beyond the level of precision 
!			we are considering here.			  
!			Z check	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(New_Z_ms(i)< 0.2d-3)then 
				flagZ = flagZ + 1
			end if
			if (flagZ== N_slave)then
!				write(*, 21) U, JH, Z_ms(:)
				write(31, 21) U, JH, New_Z_ms(:)				
				true_iters_f = iters_f !saves the true values of iters_f
				iters_f=2000 !prevent others fermion-cycle		
			end if	
!			MOTT phase - In this case it make no sense go ahead 
		end do !i-loop 
!		Z L update	  		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		if(iters_PS==2000)then 
!			rejects changes if the run didn't reach the convergence
			write(30, *)'No convergence found after 2000 PSsteps at fermion step', true_iters_f, flagZ 
			write(30,21) U, JH , New_Lambda_ms(:)
!			probably the problem arise from strange values of Lambda 			
			write(30,*)'  '		
			do i=1,N_slave
				New_Z_ms(i)= Old_Z_ms(i)
				New_Lambda_ms(i) = Old_Lambda_ms(i)
			end do 
			true_iters_f = iters_f !saves the true values of iters_f			
			iters_f=2000 !prevent others fermion-cycle	
		else 
!			standard update 		
			do i=1,N_slave	
			Z_ms(i)= mix*New_Z_ms(i)  + (1-mix)*Old_Z_ms(i) 		  
			Lambda_ms(i)= mix*New_Lambda_ms(i)  + (1-mix)*Old_Lambda_ms(i) 
			end do !i-loop 
		end if		
	end do! fermion loop 

!	If we exit form the previous loop we don't need to mix again the Lambda	
	flagCorr=1
!	Last call to compute the correlation functions	
	call pspin_corr(U,Up,JH,N_orb,N_slave,Lambda_ms,h_ms,O_ms,OT_ms,c_ms,Sz_ms,flagCorr, corr_1d1d, corr_1d1u, corr_1d2d,corr_1d2u) 
	do i=1,N_slave
		Lambda_ms(i)= New_Lambda_ms(i)
		Z_ms(i)= New_Z_ms(i)
	end do !i-loop 	
!	close(40)
!	close(41)
!	close(42)
!	close(43)
!	close(44)	
	
	11   format(I5,I5,I5, 11(1x,E12.6,2x))	
	12   format(I5, 11(1x,E12.6,2x))		
	21   format(20(1x,E12.6,2x))	
	end subroutine slave_equiv_r



