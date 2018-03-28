	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Bulk Aerosol Module (BAM): 
	!>Solves for droplet activation using Abdul-Razzak et al. (1998) or 
	!> Fountoukis and Nenes (2004)
	!> compile using the Makefile and then run using: <br>
	!> ./main.exe  namelist.in
	!> <br><br>



	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main programme reads namelist, allocates arrays and calls ctmm_activation to do
	!> calculation
	
    program main
        use nrtype1
        use random, only : random_normal
        use sub, only : n_mode, n_sv, giant_flag, method_flag, sv_flag, &
        	n_aer1, d_aer1, sig_aer1, molw_core1, density_core1, nu_core1, org_content1, &
        	molw_org1, density_org1, delta_h_vap1, nu_org1, log_c_star1, p_test, t_test, &
    		w_test, act_frac1, &
    		a_eq_7, b_eq_7, &
    		allocate_arrays, ctmm_activation, initialise_arrays
        implicit none

        real(sp) :: w1,t1,p1
        
        character (len=200) :: nmlfile = ' '
        
		! for random number:
		real(sp) :: r, mean_w, sigma_w
		real(sp), dimension(10,10) :: rs
		integer(i4b), allocatable, dimension(:) :: seed
		integer(i4b) :: l, i, n_rand
		logical :: rand_dist=.false.




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! set up flags, etc
        namelist /bulk_aerosol_setup/ n_mode, n_sv, sv_flag, method_flag, giant_flag, &
        			a_eq_7, b_eq_7      
        ! run parameters / input arrays
        namelist /bulk_aerosol_spec/ n_aer1, d_aer1, sig_aer1, molw_core1, density_core1, &
        					nu_core1, org_content1, molw_org1, density_org1, &
        					delta_h_vap1, &
        					nu_org1, log_c_star1, p_test, t_test, w_test, &
        					rand_dist, n_rand, mean_w,sigma_w
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists	and allocate arrays								   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=bulk_aerosol_setup)
        ! allocate memory / init
		call allocate_arrays(n_mode,n_sv,n_aer1,d_aer1,sig_aer1, &
			molw_core1,density_core1,nu_core1,org_content1, &
			molw_org1, density_org1,delta_h_vap1,nu_org1,log_c_star1, &
			act_frac1)
        
        read(8,nml=bulk_aerosol_spec)
        close(8)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! initialise arrays								                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call initialise_arrays(n_mode,n_sv,p_test,t_test,w_test, &
					n_aer1,d_aer1,sig_aer1, molw_org1,density_core1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! call activation code												   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(.not. rand_dist) then
        	print *,'act_frac1, n_act'
			call ctmm_activation(n_mode,n_sv,sv_flag, &
						n_aer1, d_aer1,sig_aer1,molw_core1, &
						density_core1, nu_core1, &
						org_content1,molw_org1, density_org1, delta_h_vap1, nu_org1,  &
						log_c_star1, &
						w_test, t_test,p_test, a_eq_7, b_eq_7, &
						act_frac1)
			print *,act_frac1,  sum(act_frac1*n_aer1)
		else
			! random number generator
			call random_seed(size=l)
			allocate(seed(1:l))
			seed(:)=2
			call random_seed(put=seed)
			
			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! call activation code for n_rand samples							   !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			print *,'w,act_frac1, n_act'
			do i=1,n_rand
				r=random_normal() ! from the Netlib
				r=r*sigma_w
				r=r+mean_w
				r=abs(r) ! only positive values
				w_test=r
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! initialise arrays								                       !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call initialise_arrays(n_mode,n_sv,p_test,t_test,w_test, &
							n_aer1,d_aer1,sig_aer1, molw_org1,density_core1)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				call ctmm_activation(n_mode,n_sv,sv_flag, &
							n_aer1, d_aer1,sig_aer1,molw_core1, &
							density_core1, nu_core1, &
							org_content1,molw_org1, density_org1, delta_h_vap1, nu_org1,  &
							log_c_star1, &
							w_test, t_test,p_test, a_eq_7, b_eq_7, &
							act_frac1)
				print *,w_test,act_frac1,  sum(act_frac1*n_aer1)
			enddo
			deallocate(seed)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		
    end program main
