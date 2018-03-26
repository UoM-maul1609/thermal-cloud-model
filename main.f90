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
        use sub, only : n_mode, n_sv, giant_flag, method_flag, sv_flag, &
        		act_frac2, &
        		n_aer, d_aer, sig_aer, molw_core, density_core, nu_core, org_content, &
        		molw_org, density_org, delta_h_vap, nu_org, log_c_star, p_test, t_test, &
        		w_test, act_frac,mass_initial, mass_org_condensed, &
        		a_eq_7, b_eq_7, &
        		allocate_arrays, ctmm_activation, initialise_arrays, &
        		solve_semivolatiles
        implicit none

        real(sp) :: w1,t1,p1,act_frac1
        
        character (len=200) :: nmlfile = ' '




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! set up flags, etc
        namelist /bulk_aerosol_setup/ n_mode, n_sv, sv_flag, method_flag, giant_flag, &
        			a_eq_7, b_eq_7      
        ! run parameters / input arrays
        namelist /bulk_aerosol_spec/ n_aer, d_aer, sig_aer, molw_core, density_core, &
        					nu_core, org_content, molw_org, density_org, &
        					delta_h_vap, &
        					nu_org, log_c_star, p_test, t_test, w_test
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists	and allocate arrays								   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=bulk_aerosol_setup)
        ! allocate memory / init
		call allocate_arrays(n_mode,n_sv)
        
        read(8,nml=bulk_aerosol_spec)
        close(8)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! initialise arrays								                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call initialise_arrays(n_mode,p_test,t_test,w_test)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! call activation code												   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call ctmm_activation(n_mode,n_sv,sv_flag, &
        					n_aer, d_aer,sig_aer,molw_core, &
                            density_core, nu_core, &
                            org_content,molw_org, density_org, delta_h_vap, nu_org,  &
                            log_c_star, &
                            w_test, t_test,p_test, &
                            act_frac)
        print *,'act_frac1, n_act',act_frac,  sum(act_frac*n_aer)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end program main
