	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2017
	!>@brief
	!>Thermal Cloud Model (TCM):
	!>Solves solves microphysics in an analytical thermal framework for investigating
	!> mixing processes, resolution, etc.
	!>\f$ F\left(t,x,z \right)
	!>   = initialisation,microphysics,etc \f$
	!> <br><br>
	!> compile using the Makefile (note requires netcdf) and then run using: <br>
	!> ./main.exe namelist.in
	!> <br><br>
	!> (namelist used for initialisation).



	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main programme reads in information, allocates arrays, then calls the model driver

    program main
        use nrtype
        use variables
        use initialisation
        use drivers
        implicit none

        character (len=200) :: nmlfile = ' '

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! define namelist for environment
        namelist /sounding_spec/ adiabatic_prof, adiabatic_frac, &
        					psurf, tsurf, t_cbase, t_ctop, &
        					n_levels_s, q_read, theta_read, rh_read, z_read
        ! define namelist for run parameters
        namelist /run_vars/ outputfile, runtime, ip, kp, dx, dz, dt, &
                    nq, nprec, &
        			advection_scheme, ord, halo, &
        			monotone, microphysics_flag,hm_flag, theta_flag, ice_init, &
        			num_ice, mass_ice
        namelist /run_vars2/ q_type, q_init
        			
        ! define namelist for thermal
        namelist /thermal_vars/ k, dsm_by_dz_z_eq_zc, b, del_gamma_mac, &
 					del_c_s, del_c_t, epsilon_therm, w_peak, z_offset
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        print *,'Running the Thermal Cloud Model (TCM)'

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists													   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=run_vars)
        call allocate_arrays(nq,nlevels_r,q_type,q_init,q_read) 
        read(8,nml=run_vars2)
        read(8,nml=sounding_spec)
        read(8,nml=thermal_vars)
        close(8)
        o_halo=ord+2 !ord+1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! allocate and initialise the grid                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call calc_profile_2d(nq,nprec,n_levels_s,psurf,tsurf,t_cbase,t_ctop, &
        					adiabatic_prof, adiabatic_frac, &
        					q_type,q_init, z_read,theta_read, &
                            q_read,ip,kp,o_halo,dx,dz,grid1%q, grid1%precip, &
                            grid1%theta, grid1%p,grid1%x,grid1%xn, &
                            grid1%z,grid1%zn,grid1%t,grid1%rho,grid1%u,grid1%w, &
                            ice_init, num_ice, mass_ice)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! run the model                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        io1%new_file=.true.
        call model_driver_2d(nq,nprec,ip,kp,ord,o_halo,runtime,dt, &
        					grid1%q,grid1%precip,grid1%theta, &
                            grid1%p,dx,dz,grid1%x,grid1%xn,grid1%z,grid1%zn,&
							grid1%t,grid1%rho,grid1%u,grid1%w,io1%new_file, &
                            micro_init,advection_scheme, &
                            monotone,microphysics_flag,hm_flag,theta_flag, &
                            mass_ice, &
                            k,dsm_by_dz_z_eq_zc,b,del_gamma_mac, & ! variables associated 
                            del_c_s,del_c_t,epsilon_therm,w_peak, &
                            z_offset, & ! with thermal props
                            therm_init) 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *,'TCM has finished'
    end program main



