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
        use micro_module, only : set_qnames
        use thermal, only : adjust_thermal
        use w_micro_module, only : read_in_wmm_bam_namelist
        use p_micro_module, only : read_in_pamm_bam_namelist, p_initialise_aerosol
        
        implicit none

        character (len=200) :: nmlfile = ' '

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! define namelist for environment
        namelist /sounding_spec/ adiabatic_prof, adiabatic_frac, above_cloud, &
        					psurf, tsurf, t_cbase, t_ctop, &
        					n_levels_s, q_read, theta_read, rh_read, z_read
        ! define namelist for run parameters
        namelist /run_vars/ outputfile, runtime, ip, kp, dx, dz, dt, &
                    output_interval, cvis, &
                    nq, nprec, &
        			advection_scheme, ord, halo, &
        			monotone, viscous_dissipation, microphysics_flag, ice_flag, &
        			bam_nmlfile, aero_nmlfile, aero_prof_flag, hm_flag, &
        			wr_flag, rm_flag, theta_flag, &
        			drop_num_init, num_drop, ice_init, &
        			num_ice, mass_ice
        namelist /run_vars2/ q_type, q_init
        			
        ! define namelist for thermal
        namelist /thermal_vars/ k, dsm_by_dz_z_eq_zc, b, del_gamma_mac, &
 					del_c_s, del_c_t, epsilon_therm, w_peak, z_offset, &
 					adjust_thermal_flag, offset_equal_zbase
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
        ! read in bulk aerosol namelists                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        grid1%nq=nq
        grid1%nprec=nprec
        select case(microphysics_flag)
            case(0:1) ! standard
                call set_qnames(grid1%q_name,grid1%q_type,grid1%c_s,grid1%c_e,&
                    grid1%nq,grid1%ncat,grid1%nprec, &
                    grid1%iqv, grid1%iqc, grid1%ini, grid1%iqi)
            case(2) ! wmm
                call read_in_wmm_bam_namelist(bam_nmlfile, &
                    grid1%q_name,grid1%q_type,grid1%c_s,grid1%c_e,grid1%nq,&
                    grid1%ncat, &
                    grid1%nprec, &
                    grid1%iqv, grid1%iqc, grid1%inc)
            case(3) ! pamm
                call read_in_pamm_bam_namelist(bam_nmlfile,aero_nmlfile, &
                    aero_prof_flag, &
                    ice_flag, &
                    grid1%q_name,grid1%q_type,grid1%c_s,grid1%c_e,grid1%nq,&
                    grid1%ncat, &
                    grid1%nprec, grid1%n_mode, &
                    grid1%iqv, grid1%iqc, grid1%inc, grid1%iqr,grid1%inr,&
                    grid1%iqi,grid1%ini,grid1%iai,grid1%cat_am, &
                    grid1%cat_c, grid1%cat_r,grid1%cat_i)    
                    
                    nq=grid1%nq
                    nprec=grid1%nprec                      
			case default
				print *, 'error'
				stop
		end select
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        k=2.*pi/((real(ip,sp))*dx)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! allocate and initialise the grid                                     !
        ! at the moment no extra initialisation for microphysics_flag==3       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call calc_profile_2d(grid1%nq,grid1%nprec,n_levels_s,&
                            psurf,tsurf,t_cbase,t_ctop, &
        					adiabatic_prof, adiabatic_frac, &
        					q_type,q_init, z_read,theta_read, &
                            q_read,ip,kp,o_halo,dx,dz,grid1%dx2,grid1%dz2, &
                            grid1%q, grid1%qold, &
                            grid1%iqv,grid1%iqc,grid1%inc,grid1%iqi,grid1%ini,&
                            grid1%precip, &
                            grid1%theta, grid1%th_old, grid1%p,grid1%x,grid1%xn, &
                            grid1%z,grid1%zn,grid1%t,grid1%rho,grid1%u,grid1%w, &
                            grid1%delsq,grid1%vis, &
                            drop_num_init, num_drop, ice_init, num_ice, mass_ice, &
                            grid1%zbase,grid1%ztop, &
                            microphysics_flag, above_cloud)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! adjust thermal properties?                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(adjust_thermal_flag.and.adiabatic_prof) then
        
            call adjust_thermal(k,dsm_by_dz_z_eq_zc,b,del_gamma_mac,&
                                del_c_s,del_c_t,epsilon_therm, &
                                z_offset, grid1%zbase,grid1%ztop, &
                                offset_equal_zbase)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! initialise aerosol for microphysics_flag==3                          !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(microphysics_flag .eq. 3) then
            call p_initialise_aerosol(aero_prof_flag, &
                        grid1%nq,grid1%ncat,grid1%c_s,grid1%c_e, &
                        grid1%inc, &
                        ip,kp,o_halo, &
                        grid1%x,grid1%z,grid1%rho,&
                        grid1%p, grid1%t, &
                        grid1%q,grid1%qold)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! run the model                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        io1%new_file=.true.
        call model_driver_2d(grid1%nq,grid1%nprec, grid1%ncat,grid1%n_mode, &
                            ip,kp,ord,o_halo,runtime,dt, cvis, &
                            grid1%c_s,grid1%c_e, &
                            grid1%inc,grid1%iqc, &
                            grid1%cat_am, grid1%cat_c, grid1%cat_r, &
                            grid1%q_name, &
        					grid1%q,grid1%qold, &
        					grid1%precip,grid1%theta, grid1%th_old, &
                            grid1%p,dx,dz,grid1%dx2,grid1%dz2,&
                            grid1%x,grid1%xn,grid1%z,grid1%zn,&
							grid1%t,grid1%rho,grid1%u,grid1%w,grid1%delsq,grid1%vis, &
							io1%new_file, &
                            micro_init,advection_scheme, &
                            monotone,viscous_dissipation, &
                            microphysics_flag,ice_flag, hm_flag,wr_flag, rm_flag, &
                            theta_flag, &
                            mass_ice, &
                            output_interval, &
                            k,dsm_by_dz_z_eq_zc,b,del_gamma_mac, & ! variables associated 
                            del_c_s,del_c_t,epsilon_therm,w_peak, &
                            z_offset, & ! with thermal props
                            therm_init) 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *,'TCM has finished'
    end program main



