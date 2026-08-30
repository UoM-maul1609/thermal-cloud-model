	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Bulk Aerosol Module (BAM):
	!>Solves for droplet activation using Abdul-Razzak and Ghan,
	!>Ghosh-modified ARG, or Fountoukis and Nenes parameterisations.

    program main
        use numerics_type
        use random, only : random_normal
        use bam, only : n_mode, n_sv, giant_flag, method_flag, sv_flag, &
            n_aer1, d_aer1, sig_aer1, molw_core1, density_core1, nu_core1, org_content1, &
            molw_org1, density_org1, delta_h_vap1, nu_org1, log_c_star1, p_test, t_test, &
            w_test, act_frac1, a_eq_7, b_eq_7, r, mean_w, sigma_w, seed, l, n_rand, &
            rand_dist, ctmm_activation, read_in_bam_namelist, smax1, dcrit2
        implicit none

        character(len=512) :: nmlfile
        integer(i4b) :: i
        real(wp) :: total_number, activated_number, activated_fraction

        if (command_argument_count() < 1) then
            write(*,'(a)') 'Usage: ./main.exe namelist.in'
            error stop 'BAM: namelist filename is required'
        endif
        call get_command_argument(1,nmlfile)
        call read_in_bam_namelist(trim(nmlfile))

        total_number=sum(n_aer1)

        if(.not.rand_dist) then
            print *,'w (m s-1), act_frac1, n_drop (kg-1 dry air), activated_fraction'
            call ctmm_activation(n_mode,n_sv,sv_flag, &
                        n_aer1,d_aer1,sig_aer1,molw_core1, &
                        density_core1,nu_core1, &
                        org_content1,molw_org1,density_org1,delta_h_vap1,nu_org1, &
                        log_c_star1, &
                        w_test,t_test,p_test,a_eq_7,b_eq_7, &
                        act_frac1,smax1,dcrit2)

            activated_number=sum(act_frac1*n_aer1)
            if (total_number > 0._wp) then
                activated_fraction=activated_number/total_number
            else
                activated_fraction=0._wp
            endif
            print *,w_test,act_frac1,activated_number,activated_fraction
        else
            if (n_rand <= 0) error stop 'BAM: n_rand must be > 0'
            if (sigma_w < 0._wp) error stop 'BAM: sigma_w must be >= 0'
            if (sigma_w <= tiny(1._wp) .and. mean_w < 0._wp) &
                error stop 'BAM: random updraft distribution has no non-negative values'

            call random_seed(size=l)
            allocate(seed(1:l))
            seed(:)=2
            call random_seed(put=seed)

            print *,'w (m s-1), act_frac1, n_drop (kg-1 dry air)'
            do i=1,n_rand
                if (sigma_w <= tiny(1._wp)) then
                    r=mean_w
                else
                    ! Sample a Gaussian conditioned on w >= 0.  abs(normal) would
                    ! produce a folded normal and is biased when mean_w /= 0.
                    do
                        r=mean_w+sigma_w*random_normal()
                        if (r >= 0._wp) exit
                    enddo
                endif
                w_test=r

                call ctmm_activation(n_mode,n_sv,sv_flag, &
                            n_aer1,d_aer1,sig_aer1,molw_core1, &
                            density_core1,nu_core1, &
                            org_content1,molw_org1,density_org1,delta_h_vap1,nu_org1, &
                            log_c_star1, &
                            w_test,t_test,p_test,a_eq_7,b_eq_7, &
                            act_frac1,smax1,dcrit2)
                print *,w_test,act_frac1,sum(act_frac1*n_aer1)
            enddo
            deallocate(seed)
        endif

    end program main
