	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>code to allocate arrays, and call activation 
	module bam
	use numerics_type
	implicit none
		real(wp), parameter :: grav=9.8_wp, lv=2.5e6_wp,cp=1005._wp,molw_air=29.e-3_wp,     &
									   r_air=287._wp,r_vap=461._wp,r_gas=8.314_wp,molw_vap=18.e-3_wp, &
									   eps=r_air/r_vap,kappa=r_air/cp,rhow=1000._wp,sigma=72.e-3_wp
									   ! private
		real(wp) :: rhinit,tinit,pinit,w, ndrop_test, &
							mass_dummy,density_dummy,n_dummy,sig_dummy,d_dummy, &
							tcb, pcb,xmin,a,smax,smax1, &
							alpha_sup, sigma_sup, g, chi, sd_dummy, s,c0   
		! size n_mode
		real(wp), allocatable, dimension(:) :: n_aer1, sig_aer1, d_aer1, &
										n_aer, sig_aer, d_aer, d_aer_new, sgi, &
										density_final,mass_initial, & 
									   mass_final,sd,b,sm,eta,f1,f2 
							  
		! size n_mode
		real(wp), allocatable, dimension(:) :: density_core, & 
								  molw_core, & 
								  nu_core,act_frac,act_frac1, act_frac2 , &
								  density_core1, & 
								  molw_core1, & 
								  nu_core1, dcrit2
		! size n_sv
		real(wp), allocatable, dimension(:) :: molw_org, r_org, log_c_star, cstar, &
												org_content, org_content1, &
											  density_org,nu_org,mass_org_condensed, &
											  delta_h_vap, epsilon1, c_ions, &
											  molw_org1, log_c_star1, &
											  density_org1,nu_org1,&
											  delta_h_vap1
											  
		real(wp), dimension(6) :: c1	! private
		integer(i4b) :: n_mode_s		! private
		real(wp) :: p_test, t_test, w_test, a_eq_7, b_eq_7 ! public
		real(wp) :: ghosh_sigma_acc=1.6_wp, ghosh_dmin=80.e-9_wp, &
		            ghosh_dmax=1.e-6_wp, ghosh_sigma_used=1.6_wp
		
		integer(i4b) :: n_mode, n_sv, method_flag, giant_flag, sv_flag
		integer(i4b) :: ghosh_sigma_mode=0
		! method_flag: 1=Abdul-Razzak and Ghan; 2=Fountoukis and Nenes;
		! 3=Fountoukis and Nenes with quadrature; 4=Ghosh-modified ARG.
		! ghosh_sigma_mode: 0=fixed ghosh_sigma_acc (published-style);
		! 1=effective width diagnosed from the combined PSD (experimental).
	
		! for random number:
		real(wp) :: r, mean_w, sigma_w
		real(wp), dimension(10,10) :: rs
		integer(i4b), allocatable, dimension(:) :: seed
		integer(i4b) :: l, n_rand
		logical :: rand_dist=.false.

	private 
	public :: ctmm_activation, allocate_arrays, initialise_arrays, &
		read_in_bam_namelist, find_d_and_s_crits, &
		n_mode, n_sv, method_flag, giant_flag, sv_flag, &
		p_test, t_test, w_test, a_eq_7, b_eq_7, &
		ghosh_sigma_mode, ghosh_sigma_acc, ghosh_dmin, ghosh_dmax, ghosh_sigma_used, &
		n_aer1, d_aer1, sig_aer1, &
		org_content1, molw_org1, log_c_star1, density_org1, nu_org1, delta_h_vap1, &
		molw_core1, density_core1, nu_core1, act_frac1, &
		r, mean_w, sigma_w, rs, seed, l, n_rand, rand_dist, &
		smax1, dcrit2
				
	contains
	
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>root-find to solve for smax and dcrit needed so that drop number is achieved
	!>parameterisation developed at university of manchester
	!>@param[in] p,t
	!>@param[inout] ndrop drop number mixing ratio (kg-1 dry air)
	!>@param[inout] scrit, dscrit
	subroutine find_d_and_s_crits(p,t,ndrop,w,smax,dcrit)
	    use numerics_type
	    use numerics, only : zeroin
	    implicit none
	    real(wp), intent(in) :: p,t
	    real(wp), intent(inout) :: ndrop, smax, w
	    real(wp), dimension(n_mode), intent(inout) :: dcrit
	    
	    
	    ! adjust ndrop to 99% of total aerosol number:
	    ndrop=min(ndrop,0.99_wp*sum(n_aer1))
        ndrop_test=ndrop
        
        ! for use inside the root-finder
        pinit=p
        tinit=t
        ! find the updraft speed required to activate these aerosol particles
        w=zeroin(1.e-40_wp,100._wp,find_wcbase,1.e-20_wp)
        ! call again so that dcrit1 and smax are properly set
        call initialise_arrays(n_mode,n_sv,pinit,tinit,w, &
                    n_aer1,d_aer1,sig_aer1, molw_org1,density_core1)

        call ctmm_activation(n_mode,n_sv,sv_flag, &
                    n_aer1, d_aer1,sig_aer1,molw_core1, &
                    density_core1, nu_core1, &
                    org_content1,molw_org1, density_org1, delta_h_vap1, nu_org1,  &
                    log_c_star1, &
                    w, tinit,pinit, a_eq_7, b_eq_7, &
                    act_frac1,smax,dcrit)
        
	    
	end subroutine find_d_and_s_crits
	
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>root-finding helper function
	!>@param[in] w
	function find_wcbase(w)
	    use numerics_type
	    implicit none
	    real(wp), intent(in) :: w
	    real(wp) :: find_wcbase, smax
	
        call initialise_arrays(n_mode,n_sv,pinit,tinit,w, &
                    n_aer1,d_aer1,sig_aer1, molw_org1,density_core1)

        call ctmm_activation(n_mode,n_sv,sv_flag, &
                    n_aer1, d_aer1,sig_aer1,molw_core1, &
                    density_core1, nu_core1, &
                    org_content1,molw_org1, density_org1, delta_h_vap1, nu_org1,  &
                    log_c_star1, &
                    w, tinit,pinit, a_eq_7, b_eq_7, &
                    act_frac1,smax,dcrit2)

	    find_wcbase=ndrop_test-sum(act_frac1*n_aer1)
	    
	end function find_wcbase
	
	
	
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the activated fraction of a lognormally distributed
	!>aerosol distribution including condensation of semi-volatile organics.
	!>code follows a paper by Connolly, Topping, Malavelle and Mcfiggans (in ACP 2014)
	!>parameterisation developed at university of manchester
	!>@param[in] n_modes1, n_sv1 : number of aerosol and volatility modes
	!>@param[in] sv_flag: flag for cocondensation
	!>@param[inout] n_aer, d_aer, sig_aer, molw_core, density_core, nu_core
	!> n_aer1 is number mixing ratio (kg-1 dry air); n_aer is converted internally to m-3
	!>@param[in] org_content1: amount of organic in volatility bins (microgram kg-1 air)
	!>@param[in] molw_org1: molecular weight in volatility bins
	!>@param[in] density_org1: density of organic in volatility bins
	!>@param[in] delta_h_vap1: enthalpy change in volatility bins
	!>@param[in] nu_org1: van hoff factor in volatility bins
	!>@param[in] log_c_star1: volatility bins
	!>@param[in] w1, t1, p1, a, b: vertical wind, temperature, pressure + params in ARG
	!>@param[inout] act_frac1: activated fraction in each mode
	!>@param[inout] smax1: maximum supersaturation
	!>@param[inout] dcrit1: critical diameters in each mode
	subroutine ctmm_activation(n_modes1,n_sv1,sv_flag, n_aer1,d_aer1,sig_aer1,molw_core1, &
                               density_core1, nu_core1, org_content1, &
                               molw_org1, density_org1, delta_h_vap1, nu_org1,  &
                               log_c_star1, &
                               w1, t1,p1,a_arg,b_arg, &
                               act_frac1,smax1,dcrit1)

        use numerics_type
        use numerics, only : zeroin, fmin
        implicit none
        real(wp), dimension(:), intent(in) :: n_aer1
        real(wp), dimension(:), intent(inout) :: d_aer1, sig_aer1, molw_core1, &
                                                density_core1, nu_core1
        real(wp), dimension(:), intent(in) :: org_content1, molw_org1, &
                                               density_org1, delta_h_vap1, nu_org1, log_c_star1
        real(wp), intent(in) :: w1,t1,p1, a_arg, b_arg
        integer, intent(in) :: n_modes1, n_sv1, sv_flag
        real(wp), dimension(:), intent(inout) :: act_frac1, dcrit1
        real(wp), intent(inout) :: smax1

        integer(i4b) :: i
        real(wp) :: n_total, sum_term, org_mass_total, org_volume_total, &
                    org_ion_total, mode_weight, solute_moles_per_mass, &
                    ghosh_f, ghosh_g, ghosh_p, p_exp(n_modes1)
        logical :: active(n_modes1)

        if (n_modes1 <= 0) error stop 'BAM: n_modes must be positive'
        if (n_sv1 <= 0) error stop 'BAM: n_sv must be positive'
        if (size(n_aer1) < n_modes1 .or. size(d_aer1) < n_modes1 .or. &
            size(sig_aer1) < n_modes1) error stop 'BAM: aerosol array too small'
        if (any(n_aer1(1:n_modes1) < 0._wp)) error stop 'BAM: n_aer1 must be >= 0 kg-1'
        if (any(d_aer1(1:n_modes1) <= 0._wp)) error stop 'BAM: d_aer1 must be > 0 m'
        if (any(sig_aer1(1:n_modes1) <= 0._wp)) error stop 'BAM: sig_aer1=ln(sigma_g) must be > 0'
        if (any(molw_core1(1:n_modes1) <= 0._wp)) error stop 'BAM: molw_core1 must be > 0'
        if (any(density_core1(1:n_modes1) <= 0._wp)) error stop 'BAM: density_core1 must be > 0'
        if (any(nu_core1(1:n_modes1) <= 0._wp)) error stop 'BAM: nu_core1 must be > 0'
        if (method_flag < 1 .or. method_flag > 4) error stop 'BAM: method_flag must be 1, 2, 3 or 4'
        if (sv_flag < 0 .or. sv_flag > 1) error stop 'BAM: sv_flag must be 0 or 1'

        ! Make this routine self-contained.  The public aerosol number interface is
        ! kg-1 dry air, but the published ARG/FN supersaturation balances require
        ! number per unit volume. initialise_arrays performs that internal conversion
        ! and constructs mass_initial in kg m-3.
        call initialise_arrays(n_modes1,n_sv1,p1,t1,w1,n_aer1, &
                               d_aer1,sig_aer1,molw_org1,density_core1)

        n_mode_s=n_modes1
        molw_core=molw_core1
        density_core=density_core1
        nu_core=nu_core1

        org_content=0._wp
        molw_org=molw_org1
        density_org=density_org1
        delta_h_vap=delta_h_vap1
        nu_org=nu_org1
        log_c_star=log_c_star1
        mass_org_condensed=0._wp

        active = n_aer(1:n_modes1) > 0._wp
        n_total=sum(n_aer(1:n_modes1),mask=active)
        act_frac1=0._wp
        act_frac2=0._wp
        dcrit2=huge(1._wp)
        dcrit1(1:n_modes1)=huge(1._wp)
        smax=0._wp
        smax1=0._wp

        ! No ascent means no new cloud-base activation.  This also avoids the
        ! singular w=0 limits in both activation parameterisations.
        if (w1 <= 0._wp .or. n_total <= 0._wp) return

        if(sv_flag.eq.1) then
            if (any(org_content1(1:n_sv1) < 0._wp)) error stop 'BAM: org_content1 must be >= 0 ug kg-1'
            if (any(molw_org1(1:n_sv1) <= 0._wp)) error stop 'BAM: molw_org1 must be > 0'
            if (any(density_org1(1:n_sv1) <= 0._wp)) error stop 'BAM: density_org1 must be > 0'
            if (any(nu_org1(1:n_sv1) <= 0._wp)) error stop 'BAM: nu_org1 must be > 0'

            ! External organic mass is microgram per kg dry air.  The equilibrium
            ! partitioning calculation is volumetric, so convert to kg m-3 using
            ! the same dry-air density used for n_aer.
            org_content=org_content1*1.e-9_wp*dry_air_density(t1,p1,rhinit)

            call solve_semivolatiles(n_modes1,n_sv1, &
                    org_content, log_c_star, delta_h_vap, nu_org, molw_org, &
                    mass_initial, nu_core, molw_core,rhinit, t1, &
                    mass_org_condensed)

            ! Retain the original Connolly-style multi-mode approximation: total
            ! condensed organic mass is distributed in proportion to particle
            ! number.  This is not the full Crooks multi-mode treatment.
            org_mass_total=sum(mass_org_condensed)
            org_volume_total=sum(mass_org_condensed/density_org)
            do i=1,n_modes1
                if (.not.active(i)) then
                    mass_final(i)=0._wp
                    density_final(i)=density_core(i)
                    cycle
                endif
                mode_weight=n_aer(i)/n_total
                mass_final(i)=mass_initial(i)+org_mass_total*mode_weight
                density_final(i)=mass_initial(i)/density_core(i)+org_volume_total*mode_weight
                if (density_final(i) <= 0._wp) error stop 'BAM: non-positive mixed aerosol volume'
                density_final(i)=mass_final(i)/density_final(i)
            enddo

            ! Arithmetic standard deviation of a lognormal distribution. sig_aer
            ! is ln(sigma_g), hence the 0.5*sig_aer**2 term.
            sd=0._wp
            do i=1,n_modes1
                if (active(i)) sd(i)=exp(log(d_aer(i))+0.5_wp*sig_aer(i)**2)* &
                                      sqrt(exp(sig_aer(i)**2)-1._wp)
            enddo

            ! Shift each active mode while conserving its final mass and keeping
            ! its arithmetic standard deviation fixed (Connolly et al., 2014).
            d_aer_new=d_aer
            do i=1,n_modes1
                if (.not.active(i)) cycle
                density_dummy=density_final(i)
                mass_dummy=mass_final(i)
                n_dummy=n_aer(i)
                sd_dummy=sd(i)
                d_aer_new(i)=fmin(d_aer(i),2000.e-9_wp,mass_integrate,1.e-30_wp)
                xmin=mass_integrate(d_aer_new(i))
            enddo

            d_aer=d_aer_new
            do i=1,n_modes1
                if (.not.active(i)) cycle
                d_dummy=d_aer_new(i)
                sd_dummy=sd(i)
                sig_aer(i)=zeroin(1.e-9_wp,2._wp,find_sig_aer,1.e-6_wp)
            enddo
        else
            mass_org_condensed=0._wp
            mass_final=mass_initial
            density_final=density_core
        endif

        ! Mixed-solute Kohler B.  For sv_flag=0 this reduces exactly to
        ! B = nu * M_w * rho_a / (M_a * rho_w).  For sv_flag=1 the core and
        ! condensed-organic solute mole contributions are summed explicitly.
        org_ion_total=0._wp
        if (sv_flag.eq.1) org_ion_total=sum(mass_org_condensed*nu_org/molw_org)
        b=1._wp
        do i=1,n_modes1
            if (.not.active(i)) cycle
            if (mass_final(i) <= 0._wp) error stop 'BAM: non-positive final aerosol mass'
            mode_weight=n_aer(i)/n_total
            solute_moles_per_mass=(mass_initial(i)*nu_core(i)/molw_core(i) + &
                                   org_ion_total*mode_weight)/mass_final(i)
            b(i)=molw_vap*density_final(i)/rhow*solute_moles_per_mass
            if (b(i) <= 0._wp) error stop 'BAM: non-positive Kohler B parameter'
        enddo

        if(method_flag.eq.1 .or. method_flag.eq.4) then
            a=2._wp*sigma*molw_vap/(rhow*r_gas*tcb)

            sm=huge(1._wp)
            eta=huge(1._wp)
            f1=0._wp
            f2=0._wp
            p_exp=1.5_wp
            do i=1,n_modes1
                if (.not.active(i)) cycle
                sm(i)=2._wp/sqrt(b(i))*(a/(3._wp*d_aer(i)/2._wp))**1.5_wp
            enddo

            alpha_sup=grav*molw_vap*lv/(cp*r_gas*tcb**2)- &
                      grav*molw_air/(r_gas*tcb)
            sigma_sup=r_gas*tcb/(svp(tcb)*molw_vap)+ &
                      molw_vap*lv**2/(cp*pcb*molw_air*tcb)
            g=rhow*r_gas*tcb/(svp(tcb)*dd(tcb,pcb)*molw_vap) + &
              lv*rhow/(ka(tcb)*tcb)*(lv*molw_vap/(r_gas*tcb)-1._wp)
            g=1._wp/g

            chi=(2._wp/3._wp)*(alpha_sup*w/g)**0.5_wp*a
            do i=1,n_modes1
                if (.not.active(i)) cycle
                eta(i)=(alpha_sup*w/g)**1.5_wp/ &
                       (2._wp*pi*rhow*sigma_sup*n_aer(i))
            enddo

            if (method_flag.eq.1) then
                ! Original ARG2000 form.  BAM retains a_arg and b_arg as
                ! configurable coefficients in the fitted f_i function.
                do i=1,n_modes1
                    if (.not.active(i)) cycle
                    f1(i)=a_arg*exp(b_arg*sig_aer(i)**2)
                    f2(i)=1._wp+0.25_wp*sig_aer(i)
                enddo
            else
                ! Ghosh et al. (2025) modified ARG.  The same f, g and
                ! kinetically-limited p are used for every aerosol mode and
                ! depend only on one accumulation-mode reference width.
                select case (ghosh_sigma_mode)
                case (0)
                    ghosh_sigma_used=ghosh_sigma_acc
                    if (ghosh_sigma_used < 1.4_wp .or. ghosh_sigma_used > 2.1_wp) &
                        error stop 'BAM: fixed ghosh_sigma_acc must be in calibrated range 1.4 to 2.1'
                case (1)
                    ghosh_sigma_used=effective_accum_sigma(n_modes1,n_aer,d_aer,sig_aer, &
                                                          ghosh_dmin,ghosh_dmax)
                case default
                    error stop 'BAM: ghosh_sigma_mode must be 0 (fixed) or 1 (effective)'
                end select

                if (ghosh_sigma_mode.eq.1 .and. &
                    (ghosh_sigma_used < 1.4_wp .or. ghosh_sigma_used > 2.1_wp)) then
                    write(*,'(a,f10.5,a)') 'BAM warning: diagnosed Ghosh sigma_g=', &
                        ghosh_sigma_used,' is outside published calibration range 1.4-2.1'
                endif

                ghosh_f=0.0135_wp*exp(2.367_wp*ghosh_sigma_used)
                ghosh_g=1.1058_wp-0.315_wp*ghosh_sigma_used
                ghosh_p=-0.5073_wp+1.5088_wp*ghosh_sigma_used- &
                        0.3699_wp*ghosh_sigma_used**2

                if (ghosh_f <= 0._wp .or. ghosh_g <= 0._wp .or. ghosh_p <= 0._wp) &
                    error stop 'BAM: non-positive Ghosh f, g or p; check effective/fixed sigma'

                do i=1,n_modes1
                    if (.not.active(i)) cycle
                    f1(i)=ghosh_f
                    f2(i)=ghosh_g
                    if (chi/eta(i) > 1._wp) p_exp(i)=ghosh_p
                enddo
            endif

            sum_term=0._wp
            do i=1,n_modes1
                if (.not.active(i)) cycle
                sum_term=sum_term+1._wp/sm(i)**2* &
                    (f1(i)*(chi/eta(i))**p_exp(i) + &
                     f2(i)*(sm(i)**2/(eta(i)+3._wp*chi))**0.75_wp)
            enddo
            if (sum_term <= 0._wp) return
            smax=1._wp/sqrt(sum_term)

            do i=1,n_modes1
                if (.not.active(i)) cycle
                act_frac2(i)=0.5_wp*(1._wp-erf(2._wp*log(sm(i)/smax)/ &
                               (3._wp*sqrt(2._wp)*sig_aer(i))))
                act_frac2(i)=min(1._wp,max(0._wp,act_frac2(i)))
                dcrit2(i)=2._wp*a/3._wp*(2._wp/smax/sqrt(b(i)))**(2._wp/3._wp)
            enddo

        else
            a=4._wp*sigma*molw_vap/(rhow*r_gas*tcb)
            sgi=huge(1._wp)
            do i=1,n_modes1
                if (.not.active(i)) cycle
                sgi(i)=sqrt(4._wp*a**3/(27._wp*b(i)*d_aer(i)**3))
            enddo

            smax=max(zeroin(1.e-20_wp,100.e-2_wp,fountoukis_nenes,1.e-20_wp),1.e-20_wp)

            do i=1,n_modes1
                if (.not.active(i)) cycle
                act_frac2(i)=0.5_wp*(1._wp-erf(2._wp*log(sgi(i)/smax)/ &
                               (3._wp*sqrt(2._wp)*sig_aer(i))))
                act_frac2(i)=min(1._wp,max(0._wp,act_frac2(i)))
                ! Use the smax solved in this call; smax1 is an output and may be
                ! uninitialised on entry.
                dcrit2(i)=((4._wp*a**3)/(smax**2*(27._wp*b(i))))**(1._wp/3._wp)
            enddo
        endif

        act_frac1(1:n_modes1)=act_frac2(1:n_modes1)
        smax1=smax
        dcrit1(1:n_modes1)=dcrit2(1:n_modes1)

    end subroutine ctmm_activation
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! diffusivity of water vapour in air										   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the diffusivity of water vapour in air
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@return dd: diffusivity of water vapour in air
	function dd(t,p)
	    use numerics_type
		implicit none
		real(wp), intent(in) :: t, p
		real(wp) :: dd
		dd=2.11e-5_wp*(t/273.15_wp)**1.94*(101325._wp/p)
	end function dd

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! conductivity of water vapour												   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the thermal conductivity of air
	!>@param[in] t: temperature
	!>@return ka: thermal conductivity of air
	function ka(t)
	    use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: ka
		ka=(5.69_wp+0.017_wp*(t-273.15_wp))*1.e-3_wp*4.187_wp
	end function ka

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! dry potential temperature 												   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the thermal conductivity of air
	!> Note uses tinit, pinit, rhinit from module
	!>@param[in] p: pressure (Pa)
	!>@return dry_potential temperature (K)
	function dry_potential(p)
	    use numerics_type
		implicit none
		real(wp), intent(in) :: p
		real(wp) :: dry_potential
		real(wp) :: total_water1, total_water2, tcalc 
		total_water1=rhinit*eps*svp(tinit)/(pinit-svp(tinit))

		!print *,'kappa',svp(tinit),eps
		tcalc=tinit*(p/pinit)**kappa

		total_water2=eps*svp(tcalc)/(p-svp(tcalc))

		dry_potential=total_water2-total_water1

	end function dry_potential

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! saturation vapour pressure over liquid                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the saturation vapour pressure over liquid water according to buck fit
	!>@param[in] t: temperature
	!>@return svp_liq: saturation vapour pressure over liquid water
	function svp(t)
	    use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: svp
		svp = 100._wp*6.1121_wp* &
			  exp((18.678_wp - (t-273.15_wp)/ 234.5_wp)* &
			  (t-273.15_wp)/(257.14_wp + (t-273.15_wp)))
	end function svp

    ! Dry-air density corresponding to a total pressure, temperature and RH.
    ! The BAM public number/mass mixing ratios are per kg of dry air.
    function dry_air_density(t,p,rh)
        use numerics_type
        implicit none
        real(wp), intent(in) :: t,p,rh
        real(wp) :: dry_air_density, e
        e=max(0._wp,min(rh*svp(t),0.999999_wp*p))
        dry_air_density=(p-e)/(r_air*t)
    end function dry_air_density


    !>@brief
    !>Diagnose one equivalent geometric standard deviation for the combined
    !>aerosol PSD over a specified dry-diameter interval.  This is an
    !>experimental extension for use with the Ghosh-modified ARG scheme; it is
    !>not part of the published Ghosh et al. (2025) formulation.
    !>
    !>The calculation treats ln(D) for each lognormal mode as a normal
    !>distribution, integrates its zeroth/first/second ln(D) moments between
    !>dmin and dmax analytically, combines the modes by number, and returns
    !>sigma_g,eff = exp(sqrt(var(ln D))).
    function effective_accum_sigma(nmodes,n,d,logsig,dmin,dmax)
        use numerics_type
        implicit none
        integer(i4b), intent(in) :: nmodes
        real(wp), dimension(nmodes), intent(in) :: n,d,logsig
        real(wp), intent(in) :: dmin,dmax
        real(wp) :: effective_accum_sigma

        integer(i4b) :: i
        real(wp) :: mu, ss, ya, yb, za, zb, prob, phia, phib
        real(wp) :: m0, m1, m2, mode_m1, mode_m2, mean_y, var_y
        real(wp), parameter :: sqrt2=1.4142135623730950488_wp
        real(wp), parameter :: inv_sqrt_2pi=0.39894228040143267794_wp

        if (dmin <= 0._wp .or. dmax <= dmin) &
            error stop 'BAM: require 0 < ghosh_dmin < ghosh_dmax'

        ya=log(dmin)
        yb=log(dmax)
        m0=0._wp
        m1=0._wp
        m2=0._wp

        do i=1,nmodes
            if (n(i) <= 0._wp) cycle
            if (d(i) <= 0._wp .or. logsig(i) <= 0._wp) cycle

            mu=log(d(i))
            ss=logsig(i)
            za=(ya-mu)/ss
            zb=(yb-mu)/ss

            prob=0.5_wp*(erf(zb/sqrt2)-erf(za/sqrt2))
            if (prob <= tiny(1._wp)) cycle

            phia=inv_sqrt_2pi*exp(-0.5_wp*za**2)
            phib=inv_sqrt_2pi*exp(-0.5_wp*zb**2)

            mode_m1=mu*prob + ss*(phia-phib)
            mode_m2=(mu**2+ss**2)*prob + 2._wp*mu*ss*(phia-phib) + &
                    ss**2*(za*phia-zb*phib)

            m0=m0+n(i)*prob
            m1=m1+n(i)*mode_m1
            m2=m2+n(i)*mode_m2
        enddo

        if (m0 <= tiny(1._wp)) &
            error stop 'BAM: no aerosol number lies in Ghosh effective-width interval'

        mean_y=m1/m0
        var_y=max(0._wp,m2/m0-mean_y**2)
        effective_accum_sigma=exp(sqrt(var_y))

        if (effective_accum_sigma <= 1._wp) &
            error stop 'BAM: diagnosed Ghosh effective sigma_g is not > 1'
    end function effective_accum_sigma


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 3rd moment - for integration												   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the third moment of a lognormal
	!>@param[in] d: diameter (m)
	!>@return ln3: third moment in a size interval
	function ln3(d)
	    use numerics_type
		implicit none
		real(wp), dimension(:), intent(in) :: d
		real(wp), dimension(size(d)) :: ln3

		! add all modes together
		ln3=pi*d**2/(sqrt(twopi*sig_dummy**2)*6.)* &
		exp(-log(d/d_dummy)**2/(2.*sig_dummy**2))* &
		density_dummy*n_dummy
	end function ln3

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! mass in a lognormal       												   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the mass across one lognormal distribution
	!>@param[in] d1: mode diameter (m)
	!>@return mass_integrate: total mass in the distribution
	function mass_integrate(d1)
	    use numerics_type
		use numerics, only : zeroin
		implicit none
		real(wp), intent(in) :: d1
		real(wp) :: mass_integrate

		d_dummy=d1  ! guess at d_aer, used to calculate the new standard deviation
		sig_dummy=zeroin(1.e-9_wp,2._wp,find_sig_aer,1.e-6_wp)
		!  mass_integrate=abs(qromb(ln3,0.d-10,3.d-6)-mass_dummy)
		! moment generating function
		! http://www.mlahanas.de/math/lognormal.htm
		mass_integrate=n_dummy*exp(3._wp*log(d_dummy) + &
					3._wp**2*sig_dummy**2/2._wp) &
				   *density_dummy*pi/6._wp
		mass_integrate=abs(mass_integrate-mass_dummy)

	end function mass_integrate

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! find sigma for the size distribution										   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the thermal conductivity of air
	!>@param[in] sig_aer_new: logarithmic geometric width ln(sigma_g)
	!>@return find_sig_aer: arithmetic standard deviation
	function find_sig_aer(sig_aer_new)
	    use numerics_type
		real(wp), intent(in) :: sig_aer_new
		real(wp) :: find_sig_aer
		real(wp) :: sd1
		sd1=exp(log(d_dummy)+0.5_wp*sig_aer_new**2)*sqrt(exp(sig_aer_new**2)-1._wp)
		find_sig_aer=sd1-sd_dummy
	end function find_sig_aer

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! fountoukis and nenes integrals											   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the integral in the fountoukis and nenes method
	!>@param[in] smax1: max supersaturation
	!>@return fountoukis_nenes: integral for activation
	function fountoukis_nenes(smax1)
	    use numerics_type
		use numerics, only : romb
		implicit none
		real(wp), intent(in) :: smax1
		real(wp) :: fountoukis_nenes

		real(wp) :: integral, integral1,smax2,spart,upart,umax,i1,i2,discriminant, &
				  deq,dw3
		integer(i4b):: i

		smax2=max(smax1,1e-20_wp)
		alpha_sup=grav*molw_vap*lv/(cp*r_gas*tcb**2)- &
		grav*molw_air/(r_gas*tcb)                ! eq 11 arg
		sigma_sup=r_gas*tcb/(svp(tcb)*molw_vap)+ &
				molw_vap*lv**2/(cp*pcb*molw_air*tcb)  ! eq 11?
		g=rhow*r_gas*tcb/(svp(tcb)*dd(tcb,pcb)*molw_vap) + &
		lv*rhow/(ka(tcb)*tcb)*(lv*molw_vap/(r_gas*tcb)-1_wp)  ! eq 12
		g=4_wp/g
		! do the integral 
		! f-n (2005 use approximate form for integral, here we use quadrature)
		c1(1)=2._wp*a/3._wp
		c1(2)=g/alpha_sup/w
		c1(6)=smax2


		!  discriminant=smax2**4_wp-16._wp*a**2_wp*alpha_sup*w/(9_wp*g)
		discriminant=1._wp-16._wp*a**2._wp*alpha_sup*w/(9._wp*g)/smax2**4._wp

		if(discriminant.ge.0._wp) then
		 spart=smax2*(0.5_wp*(1._wp+ &
		   (1._wp-16_wp*a**2*alpha_sup*w/(9._wp*g*smax2**4 ))**0.5_wp) &
		   )**0.5_wp
		else
		 spart=smax2*min(2.e7_wp*a/3_wp*smax2**(-0.3824_wp),1._wp) 
		endif

		integral=0._wp
		integral1=0._wp
		do i=1,n_mode_s
          if (n_aer(i) <= 0._wp) cycle
		  c1(3)=2_wp*n_aer(i)/(3._wp*sqrt(2._wp*pi)*sig_aer(i))
		  c1(4)=sgi(i)
		  c1(5)=2._wp*sig_aer(i)**2
		  ! note for multiple modes, have to calculate the integral below for each mode
		  if(method_flag.eq.3) then
			 integral1=integral1+romb(integral3_fn,0._wp,smax2)
		  endif

		! approximate method
		  upart=2._wp*log(sgi(i)/spart)/(3._wp*sqrt(2._wp)*sig_aer(i))
		  umax =2._wp*log(sgi(i)/smax2)/(3._wp*sqrt(2._wp)*sig_aer(i))

		  i1=n_aer(i)/2._wp*sqrt(g/alpha_sup/w)*smax2* &
			 (erfc(upart)-5e-1_wp*(sgi(i)/smax2)**2_wp*exp(4.5_wp*sig_aer(i)**2)* &
			 erfc(upart+3_wp*sig_aer(i)/sqrt(2._wp)))
		  i2=a*n_aer(i)/(3_wp*sgi(i))*exp(9_wp/8_wp*sig_aer(i)**2_wp)* &
			 (erf(upart-3_wp*sig_aer(i)/(2._wp*sqrt(2._wp))) - &
			 erf(umax-3._wp*sig_aer(i)/(2._wp*sqrt(2._wp) )) )
		!      print *,'i1, i2',i1,i2, upart, umax

		  if(method_flag.eq.2) then
			 integral1=integral1+(i1+i2)
		  endif   

		  ! for the giant ccn - barahona et al (2010)
		  if(giant_flag.eq.1) then
			dw3  = (sqrt(2._wp)*log(sgi(i)/spart)/3_wp/sig_aer(i))- &
				   (1.5_wp*sig_aer(i)/sqrt(2._wp))
			deq= a*2_wp/sgi(i)/3_wp/sqrt(3._wp)      
			dw3=n_aer(i)*deq*exp(9_wp/8_wp*sig_aer(i)**2._wp)*smax2* &
			   (erfc(dw3))*((g*alpha_sup*w)**0.5_wp)   
		!         
			dw3=dw3/(2._wp*g*smax2)*(g/alpha_sup/w)**0.5_wp
			integral1=integral1 +dw3
		  endif
		enddo
		!  print *,'integral',integral,integral1
		! cost function - eq 10 fountoukis and nenes

		fountoukis_nenes=(2._wp*alpha_sup*w/(pi*sigma_sup*rhow)- &
					   g*smax2*integral1)

	end function fountoukis_nenes

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! integrand 3 in FN          												   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the integrand in the FN method (for the quadrature case)
	!>@param[in] s: s parameter
	!>@return integral3_fn: integrand in FN (for use with quadrature)
	function integral3_fn(s)
	    use numerics_type
		implicit none
		real(wp), dimension(:), intent(in) :: s
		real(wp), dimension(size(s)) :: integral3_fn

		integral3_fn=((c1(1)/(s+1.e-50_wp))**2+c1(2)*(c1(6)**2-s**2))**0.5_wp* &
				   c1(3)/(s+1.e-50_wp)*exp(-log((c1(4)/(s+1.e-50_wp))**(2._wp/3._wp))**2/c1(5)) ! eq 14 f-n
	end function integral3_fn
		
		
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the condensed organics (in Crooks, GMD 2016)
	!>parameterisation developed at university of manchester
	!>@param[in] n_modes1: number of aerosol modes
	!>@param[in] n_sv1: number of volatility bins
	!>@param[in] org_content1: organic mass concentration internal to solver (kg m-3)
	!>@param[in] nu_org1: van hoff factor in volatility bins
	!>@param[in] molw_org1: molecular weight in volatility bins
	!>@param[in] mass_core: mass in aerosol modes
	!>@param[in] nu_core1: van hoff factor in core
	!>@param[in] molw_core1: molecular weight of modes (core only)
	!>@param[in] s1, t1: rh and temperature.
	!>@param[inout] mass_org_condensed1: condensed organics
	subroutine solve_semivolatiles(n_modes1,n_sv1, &
                    org_content1, log_c_star1, delta_h_vap1, &
                    nu_org1, molw_org1, &
                    mass_core1, nu_core1, molw_core1,s1, t1, &
                    mass_org_condensed1)
        use numerics_type
        use numerics, only : zeroin
        implicit none
        integer(i4b), intent(in) :: n_modes1, n_sv1
        real(wp), dimension(n_sv1), intent(in) :: org_content1, nu_org1, molw_org1, &
                            log_c_star1, delta_h_vap1
        real(wp), dimension(n_modes1), intent(in) :: mass_core1, nu_core1, molw_core1
        real(wp), intent(in) :: s1, t1
        real(wp), dimension(n_sv1), intent(inout) :: mass_org_condensed1

        real(wp) :: ct, ct_lo, ct_hi, tol
        real(wp), dimension(n_sv1) :: c_c

        nu_org=nu_org1
        molw_org=molw_org1
        nu_core=nu_core1
        molw_core=molw_core1
        s=s1
        mass_org_condensed1=0._wp
        epsilon1=0._wp

        if (s >= 1._wp) error stop 'BAM: semi-volatile equilibrium requires RH < 1'
        if (any(molw_org1 <= 0._wp)) error stop 'BAM: organic molecular weights must be > 0'

        ! log_c_star is log10(C* / (microgram m-3)), following the 2014
        ! semi-volatile activation formulation.  Convert the temperature-adjusted
        ! C* from microgram m-3 -> kg m-3 -> mol m-3 so it is consistent with
        ! the molar partitioning concentrations below.
        cstar = 10._wp**log_c_star1*(298.15_wp/t1)* &
                exp(-delta_h_vap1*1.e3_wp/r_gas*(1._wp/t1-1._wp/298.15_wp))* &
                1.e-9_wp/molw_org1

        ! org_content1 and mass_core1 arrive here in kg m-3.  Convert to
        ! effective solute molar concentrations (mol m-3).
        c_ions=org_content1*nu_org1/molw_org1
        c0=sum(mass_core1*nu_core1/molw_core1)

        if (sum(c_ions) <= 0._wp .or. c0 <= 0._wp) return

        ! Equation 5 of the absorptive partitioning treatment brackets C_OA
        ! between the no-organic-condensation and all-organic-condensation limits.
        ct_lo=c0/(1._wp-s)
        ct_hi=(c0+sum(c_ions))/(1._wp-s)
        tol=max(1.e-14_wp,1.e-8_wp*ct_hi)
        if (ct_hi-ct_lo <= tol) then
            ct=ct_lo
        else
            ct=zeroin(ct_lo,ct_hi,partition01,tol)
        endif

        epsilon1=(1._wp+cstar/max(ct,tiny(1._wp)))**(-1)
        c_c=c_ions*epsilon1
        mass_org_condensed1=c_c/nu_org1*molw_org1

    end subroutine solve_semivolatiles


    function partition01(ct)
        use numerics_type
        implicit none
        real(wp), intent(in) :: ct
        real(wp), dimension(size(epsilon1)) :: c_c
        real(wp) :: partition01
        real(wp) :: ct1, ct2

        ct2=max(abs(ct),tiny(1._wp))
        epsilon1=(1._wp+cstar/ct2)**(-1)
        c_c=c_ions*epsilon1

        ! c_ions already contains the effective van't Hoff multiplier, so do
        ! not multiply by nu_org a second time here.
        ct1=1._wp/(1._wp-s)*(sum(c_c)+c0)
        partition01=ct1-ct2

    end function partition01
	
	
		
		
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays for activation code
	!>@param[in] n_modes: number of aerosol modes
	!>@param[in] n_sv: number of organic / volatility modes
	!>@param[inout] n_aer1: number mixing ratio in modes (kg-1 dry air)
	!>@param[inout] d_aer1: diameter in modes
	!>@param[inout] sig_aer1: geo std in modes
	!>@param[inout] molw_core1:molw in core
	!>@param[inout] density_core1: solute density
	!>@param[inout] nu_core1: van hoff factor
	!>@param[inout] org_content1: organic content in vol bins
	!>@param[inout] molw_org1: molw in volatility bins
	!>@param[inout] density_org1: density in volatility bins
	!>@param[inout] delta_h_vap1: enthalpy in volatility bins
	!>@param[inout] nu_org1: van hoff factor in volatility bins
	!>@param[inout] log_c_star1: log_c_star in volatility bins
	!>@param[inout] act_frac1: activated fraction in modes
	!>@param[inout] dcrit1: critical diameter in modes
	subroutine allocate_arrays(n_mode,n_sv,n_aer1,d_aer1,sig_aer1, &
			molw_core1,density_core1,nu_core1,org_content1, &
			molw_org1, density_org1,delta_h_vap1,nu_org1,log_c_star1, act_frac1, dcrit1)
	    use numerics_type
		implicit none
		integer(i4b), intent(in) :: n_mode, n_sv
		real(wp), dimension(:), allocatable, intent(inout) :: n_aer1,d_aer1,sig_aer1, &
							molw_core1, density_core1, nu_core1, org_content1, &
							molw_org1, density_org1, delta_h_vap1, nu_org1, log_c_star1, &
							act_frac1, dcrit1
		
		integer(i4b) :: AllocateStatus
		allocate( n_aer(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( d_aer(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sig_aer(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( n_aer1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( d_aer1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sig_aer1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( d_aer_new(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sgi(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_final(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( mass_initial(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( mass_final(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sd(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( b(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sm(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( eta(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( f1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( f2(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_core(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_core1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( molw_core(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( molw_core1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nu_core(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nu_core1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( act_frac(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( act_frac1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( act_frac2(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( dcrit1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		
		allocate( molw_org(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( molw_org1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( r_org(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( log_c_star(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( log_c_star1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( cstar(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( c_ions(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( epsilon1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( org_content(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( org_content1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_org(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_org1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nu_org(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nu_org1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( mass_org_condensed(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( delta_h_vap(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( delta_h_vap1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        n_aer=0._wp; d_aer=0._wp; sig_aer=0._wp; d_aer_new=0._wp
        sgi=0._wp; density_final=0._wp; mass_initial=0._wp; mass_final=0._wp
        sd=0._wp; b=0._wp; sm=0._wp; eta=0._wp; f1=0._wp; f2=0._wp
        density_core=0._wp; molw_core=0._wp; nu_core=0._wp
        act_frac=0._wp; act_frac1=0._wp; act_frac2=0._wp; dcrit1=0._wp
        molw_org=0._wp; r_org=0._wp; log_c_star=0._wp; cstar=0._wp
        c_ions=0._wp; epsilon1=0._wp; org_content=0._wp; density_org=0._wp
        nu_org=0._wp; mass_org_condensed=0._wp; delta_h_vap=0._wp

	end subroutine allocate_arrays
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>initialise arrays for activation code
	!>@param[in] n_modes: number of aerosol modes
	!>@param[in] n_sv: number of volatility bins
	!>@param[in] p1: pressure (Pa)
	!>@param[in] t1: temperature (K)
	!>@param[in] w1: vertical wind (m/s)
	!>@param[in] n_aer1: number mixing ratio in modes (kg-1 dry air)
	!>@param[in] d_aer1: diameter in modes
	!>@param[in] sig_aer1: logarithmic geometric width ln(sigma_g) in modes
	!>@param[in] molw_org1: molecular weight in volatility bins
	!>@param[in] density_core1: density in modes
	subroutine initialise_arrays(n_modes,n_sv,p1,t1,w1,n_aer1, &
                                d_aer1,sig_aer1, molw_org1,density_core1)
        use numerics_type
        implicit none
        integer(i4b), intent(in) :: n_modes, n_sv
        real(wp), intent(in) :: p1,t1,w1
        real(wp), dimension(n_modes), intent(in) :: n_aer1,d_aer1,sig_aer1, density_core1
        real(wp), dimension(n_sv), intent(in) :: molw_org1

        integer(i4b) :: i
        real(wp) :: rho_dry_local

        rhinit=0.999_wp
        pinit=p1
        tinit=t1
        w=w1
        pcb=pinit
        tcb=tinit

        ! Public aerosol number is a dry-air number mixing ratio (# kg-1).
        ! ARG/FN and the equilibrium partitioning calculation use volumetric
        ! concentrations, so convert once at the internal boundary.
        rho_dry_local=dry_air_density(t1,p1,rhinit)
        n_aer=n_aer1*rho_dry_local
        d_aer=d_aer1
        sig_aer=sig_aer1
        molw_org=molw_org1
        density_core=density_core1

        r_org=0._wp
        where(molw_org > 0._wp) r_org=r_gas/molw_org

        mass_initial=0._wp
        do i=1,n_modes
            if (n_aer(i) <= 0._wp) cycle
            density_dummy=density_core(i)
            n_dummy=n_aer(i)
            sd_dummy=sig_aer(i)
            d_dummy=d_aer(i)
            mass_initial(i)=n_dummy*exp(3._wp*log(d_dummy) + &
                            4.5_wp*sd_dummy**2)*density_dummy*pi/6._wp
        enddo

    end subroutine initialise_arrays
	
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>read in the data from the namelists for the BAM module
	!>@param[in] nmlfile
	subroutine read_in_bam_namelist(nmlfile)
        implicit none
        character(len=*), intent(in) :: nmlfile
        integer :: iu, ios
        character(len=512) :: iomsg

        namelist /bulk_aerosol_setup/ n_mode, n_sv, sv_flag, &
                    method_flag, giant_flag, a_eq_7, b_eq_7, &
                    ghosh_sigma_mode, ghosh_sigma_acc, ghosh_dmin, ghosh_dmax
        namelist /bulk_aerosol_spec/ n_aer1, d_aer1, sig_aer1, &
                    molw_core1, density_core1, nu_core1, org_content1, &
                    molw_org1, density_org1, delta_h_vap1, nu_org1, log_c_star1, &
                    p_test, t_test, w_test, rand_dist, n_rand, mean_w, sigma_w

        n_mode=0
        n_sv=0
        sv_flag=0
        method_flag=1
        giant_flag=0
        a_eq_7=0.21_wp
        b_eq_7=1.58_wp
        ghosh_sigma_mode=0
        ghosh_sigma_acc=1.6_wp
        ghosh_dmin=80.e-9_wp
        ghosh_dmax=1.e-6_wp
        ghosh_sigma_used=ghosh_sigma_acc

        open(newunit=iu,file=trim(nmlfile),status='old',action='read',iostat=ios,iomsg=iomsg)
        if (ios /= 0) error stop 'BAM: cannot open namelist: '//trim(iomsg)
        read(iu,nml=bulk_aerosol_setup,iostat=ios,iomsg=iomsg)
        if (ios /= 0) error stop 'BAM: error reading bulk_aerosol_setup: '//trim(iomsg)

        if (n_mode <= 0) error stop 'BAM: n_mode must be > 0'
        if (n_sv <= 0) error stop 'BAM: n_sv must be > 0'
        if (sv_flag < 0 .or. sv_flag > 1) error stop 'BAM: sv_flag must be 0 or 1'
        if (method_flag < 1 .or. method_flag > 4) error stop 'BAM: method_flag must be 1, 2, 3 or 4'
        if (giant_flag < 0 .or. giant_flag > 1) error stop 'BAM: giant_flag must be 0 or 1'
        if (ghosh_sigma_mode < 0 .or. ghosh_sigma_mode > 1) &
            error stop 'BAM: ghosh_sigma_mode must be 0 or 1'
        if (ghosh_sigma_acc <= 1._wp) error stop 'BAM: ghosh_sigma_acc must be > 1'
        if (ghosh_dmin <= 0._wp .or. ghosh_dmax <= ghosh_dmin) &
            error stop 'BAM: require 0 < ghosh_dmin < ghosh_dmax'
        if (method_flag.eq.4 .and. ghosh_sigma_mode.eq.0 .and. &
            (ghosh_sigma_acc < 1.4_wp .or. ghosh_sigma_acc > 2.1_wp)) &
            error stop 'BAM: fixed Ghosh sigma must be within calibrated range 1.4 to 2.1'

        call allocate_arrays(n_mode,n_sv,n_aer1,d_aer1,sig_aer1, &
            molw_core1,density_core1,nu_core1,org_content1, &
            molw_org1,density_org1,delta_h_vap1,nu_org1,log_c_star1, &
            act_frac1,dcrit2)

        ! Sentinels make missing required namelist entries fail clearly rather
        ! than propagating uninitialised memory into the activation equations.
        n_aer1=-1._wp
        d_aer1=-1._wp
        sig_aer1=-1._wp
        molw_core1=-1._wp
        density_core1=-1._wp
        nu_core1=-1._wp
        org_content1=0._wp
        molw_org1=-1._wp
        density_org1=-1._wp
        delta_h_vap1=0._wp
        nu_org1=-1._wp
        log_c_star1=0._wp
        p_test=-1._wp
        t_test=-1._wp
        w_test=0._wp
        rand_dist=.false.
        n_rand=1
        mean_w=0._wp
        sigma_w=1._wp

        read(iu,nml=bulk_aerosol_spec,iostat=ios,iomsg=iomsg)
        close(iu)
        if (ios /= 0) error stop 'BAM: error reading bulk_aerosol_spec: '//trim(iomsg)

        if (any(n_aer1 < 0._wp)) error stop 'BAM: all n_aer1 values are required and must be >= 0 kg-1'
        if (any(d_aer1 <= 0._wp)) error stop 'BAM: all d_aer1 values are required and must be > 0 m'
        if (any(sig_aer1 <= 0._wp)) error stop 'BAM: sig_aer1 is ln(sigma_g) and must be > 0'
        if (any(molw_core1 <= 0._wp)) error stop 'BAM: all molw_core1 values must be > 0 kg mol-1'
        if (any(density_core1 <= 0._wp)) error stop 'BAM: all density_core1 values must be > 0 kg m-3'
        if (any(nu_core1 <= 0._wp)) error stop 'BAM: all nu_core1 values must be > 0'
        if (p_test <= 0._wp .or. t_test <= 0._wp) error stop 'BAM: p_test and t_test must be > 0'
        if (sv_flag.eq.1) then
            if (any(org_content1 < 0._wp)) error stop 'BAM: org_content1 must be >= 0 microgram kg-1'
            if (any(molw_org1 <= 0._wp)) error stop 'BAM: molw_org1 must be > 0 when sv_flag=1'
            if (any(density_org1 <= 0._wp)) error stop 'BAM: density_org1 must be > 0 when sv_flag=1'
            if (any(nu_org1 <= 0._wp)) error stop 'BAM: nu_org1 must be > 0 when sv_flag=1'
        else
            ! These are unused with sv_flag=0; provide benign values so callers
            ! need not specify dummy semi-volatile properties.
            where(molw_org1 <= 0._wp) molw_org1=0.2_wp
            where(density_org1 <= 0._wp) density_org1=1500._wp
            where(nu_org1 <= 0._wp) nu_org1=1._wp
        endif

    end subroutine read_in_bam_namelist
	

	
	end module bam	

