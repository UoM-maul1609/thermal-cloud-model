	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Warm Micro-Physics Module (WMM):
	!>A 2-moment warm bulk cloud microphysics module for use with other models.
	!> compile using the Makefile. Requires linking to other wrapper model for execution.

	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>microphysics code for the different cloud models
    module w_micro_module
    use nrtype
    use bam, only : n_mode, n_sv, giant_flag, method_flag, sv_flag, &
        	n_aer1, d_aer1, sig_aer1, molw_core1, density_core1, nu_core1, org_content1, &
        	molw_org1, density_org1, delta_h_vap1, nu_org1, log_c_star1, p_test, t_test, &
    		w_test, act_frac1, smax1, dcrit2,  &
    		a_eq_7, b_eq_7, &
    		ctmm_activation, initialise_arrays, read_in_bam_namelist
    
    private
    public :: w_microphysics_2d, w_microphysics_1d, read_in_wmm_bam_namelist
    
    ! physical constants
    real(sp), parameter :: rhow=1000._sp, rhoi=920._sp,lv=2.5e6_sp,ls=2.8e6_sp,lf=ls-lv, &
    					   cp=1005._sp, cw=4187._sp, cice=2093._sp, r=8.314_sp, &
    						mw=18e-3_sp, ma=29e-3_sp, ra=r/ma,rv=r/mw, eps1=ra/rv, &
    						ttr=273.15_sp, joules_in_an_erg=1.0e-7_sp, &
    						joules_in_a_cal=4.187_sp
    						
    						
    ! mass-diameter and size spectra relations
    real(sp), parameter :: cr=523.6_sp, cc=523.6_sp, &
                        cs=52.36_sp, cg=261.8_sp, ci=104._sp, &
    					dr=3_sp, dc=3_sp, ds=3._sp, dg=3._sp, di=3._sp, &
    					alpha_r=2.5_sp, alpha_c=0.0_sp, & ! note, alpha_c is for a "mass" - number distribution
    					alpha_s=2.5_sp, alpha_g=2.5_sp, alpha_i=0._sp
    					
	! terminal fall-speed relations
	real(sp), parameter :: a_r=362._sp, a_c=362._sp, &
	                        a_s=4.84_sp, a_g=253._sp, a_i=71.34_sp, &
							b_r=0.65_sp, b_c=0.65_sp, &
							b_s=0.25_sp, b_g=0.734_sp, b_i=0.6635_sp, &
							f_r=0._sp, f_c=0._sp, f_s=0._sp, f_g=0._sp, f_i=0._sp
							
	! autoconversion
	real(sp), parameter :: aw0=1e-3_sp, dwa=20e-6_sp, nl=2.4e8_sp, &
						lw0=rhow*pi/6._sp*nl*dwa**3, &
						tsaut=60._sp, dimax=0.3e-3_sp, di2s=0.33e-3_sp, &
						lambda_imin=(1._sp+di+alpha_i)/dimax, &
						tsbreak=60._sp, lambda_s_break=1000._sp
	
    ! microphysical values:
    real(sp), parameter :: hm_rate=3.5e8_sp, nar=1.1e15_sp, nbr=0._sp, &
    						rho0=1.2_sp, bbigg=100._sp, abigg=0.66_sp
    real(sp) :: mi0=1.e-14_sp

	! coalescence efficiencies
	real(sp), parameter :: erw=1._sp, erg=1._sp, ers=1._sp, eri=1._sp, esw=1._sp, &
						egw=1._sp, eiw=1._sp, egs_wet=1._sp, egi_wet=1._sp
	
	! variables used in various process rates:
	real(sp) :: gam1r,gam2r,gam1c, gam2c, gam1i,gam2i, gam1s, gam2s,gam1g,gam2g, &
				fall_q_r, fall_q_c, fall_q_s, fall_q_g, fall_n_r, fall_n_s, fall_n_g, &
				fall_q_i, fall_n_i, fall_n_c, &
				phi_r, mass_iacr,num_iacr, mass_sacw_i, mass_iacw, &
				mass_racs1,mass_racs2,mass_racs3, &
				mass_racg1,mass_racg2,mass_racg3, &
				mass_sacr1,mass_sacr2,mass_sacr3, &
				mass_sacg1,mass_sacg2,mass_sacg3, &
				mass_gacr1,mass_gacr2,mass_gacr3, &
				mass_gacs1,mass_gacs2,mass_gacs3, &
				num_racs1,num_racs2,num_racs3, num_racg1,num_racg2,num_racg3, &
				num_sacg1, num_sacg2,num_sacg3, &
				mass_gacw, mass_gaci, &
				nu_r1,nu_r2,nu_i1, nu_i2, nu_s1, nu_s2, nu_g1, nu_g2, &
				mass_imm, num_imm, q0sat, &
				chi_rain, chi_cloud, chi_ice, chi_snow, chi_graupel, &
				chi_rain1, chi_cloud1, chi_ice1, chi_snow1, chi_graupel1
				
	! Seifert and Beheng autoconversion
	real(sp) :: kc, kr, xstar
				
	real(sp), dimension(3) :: c=[1._sp,2._sp,1._sp]
	integer(i4b) :: k
	real(sp) :: isnow, iice, f1,f2,a,b, qsmall=1e-30_sp
	
	
    contains


    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>read in the data from the namelists for the BAM module
	!> and set variables for microphysics
	!>@param[in] nmlfile
	!>@param[inout] q_name, q_type, c_s, c_e
	!>@param[inout] nq,ncat, nprec, iqv, iqc, inc
	subroutine read_in_wmm_bam_namelist(nmlfile, &
                q_name,q_type,c_s,c_e,nq,ncat,nprec, &
                iqv, iqc, inc)
		use bam, only : read_in_bam_namelist, n_mode
		implicit none
        character (len=200), intent(in) :: nmlfile
        integer(i4b), intent(inout) :: nq, ncat,nprec, iqv, iqc, inc
        integer(i4b), intent(inout), dimension(:), allocatable :: q_type, c_s, c_e
        character(len=20), dimension(:), allocatable :: q_name
        
        integer(i4b) :: i
        
        call read_in_bam_namelist(nmlfile)
        
        ncat=5
        
        
        nq=5 ! vapour, cloud mass, rain mass, cloud water, rain water
        
        allocate(q_name(nq))
        allocate(q_type(ncat))
        allocate(c_s(ncat))
        allocate(c_e(ncat))
        
        q_type(1)=0 ! vapour
        c_s(1)=1
        c_e(1)=1
        q_type(2)=1 ! cloud water
        c_s(2)=2
        c_e(2)=2
        q_type(3)=1 ! rain water
        c_s(3)=3
        c_e(3)=3
        q_type(4)=2 ! cloud water number
        c_s(4)=4
        c_e(4)=4
        q_type(5)=2 ! rain water number
        c_s(5)=5
        c_e(5)=5
        
        q_name=["qv","qc","qr","nc","nr"]
        
        nprec=1
        
        iqv=1
        iqc=2
        inc=4
        
	end subroutine read_in_wmm_bam_namelist




	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>initialises variables for use with the microphysics
    subroutine initialise_microphysics_vars
    use hypergeo, only : hygfx
    implicit none

	! used to calculate intercept and slopes
	gam1r=gamma(1._sp+alpha_r)
	gam2r=gamma(1._sp+alpha_r+dr)
	gam1c=gamma(1._sp+alpha_c)
	gam2c=gamma(1._sp+alpha_c+1._sp) ! note the 1, instead of dc - drop distribution
	                                ! is a mass distribution
	gam1i=gamma(1._sp+alpha_i)
	gam2i=gamma(1._sp+alpha_i+di)
	gam1s=gamma(1._sp+alpha_s)
	gam2s=gamma(1._sp+alpha_s+ds)
	gam1g=gamma(1._sp+alpha_g)
	gam2g=gamma(1._sp+alpha_g+dg)

    ! mass weighted fall for r, c, s, g, i
    fall_q_r=a_r*gamma(1._sp+alpha_r+dr+b_r) / gamma(1._sp+alpha_r+dr)
    fall_q_c=a_c*gamma(1._sp+alpha_c+1._sp+b_c) / gamma(1._sp+alpha_c+1._sp)
    fall_q_s=a_s*gamma(1._sp+alpha_s+ds+b_s) / gamma(1._sp+alpha_s+ds)
    fall_q_g=a_g*gamma(1._sp+alpha_g+dg+b_g) / gamma(1._sp+alpha_g+dg)
    fall_q_i=a_i*gamma(1._sp+alpha_i+di+b_i) / gamma(1._sp+alpha_i+di)

    ! number weighted fall for r, c, s, g
    fall_n_r=a_r*gamma(1._sp+alpha_r+b_r) / gamma(1._sp+alpha_r)
    fall_n_c=a_c*gamma(1._sp+alpha_c+b_c) / gamma(1._sp+alpha_c)
    fall_n_s=a_s*gamma(1._sp+alpha_s+b_s) / gamma(1._sp+alpha_s)
    fall_n_g=a_g*gamma(1._sp+alpha_g+b_g) / gamma(1._sp+alpha_g)
    fall_n_i=a_i*gamma(1._sp+alpha_i+b_i) / gamma(1._sp+alpha_i)
    
    ! sweep out of rain
    phi_r=pi*a_r*gamma(3._sp+b_r+alpha_r) / 4._sp
    
    ! ice accreting rain
    mass_iacr=pi*eri*a_r*cr*gamma(3._sp+b_r+dr+alpha_r)/4._sp
    num_iacr =pi*eri*a_r*gamma(3._sp+b_r+alpha_r)/4._sp

	! collection of cloud by snow and ice
	mass_sacw_i=pi*a_s*gamma(3._sp+b_s+alpha_s)/4._sp
	mass_iacw=pi*a_i*gamma(3._sp+b_i+alpha_i)/4._sp
	
	! collisions between precipitating particles of different species
	! rain-snow
	mass_racs1=gamma(1._sp+alpha_r)*gamma(3._sp+alpha_s+ds)
	mass_racs2=2._sp*gamma(2._sp+alpha_r)*gamma(2._sp+alpha_s+ds)
	mass_racs3=gamma(3._sp+alpha_r)*gamma(1._sp+alpha_s+ds)
	
	num_racs1=gamma(1._sp+alpha_r)*gamma(3._sp+alpha_s)
	num_racs2=2._sp*gamma(2._sp+alpha_r)*gamma(2._sp+alpha_s)
	num_racs3=gamma(3._sp+alpha_r)*gamma(1._sp+alpha_s)
    ! rain-graupel
	mass_racg1=gamma(1._sp+alpha_r)*gamma(3._sp+alpha_g+dg)
	mass_racg2=2._sp*gamma(2._sp+alpha_r)*gamma(2._sp+alpha_g+dg)
	mass_racg3=gamma(3._sp+alpha_r)*gamma(1._sp+alpha_g+dg)

	num_racg1=gamma(1._sp+alpha_r)*gamma(3._sp+alpha_g)
	num_racg2=2._sp*gamma(2._sp+alpha_r)*gamma(2._sp+alpha_g)
	num_racg3=gamma(3._sp+alpha_r)*gamma(1._sp+alpha_g)
    ! snow-rain
	mass_sacr1=gamma(1._sp+alpha_s)*gamma(3._sp+alpha_r+dr)
	mass_sacr2=2._sp*gamma(2._sp+alpha_s)*gamma(2._sp+alpha_r+dr)
	mass_sacr3=gamma(3._sp+alpha_s)*gamma(1._sp+alpha_r+dr)
    ! snow-graupel
	mass_sacg1=gamma(1._sp+alpha_s)*gamma(3._sp+alpha_g+dg)
	mass_sacg2=2._sp*gamma(2._sp+alpha_s)*gamma(2._sp+alpha_g+dg)
	mass_sacg3=gamma(3._sp+alpha_s)*gamma(1._sp+alpha_g+dg)

	num_sacg1=gamma(1._sp+alpha_s)*gamma(3._sp+alpha_g)
	num_sacg2=2._sp*gamma(2._sp+alpha_s)*gamma(2._sp+alpha_g)
	num_sacg3=gamma(3._sp+alpha_s)*gamma(1._sp+alpha_g)
    ! graupel-rain
	mass_gacr1=gamma(1._sp+alpha_g)*gamma(3._sp+alpha_r+dr)
	mass_gacr2=2._sp*gamma(2._sp+alpha_g)*gamma(2._sp+alpha_r+dr)
	mass_gacr3=gamma(3._sp+alpha_g)*gamma(1._sp+alpha_r+dr)
    ! graupel-snow
	mass_gacs1=gamma(1._sp+alpha_g)*gamma(3._sp+alpha_s+ds)
	mass_gacs2=2._sp*gamma(2._sp+alpha_g)*gamma(2._sp+alpha_s+ds)
	mass_gacs3=gamma(3._sp+alpha_g)*gamma(1._sp+alpha_s+ds)
	
	
	! accretion and riming by graupel
	mass_gacw=pi*egw*a_g*gamma(3._sp+b_g+alpha_g)/4._sp
	mass_gaci=pi*a_g*gamma(3._sp+b_g+alpha_g)/4._sp
    
    ! gauss hypergeometric equations aggregation of ice with ice (and snow with snow)
    ! See Ferrier (1994, JAS part 1, equation B.21)
    ! snow:
    a=1._sp
    b=4._sp+2._sp*alpha_s+b_s
    isnow=0._sp
    do k=1,3
	    call hygfx(a, b, real(k,sp)+alpha_s+1.0_sp,0.5_sp,f1)
	    call hygfx(a, b, real(k,sp)+alpha_s+b_s+1.0_sp, 0.5_sp,f2)
	    isnow=isnow+c(k)*(f1/(real(k,sp)+alpha_s)+f2/(real(k,sp)+alpha_s+b_s))
	enddo
	isnow=a_s*pi*gamma(b)/(2._sp**(6._sp+2._sp*alpha_s+b_s)) * isnow
	
    ! ice:
    a=1._sp
    b=4._sp+2._sp*alpha_i+b_i
    iice=0._sp
    do k=1,3
	    call hygfx(a, b, real(k,sp)+alpha_i+1.0_sp,0.5_sp,f1)
	    call hygfx(a, b, real(k,sp)+alpha_i+b_i+1.0_sp, 0.5_sp,f2)
	    iice=iice+c(k)*(f1/(real(k,sp)+alpha_i)+f2/(real(k,sp)+alpha_i+b_i))
	enddo
	iice=a_i*pi*gamma(b)/(2._sp**(6._sp+2._sp*alpha_i+b_i)) * iice
	
	
	! ventilation
	! rain:
	nu_r1=0.78_sp*gamma(2._sp+alpha_r)
	nu_r2=gamma(0.5_sp*b_r+alpha_r+2.5)
	! ice:
	nu_i1=0.78_sp*gamma(2._sp+alpha_i)
	nu_i2=0.31_sp*gamma(0.5_sp*b_i+alpha_i+2.5)
	! snow:
	nu_s1=0.78_sp*gamma(2._sp+alpha_s)
	nu_s2=0.31_sp*gamma(0.5_sp*b_s+alpha_s+2.5)
	! graupel:
	nu_g1=0.78_sp*gamma(2._sp+alpha_g)
	nu_g2=0.31_sp*gamma(0.5_sp*b_g+alpha_g+2.5)
	
	! immersion freezing by bigg
	mass_imm=gamma(4._sp+dr+alpha_r)*pi*cr*bbigg/6._sp
	num_imm=gamma(4._sp+alpha_r)*pi*bbigg/6._sp
	
	
	! precipitation
	chi_rain=gamma(1._sp+alpha_r+b_r+dr)
	chi_cloud=gamma(1._sp+alpha_c+b_c+1._sp)
	chi_ice=gamma(1._sp+alpha_i+b_i+di)
	chi_snow=gamma(1._sp+alpha_s+b_s+ds)
	chi_graupel=gamma(1._sp+alpha_g+b_g+dg)
	
	chi_rain1=gamma(1._sp+alpha_r+dr)
	chi_cloud1=gamma(1._sp+alpha_c+1._sp)
	chi_ice1=gamma(1._sp+alpha_i+di)
	chi_snow1=gamma(1._sp+alpha_s+ds)
	chi_graupel1=gamma(1._sp+alpha_g+dg)
	

    ! Seifert and Beheng autoconversion:
    kc=9.44e9_sp ! m3 kg-2 s-1
    kr=5.78e0_sp ! m3 kg-2 s-1
    xstar=2.6e-10_sp ! kg
    end subroutine initialise_microphysics_vars
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>solves one time-step of the microphysics - calls w_microphysics_1d
	!>@param[in] nq: number of q-fields
	!>@param[in] ip: number of horizontal levels
	!>@param[in] kp: number of vertical levels
	!>@param[in] dt: time-step
	!>@param[in] dz: dz
	!>@param[in] o_halo: extra points for advection
	!>@param[inout] q: q-variables 
	!>@param[inout] precip: precip in rain, snow, graupel, ice cats - diagnostic
	!>@param[inout] theta: theta 
	!>@param[inout] p: pressure
	!>@param[inout] z: vertical levels 
	!>@param[inout] t: temperature 
	!>@param[inout] rho: density 
	!>@param[in] w: vertical wind 
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag: switch hm-process on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] theta_flag: whether to alter theta
    subroutine w_microphysics_2d(nq,ip,kp,o_halo,dt,dz,q,precip,theta,p, z,t,rho,w, &
    						micro_init,hm_flag, mass_ice, theta_flag)
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, ip,kp, o_halo
    real(sp), intent(in) :: dt,dz
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo,nq), intent(inout) :: q
    real(sp), dimension(1:kp,1:ip,1), intent(inout) :: precip
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: &
    					theta, p, t, rho
    real(sp), dimension(-o_halo+1:kp+o_halo) :: z
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(in) :: w
    logical, intent(in) :: hm_flag, theta_flag
    logical , intent(inout) :: micro_init
    real(sp), intent(in) :: mass_ice

	! locals
	integer(i4b) :: i
	
	do i=1,ip
		call w_microphysics_1d(nq,kp,o_halo,dt,dz,q(:,i,:),precip(:,i,:),theta(:,i),p(:,i), &
							z(:),t(:,i),rho(:,i),w(:,i), &
    						micro_init,hm_flag, mass_ice, theta_flag)
	enddo


	end subroutine w_microphysics_2d
	
	    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>solves one time-step of the microphysics
	!>@param[in] nq: number of q-fields
	!>@param[in] kp: number of vertical levels
	!>@param[in] dt: time-step
	!>@param[in] dz: dz
	!>@param[in] o_halo: extra points for advection
	!>@param[inout] q: q-variables 
	!>@param[inout] precip: precip in rain, snow, graupel, ice cats - diagnostic
	!>@param[inout] theta: theta 
	!>@param[inout] p: pressure
	!>@param[inout] z: vertical levels 
	!>@param[inout] t: temperature 
	!>@param[inout] rho: density 
	!>@param[in] u: vertical wind 
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag: switch hm-process on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] theta_flag: whether to alter theta
    subroutine w_microphysics_1d(nq,kp,o_halo,dt,dz,q,precip,theta,p, z,t,rho,u, &
    						micro_init,hm_flag, mass_ice,theta_flag)
	use advection_1d
	use nr, only : dfridr
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, kp, o_halo
    real(sp), intent(in) :: dt,dz
    real(sp), dimension(-o_halo+1:kp+o_halo,nq), intent(inout) :: q
    real(sp), dimension(1:kp,1), intent(inout) :: precip
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(inout) :: theta, p, z, t, rho
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: u
    logical, intent(in) :: hm_flag, theta_flag
    logical , intent(inout) :: micro_init
    real(sp), intent(in) :: mass_ice
    ! locals:
    integer(i4b) :: k,k1,iter, n_step
    real(sp) :: temp, qtot,qaut, a, b, ab_ice, ab_liq, ice_dep,snow_dep,graup_dep, &
    			nu_ice, nu_snow, nu_graup, diff1, ktherm1, tc, nu_vis, sc, nu_rain, rain_evap, &
    			sb_aut, sb_acr, sb_cwaut, sb_cwacr, sb_raut, sb_rsel, sb_cwsel
    real(sp), dimension(kp) :: smr, smr_i
    
    real(sp), dimension(kp) :: &
    			pgaci,  & ! accretion of cloud ice by graupel
				pgacr, & ! riming of graupel by rain
				pgacs,rgacs, & ! accretion of snow by graupel
				pgacw, & ! riming of graupel by liquid cloud
				pgaut, & ! autoconversion from snow to graupel due to heavy riming
				pgdep, & ! deposition of vapour onto graupel
				pgfr,rgfr,  & ! freezing of rain to form graupel
				pgmlt, & ! melting of graupel to form rain
				pgshd, & ! shedding of liquid from wet graupel to rain
				pgsub, & ! sublimation of graupel or evaporation from wet graupel
				riaci, & ! aggregation of ice crystals
				piacr_g,riacr_g, & ! riming of cloud ice by large rain drops to form graupel
				piacr_s,riacr_s, & ! riming of cloud ice by small rain drops to form snow
				piacw, & ! riming of ice cloud ice by liquid cloud 
				picnt, & ! nucleation of ice crystals by contact freezing
				pidep, & ! deposition of water vapour onto cloud ice
				piprm, & ! primary nucleation of ice crystals by INPs
				pifrw, & ! nucleation of ice crystals by homogeneous freezing of cloud
				pihal, & ! production of ice crystals by hm process
				pimlt, & ! cloud ice melting to form rain
				pisub, & ! sublimation of cloud ice
				praci_g, & ! accretion of cloud ice by rain to form graupel
				praci_s, & ! accretion of cloud ice by rain to form snow
				pracs, & ! accretion of snow by rain to form graupel
				pracw, & ! accretion of liquid cloud by rain
				praut, & ! autoconversion from liquid cloud to rain (coalescence)
				prevp, & ! evaporation of rain
				psacr,rsacr, & ! accretion of rain by snow to form graupel
				psaci, & ! accretion of cloud ice by snow
				rsacs, & ! aggregation of snowflakes
				psacw, & ! accretion of liquid cloud by snow
				psaut,rsaut, & ! conversion of large ice crystals to snow
				rsbrk, & ! break-up of large snowflakes
				psdep, & ! deposition of vapour onto snow
				psmlt, & ! melting of snow to form rain
				pssub, & ! sublimation of snow
				rcwaut, & ! reduction in cloud number due to autoconversion
				rcwacr, & ! cloud water number accreted onto rain
				rraut, &  ! increase in rain number due to autoconversion
    			rrsel, &     ! rain self accretion - number
    			rcwsel     ! cloud water self accretion - number
    				    
    real(sp) :: pgwet ! amount of liquid that graupel can freeze without shedding
    								

    real(sp), dimension(kp) :: n_r, lam_r, n_i, lam_i, n_s, lam_s, n_g, lam_g, lam_c, n_c
    real(sp), dimension(kp) :: rho_fac
	real(sp), dimension(-o_halo:kp+o_halo) :: vqr, vqs, vqg, vqi, vnr, vns, vng, vni, &
	                                        vqc, vnc
	! coalescence efficiencies
	real(sp), dimension(kp) :: egi_dry, egs_dry, esi, eii, ess
	real(sp) :: qold,des_dt,dqs_dt,err,cond,temp1
	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialise some variables that do not depend on prognostics                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (micro_init) then 
    	call initialise_microphysics_vars
    	mi0=mass_ice
    	micro_init=.false.
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! zero arrays                                                                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	pgaci=0._sp
	pgacr=0._sp
	pgacs=0._sp
	rgacs=0._sp
	pgacw=0._sp
	pgaut=0._sp
	pgdep=0._sp
	pgfr=0._sp
	rgfr=0._sp
	pgmlt=0._sp
	pgshd=0._sp
	pgsub=0._sp
	riaci=0._sp
	piacr_g=0._sp
	riacr_g=0._sp
	piacr_s=0._sp
	riacr_s=0._sp
	piacw=0._sp
	picnt=0._sp
	pidep=0._sp
	piprm=0._sp
	pifrw=0._sp
	pihal=0._sp
	pimlt=0._sp
	pisub=0._sp
	praci_g=0._sp
	praci_s=0._sp
	pracs=0._sp
	pracw=0._sp
	praut=0._sp
	prevp=0._sp
	psacr=0._sp
	rsacr=0._sp
	psaci=0._sp
	rsacs=0._sp
	psacw=0._sp
	psaut=0._sp
	rsaut=0._sp
	rsbrk=0._sp
	psdep=0._sp
	psmlt=0._sp
	pssub=0._sp
	rcwaut=0._sp
    rcwacr=0._sp
    rraut=0._sp
    rrsel=0._sp
    rcwsel=0._sp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! some commonly used variables that depend on prognostics                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t=theta*(p/1.e5_sp)**(ra/cp) ! temperature
    rho=p / (ra*t) ! air density    
    rho_fac=(rho0/rho(1:kp))**0.5_sp
    ! rain n0, lambda
    lam_r=(max(q(1:kp,5),1._sp)*cr*gam2r / (max(q(1:kp,3),1.e-10_sp)*gam1r))**(1._sp/dr)
    n_r=rho(1:kp)*max(q(1:kp,5),0._sp)*lam_r**(1._sp+alpha_r) / gam1r
    ! cloud n0, lambda    
    lam_c=(max(q(1:kp,4),1._sp)*cc*gam2c / (max(q(1:kp,2),1.e-10_sp)*gam1c))**(1._sp/1._sp)
    n_c=rho(1:kp)*max(q(1:kp,4),0._sp)*lam_c**(1._sp+alpha_c) / gam1c

    
    ! precipitation
	precip(1:kp,1)=cr*n_r*(a_r*chi_rain/(lam_r**(alpha_r+b_r+dr+1._sp)) - &
					u(1:kp)*chi_rain1/(lam_r**(alpha_r+dr+1._sp))) &
					/rho(1:kp) *3600._sp
    
    ! fall speeds
    ! rain
    vqr(1:kp)=max(fall_q_r*rho_fac * lam_r**(1._sp+alpha_r+dr) / &
    	(lam_r+f_r)**(1._sp+alpha_r+dr+b_r), 0._sp)
    
    vnr(1:kp)=max(fall_n_r*rho_fac * lam_r**(1._sp+alpha_r) / &
    	(lam_r+f_r)**(1._sp+alpha_r+b_r), 0._sp)
    
    ! cloud
    vqc(1:kp)=max(fall_q_c*rho_fac * lam_c**(1._sp+alpha_c+1._sp) / &
    	(lam_c+f_c)**(1._sp+alpha_c+1._sp+b_c), 1.e-3_sp)
    
    vnc(1:kp)=max(fall_n_c*rho_fac * lam_c**(1._sp+alpha_c) / &
    	(lam_c+f_c)**(1._sp+alpha_c+b_c), 1.e-3_sp)
    ! coalescence efficiencies
    egi_dry=0.2_sp*exp(0.08*(t(1:kp)-ttr))
    egs_dry=0.2_sp*exp(0.08*(t(1:kp)-ttr))
    esi=0.2_sp*exp(0.08*(t(1:kp)-ttr))
    eii=0.2_sp*exp(0.08*(t(1:kp)-ttr))
    ess=0.2_sp*exp(0.08*(t(1:kp)-ttr))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
     
    
   
    
    ! loop over all levels
    do k=1,kp
 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! condensation of liquid water                                                   !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! smr at 0c
		!q0sat=eps1*svp_liq(ttr)/(p(k)-svp_liq(ttr))
		q0sat=eps1*svp_liq(ttr)/(p(k)-svp_liq(ttr))
    	smr(k)=eps1*svp_liq(t(k))/(p(k)-svp_liq(t(k))) ! saturation mixing ratio

        des_dt=dfridr(svp_liq,t(k),1.e0_sp,err)
        dqs_dt=eps1*p(k)*des_dt/(p(k)-svp_liq(t(k)))**2
        qold=q(k,2)
        qtot=q(k,1)+q(k,2)
		
        q(k,2)=q(k,1)+q(k,2)-smr(k)
        if (theta_flag) q(k,2)=(q(k,2)+(lv/cp*qold)*dqs_dt) / (1._sp+lv/cp*dqs_dt)
        q(k,2)=max(q(k,2),0._sp)
        t(k)=t(k)
        if(theta_flag) t(k)=t(k)+lv/cp*(q(k,2)-qold)
		
		tc=t(k)-ttr
    	smr(k)=eps1*svp_liq(t(k))/(p(k)-svp_liq(t(k))) ! saturation mixing ratio
    	q0sat=smr(k)	
    	smr_i(k)=eps1*svp_ice(t(k))/(p(k)-svp_ice(t(k))) ! saturation mixing ratio - ice	
    	
    	cond=(q(k,2)-qold)
    	q(k,1)=q(k,1)-cond
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! inhomogeneous mixing -https://journals.ametsoc.org/doi/pdf/10.1175/2007JAS2374.1
 !        if(q(2,k)<qold) then
!             q(4,k)=q(4,k)*(q(2,k)/qold)**1._sp
!         endif
        




        k1=max(k-1,1)
	    if((q(k,2) .gt. qsmall) .and. (q(k1,2) .lt. qsmall)) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Bulk Aerosol Activation - number of drops
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            p_test=p(k)
            t_test=t(k)
            w_test=max(u(k),0.001_sp)
            call initialise_arrays(n_mode,n_sv,p_test,t_test,w_test, &
                        n_aer1,d_aer1,sig_aer1, molw_org1,density_core1)
        
            call ctmm_activation(n_mode,n_sv,sv_flag, &
                        n_aer1, d_aer1,sig_aer1,molw_core1, &
                        density_core1, nu_core1, &
                        org_content1,molw_org1, density_org1, delta_h_vap1, nu_org1,  &
                        log_c_star1, &
                        w_test, t_test,p_test, a_eq_7, b_eq_7, &
                        act_frac1,smax1,dcrit2)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            temp1=sum(n_aer1*act_frac1)
            !temp1=10.e6_sp
!             q(k-1,4)=temp1
            q(k,4)=temp1
!             q(k+1,4)=temp1
!             q(k+2,4)=temp1
        endif      
        
    





		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! evaporation of rain                                                            !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		diff1=dd(t(k),p(k))
		ktherm1=ka(t(k))
		ab_liq=lv**2 / (ktherm1*rv*t(k)**2) + 1._sp/(rho(k)*smr(k)*diff1)
		nu_vis=viscosity_air(t(k)) / rho(k)
		sc=nu_vis / diff1
		nu_rain=2._sp*pi*n_r(k) / rho(k) * &
				(nu_r1 / lam_r(k)**(2._sp+alpha_r) + &
				(a_r/nu_vis)**0.5_sp*sc**(1._sp/3._sp)* &
				(rho(k)*rho0)**0.25_sp*nu_r2 / &
				(lam_r(k)+0.5_sp*f_r)**(0.5_sp*b_r+alpha_r+2.5_sp))

	
		rain_evap=(q(k,1)/smr(k)-1._sp) / (rho(k)*ab_liq)*nu_rain
		if(q(k,1).gt.smr(k)) then
			prevp(k)=0._sp
		else
			prevp(k)=-min(rain_evap,0._sp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! warm rain autoconversion based on Seifert and Beheng (2006)                    !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call seifert_beheng(sb_aut,sb_acr, sb_cwaut, sb_cwacr, sb_raut, &
		                    sb_rsel, sb_cwsel, q(k,2),q(k,4),q(k,3),q(k,5),rho(k),dt)
		praut(k)=sb_aut
		pracw(k)=sb_acr
		rcwaut(k)=sb_cwaut
		rcwacr(k)=sb_cwacr
		rraut(k)=sb_raut
		rrsel(k)=sb_rsel
		rcwsel(k)=sb_cwsel
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    enddo
    
    
    
    
    ! update variables
    ! vapour mass
    q(1:kp,1)=q(1:kp,1)+(pgsub+pssub+pisub-(psdep+pidep+piprm+pgdep))*dt
    ! liquid mass
    q(1:kp,2)=q(1:kp,2)-((pgacw+praut+psacw+pracw+piacw+pihal+picnt+pifrw))*dt
    ! rain mass
    q(1:kp,3)=q(1:kp,3)+(pgmlt+praut+pgshd+pracw+psmlt+pimlt- &
    			(pgacr+pgfr+psacr+piacr_g+piacr_s))*dt
    prevp=min(prevp,q(1:kp,3)/dt)
    t(1:kp)=t(1:kp)-lv/cp*prevp*dt
    q(1:kp,3)=q(1:kp,3)-prevp*dt
    q(1:kp,1)=q(1:kp,1)+prevp*dt
    
    ! liquid number
    q(1:kp,4)=q(1:kp,4)+(rcwaut+rcwacr+rcwsel)*dt
    ! rain number
    q(1:kp,5)=q(1:kp,5)+(rraut+rrsel-prevp*(q(1:kp,5)/(q(1:kp,3)+qsmall)))*dt
    
    where(q(:,2) .lt. qsmall)
        q(:,4) = 0.0_sp
    end where
    where(q(:,3) .lt. qsmall)
        q(:,5) = 0.0_sp
    end where

    q=max(q,0._sp)	 
    if (theta_flag) theta=t*(1.e5_sp/p)**(ra/cp)



    
 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! advection rain 0th order Bott, a.k.a. upstream advection                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! rain 
    if(sum(q(1:kp,3)).gt.qsmall) then
		where(isnan(vqr))
			vqr=0_sp
		end where
		vqr(-o_halo:0)=vqr(1)
		vqr(kp+1:kp+o_halo)=vqr(kp)
		n_step=max(ceiling(maxval(vqr)*dt/dz*2_sp),1)
		vqr(-o_halo:kp+o_halo-1)=-vqr(-o_halo+1:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vqr,q(:,3),.false.)
		enddo
		where(isnan(vnr))
			vnr=0_sp
		end where
		vnr(-o_halo:0)=vnr(1)
		vnr(kp+1:kp+o_halo)=vnr(kp)
		n_step=max(ceiling(maxval(vnr)*dt/dz*2_sp),1)
		vnr(-o_halo:kp+o_halo-1)=-vnr(-o_halo+1:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vnr,q(:,5),.false.)
		enddo
	endif
    ! cloud 
    if(sum(q(1:kp,2)).gt.qsmall) then
		where(isnan(vqc))
			vqc=0_sp
		end where
		vqc(-o_halo:0)=vqc(1)
		vqc(kp+1:kp+o_halo)=vqc(kp)
		n_step=max(ceiling(maxval(vqc)*dt/dz*2_sp),1)
		vqc(-o_halo:kp+o_halo-1)=-vqc(-o_halo+1:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vqc,q(:,2),.false.)
		enddo
		where(isnan(vnc))
			vnc=0_sp
		end where
		vnc(-o_halo:0)=vnc(1)
		vnc(kp+1:kp+o_halo)=vnc(kp)
		n_step=max(ceiling(maxval(vnc)*dt/dz*2_sp),1)
		vnc(-o_halo:kp+o_halo-1)=-vnc(-o_halo+1:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vnc,q(:,4),.false.)
		enddo
	endif
 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
    
    

    end subroutine w_microphysics_1d
    
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates variables for Seifert and Beheng (2001) autoconversion scheme
	!>@param[inout] praut, pracw, rcwaut, rcwacr, rraut, rrsel, rwsel: variables for SB 
	!>@param[in] qc: qcloud
	!>@param[in] nc: ncloud
	!>@param[in] qr: qrain
	!>@param[in] nr: nrain
	!>@param[in] rho: density of air
	!>@param[in] dt: timestep
    subroutine seifert_beheng(praut,pracw, rcwaut, rcwacr, rraut, &
		                    rrsel, rcwsel, qc,nc,qr,nr,rho,dt)
	use advection_1d
    implicit none
    ! arguments:
    real(sp), intent(inout) :: praut,pracw, rcwaut, rcwacr, rraut, &
		                    rrsel, rcwsel
	real(sp), intent(in) :: qc,nc,qr,nr,rho, dt
	real(sp) :: lc, lr, nc1,nr1,xc_bar, phi_au, phi_ac, b_slope, tau, factor1, factor2, &
	            test
	
	
	    praut=0._sp
	    pracw=0._sp
	    rcwaut=0._sp
	    rcwacr=0._sp
	    rraut=0._sp
	    rrsel=0._sp
	    rcwsel=0._sp
	    
	    lc=qc*rho
	    lr=qr*rho
	    
	    nc1=max(nc*rho,lc/xstar)
	    nr1=max(nr*rho,1._sp)
	    
		b_slope=((nc1+1.e-20_sp)/(lc+1.e-20_sp))*(gam2c+1.e-20_sp)/(gam1c +1.e-20_sp)
		xc_bar=gam2c/(gam1c*b_slope+1.e-20_sp)
		
		tau=1._sp-lc/(lc+lr+1.e-20_sp)
		tau=max(tau,1.e-6_sp)
		phi_au=600._sp*tau**0.68*(1._sp-tau**0.68)**3
		phi_ac=(tau/(tau+5.e-4_sp))**4

		if (lc .gt. qsmall) then
            ! autoconversion: equation a1 in Seifert and Beheng (2001, atmos res)
            praut = kc/(20._sp*xstar)*(alpha_c+2._sp)*(alpha_c+4._sp)/(alpha_c+1._sp)**2 * &
                    (lc*xc_bar)**2*(1._sp+phi_au/(1._sp-tau+1.e-20_sp)**2)
            
            ! accretion: equation a2 in Seifert and Beheng (2001, atmos res)
            pracw=kr*lc*lr*phi_ac
        
            ! cloud number autoconversion: equation a5 in Seifert and Beheng (2001, atmos res)
            rcwaut=-2._sp/xstar*praut
            ! cloud num accretion: equation a6 in Seifert and Beheng (2001, atmos res)
            rcwacr=-1._sp/xc_bar*pracw
            ! rain num autoconversion: equation a7 in Seifert and Beheng (2001, atmos res)
            rraut=-1._sp/2._sp*rcwaut
            ! rain num self collection: equation a8 in Seifert and Beheng (2001, atmos res)
            rrsel=-kr*nr*lr
            ! cloud num self collection: equation a9 in Seifert and Beheng (2001, atmos res)
            rcwsel=-kr*(alpha_c+2._sp)/(alpha_c+1._sp)*lc**2-rcwaut

            factor1=min(lc/dt,praut+pracw)/(praut+pracw)
            factor2=min(nc1/dt,-(rcwaut+rcwacr+rcwsel))/(-(rcwaut+rcwacr+rcwsel))

            praut=praut*factor1
            pracw=pracw*factor1
            rcwaut=rcwaut*factor2
            rcwacr=rcwacr*factor2
            rcwsel=rcwsel*factor2


            
            praut=praut/rho
            pracw=pracw/rho
            rcwaut=rcwaut/rho
            rcwacr=rcwacr/rho
            rraut=rraut/rho   
            rrsel=rrsel/rho
            rcwsel=rcwsel/rho
        endif
        

    end subroutine seifert_beheng  
    
    
      
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! saturation vapour pressure over liquid                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the saturation vapour pressure over liquid water according to buck fit
	!>@param[in] t: temperature
	!>@return svp_liq: saturation vapour pressure over liquid water
	function svp_liq(t)
		use nrtype
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: svp_liq
		svp_liq = 100._sp*6.1121_sp* &
			  exp((18.678_sp - (t-ttr)/ 234.5_sp)* &
			  (t-ttr)/(257.14_sp + (t-ttr)))
	end function svp_liq


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! saturation vapour pressure over ice										   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the saturation vapour pressure over ice water according to buck fit
	!>@param[in] t: temperature
	!>@return svp_ice: saturation vapour pressure over ice water
	function svp_ice(t)
		use nrtype
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: svp_ice
		svp_ice = 100._sp*6.1115_sp* &
			  exp((23.036_sp - (t-ttr)/ 333.7_sp)* &
			  (t-ttr)/(279.82_sp + (t-ttr)))
	end function svp_ice

    
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Viscosity of air - Page 417 Pruppacher and Klett							   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the viscosity of air vs temperature
	!>@param[in] t: temperature
	!>@return viscosity_air: viscosity of air
	function viscosity_air(t)
		use nrtype
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: viscosity_air
		real(sp) :: tc

		tc = t-ttr
		tc = max(tc,-200._sp)

		if( tc.ge.0._sp) then
			viscosity_air = (1.718_sp+0.0049_sp*tc) * 1E-5_sp ! the 1d-5 converts from poise to si units
		else
			viscosity_air = (1.718_sp+0.0049_sp*tc-1.2e-5_sp*tc**2) * 1e-5_sp
		end if
	end function viscosity_air

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
	  use nrtype
	  implicit none
	  real(sp), intent(in) :: t, p
	  real(sp) :: dd, t1
	  t1=max(t,200._sp)
	  dd=2.11e-5_sp*(t1/ttr)**1.94_sp*(101325_sp/p)
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
	  use nrtype
	  implicit none
	  real(sp), intent(in) :: t
	  real(sp) :: ka, t1
	  t1=max(t,200._sp)
	  ka=(5.69_sp+0.017_sp*(t1-ttr))*1e-3_sp*joules_in_a_cal
	end function ka


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! hm factor                                                                    !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the factor to multiply the peak production rate by
	!>@param[in] t: temperature
	!>@return hm_func: factor to multiple hm by
	function hm_func(t)
		use nrtype
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: hm_func
		hm_func=(min(max((t-265.65) / 2.5_sp,0._sp),1._sp) + &
		        min(max((270.65-t) / 2.5_sp,0._sp),1._sp)) -1.0_sp
	end function hm_func


    end module w_micro_module
    
