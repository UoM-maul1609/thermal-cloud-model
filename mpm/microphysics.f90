	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2017
	!>@brief
	!>Micro-Physics Module (MPM):
	!>A simple bulk cloud microphysics module for use with other models.
	!> compile using the Makefile. Requires linking to other wrapper model for execution.

	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>microphysics code for the simple cloud model
    module micro_module
    use nrtype
    
    private
    public :: microphysics_2d, microphysics_1d
    
    ! physical constants
    real(sp), parameter :: rhow=1000._sp, rhoi=920._sp,lv=2.5e6_sp,ls=2.8e6_sp,lf=ls-lv, &
    					   cp=1005._sp, cw=4187._sp, cice=2093._sp, r=8.314_sp, &
    						mw=18e-3_sp, ma=29e-3_sp, ra=r/ma,rv=r/mw, eps1=ra/rv, &
    						ttr=273.15_sp, joules_in_an_erg=1.0e-7_sp, &
    						joules_in_a_cal=4.187_sp
    						
    						
    ! mass-diameter and size spectra relations
    real(sp), parameter :: cr=523.6_sp, cs=52.36_sp, cg=261.8_sp, ci=104._sp, &
    					dr=3_sp, ds=3._sp, dg=3._sp, di=3._sp, &
    					alpha_r=2.5_sp, alpha_s=2.5_sp, alpha_g=2.5_sp, alpha_i=0._sp
    					
	! terminal fall-speed relations
	real(sp), parameter :: a_r=362._sp, a_s=4.84_sp, a_g=253._sp, a_i=71.34_sp, &
							b_r=0.65_sp, b_s=0.25_sp, b_g=0.734_sp, b_i=0.6635_sp, &
							f_r=0._sp, f_s=0._sp, f_g=0._sp, f_i=0._sp
							
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
	real(sp) :: gam1r,gam2r,gam1i,gam2i, gam1s, gam2s,gam1g,gam2g, &
				fall_q_r, fall_q_s, fall_q_g, fall_n_r, fall_n_s, fall_n_g, &
				fall_q_i, fall_n_i,  &
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
				chi_rain, chi_ice, chi_snow, chi_graupel, &
				chi_rain1, chi_ice1, chi_snow1, chi_graupel1
				
	real(sp), dimension(3) :: c=[1._sp,2._sp,1._sp]
	integer(i4b) :: k
	real(sp) :: isnow, iice, f1,f2,a,b, qsmall=1e-30_sp
	
	
    contains
    
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
	gam1i=gamma(1._sp+alpha_i)
	gam2i=gamma(1._sp+alpha_i+di)
	gam1s=gamma(1._sp+alpha_s)
	gam2s=gamma(1._sp+alpha_s+ds)
	gam1g=gamma(1._sp+alpha_g)
	gam2g=gamma(1._sp+alpha_g+dg)

    ! mass weighted fall for r, s, g, i
    fall_q_r=a_r*gamma(1._sp+alpha_r+dr+b_r) / gamma(1._sp+alpha_r+dr)
    fall_q_s=a_s*gamma(1._sp+alpha_s+ds+b_s) / gamma(1._sp+alpha_s+ds)
    fall_q_g=a_g*gamma(1._sp+alpha_g+dg+b_g) / gamma(1._sp+alpha_g+dg)
    fall_q_i=a_i*gamma(1._sp+alpha_i+di+b_i) / gamma(1._sp+alpha_i+di)

    ! number weighted fall for r, s, g
    fall_n_r=a_r*gamma(1._sp+alpha_r+b_r) / gamma(1._sp+alpha_r)
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
	chi_ice=gamma(1._sp+alpha_i+b_i+di)
	chi_snow=gamma(1._sp+alpha_s+b_s+ds)
	chi_graupel=gamma(1._sp+alpha_g+b_g+dg)
	
	chi_rain1=gamma(1._sp+alpha_r+dr)
	chi_ice1=gamma(1._sp+alpha_i+di)
	chi_snow1=gamma(1._sp+alpha_s+ds)
	chi_graupel1=gamma(1._sp+alpha_g+dg)
	
    end subroutine initialise_microphysics_vars
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>solves one time-step of the microphysics - calls microphysics_1d
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
    subroutine microphysics_2d(nq,ip,kp,o_halo,dt,dz,q,precip,theta,p, z,t,rho,w, &
    						micro_init,hm_flag, mass_ice,theta_flag)
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, ip,kp, o_halo
    real(sp), intent(in) :: dt,dz
    real(sp), dimension(nq,-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: q
    real(sp), dimension(4,1:kp,1:ip), intent(inout) :: precip
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
		call microphysics_1d(nq,kp,o_halo,dt,dz,q(:,:,i),precip(:,:,i),theta(:,i),p(:,i), &
							z(:),t(:,i),rho(:,i),w(:,i), &
    						micro_init,hm_flag, mass_ice,theta_flag)
	enddo


	end subroutine microphysics_2d
	
	    
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
    subroutine microphysics_1d(nq,kp,o_halo,dt,dz,q,precip,theta,p, z,t,rho,u, &
    						micro_init,hm_flag, mass_ice, theta_flag)
	use advection_1d
	use nr, only : dfridr
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, kp, o_halo
    real(sp), intent(in) :: dt,dz
    real(sp), dimension(nq,-o_halo+1:kp+o_halo), intent(inout) :: q
    real(sp), dimension(4,1:kp), intent(inout) :: precip
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(inout) :: theta, p, z, t, rho
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: u
    logical, intent(in) :: hm_flag, theta_flag
    logical , intent(inout) :: micro_init
    real(sp), intent(in) :: mass_ice
    ! locals:
    integer(i4b) :: k,iter, n_step
    real(sp) :: temp, qtot,qaut, a, b, ab_ice, ab_liq, ice_dep,snow_dep,graup_dep, &
    			nu_ice, nu_snow, nu_graup, diff1, ktherm1, tc, nu_vis, sc, nu_rain, rain_evap
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
				pssub    ! sublimation of snow
    				    
    real(sp) :: pgwet ! amount of liquid that graupel can freeze without shedding
    								

    real(sp), dimension(kp) :: n_r, lam_r, n_i, lam_i, n_s, lam_s, n_g, lam_g
    real(sp), dimension(kp) :: rho_fac
	real(sp), dimension(-o_halo:kp+o_halo) :: vqr, vqs, vqg, vqi, vnr, vns, vng, vni
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! some commonly used variables that depend on prognostics                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t=theta*(p/1e5_sp)**(ra/cp) ! temperature
    rho=p / (ra*t) ! air density    
    rho_fac=(rho0/rho(1:kp))**0.5_sp
    ! rain n0, lambda
    lam_r=(nar*cr*gam2r / (rho(1:kp)*max(q(3,1:kp),1.e-10_sp)))**(1._sp/(1._sp+alpha_r+dr-nbr))
    n_r=nar*lam_r**nbr
    ! ice n0, lambda
    lam_i=(max(q(7,1:kp),1._sp)*ci*gam2i / (max(q(6,1:kp),1.e-10_sp)*gam1i))**(1._sp/di)
    n_i=rho(1:kp)*max(q(7,1:kp),0._sp)*lam_i**(1._sp+alpha_i) / gam1i
    ! snow n0, lambda
    lam_s=(max(q(8,1:kp),1._sp)*cs*gam2s / (q(4,1:kp)*gam1s))**(1._sp/ds)
    n_s=rho(1:kp)*max(q(8,1:kp),0._sp)*lam_s**(1._sp+alpha_s) / gam1s
    ! graupel n0, lambda
    lam_g=(max(q(9,1:kp),1._sp)*cg*gam2g / (max(q(5,1:kp),1.e-10)*gam1g))**(1._sp/dg)
    n_g=rho(1:kp)*max(q(9,1:kp),0._sp)*lam_g**(1._sp+alpha_g) / gam1g
    
    
    ! precipitation
	precip(1,1:kp)=cr*n_r*(a_r*chi_rain/(lam_r**(alpha_r+b_r+dr+1._sp)) - &
					u(1:kp)*chi_rain1/(lam_r**(alpha_r+dr+1._sp))) &
					/rho(1:kp) *3600._sp
	precip(2,1:kp)=cs*n_s*(a_s*chi_snow/(lam_s**(alpha_s+b_s+ds+1._sp)) - &
					u(1:kp)*chi_snow1/(lam_s**(alpha_s+ds+1._sp))) &
					/rho(1:kp)*3600._sp
	precip(3,1:kp)=cg*n_g*(a_g*chi_graupel/(lam_g**(alpha_g+b_g+dg+1._sp)) - &
					u(1:kp)*chi_graupel1/(lam_g**(alpha_g+dg+1._sp))) &
					/rho(1:kp)*3600._sp
	precip(4,1:kp)=ci*n_i*(a_i*chi_ice/(lam_i**(alpha_i+b_i+di+1._sp)) - &
					u(1:kp)*chi_ice1/(lam_i**(alpha_i+di+1._sp))) &
					/rho(1:kp)*3600._sp
    
    ! fall speeds
    vqr(1:kp)=max(fall_q_r*rho_fac * lam_r**(1._sp+alpha_r+dr) / &
    	(lam_r+f_r)**(1._sp+alpha_r+dr+b_r), 0._sp)
    vqs(1:kp)=max(fall_q_s*rho_fac * lam_s**(1._sp+alpha_s+ds) / &
    	(lam_s+f_s)**(1._sp+alpha_s+ds+b_s), 0._sp)
    vqg(1:kp)=max(fall_q_g*rho_fac * lam_g**(1._sp+alpha_g+dg) / &
    	(lam_g+f_g)**(1._sp+alpha_g+dg+b_g), 0._sp)
    vqi(1:kp)=max(fall_q_i*rho_fac * lam_i**(1._sp+alpha_i+di) / &
    	(lam_i+f_i)**(1._sp+alpha_i+di+b_i), 0._sp)
    
    vnr(1:kp)=max(fall_n_r*rho_fac * lam_r**(1._sp+alpha_r) / &
    	(lam_r+f_r)**(1._sp+alpha_r+b_r), 0._sp)
    vns(1:kp)=max(fall_n_s*rho_fac * lam_s**(1._sp+alpha_s) / &
    	(lam_s+f_s)**(1._sp+alpha_s+b_s), 0._sp)
    vng(1:kp)=max(fall_n_g*rho_fac * lam_g**(1._sp+alpha_g) / &
    	(lam_g+f_g)**(1._sp+alpha_g+b_g), 0._sp)
    vni(1:kp)=max(fall_n_i*rho_fac * lam_i**(1._sp+alpha_i) / &
    	(lam_i+f_i)**(1._sp+alpha_i+b_i), 0._sp)
    	
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
        qold=q(2,k)
        qtot=q(1,k)+q(2,k)
        q(2,k)=q(1,k)+q(2,k)-smr(k)
        if (theta_flag) q(2,k)=q(2,k)+(lv/cp*qold)*dqs_dt / (1._sp+lv/cp*dqs_dt)
        q(2,k)=max(q(2,k),0._sp)
        t(k)=t(k)
        if(theta_flag) t(k)=t(k)+lv/cp*(q(2,k)-qold)
		
		tc=t(k)-ttr
    	smr(k)=eps1*svp_liq(t(k))/(p(k)-svp_liq(t(k))) ! saturation mixing ratio
    	q0sat=smr(k)	
    	smr_i(k)=eps1*svp_ice(t(k))/(p(k)-svp_ice(t(k))) ! saturation mixing ratio - ice	
    	
    	cond=(q(2,k)-qold)
    	q(1,k)=q(1,k)-cond
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! inhomogeneous mixing -https://journals.ametsoc.org/doi/pdf/10.1175/2007JAS2374.1
 !        if(q(2,k)<qold) then
!             q(4,k)=q(4,k)*(q(2,k)/qold)**1._sp
!         endif


	
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

	
		rain_evap=(q(1,k)/smr(k)-1._sp) / (rho(k)*ab_liq)*nu_rain
		if(q(1,k).gt.smr(k)) then
			prevp(k)=0._sp
		else
			prevp(k)=-min(rain_evap,0._sp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		


	
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! deposition & sublimation onto ice, snow, graupel                               !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		nu_snow=2._sp*pi*n_s(k) / rho(k) * &
				(nu_s1 / lam_s(k)**(2._sp+alpha_s) + &
				(a_s/nu_vis)**0.5_sp*sc**(1._sp/3._sp)* &
				(rho(k)*rho0)**0.25_sp*nu_s2 / &
				(lam_s(k)+0.5_sp*f_s)**(0.5_sp*b_s+alpha_s+2.5_sp))

		nu_graup=2._sp*pi*n_g(k) / rho(k) * &
				(nu_g1 / lam_g(k)**(2._sp+alpha_g) + &
				(a_g/nu_vis)**0.5_sp*sc**(1._sp/3._sp)* &
				(rho(k)*rho0)**0.25_sp*nu_g2 / &
				(lam_g(k)+0.5_sp*f_g)**(0.5_sp*b_g+alpha_g+2.5_sp))
	
		if (t(k).le.ttr) then
			nu_ice=2._sp*pi*n_i(k) / rho(k) * &
					(nu_i1 / lam_i(k)**(2._sp+alpha_i) + &
					(a_i/nu_vis)**0.5_sp*sc**(1._sp/3._sp)* &
					(rho(k)*rho0)**0.25_sp*nu_i2 / &
					(lam_i(k)+0.5_sp*f_i)**(0.5_sp*b_i+alpha_i+2.5_sp))

			ab_ice=ls**2 / (ktherm1*rv*t(k)**2) + 1._sp/(rho(k)*smr_i(k)*diff1)
		
			ice_dep=(q(1,k)/smr_i(k)-1._sp) / (rho(k)*ab_ice)*nu_ice
			snow_dep=(q(1,k)/smr_i(k)-1._sp) / (rho(k)*ab_ice)*nu_snow
			graup_dep=(q(1,k)/smr_i(k)-1._sp) / (rho(k)*ab_ice)*nu_graup
			if(q(1,k).gt.smr_i(k)) then
				pisub(k)=0._sp
				pssub(k)=0._sp
				pgsub(k)=0._sp
				pidep(k)=max(ice_dep,0._sp)
				psdep(k)=max(snow_dep,0._sp)
				pgdep(k)=max(graup_dep,0._sp)
			else
				pidep(k)=0._sp
				psdep(k)=0._sp
				pgdep(k)=0._sp			
				pisub(k)=-min(ice_dep,0._sp)
				pssub(k)=-min(snow_dep,0._sp)
				pgsub(k)=-min(graup_dep,0._sp)
			endif
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! warm rain autoconversion                                                       !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		praut(k)=aw0*max(q(2,k)-lw0/rho(k), 0._sp) ! kessler scheme
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! collection of cloud by rain, cloud by snow and cloud by ice                    !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		pracw(k)=max(phi_r* n_r(k)* erw *q(2,k)*rho_fac(k) / &
				(lam_r(k)+f_r)**(3._sp+b_r+alpha_r),0._sp)
		! snow-cloud water
		psacw(k)=max(mass_sacw_i * n_s(k)* esw *q(2,k)*rho_fac(k) / &
				(lam_s(k)+f_s)**(3._sp+b_s+alpha_s),0._sp)
		psaci(k)=max(mass_sacw_i * n_s(k)* esi(k) *q(6,k)*rho_fac(k) / &
				(lam_s(k)+f_s)**(3._sp+b_s+alpha_s),0._sp)
		if (t(k).le.ttr) then
			piacw(k)=max(mass_iacw * n_i(k)* eiw *q(2,k)*rho_fac(k) / &
					(lam_i(k)+f_i)**(3._sp+b_i+alpha_i),0._sp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! collection of rain by ice to make snow or graupel                              !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (t(k).le.ttr) then
			if(q(3,k).lt.1e-4_sp) then
				praci_g(k)=0._sp
				praci_s(k)=max(phi_r* n_r(k)* eri *q(6,k)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_sp+b_r+alpha_r),0._sp)
				piacr_g(k)=0._sp
				piacr_s(k)=max(mass_iacr*n_r(k)*q(7,k)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_sp+b_r+dr+alpha_r),0._sp)
				riacr_g(k)=0._sp
				riacr_s(k)=max(num_iacr*n_r(k)*q(7,k)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_sp+b_r+alpha_r),0._sp)
			else
				praci_s(k)=0._sp
				praci_g(k)=max(phi_r* n_r(k)* eri *q(6,k)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_sp+b_r+alpha_r),0._sp)
				piacr_s(k)=0._sp
				piacr_g(k)=max(mass_iacr*n_r(k)*q(7,k)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_sp+b_r+dr+alpha_r),0._sp)
				riacr_s(k)=0._sp
				riacr_g(k)=max(num_iacr*n_r(k)*q(7,k)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_sp+b_r+alpha_r),0._sp)
		
			endif
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! snow autoconversion                                                            !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (q(6,k).gt.1e-10_sp.and.(q(7,k).gt.1e-5_sp)) then
			psaut(k)=max(q(6,k)* &
				(min(lambda_imin /lam_i(k),2._sp)**di - 1._sp) / tsaut,0._sp)
			rsaut(k)=max(psaut(k)/(ci*di2s**di),0._sp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! graupel autoconversion                                                         !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (rho(k)*q(4,k).gt.3e-4_sp.and.(t(k).lt.269.15_sp)) then
			pgaut(k)=max(0.5_sp*max(0._sp,psacw(k)-psdep(k)-psaci(k)), 0._sp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! snow and ice aggregation                                                       !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		rsacs(k)=max(isnow*n_s(k)**2._sp*ess(k) *rho_fac(k) / &
 				lam_s(k)**(4._sp+2.*sp*alpha_s+b_s),0._sp)
		riaci(k)=max(iice*n_i(k)**2._sp*eii(k) *rho_fac(k) / &
				lam_i(k)**(4._sp+2.*sp*alpha_i+b_i),0._sp)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! snow break-up                                                                  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		rsbrk(k)=max((lambda_s_break/lam_s(k)-1._sp)**ds*q(8,k)/(tsbreak),0._sp)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! immersion freezing of rain                                                     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		if(t(k).le.268.15_sp) then
! 			pgfr(k)=max( mass_imm*(exp(-abigg*(t(k)-ttr))-1._sp)*n_r(k) / &
! 						(lam_r(k))**(4_sp+dr+alpha_r),0._sp)
! 			rgfr(k)=num_imm*(exp(-abigg*(t(k)-ttr))-1._sp)*n_r(k)/ &
! 						(lam_r(k))**(4_sp+alpha_r)
! 		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! collisions between precipitating particles of different species                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(t(k).lt.268._sp) then
			pracs(k)=max(n_r(k)*n_s(k)*pi/(4._sp*rho(k))*ers*cs*max((vqs(k)+vqr(k))/8._sp,abs(vqs(k)-vqr(k))) * &
					( &
					mass_racs1/(lam_r(k)**(1._sp+alpha_r) *lam_s(k)**(3._sp+alpha_s+ds)) + &
					mass_racs2/(lam_r(k)**(2._sp+alpha_r) *lam_s(k)**(2._sp+alpha_s+ds)) + &
					mass_racs3/(lam_r(k)**(3._sp+alpha_r) *lam_s(k)**(1._sp+alpha_s+ds))  &
					) , 0._sp)  
				
	!         pracg(k)=max(n_r(k)*n_g(k)*pi/(4._sp*rho(k))*erg*cg*max((vqg(k)+vqr(k))/8._sp,abs(vqg(k)-vqr(k))) * &
	!         		( &
	!         		mass_racg1/(lam_r(k)**(1._sp+alpha_r) *lam_g(k)**(3._sp+alpha_g+dg)) + &
	!         		mass_racg2/(lam_r(k)**(2._sp+alpha_r) *lam_g(k)**(2._sp+alpha_g+dg)) + &
	!         		mass_racg3/(lam_r(k)**(3._sp+alpha_r) *lam_g(k)**(1._sp+alpha_g+dg))  &
	!         		)   , 0._sp)        		
			psacr(k)=max(n_s(k)*n_r(k)*pi/(4._sp*rho(k))*ers*cr*max((vqr(k)+vqs(k))/8._sp,abs(vqr(k)-vqs(k))) * &
					( &
					mass_sacr1/(lam_s(k)**(1._sp+alpha_s) *lam_r(k)**(3._sp+alpha_r+dr)) + &
					mass_sacr2/(lam_s(k)**(2._sp+alpha_s) *lam_r(k)**(2._sp+alpha_r+dr)) + &
					mass_sacr3/(lam_s(k)**(3._sp+alpha_s) *lam_r(k)**(1._sp+alpha_r+dr))  &
					)    , 0._sp)       		
			rsacr(k)=max(n_s(k)*n_r(k)*pi/(4._sp*rho(k))*ers*max((vqr(k)+vqs(k))/8._sp,abs(vqr(k)-vqs(k))) * &
					( &
					num_racs1/(lam_s(k)**(1._sp+alpha_s) *lam_r(k)**(3._sp+alpha_r)) + &
					num_racs2/(lam_s(k)**(2._sp+alpha_s) *lam_r(k)**(2._sp+alpha_r)) + &
					num_racs3/(lam_s(k)**(3._sp+alpha_s) *lam_r(k)**(1._sp+alpha_r))  &
					)    , 0._sp)    
		
	!         psacg(k)=max(n_s(k)*n_g(k)*pi/(4._sp*rho(k))*egs_dry(k)*cg*max((vqg(k)+vqs(k))/8._sp,abs(vqg(k)-vqs(k))) * &
	!         		( &
	!         		mass_sacg1/(lam_s(k)**(1._sp+alpha_s) *lam_g(k)**(3._sp+alpha_g+dg)) + &
	!         		mass_sacg2/(lam_s(k)**(2._sp+alpha_s) *lam_g(k)**(2._sp+alpha_g+dg)) + &
	!         		mass_sacg3/(lam_s(k)**(3._sp+alpha_s) *lam_g(k)**(1._sp+alpha_g+dg))  &
	!         		)    , 0._sp)       		
			pgacr(k)=max(n_g(k)*n_r(k)*pi/(4._sp*rho(k))*erg*cr*max((vqr(k)+vqg(k))/8._sp,abs(vqr(k)-vqg(k))) * &
					( &
					mass_gacr1/(lam_g(k)**(1._sp+alpha_g) *lam_r(k)**(3._sp+alpha_r+dr)) + &
					mass_gacr2/(lam_g(k)**(2._sp+alpha_g) *lam_r(k)**(2._sp+alpha_r+dr)) + &
					mass_gacr3/(lam_g(k)**(3._sp+alpha_g) *lam_r(k)**(1._sp+alpha_r+dr))  &
					)    , 0._sp)      
			! set pgacs and rgacs to wet first 		
			pgacs(k)=max(n_g(k)*n_s(k)*pi/(4._sp*rho(k))*egs_wet*cs*max((vqs(k)+vqg(k))/8._sp,abs(vqs(k)-vqg(k))) * &
					( &
					mass_gacs1/(lam_g(k)**(1._sp+alpha_g) *lam_s(k)**(3._sp+alpha_s+ds)) + &
					mass_gacs2/(lam_g(k)**(2._sp+alpha_g) *lam_s(k)**(2._sp+alpha_s+ds)) + &
					mass_gacs3/(lam_g(k)**(3._sp+alpha_g) *lam_s(k)**(1._sp+alpha_s+ds))  &
					)    , 0._sp)       		
				
			rgacs(k)=max(n_g(k)*n_s(k)*pi/(4._sp*rho(k))*egs_wet*max((vqs(k)+vqg(k))/8._sp,abs(vqs(k)-vqg(k))) * &
					( &
					num_sacg1/(lam_g(k)**(1._sp+alpha_g) *lam_s(k)**(3._sp+alpha_s)) + &
					num_sacg2/(lam_g(k)**(2._sp+alpha_g) *lam_s(k)**(2._sp+alpha_s)) + &
					num_sacg3/(lam_g(k)**(3._sp+alpha_g) *lam_s(k)**(1._sp+alpha_s))  &
					)  , 0._sp)   
		endif
        		
        !rracg?, psacg? pracg? rgacr?   		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! riming of graupel and accretion of ice                                         !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(q(5,k).gt.0._sp) then
			pgacw(k)=max(mass_gacw*rho_fac(k)*n_g(k) / &
					((lam_g(k)+f_g)**(3._sp+b_g+alpha_g))*q(2,k), 0._sp)
		endif
		
		if((q(5,k).gt.0._sp).and.(t(k).lt.ttr)) then
		
			! below assuming wet-growth first
			pgaci(k)=max(mass_gaci*egi_wet*rho_fac(k)*n_g(k) / &
					((lam_g(k)+f_g)**(3._sp+b_g+alpha_g))*q(6,k), 0._sp)
					
			pgwet=(910._sp/(cg*6._sp/pi))**0.625_sp*  &
				(rho(k)*lv*(q0sat-q(1,k))-ktherm1*tc) / (rho(k)*(lf-cw*tc))*nu_graup + &
				(pgaci(k)+pgacs(k))*(1._sp-ci*tc/(lf+cw*tc))
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! redefine as dry growth, if water can be frozen without shedding            !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if((pgacr(k)+pgacw(k)+pgaci(k)+pgacs(k)).gt.pgwet) then
				! redefine graupel accreting rain so that it is
				pgshd(k)=max(0._sp,pgacw(k)+pgaci(k)+pgacs(k)-pgwet)
				pgacr(k)=max(0._sp,pgwet-pgacw(k)-pgaci(k)-pgacs(k))
			else
				! dry growth: graupel-ice
				pgaci(k)=max(mass_gaci*egi_dry(k)*rho_fac(k)*n_g(k) / &
						((lam_g(k)+f_g)**(3._sp+b_g+alpha_g))*q(6,k), 0._sp)
									
				! dry growth: graupel-snow -- mass
				pgacs(k)=max(n_g(k)*n_s(k)*pi/(4._sp*rho(k))*egs_dry(k)* &
						cs*max((vqs(k)+vqg(k))/8._sp,abs(vqs(k)-vqg(k))) * &
						( &
				 mass_gacs1/(lam_g(k)**(1._sp+alpha_g) *lam_s(k)**(3._sp+alpha_s+ds)) + &
				 mass_gacs2/(lam_g(k)**(2._sp+alpha_g) *lam_s(k)**(2._sp+alpha_s+ds)) + &
				 mass_gacs3/(lam_g(k)**(3._sp+alpha_g) *lam_s(k)**(1._sp+alpha_s+ds))  &
						)    , 0._sp)       		
				
				! dry growth: graupel-snow -- number
				rgacs(k)=max(n_g(k)*n_s(k)*pi/(4._sp*rho(k))*egs_dry(k)* &
						max((vqs(k)+vqg(k))/8._sp,abs(vqs(k)-vqg(k))) * &
						( &
					num_sacg1/(lam_g(k)**(1._sp+alpha_g) *lam_s(k)**(3._sp+alpha_s)) + &
					num_sacg2/(lam_g(k)**(2._sp+alpha_g) *lam_s(k)**(2._sp+alpha_s)) + &
					num_sacg3/(lam_g(k)**(3._sp+alpha_g) *lam_s(k)**(1._sp+alpha_s))  &
						)  , 0._sp)   
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! h-m process                                                                    !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(hm_flag) &
			pihal(k)=max(hm_rate*mi0*(pgacw(k)+psacw(k))*hm_func(t(k)),0._sp)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! melting of ice, snow, and graupel                                              !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(t(k).gt.ttr) then
			pimlt(k)=q(6,k)/dt ! ice melts instantaneously
			q(7,k)=0._sp
			if(q(5,k).gt.0._sp) then
				pgmlt(k)=max(1._sp/(rho(k)*lf) * &
					(ktherm1*tc-lv*diff1*rho(k)*(q(1,k)-q0sat))*nu_graup &
					+ cw*tc/lf*(pgacw(k)+pgacr(k)-pgshd(k)),0._sp)
			endif
				
			psmlt(k)=max(1._sp/(rho(k)*lf) * &
				(ktherm1*tc-lv*diff1*rho(k)*(q(1,k)-q0sat))*nu_snow &
				+ cw*tc/lf*(psacw(k)+psacr(k)),0._sp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    enddo
    
    
    
    
    
    
    ! update variables
    ! ice number
    
    
    q(7,1:kp)=q(7,1:kp)+((piprm+pihal+picnt+pifrw)/mi0 - riaci - &
    		(pgaci+psaci+praci_g+praci_s)*q(7,1:kp)/(q(6,1:kp)+qsmall)  &
    		-riaci - rsaut)*dt
    where(pimlt.gt.0._sp)
    	q(7,1:kp)=0._sp
    end where
    ! snow number
    q(8,1:kp)=q(8,1:kp)+(riacr_s+rsaut+rsbrk &
    			-(pssub+psmlt+pgaut)*q(8,1:kp)/(q(4,1:kp)+qsmall) &
    			-(rgacs+rsacr+rsacs))*dt
    ! graupel number
    q(9,1:kp)=q(9,1:kp)+(pgaut*q(9,1:kp)/(q(5,1:kp)+qsmall) + rsacr+riacr_g+rgfr &
    			-(pgsub+pgmlt)*q(9,1:kp)/(q(5,1:kp)+qsmall))*dt
    			

    ! vapour mass
    q(1,1:kp)=q(1,1:kp)+(pgsub+pssub+prevp+pisub-(psdep+pidep+piprm+pgdep))*dt
    ! liquid mass
    q(2,1:kp)=q(2,1:kp)-((pgacw+praut+psacw+pracw+piacw+pihal+picnt+pifrw))*dt
    ! rain mass
    q(3,1:kp)=q(3,1:kp)+(pgmlt+praut+pgshd+pracw+psmlt+pimlt- &
    			(pgacr+prevp+pgfr+psacr+piacr_g+piacr_s))*dt
    ! snow mass
    q(4,1:kp)=q(4,1:kp)+(psaut+psdep+psaci+praci_s+piacr_s+psacw-(pssub+pgacs+pracs+pgaut+psmlt))*dt
    ! graupel mass
    q(5,1:kp)=q(5,1:kp)+(pgaci+pgacw+pgacs+pgacr+psacr+pracs+pgaut+pgfr+praci_g+piacr_g+pgdep- &
    				(pgsub+pgmlt+pgshd))*dt
!     if(sum(piacw).gt.0._sp) stop
    ! ice mass
    q(6,1:kp)=q(6,1:kp)+(pidep+piprm+pihal+picnt+piacw- &
    			(psaut+pgaci+psaci+pisub+pifrw+pimlt+praci_g+praci_s))*dt
    			
    q=max(q,0._sp)	    
    if (theta_flag) theta=t*(1.e5_sp/p)**(ra/cp)
    
 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! advection rain 0th order Bott, a.k.a. upstream advection                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! rain 
    if(sum(q(3,1:kp)).gt.qsmall) then
		where(isnan(vqr))
			vqr=0_sp
		end where
		vqr(-o_halo:0)=vqr(1)
		vqr(kp+1:kp+o_halo)=vqr(kp)
		n_step=max(ceiling(maxval(vqr)*dt/dz*2_sp),1)
		vqr(-o_halo:kp+o_halo-1)=-vqr(-o_halo+1:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vqr,q(3,:),.false.)
		enddo
	endif
	
    ! graupel 
    if(sum(q(5,1:kp)).gt.qsmall) then
		where(isnan(vqg))
			vqg=0_sp
		end where
		vqg(-o_halo:0)=vqg(1)
		vqg(kp+1:kp+o_halo)=vqg(kp)
		n_step=max(ceiling(maxval(vqg)*dt/dz*2_sp),1)
		vqg(-o_halo:kp+o_halo-1)=-vqg(-o_halo+1:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vqg,q(5,:),.false.)
		enddo
		where(isnan(vng))
			vng=0_sp
		end where
		vng(-o_halo:0)=vng(1)
		vng(kp+1:kp+o_halo)=vng(kp)
		n_step=max(ceiling(maxval(vng)*dt/dz*2_sp),1)
		vng(-o_halo:kp+o_halo-1)=-vng(-o_halo+1:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vng,q(9,:),.false.)
		enddo
	endif
	
    ! snow 
    if(sum(q(4,1:kp)).gt.qsmall) then
		where(isnan(vqs))
			vqs=0_sp
		end where
		vqs(-o_halo:0)=vqs(1)
		vqs(kp+1:kp+o_halo)=vqs(kp)
		n_step=max(ceiling(maxval(vqs)*dt/dz*2_sp),1)
		vqs(-o_halo:kp+o_halo-1)=-vqs(-o_halo+1:kp+o_halo)
		do iter=1,n_step		
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vqs,q(4,:),.false.)
		enddo
		where(isnan(vns))
			vns=0_sp
		end where
		vns(-o_halo:0)=vns(1)
		vns(kp+1:kp+o_halo)=vns(kp)
		n_step=max(ceiling(maxval(vns)*dt/dz*2_sp),1)
		vns(-o_halo:kp+o_halo-1)=-vns(-o_halo+1:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vns,q(8,:),.false.)
		enddo
	endif
    ! ice
    if(sum(q(6,1:kp)).gt.qsmall) then     
		where(isnan(vqi))
			vqi=0_sp
		end where
		vqi(-o_halo:0)=vqi(1)
		vqi(kp+1:kp+o_halo)=vqi(kp)
		n_step=max(ceiling(maxval(vqi)*dt/dz*2_sp),1)
		vqi(-o_halo:kp+o_halo-1)=-vqi(-o_halo+1:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vqg,q(6,:),.false.)
		enddo
		where(isnan(vni))
			vni=0_sp
		end where
		vni(-o_halo:0)=vni(1)
		vni(kp+1:kp+o_halo)=vni(kp)
		n_step=max(ceiling(maxval(vni)*dt/dz*2_sp),1)
		vni(-o_halo:kp+o_halo-1)=-vni(-o_halo+1:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,sp),dz,z,vng,q(7,:),.false.)
		enddo
	endif
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
    
    

    end subroutine microphysics_1d
    
    
    
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


    end module micro_module
    
