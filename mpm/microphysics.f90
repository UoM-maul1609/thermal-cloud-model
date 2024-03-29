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
    use numerics_type
    
    private
    public :: microphysics_2d, microphysics_1d,set_qnames
    
    ! physical constants
    real(wp), parameter :: rhow=1000._wp, rhoi=920._wp,lv=2.5e6_wp,ls=2.8e6_wp,lf=ls-lv, &
    					   cp=1005._wp, cw=4187._wp, cice=2093._wp, r=8.314_wp, &
    						mw=18e-3_wp, ma=29e-3_wp, ra=r/ma,rv=r/mw, eps1=ra/rv, &
    						ttr=273.15_wp, joules_in_an_erg=1.0e-7_wp, &
    						joules_in_a_cal=4.187_wp
    						
    						
    ! mass-diameter and size spectra relations
    real(wp), parameter :: cr=523.6_wp, cs=52.36_wp, cg=261.8_wp, ci=104._wp, &
    					dr=3_wp, ds=3._wp, dg=3._wp, di=3._wp, &
    					alpha_r=2.5_wp, alpha_s=2.5_wp, alpha_g=2.5_wp, alpha_i=0._wp
    					
	! terminal fall-speed relations
	real(wp), parameter :: a_r=362._wp, a_s=4.84_wp, a_g=253._wp, a_i=71.34_wp, &
							b_r=0.65_wp, b_s=0.25_wp, b_g=0.734_wp, b_i=0.6635_wp, &
							f_r=0._wp, f_s=0._wp, f_g=0._wp, f_i=0._wp
							
	! autoconversion
	real(wp), parameter :: aw0=1e-3_wp, dwa=20e-6_wp, nl=2.4e8_wp, &
						lw0=rhow*pi/6._wp*nl*dwa**3, &
						tsaut=60._wp, dimax=0.3e-3_wp, di2s=0.33e-3_wp, &
						lambda_imin=(1._wp+di+alpha_i)/dimax, &
						tsbreak=60._wp, lambda_s_break=1000._wp
	
    ! microphysical values:
    real(wp), parameter :: hm_rate=3.5e8_wp, nar=1.1e15_wp, nbr=0._wp, &
    						rho0=1.2_wp, bbigg=100._wp, abigg=0.66_wp
    real(wp) :: mi0=1.e-14_wp

	! coalescence efficiencies
	real(wp), parameter :: erw=1._wp, erg=1._wp, ers=1._wp, eri=1._wp, esw=1._wp, &
						egw=1._wp, eiw=1._wp, egs_wet=1._wp, egi_wet=1._wp
	
	! variables used in various process rates:
	real(wp) :: gam1r,gam2r,gam1i,gam2i, gam1s, gam2s,gam1g,gam2g, &
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
				
	real(wp), dimension(3) :: c=[1._wp,2._wp,1._wp]
	integer(i4b) :: k
	real(wp) :: isnow, iice, f1,f2,a,b, qsmall=1e-30_wp
	
	
    contains
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>read in the data from the namelists for the BAM module
	!> and set variables for microphysics
	!>@param[inout] q_name, q_type, c_s, c_e
	!>@param[inout] nq,ncat, nprec, iqv, iqc, ini, iqi
	subroutine set_qnames(q_name,q_type,c_s,c_e,nq,ncat,nprec, &
	            iqv, iqc, ini, iqi)
		implicit none
        integer(i4b), intent(inout) :: nq, ncat, nprec, iqv, iqc, ini, iqi
        integer(i4b), intent(inout), dimension(:), allocatable :: q_type, c_s, c_e
        character(len=20), dimension(:), allocatable :: q_name
        
        integer(i4b) :: i
        
        
        ncat=9
        
        
        nq=9 ! vapour, cloud mass, rain mass, cloud water, rain water
        
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
        q_type(4)=1 ! snow water
        c_s(4)=4
        c_e(4)=4
        q_type(5)=2 ! graupel water
        c_s(5)=5
        c_e(5)=5
        q_type(6)=2 ! ice water
        c_s(6)=6
        c_e(6)=6
        q_type(7)=2 ! ice water number
        c_s(7)=7
        c_e(7)=7
        q_type(8)=2 ! snow water number
        c_s(8)=8
        c_e(8)=8
        q_type(9)=2 ! graupel water number
        c_s(9)=9
        c_e(9)=9
        
        q_name=["qv","qc","qr","qs","qg","qi","ni","ns","ng"]
        
        nprec=4
        
        iqv=1
        iqc=2
        ini=7
        iqi=6
        
	end subroutine set_qnames


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>initialises variables for use with the microphysics
    subroutine initialise_microphysics_vars
    use hypergeo, only : hygfx
    implicit none

	! used to calculate intercept and slopes
	gam1r=gamma(1._wp+alpha_r)
	gam2r=gamma(1._wp+alpha_r+dr)
	gam1i=gamma(1._wp+alpha_i)
	gam2i=gamma(1._wp+alpha_i+di)
	gam1s=gamma(1._wp+alpha_s)
	gam2s=gamma(1._wp+alpha_s+ds)
	gam1g=gamma(1._wp+alpha_g)
	gam2g=gamma(1._wp+alpha_g+dg)

    ! mass weighted fall for r, s, g, i
    fall_q_r=a_r*gamma(1._wp+alpha_r+dr+b_r) / gamma(1._wp+alpha_r+dr)
    fall_q_s=a_s*gamma(1._wp+alpha_s+ds+b_s) / gamma(1._wp+alpha_s+ds)
    fall_q_g=a_g*gamma(1._wp+alpha_g+dg+b_g) / gamma(1._wp+alpha_g+dg)
    fall_q_i=a_i*gamma(1._wp+alpha_i+di+b_i) / gamma(1._wp+alpha_i+di)

    ! number weighted fall for r, s, g
    fall_n_r=a_r*gamma(1._wp+alpha_r+b_r) / gamma(1._wp+alpha_r)
    fall_n_s=a_s*gamma(1._wp+alpha_s+b_s) / gamma(1._wp+alpha_s)
    fall_n_g=a_g*gamma(1._wp+alpha_g+b_g) / gamma(1._wp+alpha_g)
    fall_n_i=a_i*gamma(1._wp+alpha_i+b_i) / gamma(1._wp+alpha_i)
    
    ! sweep out of rain
    phi_r=pi*a_r*gamma(3._wp+b_r+alpha_r) / 4._wp
    
    ! ice accreting rain
    mass_iacr=pi*eri*a_r*cr*gamma(3._wp+b_r+dr+alpha_r)/4._wp
    num_iacr =pi*eri*a_r*gamma(3._wp+b_r+alpha_r)/4._wp

	! collection of cloud by snow and ice
	mass_sacw_i=pi*a_s*gamma(3._wp+b_s+alpha_s)/4._wp
	mass_iacw=pi*a_i*gamma(3._wp+b_i+alpha_i)/4._wp
	
	! collisions between precipitating particles of different species
	! rain-snow
	mass_racs1=gamma(1._wp+alpha_r)*gamma(3._wp+alpha_s+ds)
	mass_racs2=2._wp*gamma(2._wp+alpha_r)*gamma(2._wp+alpha_s+ds)
	mass_racs3=gamma(3._wp+alpha_r)*gamma(1._wp+alpha_s+ds)
	
	num_racs1=gamma(1._wp+alpha_r)*gamma(3._wp+alpha_s)
	num_racs2=2._wp*gamma(2._wp+alpha_r)*gamma(2._wp+alpha_s)
	num_racs3=gamma(3._wp+alpha_r)*gamma(1._wp+alpha_s)
    ! rain-graupel
	mass_racg1=gamma(1._wp+alpha_r)*gamma(3._wp+alpha_g+dg)
	mass_racg2=2._wp*gamma(2._wp+alpha_r)*gamma(2._wp+alpha_g+dg)
	mass_racg3=gamma(3._wp+alpha_r)*gamma(1._wp+alpha_g+dg)

	num_racg1=gamma(1._wp+alpha_r)*gamma(3._wp+alpha_g)
	num_racg2=2._wp*gamma(2._wp+alpha_r)*gamma(2._wp+alpha_g)
	num_racg3=gamma(3._wp+alpha_r)*gamma(1._wp+alpha_g)
    ! snow-rain
	mass_sacr1=gamma(1._wp+alpha_s)*gamma(3._wp+alpha_r+dr)
	mass_sacr2=2._wp*gamma(2._wp+alpha_s)*gamma(2._wp+alpha_r+dr)
	mass_sacr3=gamma(3._wp+alpha_s)*gamma(1._wp+alpha_r+dr)
    ! snow-graupel
	mass_sacg1=gamma(1._wp+alpha_s)*gamma(3._wp+alpha_g+dg)
	mass_sacg2=2._wp*gamma(2._wp+alpha_s)*gamma(2._wp+alpha_g+dg)
	mass_sacg3=gamma(3._wp+alpha_s)*gamma(1._wp+alpha_g+dg)

	num_sacg1=gamma(1._wp+alpha_s)*gamma(3._wp+alpha_g)
	num_sacg2=2._wp*gamma(2._wp+alpha_s)*gamma(2._wp+alpha_g)
	num_sacg3=gamma(3._wp+alpha_s)*gamma(1._wp+alpha_g)
    ! graupel-rain
	mass_gacr1=gamma(1._wp+alpha_g)*gamma(3._wp+alpha_r+dr)
	mass_gacr2=2._wp*gamma(2._wp+alpha_g)*gamma(2._wp+alpha_r+dr)
	mass_gacr3=gamma(3._wp+alpha_g)*gamma(1._wp+alpha_r+dr)
    ! graupel-snow
	mass_gacs1=gamma(1._wp+alpha_g)*gamma(3._wp+alpha_s+ds)
	mass_gacs2=2._wp*gamma(2._wp+alpha_g)*gamma(2._wp+alpha_s+ds)
	mass_gacs3=gamma(3._wp+alpha_g)*gamma(1._wp+alpha_s+ds)
	
	
	! accretion and riming by graupel
	mass_gacw=pi*egw*a_g*gamma(3._wp+b_g+alpha_g)/4._wp
	mass_gaci=pi*a_g*gamma(3._wp+b_g+alpha_g)/4._wp
    
    ! gauss hypergeometric equations aggregation of ice with ice (and snow with snow)
    ! See Ferrier (1994, JAS part 1, equation B.21)
    ! snow:
    a=1._wp
    b=4._wp+2._wp*alpha_s+b_s
    isnow=0._wp
    do k=1,3
	    call hygfx(a, b, real(k,wp)+alpha_s+1.0_wp,0.5_wp,f1)
	    call hygfx(a, b, real(k,wp)+alpha_s+b_s+1.0_wp, 0.5_wp,f2)
	    isnow=isnow+c(k)*(f1/(real(k,wp)+alpha_s)-f2/(real(k,wp)+alpha_s+b_s))
	enddo
	isnow=a_s*pi*gamma(b)/(2._wp**(6._wp+2._wp*alpha_s+b_s)) * isnow
	
    ! ice:
    a=1._wp
    b=4._wp+2._wp*alpha_i+b_i
    iice=0._wp
    do k=1,3
	    call hygfx(a, b, real(k,wp)+alpha_i+1.0_wp,0.5_wp,f1)
	    call hygfx(a, b, real(k,wp)+alpha_i+b_i+1.0_wp, 0.5_wp,f2)
	    iice=iice+c(k)*(f1/(real(k,wp)+alpha_i)-f2/(real(k,wp)+alpha_i+b_i))
	enddo
	iice=a_i*pi*gamma(b)/(2._wp**(6._wp+2._wp*alpha_i+b_i)) * iice
	
	
	! ventilation
	! rain:
	nu_r1=0.78_wp*gamma(2._wp+alpha_r)
	nu_r2=gamma(0.5_wp*b_r+alpha_r+2.5)
	! ice:
	nu_i1=0.78_wp*gamma(2._wp+alpha_i)
	nu_i2=0.31_wp*gamma(0.5_wp*b_i+alpha_i+2.5)
	! snow:
	nu_s1=0.78_wp*gamma(2._wp+alpha_s)
	nu_s2=0.31_wp*gamma(0.5_wp*b_s+alpha_s+2.5)
	! graupel:
	nu_g1=0.78_wp*gamma(2._wp+alpha_g)
	nu_g2=0.31_wp*gamma(0.5_wp*b_g+alpha_g+2.5)
	
	! immersion freezing by bigg
	mass_imm=gamma(4._wp+dr+alpha_r)*pi*cr*bbigg/6._wp
	num_imm=gamma(4._wp+alpha_r)*pi*bbigg/6._wp
	
	
	! precipitation
	chi_rain=gamma(1._wp+alpha_r+b_r+dr)
	chi_ice=gamma(1._wp+alpha_i+b_i+di)
	chi_snow=gamma(1._wp+alpha_s+b_s+ds)
	chi_graupel=gamma(1._wp+alpha_g+b_g+dg)
	
	chi_rain1=gamma(1._wp+alpha_r+dr)
	chi_ice1=gamma(1._wp+alpha_i+di)
	chi_snow1=gamma(1._wp+alpha_s+ds)
	chi_graupel1=gamma(1._wp+alpha_g+dg)
	
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
	!>@param[in] theta_ref: reference potential temperature 
	!>@param[inout] rho: density 
	!>@param[in] w: vertical wind 
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag: switch hm-process on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] theta_flag: whether to alter theta
    subroutine microphysics_2d(nq,ip,kp,o_halo,dt,dz,q,precip,theta,p, z,theta_ref,rho,w, &
    						micro_init,hm_flag, mass_ice,theta_flag)
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, ip,kp, o_halo
    real(wp), intent(in) :: dt,dz
    real(wp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo,nq), intent(inout) :: q
    real(wp), dimension(1:kp,1:ip,4), intent(inout) :: precip
    real(wp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: &
    					theta, p, rho
    real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z, theta_ref
    real(wp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(in) :: w
    logical, intent(in) :: hm_flag, theta_flag
    logical , intent(inout) :: micro_init
    real(wp), intent(in) :: mass_ice

	! locals
	integer(i4b) :: i
	
	do i=1,ip
		call microphysics_1d(nq,kp,o_halo,dt,dz,q(:,i,:),precip(:,i,:),theta(:,i),p(:,i), &
							z(:),theta_ref,rho(:,i),w(:,i), &
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
	!>@param[inout] th: theta perturbation
	!>@param[inout] p: pressure
	!>@param[in] z: vertical levels 
	!>@param[inout] theta: potential temperature 
	!>@param[inout] rho: density 
	!>@param[in] u: vertical wind 
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag: switch hm-process on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] theta_flag: whether to alter theta
    subroutine microphysics_1d(nq,kp,o_halo,dt,dz,q,precip,th,p, z,theta,rho,u, &
    						micro_init,hm_flag, mass_ice, theta_flag)
	use advection_1d
	use numerics, only : dfsid1
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, kp, o_halo
    real(wp), intent(in) :: dt,dz
    real(wp), dimension(-o_halo+1:kp+o_halo,nq), intent(inout) :: q
    real(wp), dimension(1:kp,4), intent(inout) :: precip
    real(wp), dimension(-o_halo+1:kp+o_halo), intent(inout) :: p, th, rho
    real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: u, z, theta
    logical, intent(in) :: hm_flag, theta_flag
    logical , intent(inout) :: micro_init
    real(wp), intent(in) :: mass_ice
    ! locals:
    integer(i4b) :: k,iter, n_step
    real(wp) :: temp, qtot,qaut, a, b, ab_ice, ab_liq, ice_dep,snow_dep,graup_dep, &
    			nu_ice, nu_snow, nu_graup, diff1, ktherm1, tc, nu_vis, sc, nu_rain, rain_evap
    real(wp), dimension(kp) :: smr, smr_i
    
    real(wp), dimension(kp) :: &
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
    				    
    real(wp) :: pgwet ! amount of liquid that graupel can freeze without shedding
    								

    real(wp), dimension(kp) :: n_r, lam_r, n_i, lam_i, n_s, lam_s, n_g, lam_g
    real(wp), dimension(kp) :: rho_fac
	real(wp), dimension(1-o_halo:kp+o_halo) :: vqr, vqs, vqg, vqi, vnr, vns, vng, vni
	real(wp), dimension(1-o_halo:kp+o_halo) :: t
	! coalescence efficiencies
	real(wp), dimension(kp) :: egi_dry, egs_dry, esi, eii, ess
	real(wp) :: qold,des_dt,dqs_dt,err,cond,temp1
	
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
	pgaci=0._wp
	pgacr=0._wp
	pgacs=0._wp
	rgacs=0._wp
	pgacw=0._wp
	pgaut=0._wp
	pgdep=0._wp
	pgfr=0._wp
	rgfr=0._wp
	pgmlt=0._wp
	pgshd=0._wp
	pgsub=0._wp
	riaci=0._wp
	piacr_g=0._wp
	riacr_g=0._wp
	piacr_s=0._wp
	riacr_s=0._wp
	piacw=0._wp
	picnt=0._wp
	pidep=0._wp
	piprm=0._wp
	pifrw=0._wp
	pihal=0._wp
	pimlt=0._wp
	pisub=0._wp
	praci_g=0._wp
	praci_s=0._wp
	pracs=0._wp
	pracw=0._wp
	praut=0._wp
	prevp=0._wp
	psacr=0._wp
	rsacr=0._wp
	psaci=0._wp
	rsacs=0._wp
	psacw=0._wp
	psaut=0._wp
	rsaut=0._wp
	rsbrk=0._wp
	psdep=0._wp
	psmlt=0._wp
	pssub=0._wp

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! some commonly used variables that depend on prognostics                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t=(theta+th)*(p/1e5_wp)**(ra/cp) ! temperature
    rho=p / (ra*t) ! air density    
    rho_fac=(rho0/rho(1:kp))**0.5_wp
    ! rain n0, lambda
    lam_r=(nar*cr*gam2r / (rho(1:kp)*max(q(1:kp,3),1.e-10_wp)))**(1._wp/(1._wp+alpha_r+dr-nbr))
    n_r=nar*lam_r**nbr
    ! ice n0, lambda
    lam_i=(max(q(1:kp,7),1._wp)*ci*gam2i / (max(q(1:kp,6),1.e-10_wp)*gam1i))**(1._wp/di)
    n_i=rho(1:kp)*max(q(1:kp,7),0._wp)*lam_i**(1._wp+alpha_i) / gam1i
    ! snow n0, lambda
    lam_s=(max(q(1:kp,8),1._wp)*cs*gam2s / (q(1:kp,4)*gam1s))**(1._wp/ds)
    n_s=rho(1:kp)*max(q(1:kp,8),0._wp)*lam_s**(1._wp+alpha_s) / gam1s
    ! graupel n0, lambda
    lam_g=(max(q(1:kp,9),1._wp)*cg*gam2g / (max(q(1:kp,5),1.e-10)*gam1g))**(1._wp/dg)
    n_g=rho(1:kp)*max(q(1:kp,9),0._wp)*lam_g**(1._wp+alpha_g) / gam1g
    
    
    ! precipitation
	precip(1:kp,1)=cr*n_r*(a_r*chi_rain/(lam_r**(alpha_r+b_r+dr+1._wp)) - &
					u(1:kp)*chi_rain1/(lam_r**(alpha_r+dr+1._wp))) &
					/rho(1:kp) *3600._wp
	precip(1:kp,2)=cs*n_s*(a_s*chi_snow/(lam_s**(alpha_s+b_s+ds+1._wp)) - &
					u(1:kp)*chi_snow1/(lam_s**(alpha_s+ds+1._wp))) &
					/rho(1:kp)*3600._wp
	precip(1:kp,3)=cg*n_g*(a_g*chi_graupel/(lam_g**(alpha_g+b_g+dg+1._wp)) - &
					u(1:kp)*chi_graupel1/(lam_g**(alpha_g+dg+1._wp))) &
					/rho(1:kp)*3600._wp
	precip(1:kp,4)=ci*n_i*(a_i*chi_ice/(lam_i**(alpha_i+b_i+di+1._wp)) - &
					u(1:kp)*chi_ice1/(lam_i**(alpha_i+di+1._wp))) &
					/rho(1:kp)*3600._wp
    
    ! fall speeds
    vqr(1:kp)=max(fall_q_r*rho_fac * lam_r**(1._wp+alpha_r+dr) / &
    	(lam_r+f_r)**(1._wp+alpha_r+dr+b_r), 0._wp)
    vqs(1:kp)=max(fall_q_s*rho_fac * lam_s**(1._wp+alpha_s+ds) / &
    	(lam_s+f_s)**(1._wp+alpha_s+ds+b_s), 0._wp)
    vqg(1:kp)=max(fall_q_g*rho_fac * lam_g**(1._wp+alpha_g+dg) / &
    	(lam_g+f_g)**(1._wp+alpha_g+dg+b_g), 0._wp)
    vqi(1:kp)=max(fall_q_i*rho_fac * lam_i**(1._wp+alpha_i+di) / &
    	(lam_i+f_i)**(1._wp+alpha_i+di+b_i), 0._wp)
    
    vnr(1:kp)=max(fall_n_r*rho_fac * lam_r**(1._wp+alpha_r) / &
    	(lam_r+f_r)**(1._wp+alpha_r+b_r), 0._wp)
    vns(1:kp)=max(fall_n_s*rho_fac * lam_s**(1._wp+alpha_s) / &
    	(lam_s+f_s)**(1._wp+alpha_s+b_s), 0._wp)
    vng(1:kp)=max(fall_n_g*rho_fac * lam_g**(1._wp+alpha_g) / &
    	(lam_g+f_g)**(1._wp+alpha_g+b_g), 0._wp)
    vni(1:kp)=max(fall_n_i*rho_fac * lam_i**(1._wp+alpha_i) / &
    	(lam_i+f_i)**(1._wp+alpha_i+b_i), 0._wp)
    	
    ! coalescence efficiencies
    egi_dry=0.2_wp*exp(0.08*(t(1:kp)-ttr))
    egs_dry=0.2_wp*exp(0.08*(t(1:kp)-ttr))
    esi=0.2_wp*exp(0.08*(t(1:kp)-ttr))
    eii=0.2_wp*exp(0.08*(t(1:kp)-ttr))
    ess=0.2_wp*exp(0.08*(t(1:kp)-ttr))
    
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

        des_dt=dfsid1(svp_liq,t(k),1.e0_wp,1.e-8_wp,err)
        dqs_dt=eps1*p(k)*des_dt/(p(k)-svp_liq(t(k)))**2
        qold=q(k,2)
        qtot=q(k,1)+q(k,2)
        q(k,2)=q(k,1)+q(k,2)-smr(k)
        if (theta_flag) q(k,2)=(q(k,2)+(lv/cp*qold)*dqs_dt) / (1._wp+lv/cp*dqs_dt)
        q(k,2)=max(q(k,2),0._wp)
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
!             q(4,k)=q(4,k)*(q(2,k)/qold)**1._wp
!         endif


	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! evaporation of rain                                                            !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		diff1=dd(t(k),p(k))
		ktherm1=ka(t(k))
		ab_liq=lv**2 / (ktherm1*rv*t(k)**2) + 1._wp/(rho(k)*smr(k)*diff1)
		nu_vis=viscosity_air(t(k)) / rho(k)
		sc=nu_vis / diff1
		nu_rain=2._wp*pi*n_r(k) / rho(k) * &
				(nu_r1 / lam_r(k)**(2._wp+alpha_r) + &
				(a_r/nu_vis)**0.5_wp*sc**(1._wp/3._wp)* &
				(rho(k)*rho0)**0.25_wp*nu_r2 / &
				(lam_r(k)+0.5_wp*f_r)**(0.5_wp*b_r+alpha_r+2.5_wp))

	
		rain_evap=(q(k,1)/smr(k)-1._wp) / (rho(k)*ab_liq)*nu_rain
		if(q(k,1).gt.smr(k)) then
			prevp(k)=0._wp
		else
			prevp(k)=-min(rain_evap,0._wp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		


	
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! deposition & sublimation onto ice, snow, graupel                               !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		nu_snow=2._wp*pi*n_s(k) / rho(k) * &
				(nu_s1 / lam_s(k)**(2._wp+alpha_s) + &
				(a_s/nu_vis)**0.5_wp*sc**(1._wp/3._wp)* &
				(rho(k)*rho0)**0.25_wp*nu_s2 / &
				(lam_s(k)+0.5_wp*f_s)**(0.5_wp*b_s+alpha_s+2.5_wp))

		nu_graup=2._wp*pi*n_g(k) / rho(k) * &
				(nu_g1 / lam_g(k)**(2._wp+alpha_g) + &
				(a_g/nu_vis)**0.5_wp*sc**(1._wp/3._wp)* &
				(rho(k)*rho0)**0.25_wp*nu_g2 / &
				(lam_g(k)+0.5_wp*f_g)**(0.5_wp*b_g+alpha_g+2.5_wp))
	
		if (t(k).le.ttr) then
			nu_ice=2._wp*pi*n_i(k) / rho(k) * &
					(nu_i1 / lam_i(k)**(2._wp+alpha_i) + &
					(a_i/nu_vis)**0.5_wp*sc**(1._wp/3._wp)* &
					(rho(k)*rho0)**0.25_wp*nu_i2 / &
					(lam_i(k)+0.5_wp*f_i)**(0.5_wp*b_i+alpha_i+2.5_wp))

			ab_ice=ls**2 / (ktherm1*rv*t(k)**2) + 1._wp/(rho(k)*smr_i(k)*diff1)
		
			ice_dep=(q(k,1)/smr_i(k)-1._wp) / (rho(k)*ab_ice)*nu_ice
			snow_dep=(q(k,1)/smr_i(k)-1._wp) / (rho(k)*ab_ice)*nu_snow
			graup_dep=(q(k,1)/smr_i(k)-1._wp) / (rho(k)*ab_ice)*nu_graup
			if(q(k,1).gt.smr_i(k)) then
				pisub(k)=0._wp
				pssub(k)=0._wp
				pgsub(k)=0._wp
				pidep(k)=max(ice_dep,0._wp)
				psdep(k)=max(snow_dep,0._wp)
				pgdep(k)=max(graup_dep,0._wp)
			else
				pidep(k)=0._wp
				psdep(k)=0._wp
				pgdep(k)=0._wp			
				pisub(k)=-min(ice_dep,0._wp)
				pssub(k)=-min(snow_dep,0._wp)
				pgsub(k)=-min(graup_dep,0._wp)
			endif
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! warm rain autoconversion                                                       !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		praut(k)=aw0*max(q(k,2)-lw0/rho(k), 0._wp) ! kessler scheme
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! collection of cloud by rain, cloud by snow and cloud by ice                    !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		pracw(k)=max(phi_r* n_r(k)* erw *q(k,2)*rho_fac(k) / &
				(lam_r(k)+f_r)**(3._wp+b_r+alpha_r),0._wp)
		! snow-cloud water
		psacw(k)=max(mass_sacw_i * n_s(k)* esw *q(k,2)*rho_fac(k) / &
				(lam_s(k)+f_s)**(3._wp+b_s+alpha_s),0._wp)
		psaci(k)=max(mass_sacw_i * n_s(k)* esi(k) *q(k,6)*rho_fac(k) / &
				(lam_s(k)+f_s)**(3._wp+b_s+alpha_s),0._wp)
		if (t(k).le.ttr) then
			piacw(k)=max(mass_iacw * n_i(k)* eiw *q(k,2)*rho_fac(k) / &
					(lam_i(k)+f_i)**(3._wp+b_i+alpha_i),0._wp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! collection of rain by ice to make snow or graupel                              !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (t(k).le.ttr) then
			if(q(k,3).lt.1e-4_wp) then
				praci_g(k)=0._wp
				praci_s(k)=max(phi_r* n_r(k)* eri *q(k,6)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_wp+b_r+alpha_r),0._wp)
				piacr_g(k)=0._wp
				piacr_s(k)=max(mass_iacr*n_r(k)*q(k,7)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_wp+b_r+dr+alpha_r),0._wp)
				riacr_g(k)=0._wp
				riacr_s(k)=max(num_iacr*n_r(k)*q(k,7)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_wp+b_r+alpha_r),0._wp)
			else
				praci_s(k)=0._wp
				praci_g(k)=max(phi_r* n_r(k)* eri *q(k,6)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_wp+b_r+alpha_r),0._wp)
				piacr_s(k)=0._wp
				piacr_g(k)=max(mass_iacr*n_r(k)*q(k,7)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_wp+b_r+dr+alpha_r),0._wp)
				riacr_s(k)=0._wp
				riacr_g(k)=max(num_iacr*n_r(k)*q(k,7)*rho_fac(k) / &
						(lam_r(k)+f_r)**(3_wp+b_r+alpha_r),0._wp)
		
			endif
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! snow autoconversion                                                            !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (q(k,6).gt.1e-10_wp.and.(q(k,7).gt.1e-5_wp)) then
			psaut(k)=max(q(k,6)* &
				(min(lambda_imin /lam_i(k),2._wp)**di - 1._wp) / tsaut,0._wp)
			rsaut(k)=max(psaut(k)/(ci*di2s**di),0._wp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! graupel autoconversion                                                         !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (rho(k)*q(k,4).gt.3e-4_wp.and.(t(k).lt.269.15_wp)) then
			pgaut(k)=max(0.5_wp*max(0._wp,psacw(k)-psdep(k)-psaci(k)), 0._wp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! snow and ice aggregation                                                       !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		rsacs(k)=max(isnow*n_s(k)**2._wp*ess(k) *rho_fac(k) / &
 				lam_s(k)**(4._wp+2._wp*alpha_s+b_s),0._wp)
		riaci(k)=max(iice*n_i(k)**2._wp*eii(k) *rho_fac(k) / &
				lam_i(k)**(4._wp+2._wp*alpha_i+b_i),0._wp)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! snow break-up                                                                  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		rsbrk(k)=max((lambda_s_break/lam_s(k)-1._wp)**ds*q(k,8)/(tsbreak),0._wp)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! immersion freezing of rain                                                     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		if(t(k).le.268.15_wp) then
! 			pgfr(k)=max( mass_imm*(exp(-abigg*(t(k)-ttr))-1._wp)*n_r(k) / &
! 						(lam_r(k))**(4_wp+dr+alpha_r),0._wp)
! 			rgfr(k)=num_imm*(exp(-abigg*(t(k)-ttr))-1._wp)*n_r(k)/ &
! 						(lam_r(k))**(4_wp+alpha_r)
! 		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! collisions between precipitating particles of different species                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(t(k).lt.268._wp) then
			pracs(k)=max(n_r(k)*n_s(k)*pi/(4._wp*rho(k))*ers*cs*max((vqs(k)+vqr(k))/8._wp,abs(vqs(k)-vqr(k))) * &
					( &
					mass_racs1/(lam_r(k)**(1._wp+alpha_r) *lam_s(k)**(3._wp+alpha_s+ds)) + &
					mass_racs2/(lam_r(k)**(2._wp+alpha_r) *lam_s(k)**(2._wp+alpha_s+ds)) + &
					mass_racs3/(lam_r(k)**(3._wp+alpha_r) *lam_s(k)**(1._wp+alpha_s+ds))  &
					) , 0._wp)  
				
	!         pracg(k)=max(n_r(k)*n_g(k)*pi/(4._wp*rho(k))*erg*cg*max((vqg(k)+vqr(k))/8._wp,abs(vqg(k)-vqr(k))) * &
	!         		( &
	!         		mass_racg1/(lam_r(k)**(1._wp+alpha_r) *lam_g(k)**(3._wp+alpha_g+dg)) + &
	!         		mass_racg2/(lam_r(k)**(2._wp+alpha_r) *lam_g(k)**(2._wp+alpha_g+dg)) + &
	!         		mass_racg3/(lam_r(k)**(3._wp+alpha_r) *lam_g(k)**(1._wp+alpha_g+dg))  &
	!         		)   , 0._wp)        		
			psacr(k)=max(n_s(k)*n_r(k)*pi/(4._wp*rho(k))*ers*cr*max((vqr(k)+vqs(k))/8._wp,abs(vqr(k)-vqs(k))) * &
					( &
					mass_sacr1/(lam_s(k)**(1._wp+alpha_s) *lam_r(k)**(3._wp+alpha_r+dr)) + &
					mass_sacr2/(lam_s(k)**(2._wp+alpha_s) *lam_r(k)**(2._wp+alpha_r+dr)) + &
					mass_sacr3/(lam_s(k)**(3._wp+alpha_s) *lam_r(k)**(1._wp+alpha_r+dr))  &
					)    , 0._wp)       		
			rsacr(k)=max(n_s(k)*n_r(k)*pi/(4._wp*rho(k))*ers*max((vqr(k)+vqs(k))/8._wp,abs(vqr(k)-vqs(k))) * &
					( &
					num_racs1/(lam_s(k)**(1._wp+alpha_s) *lam_r(k)**(3._wp+alpha_r)) + &
					num_racs2/(lam_s(k)**(2._wp+alpha_s) *lam_r(k)**(2._wp+alpha_r)) + &
					num_racs3/(lam_s(k)**(3._wp+alpha_s) *lam_r(k)**(1._wp+alpha_r))  &
					)    , 0._wp)    
		
	!         psacg(k)=max(n_s(k)*n_g(k)*pi/(4._wp*rho(k))*egs_dry(k)*cg*max((vqg(k)+vqs(k))/8._wp,abs(vqg(k)-vqs(k))) * &
	!         		( &
	!         		mass_sacg1/(lam_s(k)**(1._wp+alpha_s) *lam_g(k)**(3._wp+alpha_g+dg)) + &
	!         		mass_sacg2/(lam_s(k)**(2._wp+alpha_s) *lam_g(k)**(2._wp+alpha_g+dg)) + &
	!         		mass_sacg3/(lam_s(k)**(3._wp+alpha_s) *lam_g(k)**(1._wp+alpha_g+dg))  &
	!         		)    , 0._wp)       		
			pgacr(k)=max(n_g(k)*n_r(k)*pi/(4._wp*rho(k))*erg*cr*max((vqr(k)+vqg(k))/8._wp,abs(vqr(k)-vqg(k))) * &
					( &
					mass_gacr1/(lam_g(k)**(1._wp+alpha_g) *lam_r(k)**(3._wp+alpha_r+dr)) + &
					mass_gacr2/(lam_g(k)**(2._wp+alpha_g) *lam_r(k)**(2._wp+alpha_r+dr)) + &
					mass_gacr3/(lam_g(k)**(3._wp+alpha_g) *lam_r(k)**(1._wp+alpha_r+dr))  &
					)    , 0._wp)      
			! set pgacs and rgacs to wet first 		
			pgacs(k)=max(n_g(k)*n_s(k)*pi/(4._wp*rho(k))*egs_wet*cs*max((vqs(k)+vqg(k))/8._wp,abs(vqs(k)-vqg(k))) * &
					( &
					mass_gacs1/(lam_g(k)**(1._wp+alpha_g) *lam_s(k)**(3._wp+alpha_s+ds)) + &
					mass_gacs2/(lam_g(k)**(2._wp+alpha_g) *lam_s(k)**(2._wp+alpha_s+ds)) + &
					mass_gacs3/(lam_g(k)**(3._wp+alpha_g) *lam_s(k)**(1._wp+alpha_s+ds))  &
					)    , 0._wp)       		
				
			rgacs(k)=max(n_g(k)*n_s(k)*pi/(4._wp*rho(k))*egs_wet*max((vqs(k)+vqg(k))/8._wp,abs(vqs(k)-vqg(k))) * &
					( &
					num_sacg1/(lam_g(k)**(1._wp+alpha_g) *lam_s(k)**(3._wp+alpha_s)) + &
					num_sacg2/(lam_g(k)**(2._wp+alpha_g) *lam_s(k)**(2._wp+alpha_s)) + &
					num_sacg3/(lam_g(k)**(3._wp+alpha_g) *lam_s(k)**(1._wp+alpha_s))  &
					)  , 0._wp)   
		endif
        		
        !rracg?, psacg? pracg? rgacr?   		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! riming of graupel and accretion of ice                                         !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(q(k,5).gt.0._wp) then
			pgacw(k)=max(mass_gacw*rho_fac(k)*n_g(k) / &
					((lam_g(k)+f_g)**(3._wp+b_g+alpha_g))*q(k,2), 0._wp)
		endif
		
		if((q(k,5).gt.0._wp).and.(t(k).lt.ttr)) then
		
			! below assuming wet-growth first
			pgaci(k)=max(mass_gaci*egi_wet*rho_fac(k)*n_g(k) / &
					((lam_g(k)+f_g)**(3._wp+b_g+alpha_g))*q(k,6), 0._wp)
					
			pgwet=(910._wp/(cg*6._wp/pi))**0.625_wp*  &
				(rho(k)*lv*(q0sat-q(k,1))-ktherm1*tc) / (rho(k)*(lf-cw*tc))*nu_graup + &
				(pgaci(k)+pgacs(k))*(1._wp-ci*tc/(lf+cw*tc))
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! redefine as dry growth, if water can be frozen without shedding            !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if((pgacr(k)+pgacw(k)+pgaci(k)+pgacs(k)).gt.pgwet) then
				! redefine graupel accreting rain so that it is
				pgshd(k)=max(0._wp,pgacw(k)+pgaci(k)+pgacs(k)-pgwet)
				pgacr(k)=max(0._wp,pgwet-pgacw(k)-pgaci(k)-pgacs(k))
			else
				! dry growth: graupel-ice
				pgaci(k)=max(mass_gaci*egi_dry(k)*rho_fac(k)*n_g(k) / &
						((lam_g(k)+f_g)**(3._wp+b_g+alpha_g))*q(k,6), 0._wp)
									
				! dry growth: graupel-snow -- mass
				pgacs(k)=max(n_g(k)*n_s(k)*pi/(4._wp*rho(k))*egs_dry(k)* &
						cs*max((vqs(k)+vqg(k))/8._wp,abs(vqs(k)-vqg(k))) * &
						( &
				 mass_gacs1/(lam_g(k)**(1._wp+alpha_g) *lam_s(k)**(3._wp+alpha_s+ds)) + &
				 mass_gacs2/(lam_g(k)**(2._wp+alpha_g) *lam_s(k)**(2._wp+alpha_s+ds)) + &
				 mass_gacs3/(lam_g(k)**(3._wp+alpha_g) *lam_s(k)**(1._wp+alpha_s+ds))  &
						)    , 0._wp)       		
				
				! dry growth: graupel-snow -- number
				rgacs(k)=max(n_g(k)*n_s(k)*pi/(4._wp*rho(k))*egs_dry(k)* &
						max((vqs(k)+vqg(k))/8._wp,abs(vqs(k)-vqg(k))) * &
						( &
					num_sacg1/(lam_g(k)**(1._wp+alpha_g) *lam_s(k)**(3._wp+alpha_s)) + &
					num_sacg2/(lam_g(k)**(2._wp+alpha_g) *lam_s(k)**(2._wp+alpha_s)) + &
					num_sacg3/(lam_g(k)**(3._wp+alpha_g) *lam_s(k)**(1._wp+alpha_s))  &
						)  , 0._wp)   
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! h-m process                                                                    !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(hm_flag) &
			pihal(k)=max(hm_rate*mi0*(pgacw(k)+psacw(k))*hm_func(t(k)),0._wp)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! melting of ice, snow, and graupel                                              !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(t(k).gt.ttr) then
			pimlt(k)=q(k,6)/dt ! ice melts instantaneously
			q(k,7)=0._wp
			if(q(k,5).gt.0._wp) then
				pgmlt(k)=max(1._wp/(rho(k)*lf) * &
					(ktherm1*tc-lv*diff1*rho(k)*(q(k,1)-q0sat))*nu_graup &
					+ cw*tc/lf*(pgacw(k)+pgacr(k)-pgshd(k)),0._wp)
			endif
				
			psmlt(k)=max(1._wp/(rho(k)*lf) * &
				(ktherm1*tc-lv*diff1*rho(k)*(q(k,1)-q0sat))*nu_snow &
				+ cw*tc/lf*(psacw(k)+psacr(k)),0._wp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    enddo
    
    
    
    
    
    
    ! update variables
    ! ice number
    
    
    q(1:kp,7)=q(1:kp,7)+((piprm+pihal+picnt+pifrw)/mi0 - riaci - &
    		(pgaci+psaci+praci_g+praci_s)*q(1:kp,7)/(q(1:kp,6)+qsmall)  &
    		-riaci - rsaut)*dt
    where(pimlt.gt.0._wp)
    	q(1:kp,7)=0._wp
    end where
    ! snow number
    q(1:kp,8)=q(1:kp,8)+(riacr_s+rsaut+rsbrk &
    			-(pssub+psmlt+pgaut)*q(1:kp,8)/(q(1:kp,4)+qsmall) &
    			-(rgacs+rsacr+rsacs))*dt
    ! graupel number
    q(1:kp,9)=q(1:kp,9)+(pgaut*q(1:kp,9)/(q(1:kp,5)+qsmall) + rsacr+riacr_g+rgfr &
    			-(pgsub+pgmlt)*q(1:kp,9)/(q(1:kp,5)+qsmall))*dt
    			

    ! vapour mass
    q(1:kp,1)=q(1:kp,1)+(pgsub+pssub+prevp+pisub-(psdep+pidep+piprm+pgdep))*dt
    ! liquid mass
    q(1:kp,2)=q(1:kp,2)-((pgacw+praut+psacw+pracw+piacw+pihal+picnt+pifrw))*dt
    ! rain mass
    q(1:kp,3)=q(1:kp,3)+(pgmlt+praut+pgshd+pracw+psmlt+pimlt- &
    			(pgacr+prevp+pgfr+psacr+piacr_g+piacr_s))*dt
    ! snow mass
    q(1:kp,4)=q(1:kp,4)+(psaut+psdep+psaci+praci_s+piacr_s+psacw-(pssub+pgacs+pracs+pgaut+psmlt))*dt
    ! graupel mass
    q(1:kp,5)=q(1:kp,5)+(pgaci+pgacw+pgacs+pgacr+psacr+pracs+pgaut+pgfr+praci_g+piacr_g+pgdep- &
    				(pgsub+pgmlt+pgshd))*dt
!     if(sum(piacw).gt.0._wp) stop
    ! ice mass
    q(1:kp,6)=q(1:kp,6)+(pidep+piprm+pihal+picnt+piacw- &
    			(psaut+pgaci+psaci+pisub+pifrw+pimlt+praci_g+praci_s))*dt
    			
    q=max(q,0._wp)	    
    if (theta_flag) th=t*(1.e5_wp/p)**(ra/cp)-theta
    
 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! advection rain 0th order Bott, a.k.a. upstream advection                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! rain 
    if(sum(q(1:kp,3)).gt.qsmall) then
		where(isnan(vqr))
			vqr=0._wp
		end where
		vqr(1-o_halo:0)=vqr(1)
		vqr(kp+1:kp+o_halo)=vqr(kp)
		n_step=max(ceiling(maxval(vqr)*dt/dz*2_wp),1)
		vqr(1-o_halo:kp+o_halo-1)=-vqr(-o_halo+2:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,wp),dz,z,vqr,q(:,3),.false.)
		enddo
	endif
	
    ! graupel 
    if(sum(q(1:kp,5)).gt.qsmall) then
		where(isnan(vqg))
			vqg=0._wp
		end where
		vqg(1-o_halo:0)=vqg(1)
		vqg(kp+1:kp+o_halo)=vqg(kp)
		n_step=max(ceiling(maxval(vqg)*dt/dz*2_wp),1)
		vqg(1-o_halo:kp+o_halo-1)=-vqg(-o_halo+2:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,wp),dz,z,vqg,q(:,5),.false.)
		enddo
		where(isnan(vng))
			vng=0._wp
		end where
		vng(1-o_halo:0)=vng(1)
		vng(kp+1:kp+o_halo)=vng(kp)
		n_step=max(ceiling(maxval(vng)*dt/dz*2_wp),1)
		vng(1-o_halo:kp+o_halo-1)=-vng(-o_halo+2:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,wp),dz,z,vng,q(:,9),.false.)
		enddo
	endif
	
    ! snow 
    if(sum(q(1:kp,4)).gt.qsmall) then
		where(isnan(vqs))
			vqs=0._wp
		end where
		vqs(1-o_halo:0)=vqs(1)
		vqs(kp+1:kp+o_halo)=vqs(kp)
		n_step=max(ceiling(maxval(vqs)*dt/dz*2_wp),1)
		vqs(1-o_halo:kp+o_halo-1)=-vqs(-o_halo+2:kp+o_halo)
		do iter=1,n_step		
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,wp),dz,z,vqs,q(:,4),.false.)
		enddo
		where(isnan(vns))
			vns=0._wp
		end where
		vns(1-o_halo:0)=vns(1)
		vns(kp+1:kp+o_halo)=vns(kp)
		n_step=max(ceiling(maxval(vns)*dt/dz*2_wp),1)
		vns(1-o_halo:kp+o_halo-1)=-vns(-o_halo+2:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,wp),dz,z,vns,q(:,8),.false.)
		enddo
	endif
    ! ice
    if(sum(q(1:kp,6)).gt.qsmall) then     
		where(isnan(vqi))
			vqi=0._wp
		end where
		vqi(1-o_halo:0)=vqi(1)
		vqi(kp+1:kp+o_halo)=vqi(kp)
		n_step=max(ceiling(maxval(vqi)*dt/dz*2_wp),1)
		vqi(1-o_halo:kp+o_halo-1)=-vqi(-o_halo+2:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,wp),dz,z,vqi,q(:,6),.false.)
		enddo
		where(isnan(vni))
			vni=0._wp
		end where
		vni(1-o_halo:0)=vni(1)
		vni(kp+1:kp+o_halo)=vni(kp)
		n_step=max(ceiling(maxval(vni)*dt/dz*2_wp),1)
		vni(1-o_halo:kp+o_halo-1)=-vni(-o_halo+2:kp+o_halo)
		do iter=1,n_step
			call bott_scheme_1d(kp,0,o_halo,dt/real(n_step,wp),dz,z,vni,q(:,7),.false.)
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
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: svp_liq
		svp_liq = 100._wp*6.1121_wp* &
			  exp((18.678_wp - (t-ttr)/ 234.5_wp)* &
			  (t-ttr)/(257.14_wp + (t-ttr)))
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
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: svp_ice
		svp_ice = 100._wp*6.1115_wp* &
			  exp((23.036_wp - (t-ttr)/ 333.7_wp)* &
			  (t-ttr)/(279.82_wp + (t-ttr)))
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
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: viscosity_air
		real(wp) :: tc

		tc = t-ttr
		tc = max(tc,-200._wp)

		if( tc.ge.0._wp) then
			viscosity_air = (1.718_wp+0.0049_wp*tc) * 1E-5_wp ! the 1d-5 converts from poise to si units
		else
			viscosity_air = (1.718_wp+0.0049_wp*tc-1.2e-5_wp*tc**2) * 1e-5_wp
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
	  use numerics_type
	  implicit none
	  real(wp), intent(in) :: t, p
	  real(wp) :: dd, t1
	  t1=max(t,200._wp)
	  dd=2.11e-5_wp*(t1/ttr)**1.94_wp*(101325_wp/p)
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
	  real(wp) :: ka, t1
	  t1=max(t,200._wp)
	  ka=(5.69_wp+0.017_wp*(t1-ttr))*1e-3_wp*joules_in_a_cal
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
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: hm_func
		hm_func=(min(max((t-265.65) / 2.5_wp,0._wp),1._wp) + &
		        min(max((270.65-t) / 2.5_wp,0._wp),1._wp)) -1.0_wp
	end function hm_func


    end module micro_module
    
