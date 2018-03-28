	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>initialisation for the simple cloud model and solves the 
	!>hydrostatic equation for pressure:
	!>\f$ \frac{\partial P}{\partial z} = - \frac{P}{R_aT}g \f$
    module initialisation
    use nrtype
!    use variables
    private
    public :: calc_profile_2d, allocate_arrays
    contains
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays for q_type and q_init
	!>@param[in] nq number of q fields
	!>@param[in] n_levels number of levels for reading in sounding
	!>@param[inout] q_type: integer array
	!>@param[inout] q_init: logical array
	!>@param[inout] q_read: real array
    subroutine allocate_arrays(nq,n_levels,q_type,q_init,q_read)
    use nrtype
    implicit none
    integer(i4b), intent(in) :: nq, n_levels
    integer(i4b), dimension(:), allocatable, intent(inout) :: q_type
    logical, dimension(:), allocatable, intent(inout) :: q_init
    real(sp), dimension(:,:), allocatable, intent(inout) :: q_read
    ! local variables:
    integer(i4b) :: AllocateStatus
    
    ! allocate arrays
    allocate( q_type(1:nq), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( q_init(1:nq), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( q_read(1:nq,1:n_levels), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    
    
    end subroutine allocate_arrays
    
    
    
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>interpolates the sounding to the grid
	!>@param[in] nq number of q fields
	!>@param[in] nprec number of precipitation arrays
	!>@param[in] n_levels number levels for sounding
	!>@param[in] psurf surface pressure
	!>@param[in] tsurf surface temperature
	!>@param[in] t_cbase cloud base temperature
	!>@param[in] t_ctop cloud top temperature
	!>@param[in] adiabatic_prof: flag if we want an adiabatic profile
	!>@param[in] adiabatic_frac: fraction of adiabatic liquid water content in cloud
	!>@param[in] q_type flag for type of q-field
	!>@param[in] q_init flag for whether q-field initialised
	!>@param[in] z_read vertical levels for sounding
	!>@param[in] theta_read theta on vertical levels for sounding
	!>@param[in] q_read q fields on vertical levels for sounding
	!>@param[in] ip number of vertical levels of grid
	!>@param[in] kp number of vertical levels of grid
	!>@param[in] o_halo number of extra grid levels required for advection
	!>@param[in] dx horizontal resolution of grid
	!>@param[in] dz vertical resolution of grid
	!>@param[inout] q, precip, theta, pressure, x,xn,z,zn, temperature, rho,u,w
	!>@param[in] ice_init: flag to initialise ice crystals in model
	!>@param[in] number conc of ice crystals #/kg
	!>@param[in] mass of a single ice crystal kg.
    subroutine calc_profile_2d(nq,nprec,n_levels,psurf,tsurf,t_cbase, &
    						t_ctop, adiabatic_prof, adiabatic_frac,q_type,q_init, &
                             z_read,theta_read,q_read, &
                             ip,kp,o_halo,dx,dz,q,precip,theta,p,x,xn,z,zn,t,rho,u,w, &
                             ice_init, num_ice, mass_ice)
    use nrtype
    use nr, only : locate, polint, rkqs, odeint, zbrent
    use constants
    use variables, only : theta_surf, theta_q_sat, w_cb, t1old, p111
    !use micro_module, only : svp_liq

    implicit none
    ! inputs
    integer(i4b), intent(in) :: n_levels, nq,o_halo, nprec
    real(sp), dimension(n_levels), intent(in) :: z_read, theta_read
    real(sp), dimension(nq,n_levels), intent(in) :: q_read
    integer(i4b), dimension(nq), intent(in) :: q_type
    logical, dimension(nq), intent(in) :: q_init
    integer(i4b), intent(in) :: ip, kp
    real(sp), intent(in) :: dx, dz, psurf, tsurf, t_cbase, t_ctop
    logical, intent(in) :: adiabatic_prof, ice_init
    real(sp), intent(in) :: adiabatic_frac
    real(sp), intent(in) :: num_ice, mass_ice
    ! inouts
    real(sp), dimension(:,:), allocatable, intent(inout) :: theta, p, t, rho,u, w
    real(sp), dimension(:), allocatable, intent(inout) :: x, z,xn,zn
    real(sp), dimension(:,:,:), allocatable, intent(inout) :: q, precip
    ! local variables:
    integer(i4b) :: i,j, iloc, AllocateStatus, istore,istore2
    real(sp) :: var, dummy
	! variables for odesolver:
	real(sp), dimension(1) :: z1, p11,p22
	real(sp) :: htry,hmin,eps2,p1,p2,p_ctop,z11,z22,theta1

    ! allocate arrays
    allocate( precip(1:nprec,1:kp,1:ip), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( q(1:nq,-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( theta(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( p(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( x(-o_halo+1:ip+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( z(-o_halo+1:kp+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( xn(-o_halo+1:ip+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( zn(-o_halo+1:kp+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( t(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( rho(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( u(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( w(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"


 	precip=0._sp
    ! set up vertical level array
    z=dz*(/(i,i=-o_halo,kp+o_halo-1)/)
    zn=z+0.5_sp*dz
    ! set up horizontal level array
    x=dx*(/(i,i=-o_halo,ip+o_halo-1)/)!-0.5_sp*dx
    x=x+0.5_sp*dx
    xn=x+0.5_sp*dx

	q=0._sp
	if (adiabatic_prof) then
		! calculate the dry adiabat:
		theta_surf=tsurf*(1e5_sp/psurf)**(ra/cp)
		p1=psurf
		z1=0._sp
		p2=1.e5_sp*(t_cbase/theta_surf)**(cp/ra)
		htry=p2-psurf
		eps2=1e-5_sp
		call odeint(z1,p1,p2,eps2,htry,hmin,hydrostatic1,rkqs)
		p1=p2
		
		! integrate going downwards - dry adiabatic layer
		p(1,:)=psurf
		do i=1,-o_halo+2,-1
			z11=z(i)
			z22=z(i-1)
			p11=p(i,1)
			htry=-dz
			hmin=-1.e-2_sp
			call odeint(p11,z11,z22,eps2,htry,hmin,hydrostatic1a,rkqs)
			p(i-1,:)=p11(1)
		enddo
		! integrate going upwards - dry adiabatic layer
		p(1,:)=psurf
		do i=1,kp+o_halo-1
			z11=z(i)
			z22=z(i+1)
			if(z22.gt.z1(1)) exit
			p11=p(i,1)
			htry=dz
			hmin=1.e-2_sp
			call odeint(p11,z11,z22,eps2,htry,hmin,hydrostatic1a,rkqs)
			p(i+1,:)=p11(1)
		enddo
		istore=i-1
		! adiabatic temperature
		t(-o_halo+1:istore,:)=theta_surf*(p(-o_halo+1:istore,:)/1.e5_sp)**(ra/cp)
		! adiabatic vapour mixing ratio
		q(1,-o_halo+1:istore,:)=eps1*svp_liq(t_cbase)/(p1-svp_liq(t_cbase))


		! now calculate the moist adiabat
		w_cb=eps1*svp_liq(t_cbase)/(p1-svp_liq(t_cbase))
		theta_q_sat=t_cbase*(1.e5_sp/p1)**(ra/cp)*exp(lv*w_cb/cp/t_cbase)	
		! theta_q_sat is conserved. Use it to calculate the new temperature
		t1old=t_ctop

		p_ctop=zbrent(calc_theta_q2,p1,3000._sp,1.e-5_sp)

!		stop
		! integrate going upwards - moist adiabatic layer
		t1old=t_ctop
		istore=istore+1
		! now find the temperature
		p111=p(istore,1)
		t(istore,:)=theta_surf*( &
			(p(istore,:)-dz*p(istore,:)/ra/t_ctop)/1.e5_sp)**(ra/cp) ! a temperature colder than
		                                                   ! next level

		t(istore,:)=zbrent(calc_theta_q,1.01*t(istore,1),t_ctop,1.e-5_sp)

		t1old=t(istore,1)
		do i=istore,kp+o_halo-1
			
			z22=z(i+1)
			if (i.eq.istore) then
				z11=z1(1)
				p11=p1
			else
				z11=z(i)
				p11=p(i,1)
			endif
			htry=dz
			hmin=1.e-2_sp
			!print *,t1old,p11,z11,z22
			call odeint(p11,z11,z22,eps2,htry,hmin,hydrostatic2a,rkqs)
			p(i+1,:)=p11(1)
			t(i+1,:)=theta_surf*(p(i+1,:)/1e5_sp)**(ra/cp)
			t1old=t(i,1)
			p111=p(i+1,1)
			t(i+1,:)=zbrent(calc_theta_q,t(i+1,1),t1old*1.01_sp,1.e-5_sp)
			if(t(i+1,1).lt.t_ctop) exit
		enddo
		istore2=i-1
		do i=istore,istore2
			q(1,i,:)=eps1*svp_liq(t(i,1))/ &
								(p(i,1)-svp_liq(t(i,1)))
			q(2,i,:)=adiabatic_frac* &
					max(eps1*svp_liq(t_cbase)/(p1-svp_liq(t_cbase)) - q(1,i,1),0._sp)
		enddo

		! integrate going upwards - dry adiabatic layer
		theta1=t_ctop*(1.e5_sp/p_ctop)**(ra/cp)
		istore2=istore2+1
		do i=istore2,kp+o_halo-1
			z11=z(i)
			z22=z(i+1)
			p11=p(i,1)
			htry=dz
			hmin=1.e-2_sp
			call odeint(p11,z11,z22,eps2,htry,hmin,hydrostatic1a,rkqs)
			p(i+1,:)=p11(1)
		enddo
		t(istore2:kp+o_halo,:)=theta1*(p(istore2:kp+o_halo,:)/1.e5_sp)**(ra/cp)

		! initialise ice crystals
		if(ice_init) then
            where(t(istore:istore2,:).lt.ttr)
                q(6,istore:istore2,:)=num_ice*mass_ice
                q(7,istore:istore2,:)=num_ice
            end where
        endif
	else
		! use linear interpolation to put sounding on grid:
		do i=-o_halo+1,kp+o_halo
			iloc=locate(z_read(1:n_levels),z(i))
			iloc=min(n_levels-1,iloc)
			iloc=max(1,iloc)
			! linear interp theta
			call polint(z_read(iloc:iloc+1), theta_read(iloc:iloc+1), &
						min(z(i),z_read(n_levels)), var,dummy)
			theta(i,:)=var


			! q-fields
			do j=1,nq
				if(q_init(j)) then
					! linear interp q fields
					call polint(z_read(iloc:iloc+1), q_read(j,iloc:iloc+1), &
								min(z(i),z_read(n_levels)), var, dummy)
					q(j,i,:)=var
				else
					if(q_type(j).eq.2) then
						q(j,i,:) = .0_sp
					else
						q(j,i,:) = 0.0_sp
					endif
				endif
			enddo

		enddo
	
		p(1,:)=psurf
		do i=1,kp+o_halo
			! solve the hydrostatic equation to get pressure and temperature
			t(i,:)=theta(i,:)*(p(i,:)/psurf)**(ra/cp)
			if(i.lt.(kp+o_halo)) then
				p(i+1,:)=p(i,:)-dz*p(i,:)/(ra*t(i,:))*grav
			endif
		enddo

		do i=0,-o_halo+1,-1
			! solve the hydrostatic equation to get pressure and temperature
			p(i,:)=p(i+1,:)+dz*p(i+1,:)/(ra*t(i+1,:))*grav
			t(i,:)=theta(i,:)*(p(i,:)/psurf)**(ra/cp)
		enddo
	endif

    u(:,:)=0._sp
    w(:,:)=0._sp
    theta=t*(1.e5_sp/p)**(ra/cp)
    rho=p/(ra*t)
    end subroutine calc_profile_2d


	subroutine hydrostatic1(p,z,dzdp)
	use nrtype
	use constants
	use variables, only : theta_surf
	implicit none
	real(sp), intent(in) :: p
	real(sp), dimension(:), intent(in) :: z
	real(sp), dimension(:), intent(out) :: dzdp
	real(sp) :: t
	
	t=theta_surf*(p/1.e5_sp)**(ra/cp)
	dzdp(1)=-(ra*t) / (grav*p)
	
	end subroutine hydrostatic1

	subroutine hydrostatic1a(z,p,dpdz)
	use nrtype
	use constants
	use variables, only : theta_surf
	implicit none
	real(sp), intent(in) :: z
	real(sp), dimension(:), intent(in) :: p
	real(sp), dimension(:), intent(out) :: dpdz
	real(sp) :: t
	
	t=theta_surf*(p(1)/1.e5_sp)**(ra/cp)
	dpdz(1)=-(grav*p(1)) / (ra*t) 
	
	end subroutine hydrostatic1a

	subroutine hydrostatic2(p,z,dzdp)
	use nrtype
	use nr, only : zbrent
	use constants
	use variables, only : theta_surf,theta_q_sat, w_cb, t1old, p111
	implicit none
	real(sp), intent(in) :: p
	real(sp), dimension(:), intent(in) :: z
	real(sp), dimension(:), intent(out) :: dzdp
	real(sp) :: t
	
	p111=p
	t=theta_surf*(p111/1.e5_sp)**(ra/cp)
	t=zbrent(calc_theta_q,t,t1old*1.01_sp,1.e-5_sp)
!	print *,'hi',t,calc_theta_q(t)
	! find the temperature by iteration
	dzdp(1)=-(ra*t) / (grav*p)
	
	end subroutine hydrostatic2

	subroutine hydrostatic2a(z,p,dpdz)
	use nrtype
	use nr, only : zbrent
	use constants
	use variables, only : theta_surf,theta_q_sat, w_cb, t1old, p111
	implicit none
	real(sp), intent(in) :: z
	real(sp), dimension(:), intent(in) :: p
	real(sp), dimension(:), intent(out) :: dpdz
	real(sp) :: t
	
	p111=p(1)
	t=theta_surf*(p111/1.e5_sp)**(ra/cp)
	t=zbrent(calc_theta_q,t,t1old*1.01_sp,1.e-5_sp)
!	print *,'hi',t,calc_theta_q(t)
	! find the temperature by iteration
	dpdz(1)=-(grav*p(1))/(ra*t)
	
	end subroutine hydrostatic2a
	

	function calc_theta_q(t111)
	use nrtype
	use constants
	use variables, only : theta_q_sat, p111
	!use micro_module, only : svp_liq
	implicit none
	real(sp), intent(in) :: t111
	real(sp) :: calc_theta_q
	real(sp) :: ws
	ws=eps1*svp_liq(t111)/(p111-svp_liq(t111))
	calc_theta_q=t111*(1.e5_sp/p111)**(ra/cp)*exp(lv*ws/cp/t111)-theta_q_sat

	end function calc_theta_q     

	function calc_theta_q2(p)
	use nrtype
	use constants
	use variables, only : theta_q_sat, t1old
	!use micro_module, only : svp_liq
	implicit none
	real(sp), intent(in) :: p
	real(sp) :: calc_theta_q2
	real(sp) :: ws
	ws=eps1*svp_liq(t1old)/(p-svp_liq(t1old))
	calc_theta_q2=t1old*(1e5_sp/p)**(ra/cp)*exp(lv*ws/cp/t1old)-theta_q_sat

	end function calc_theta_q2    


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
		use constants, only : ttr
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: svp_liq
		svp_liq = 100._sp*6.1121_sp* &
			  exp((18.678_sp - (t-ttr)/ 234.5_sp)* &
			  (t-ttr)/(257.14_sp + (t-ttr)))
	end function svp_liq



    end module initialisation

