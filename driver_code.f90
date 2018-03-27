	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>drivers and physics code for the thermal cloud model
    module drivers
    use nrtype
    use variables
    private
    public :: model_driver_2d
    contains
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calls IO and runs model microphysics and advection for a number of steps
	!>writes to output file, solves for the microphysics over one time-step, 
	!> then advects particles
	!>@param[in] nq: number of q fields
	!>@param[in] nprec: number of precipitation arrays
	!>@param[in] ip: number of horizontal levels
	!>@param[in] kp: number of vertical levels
	!>@param[in] ord: order of advection scheme
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] runtime
	!>@param[in] dt - timestep
	!>@param[in] dx,dz - grid spacing
	!>@param[inout] q, theta, pressure, x,xn,z,zn, temperature, rho, u,w
	!>@param[inout] precip
	!>@param[inout] new_file
	!>@param[inout] micro_init - flag to initialise microphysics
	!>@param[in] advection_scheme
	!>@param[in] monotone - flag for monotonic advection
	!>@param[in] microphysics_flag - flag for calling microphysics
	!>@param[in] hm_flag - flag for switching on / off hm process
	!>@param[in] theta_flag - flag for advecting theta dry
	!>@param[in] mass_ice - mass of a single ice crystal (override)
	!>@param[inout] k
	!>@param[inout] dsm_by_dz_z_eq_zc
	!>@param[inout] b
	!>@param[inout] del_gamma_mac
	!>@param[inout] del_c_s
	!>@param[inout] del_c_t
	!>@param[inout] epsilon_therm
	!>@param[in] w_peak
	!>@param[in] z_offset
	!>@param[inout] therm_init
    subroutine model_driver_2d(nq,nprec,ip,kp,ord,o_halo,runtime, &
                               dt, &
                               q,precip,theta,p,dx,dz,x,xn,z,zn,t,rho,&
                               u,w,new_file,micro_init,advection_scheme, monotone, &
                               microphysics_flag,hm_flag,theta_flag,mass_ice, &
                               ! variables associated with thermal properties
                               k,dsm_by_dz_z_eq_zc,b,del_gamma_mac, & 
                               del_c_s,del_c_t,epsilon_therm,w_peak,z_offset, therm_init)

    use nrtype
    use thermal
    use io_module
    use advection_2d
    use micro_module

    implicit none
    integer(i4b), intent(in) :: nq,nprec,ip,kp, ord, o_halo, advection_scheme
    real(sp), intent(in) :: runtime, dt, dx,dz
    real(sp), dimension(nq,-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: q
    real(sp), dimension(nprec,1:kp,1:ip), intent(inout) :: precip
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(inout) :: z,zn
    real(sp), dimension(-o_halo+1:ip+o_halo), intent(inout) :: x,xn
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: theta, p, &
    																			t,rho,u,w
    logical, intent(inout) :: new_file, micro_init,therm_init
    logical, intent(in) :: monotone,hm_flag,theta_flag
    integer(i4b), intent(in) :: microphysics_flag
    real(sp), intent(in) :: mass_ice
    
    ! variables associated with thermals:
    real(sp), intent(inout) :: k,dsm_by_dz_z_eq_zc,b,del_gamma_mac,del_c_s,del_c_t, &
    						epsilon_therm 
    real(sp), intent(in) :: w_peak, z_offset

    ! local variables
    integer(i4b) :: nt, i, j, nsteps, iter
    real(sp) :: time


    nt=ceiling(runtime / real(dt,kind=sp) )
    do i=1,nt
    	time=real(i-1,sp)*dt
        print *,'time-step ',i,' of ',nt, ' time=',time

		
    	! wind / thermal properties:
		call thermal_2d(time,ip,kp,o_halo,k,dsm_by_dz_z_eq_zc, &
						b,del_gamma_mac,del_c_s,del_c_t, &
    						epsilon_therm,x,xn,z,zn,dx,dz,u,w,w_peak,z_offset, therm_init)   	
    	
        ! output:
        call output_2d(time,nq,nprec,ip,kp,q(:,1:kp,1:ip),precip(:,1:kp,1:ip), &
						theta(1:kp,1:ip),p(1:kp,1:ip), &
						x(1:ip),xn(1:ip),z(1:kp),zn(1:kp), &
						t(1:kp,1:ip),u(1:kp,1:ip),w(1:kp,1:ip),new_file)



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! advection                                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		nsteps=ceiling(dt/(0.7_sp*dz)*maxval(sqrt(u**2+w**2)))
		select case (advection_scheme)
			case (0) ! upstream
				do iter=1,nsteps
					do j=1,nq
						! set halos
						call set_halos_2d(ip,kp,o_halo,q(j,:,:))            
						! advection
						call first_order_upstream_2d(dt/real(nsteps,sp), &
								dx,dz,ip,kp,1,u(0:kp+1,0:ip+1),w(0:kp+1,0:ip+1), &
								q(j,0:kp+1,0:ip+1))
					enddo
					if(theta_flag) then
						! set halos
						call set_halos_2d(ip,kp,o_halo,theta)        
						! advection
						call first_order_upstream_2d(dt/real(nsteps,sp), &
								dx,dz,ip,kp,1,u(0:kp+1,0:ip+1),w(0:kp+1,0:ip+1), &
								theta(0:kp+1,0:ip+1))
					endif
				enddo
			case(1) ! mpdata
				do iter=1,nsteps
					do j=1,nq
						! set halos
						call set_halos_2d(ip,kp,o_halo,q(j,:,:))            
						! advection
						call mpdata(4,ip,kp,o_halo,xn,zn,dx,dz,dt/real(nsteps,sp),u(:,:), &
							 w(:,:),q(j,:,:),monotone)
					enddo
					if(theta_flag) then
						! set halos
						call set_halos_2d(ip,kp,o_halo,theta)        
						! advection
						call mpdata(4,ip,kp,o_halo,xn,zn,dx,dz,dt/real(nsteps,sp),u(:,:), &
						     w(:,:),theta(:,:),monotone)
					endif
				enddo
			
			case default
				print *, 'error'
				stop
		end select
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! solve microphysics. initialise constants that are commonly used, if needed
        if (microphysics_flag .eq. 1) then
			call microphysics_2d(nq,ip,kp,o_halo,dt,dz,q(:,:,:),precip(:,:,:),&
							theta(:,:),p(:,:), &
						   zn(:),t,rho(:,:),w(:,:),micro_init,hm_flag,mass_ice)			
		endif       






    enddo



    end subroutine model_driver_2d


    
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>set halos 2-d 
	!>@param[in] ip, kp:  number of grid points
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[inout] psi: field to set
    subroutine set_halos_2d(ip,kp,o_halo,psi)  

    use nrtype
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: ip, kp, o_halo
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: psi

	integer(i4b) :: i,k
	do i=1,ip
		! top
		psi(kp+1:kp+o_halo,i)=psi(kp,i)
		! bottom
		psi(-o_halo+1:0,i)=psi(1,i)
	enddo
		
	do k=1,kp
		! left
		psi(k,-o_halo+1:0)=psi(k,ip-o_halo+1:ip)
		! right
		psi(k,ip+1:ip+o_halo)=psi(k,1:o_halo)
	enddo 
		
	end subroutine set_halos_2d

    end module drivers
