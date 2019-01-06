	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>drivers and physics code for the thermal cloud model
    module thermal
    use nrtype
    use nr, only : zbrent
    private
    public :: thermal_2d, fd_thermal_2d, adjust_thermal
    
    
	! variables for thermal properties
	
    real(sp), parameter :: Md=29.e-3_sp,Mv=18.e-3_sp,grav=9.81_sp,cp=1005._sp, &
    					lv=2.5e6_sp,ttr=273.15_sp, small1=1.e-30_sp
	real(sp) :: zc,beta1=Md/Mv-1._sp, alpha1=1./ttr, klarge, z_bar, k1,cell_size
	real(sp) :: ztop1, del_gamma_mac1,dsm_by_dz_z_eq_zc1,b1,del_c_s1,del_c_t1
	complex(sp) :: n_bar_mac
    contains
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>adjusts the thermal so it is consistent with cloud base and top temperatures
	!>@param[in] k
	!>@param[in] dsm_by_dz_z_eq_zc
	!>@param[in] b
	!>@param[in] del_gamma_mac
	!>@param[in] del_c_s
	!>@param[in] del_c_t
	!>@param[in] epsilon_therm
	!>@param[in] z_offset
	!>@param[in]  zbase,ztop
	!>@param[in] offset_equal_zbase
    subroutine adjust_thermal(k,dsm_by_dz_z_eq_zc,b,del_gamma_mac,del_c_s,del_c_t, &
    						epsilon_therm,z_offset, zbase,ztop,offset_equal_zbase)
        use nrtype
    
        implicit none
        real(sp), intent(inout) :: k,dsm_by_dz_z_eq_zc,b,del_gamma_mac,del_c_s,del_c_t, &
                                epsilon_therm,z_offset
        real(sp), intent(in) :: ztop, zbase
        logical, intent(in) :: offset_equal_zbase
    
        del_gamma_mac1=del_gamma_mac
        dsm_by_dz_z_eq_zc1=dsm_by_dz_z_eq_zc
        b1=b
        del_c_t1=del_c_t
        del_c_s1=del_c_s
        dsm_by_dz_z_eq_zc1=dsm_by_dz_z_eq_zc
        if(offset_equal_zbase) z_offset=zbase
        ztop1=ztop-z_offset
        ! equation 32:
        n_bar_mac=grav*(alpha1*del_gamma_mac1+beta1*(dsm_by_dz_z_eq_zc1+b1))
        n_bar_mac=sqrt(n_bar_mac)
        ! equation 33:
        z_bar=(alpha1*del_c_t1+beta1*del_c_s1)/ &
            (alpha1*del_gamma_mac1+beta1*(dsm_by_dz_z_eq_zc1+b1))

        epsilon_therm=zbrent(calc_therm_height,1.e-30_sp,1._sp,1.e-20_sp)

        
    end subroutine adjust_thermal
    
    function calc_therm_height(epsilon1)
        use nrtype
        implicit none
        real(sp), intent(in) :: epsilon1
        real(sp) :: calc_therm_height,delta_z
        
        klarge=epsilon1*cp/lv

        k1=0.5_sp*(alpha1*epsilon1-beta1*klarge)/ &
            (alpha1*del_gamma_mac1+beta1*(dsm_by_dz_z_eq_zc1+b1))
        
        
        
        delta_z=3._sp/(4._sp*k1)*(1._sp+sqrt(1._sp+16._sp/3._sp*k1*z_bar))

        calc_therm_height=delta_z-ztop1
	end function calc_therm_height 
    
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>sets the vertical wind speed
	!>@param[in] time
	!>@param[in] ip
	!>@param[in] kp
	!>@param[in] o_halo
	!>@param[in] k
	!>@param[in] dsm_by_dz_z_eq_zc
	!>@param[in] b
	!>@param[in] del_gamma_mac
	!>@param[in] del_c_s
	!>@param[in] del_c_t
	!>@param[in] epsilon_therm
	!>@param[in] x,xn,z,zn,dx,dz
	!>@param[inout]  u, w,
	!>@param[in] w_peak , z_offset
	!>@param[inout]  therm_init
    subroutine thermal_2d(time,ip,kp,o_halo,k,dsm_by_dz_z_eq_zc,b,del_gamma_mac,del_c_s,del_c_t, &
    						epsilon_therm,x,xn,z,zn,dx,dz,u,w,w_peak,z_offset, therm_init)

    use nrtype

    implicit none
    integer(i4b), intent(in) :: ip,kp, o_halo
    real(sp), intent(in) :: time
    real(sp), intent(inout) :: k,dsm_by_dz_z_eq_zc,b,del_gamma_mac,del_c_s,del_c_t, &
    						epsilon_therm
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: &
    					u,w
    real(sp), dimension(-o_halo+1:ip+o_halo), intent(in) :: x,xn
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z,zn
    real(sp), intent(in) :: dx,dz, w_peak, z_offset
    logical, intent(inout) :: therm_init
    
    ! local variable
    real(sp), dimension(-o_halo:kp+o_halo,-o_halo:ip+o_halo) :: xx,zz,phi
    real(sp), dimension(-o_halo:kp+o_halo,-o_halo:ip+o_halo) :: wcen
	integer(i4b) :: i,j
	real(sp) :: test3
	complex(sp) :: test1,test2 

    ! calculate the thermal properties
	if(therm_init) then
		zc=zn(2)-z_offset !0._sp !500._sp
		alpha1=1._sp/ttr
		!k=2.e-3_sp                      ! changes the width
		!dsm_by_dz_z_eq_zc=-1.6e-6_sp    ! 
		!b=1.e-6_sp                      ! range from 0 to 5e-6, default 1e-6
		!del_gamma_mac=5e-4_sp
		!del_c_s=0._sp
		!del_c_t=0.5_sp                  ! changing alters height
		!epsilon_therm=3e-7_sp			! changing this alters height too
		klarge=epsilon_therm*cp/lv
		! equation 32:
		n_bar_mac=grav*(alpha1*del_gamma_mac+beta1*(dsm_by_dz_z_eq_zc+b))
		n_bar_mac=sqrt(n_bar_mac)
		! equation 33:
		z_bar=(alpha1*del_c_t+beta1*del_c_s)/ &
    		(alpha1*del_gamma_mac+beta1*(dsm_by_dz_z_eq_zc+b))
    	!z_bar=0._sp
    	k1=0.5_sp*(alpha1*epsilon_therm-beta1*klarge)/ &
    		(alpha1*del_gamma_mac+beta1*(dsm_by_dz_z_eq_zc+b))
    	cell_size=3._sp/4._sp/k1*(1._sp+sqrt(1._sp+16._sp/3._sp*k1*z_bar))
		therm_init=.false.
	endif

	wcen=0._sp	
	! order of nested do loop is 2nd dimension first
	do i=0,ip
		do j=-1,kp
			! equation 39 of ZZRA:
! 			zz(j,i)=k*n_bar_mac* &
! 				sqrt(2._sp*(zn(j)-zc)*(z_bar+(zn(j)-zc)/2._sp-k1/3._sp*(zn(j)-zc)**2))

			! equation 41 of ZZRA:
! 			xx(j,i)=1._sp/k**2._sp*cos(k*xn(i))

			! equation 42 of ZZRA:
! 			phi=zz*xx
			! u and w winds
				zc=zn(2) !-z_offset
				test3=((z(j)-z_offset)-zc)
				test1=small1+2._sp*test3*(z_bar+test3/2._sp- &
							 k1/3._sp*test3**2._sp)
				test1=-n_bar_mac/k*(z_bar+test3-k1*test3**2)/sqrt(test1) &
						*sin(k*xn(i+1))
									
				zc=zn(2) !-z_offset
				test3=((zn(j+1)-z_offset)-zc)
				test2=small1+2._sp*test3*(z_bar+test3/2._sp- &
							 k1/3._sp*test3**2._sp)

				test2=n_bar_mac* sqrt(test2) *cos(k*x(i))

				u(j,i)=real(test1)!+imag(test1)
! 				if(j.eq.122) u(j,i)=u(j,i)-imag(test1)
				w(j,i)=real(test2)!+imag(test2)
! 				if(j.eq.122) w(j,i)=w(j,i)+imag(test2)
! 				u(j,i)=u(j,i)-imag(test1)
! 				w(j,i)=w(j,i)-imag(test2)
				! w on centred points
				zc=z(2) !-z_offset
				test3=((z(j+1)-z_offset)-zc)
				test2=small1+2._sp*test3*(z_bar+test3/2._sp- &
							 k1/3._sp*test3**2._sp)
				test2=n_bar_mac* sqrt(test2) *cos(k*(xn(i)))
				wcen(j,i)=real(test2)
				
				
! 				if(abs((z(j)-zc)-(3._sp/2._sp/k1)).lt.dz) then 
! 					u(j,i)=0._sp
! 					w(j,i)=0._sp
! 				endif

		enddo
	enddo	
	! calculate u via finite difference
	wcen(0,:)=wcen(1,:)



	u=0._sp
	w=0._sp
	! centred differences:
	u(1:kp,0)=(wcen(1:kp,1)-wcen(0:kp-1,1))*dx/dz
	do i=1,ip
		u(1:kp,i)=u(1:kp,i-1)-(wcen(1:kp,i)-wcen(0:kp-1,i))*dx/dz
	enddo
	w(1:kp,1:ip)=(wcen(1:kp,1:ip))
    u(1,:)=0._sp
    w(1,:)=0._sp
	! halos
	u(1:kp,-o_halo+1:0)=u(1:kp,ip-o_halo+1:ip)
 	u(1:kp,ip+1:ip+o_halo)=u(1:kp,1:o_halo)
	w(1:kp,-o_halo+1:0)=w(1:kp,ip-o_halo+1:ip)
 	w(1:kp,ip+1:ip+o_halo)=w(1:kp,1:o_halo)
	do j=-o_halo+1,0
		u(j,:)=u(1,:)
		w(j,:)=w(1,:)
	enddo

	do j=1,kp+o_halo
		if(z(j+1) >= z_offset) exit

		w(j,:)=0._sp
		u(j,:)=0._sp
	enddo

	u=u*w_peak/maxval(w(1:kp,1:ip))
	w=w*w_peak/maxval(w(1:kp,1:ip))
	
	w(-o_halo+1:0,-o_halo+1:ip+o_halo)=0._sp
	u(-o_halo+1:0,-o_halo+1:ip+o_halo)=0._sp
	

    end subroutine thermal_2d

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>sets the vertical wind speed
	!>@param[in] time
	!>@param[in] ip
	!>@param[in] kp
	!>@param[in] o_halo
	!>@param[in] k
	!>@param[in] dsm_by_dz_z_eq_zc
	!>@param[in] b
	!>@param[in] del_gamma_mac
	!>@param[in] del_c_s
	!>@param[in] del_c_t
	!>@param[in] epsilon_therm
	!>@param[in] x,xn,z,zn,dx,dz
	!>@param[inout]  u, w,
	!>@param[in]  w_peak, z_offset
	!>@param[inout] therm_init
    subroutine fd_thermal_2d(time,ip,kp,o_halo,k,dsm_by_dz_z_eq_zc,b,del_gamma_mac,del_c_s,del_c_t, &
    						epsilon_therm,x,xn,z,zn,dx,dz,u,w,w_peak,z_offset, therm_init)

    use nrtype

    implicit none
    integer(i4b), intent(in) :: ip,kp, o_halo
    real(sp), intent(in) :: time
    real(sp), intent(inout) :: k,dsm_by_dz_z_eq_zc,b,del_gamma_mac,del_c_s,del_c_t, &
    						epsilon_therm
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: &
    					u,w
    real(sp), dimension(-o_halo+1:ip+o_halo), intent(in) :: x,xn
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z,zn
    real(sp), intent(in) :: dx,dz
    real(sp), intent(in) :: w_peak, z_offset
    logical, intent(inout) :: therm_init
    
    ! local variable
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo) :: xx,zz
    real(sp), dimension(-o_halo+1:kp+o_halo+1,-o_halo+1:ip+o_halo+1) :: phi
	integer(i4b) :: i,j
	complex(sp) :: test1,test2 

    ! calculate the thermal properties
	if(therm_init) then
		zc=zn(1)-z_offset !0._sp !500._sp
		alpha1=1._sp/ttr
		klarge=epsilon_therm*cp/lv
		! equation 32:
		n_bar_mac=grav*(alpha1*del_gamma_mac+beta1*(dsm_by_dz_z_eq_zc+b))
		n_bar_mac=sqrt(n_bar_mac)
		! equation 33:
		z_bar=(alpha1*del_c_t+beta1*del_c_s)/ &
    		(alpha1*del_gamma_mac+beta1*(dsm_by_dz_z_eq_zc+b))
    	!z_bar=0._sp
    	k1=0.5_sp*(alpha1*epsilon_therm-beta1*klarge)/ &
    		(alpha1*del_gamma_mac+beta1*(dsm_by_dz_z_eq_zc+b))
    	cell_size=3._sp/4._sp/k1*(1._sp+sqrt(1._sp+16._sp/3._sp*k1*z_bar))
		therm_init=.false.
	endif
	
	! order of nested do loop is 2nd dimension first
	do i=0,ip
		do j=1,kp
			! equation 39 of ZZRA:
			zz(j,i)=k*n_bar_mac* &
				sqrt(2._sp*((zn(j)-z_offset)-zc)* &
				(z_bar+((zn(j)-z_offset)-zc)/2._sp-k1/3._sp*((zn(j)-z_offset)-zc)**2))

			! equation 41 of ZZRA:
			xx(j,i)=1._sp/k**2._sp*cos(k*xn(i))

			! equation 42 of ZZRA:
			phi(j,i)=zz(j,i)*xx(j,i)
			if((zn(j)-z_offset).ge.(zc+cell_size-0._sp*dz)) then
				phi(j,i)=0._sp
			endif
		enddo
	enddo
	
	! some halo stuff
	do i=1,ip
		phi(-o_halo+1:0,i)=phi(1,i) ! bottom
		phi(kp+1:kp+o_halo+1,i)=phi(kp,i) ! top
	enddo
	phi(1:kp,-o_halo+1:0)=phi(1:kp,ip-o_halo+1:ip) !left
	phi(1:kp,ip+1:ip+o_halo+1)=phi(1:kp,1:o_halo+1) ! right
	
	! finite difference:
	u(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo)= &
		-(phi(-o_halo+2:kp+o_halo+1,-o_halo+1:ip+o_halo)- &
		  phi(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo)) / dx
	w(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo)= &
	 	(phi(-o_halo+1:kp+o_halo,-o_halo+2:ip+o_halo+1)- &
	 	 phi(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo)) / dz
	
! 			u and w winds
! 			if(z(j).ge.(zc+cell_size-0._sp*dz)) then
! 				u(j,i)=0._sp
! 			else
! 				
! 				test1=small1+2._sp*(z(j)-zc)*(z_bar+(z(j)-zc)/2._sp- &
! 							 k1/3._sp*(z(j)-zc)**2)
! 				test1=real(sqrt(test1))
! 						 
! 						 
! 				u(j,i)=-n_bar_mac/k*(z_bar+(z(j)-zc)-k1*(z(j)-zc)**2) / &
! 					real(test1,sp) &
! 						*cos(k*xn(i))
! 			endif
! 			
! 			if(zn(j).ge.(zc+cell_size-0._sp*dz)) then
! 				w(j,i)=0._sp !w(j-1,i) !*exp(1.e-3_sp*((zc+cell_size-3._sp*dz)-zn(j)))
! 			else
! 				test2=small1+2._sp*(zn(j)-zc)*(z_bar+(zn(j)-zc)/2._sp- &
! 							 k1/3._sp*(zn(j)-zc)**2)
! 				test2=real(sqrt(test2))
! 				
! 				w(j,i)=-n_bar_mac* &
! 					real(test2,sp) &
! 						*sin(k*x(i))
! 			endif
! 
! 		enddo
! 	enddo	
	! halos
	u(1:kp,-o_halo+1:0)=u(1:kp,ip-o_halo+1:ip)
 	u(1:kp,ip+1:ip+o_halo)=u(1:kp,1:o_halo)

	u(1,:)=0._sp
	w(1,:)=0._sp
	do j=-o_halo+1,0
		u(j,:)=u(1,:)
		w(j,:)=w(1,:)
	enddo
	do j=kp+1,kp+o_halo
		u(j,:)=u(j,:)
		w(j,:)=w(j,:)
	enddo

	u=u*w_peak/maxval(w)
	w=w*w_peak/maxval(w)
	
	do j=1,kp+o_halo
		if(zn(j) >= z_offset) exit

		w(j,:)=0._sp
		u(j,:)=0._sp
	enddo

	
    end subroutine fd_thermal_2d


  
      end module thermal
    
    
