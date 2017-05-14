	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>advection code for the thermal cloud model
    module advection_2d
    use nrtype
    
    private
    public :: mpdata, first_order_upstream_2d
    
	contains
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Simple first order upstream scheme                                                 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>perform 1 time-step of 2-d first order upstream method 
	!>@param[in] dt
	!>@param[in] dx,dy
	!>@param[in] ip,kp,halo
	!>@param[in] u
	!>@param[in] w
	!>@param[inout] psi
	subroutine first_order_upstream_2d(dt,dx,dy,ip,kp,halo,u,w,psi)
	use nrtype
	implicit none
	real(sp), intent(in) :: dt, dx, dy
	integer(i4b), intent(in) :: ip, kp, halo
	real(sp), dimension(-halo+1:kp+halo,-halo+1:ip+halo), intent(in) :: u, w
	real(sp), dimension(-halo+1:kp+halo,-halo+1:ip+halo), intent(inout) :: psi
	
	! locals
	real(sp), dimension(kp,ip) :: fx_r, fx_l, fy_r, fy_l
	
	
	fx_r=( (u(1:kp,1:ip)+abs(u(1:kp,1:ip)))*psi(1:kp,1:ip)+ &
		(u(1:kp,1:ip)-abs(u(1:kp,1:ip)))*psi(1:kp,2:ip+1) )*dt/(2._sp*dx)
		
	fx_l=( (u(1:kp,0:ip-1)+abs(u(1:kp,0:ip-1)))*psi(1:kp,0:ip-1)+ &
		(u(1:kp,0:ip-1)-abs(u(1:kp,0:ip-1)))*psi(1:kp,1:ip) )*dt/(2._sp*dx)
		
	fy_r=( (w(1:kp,1:ip)+abs(w(1:kp,1:ip)))*psi(1:kp,1:ip)+ &
		(w(1:kp,1:ip)-abs(w(1:kp,1:ip)))*psi(2:kp+1,1:ip) )*dt/(2._sp*dy)
		
	fy_l=( (w(0:kp-1,1:ip)+abs(w(0:kp-1,1:ip)))*psi(0:kp-1,1:ip)+ &
		(w(0:kp-1,1:ip)-abs(w(0:kp-1,1:ip)))*psi(1:kp,1:ip) )*dt/(2._sp*dy)
		
	
	psi(1:kp,1:ip)=psi(1:kp,1:ip)-(fx_r-fx_l)-(fy_r-fy_l)
	end subroutine first_order_upstream_2d
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! multi-dimensional advection using the smolarkiewicz scheme
	! advect using 1st order upstream
	! then re-advect using 1st order upstream with antidiffusive velocities to 
	!   correct diffusiveness of the 1st order upstream and iterate
	! nft option is also coded
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine mpdata(kord,ip,kp,o_halo,x,z,dx,dz,dt,u,w,q_k,monotone)
	use nrtype
	implicit none
	integer(i4b), intent(in) :: kord,ip, kp, o_halo
	real(sp), dimension(-o_halo+1:ip+o_halo), intent(in) :: x
	real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z
	real(sp), intent(in) :: dx,dz,dt
	real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(in) :: u,w
	real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: q_k

	real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo) :: ut,vt,wt,q_kp1,q_k_old, &
											 ut_sav,vt_sav,wt_sav
	real(sp) :: fip,fim,small=1e-15_sp, &
			u_j_bar1, u_div1, u_j_bar2, u_div2, u_j_bar3, u_div3, &
			psi_i_max, psi_i_min, psi_ip_max,psi_ip_min, beta_i_down, beta_i_up, &
			beta_ip_down, beta_ip_up
	integer(i4b) i,j,k,it
	logical :: monotone
	
	if(sum(q_k).lt.small) return
	
	! zero arrays
	ut=0._sp
	wt=0._sp
	ut_sav=0._sp
	wt_sav=0._sp
	q_kp1=0._sp
	q_k_old=0._sp
	! save old data
	q_k_old(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo)=q_k

	do it=1,kord
		if(it.eq.1) then
		  ut(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo)=u
		  wt(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo)=w
		else 
   		  do i=-o_halo+2,ip+o_halo-1
			 do k=-o_halo+2,kp+o_halo-1
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! now calculate the anti-diffusive velocities
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				u_j_bar1=0._sp
				u_div1=0._sp
	
				! calculate the part of the anti-diffusive velocities associated 
				! with cross derivatives 
				! equation 13 page 330 of smolark (1984) 
				! journal of computational physics
	
				u_j_bar1=u_j_bar1+0.5_sp*dt*ut(k,i)* &
						0.25_sp*(wt(k,i+1)+wt(k,i)+wt(k-1,i+1)+wt(k-1,i)) &
				  * (q_k_old(k+1,i+1)+q_k_old(k+1,i)-q_k_old(k-1,i+1)-q_k_old(k-1,i))/ &
				 (q_k_old(k+1,i+1)+q_k_old(k+1,i)+q_k_old(k-1,i+1)+q_k_old(k-1,i)+small )/dz
				! for divergent flow: eq 38 smolarkiewicz 1984 
				u_div1=u_div1+(wt(k,i)+wt(k,i+1)-wt(k-1,i)-wt(k-1,i+1))/dz
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				u_j_bar3=0._sp
				u_div3=0._sp
				u_j_bar3=u_j_bar3+0.5_sp*dt*wt(k,i)* &
						0.25_sp*(ut(k+1,i)+ut(k,i)+ut(k+1,i-1)+ut(k,i-1)) &
				  * (q_k_old(k+1,i+1)+q_k_old(k,i+1)-q_k_old(k+1,i-1)-q_k_old(k,i-1))/ &
				  (q_k_old(k+1,i+1)+q_k_old(k,i+1)+q_k_old(k+1,i-1)+q_k_old(k,i-1)+small )/dx
				! for divergent flow: eq 38 smolarkiewicz 1984
				u_div3=u_div3+(ut(k,i)+ut(k+1,i)-ut(k,i-1)-ut(k+1,i-1))/dx
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				!!!!
				!!!!
				ut_sav(k,i)=(abs(ut(k,i))*dx-dt*ut(k,i)*ut(k,i) ) * &
					(q_k_old(k,i+1)-q_k_old(k,i) ) / &
					(q_k_old(k,i+1)+q_k_old(k,i)+small) /dx - u_j_bar1
				! divergent flow: eq 38 smolarkiewicz 1984
! 				ut_sav(k,i)=ut_sav(k,i) - 0.25_sp*dt*ut(k,i)* &
! 					( (ut(k,i+1)-ut(k,i-1))/dx-u_div1 )
				!!!!
				wt_sav(k,i)=(abs(wt(k,i))*dz-dt*wt(k,i)*wt(k,i) ) * &
					(q_k_old(k+1,i)-q_k_old(k,i) ) / &
					(q_k_old(k+1,i)+q_k_old(k,i)+small) /dz - u_j_bar3
				! divergent flow: eq 38 smolarkiewicz 1984
! 				wt_sav(k,i)=wt_sav(k,i) - 0.25_sp*dt*wt(k,i)* &
! 					( (wt(k+1,i)-wt(k-1,i))/dz-u_div3 )
				!!!!          
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			  enddo
		  enddo
		  ut=ut_sav
		  wt=wt_sav
		endif


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Non oscillatory forward in time (NFT) flux limiter -                           !
        ! Smolarkiewicz and Grabowski (1990)                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(it.ge.2.and.monotone) then
		do i=-o_halo+3,ip+o_halo-2 
			do k=-o_halo+3,kp+o_halo-2
	
	
			
		! x direction - note: should the last q in the max/min be q+1?
		psi_i_max=max(q_k(k,i-1),q_k(k,i),q_k(k,i+1), &
					q_k_old(k,i-1),q_k_old(k,i),q_k_old(k,i-1))
					
		psi_i_min=min(q_k(k,i-1),q_k(k,i),q_k(k,i+1), &
					q_k_old(k,i-1),q_k_old(k,i),q_k_old(k,i-1))
					
		psi_ip_max=max(q_k(k,i),q_k(k,i+1),q_k(k,i+2), &
					q_k_old(k,i),q_k_old(k,i+1),q_k_old(k,i))
					
		psi_ip_min=min(q_k(k,i),q_k(k,i+1),q_k(k,i+2), &
					q_k_old(k,i),q_k_old(k,i+1),q_k_old(k,i))





		! for 3-d add another term to the betas
		beta_i_up=(psi_i_max-q_k_old(k,i)) / &
			(dt/dx*(max(ut(k,i-1),0._sp)*q_k_old(k,i-1)-min(ut(k,i),0._sp)*q_k_old(k,i+1))+ &
			dt/dz*(max(wt(k-1,i),0._sp)*q_k_old(k-1,i)-min(wt(k,i),0._sp)*q_k_old(k+1,i)) &
				+small)
		beta_i_down=(q_k_old(k,i)-psi_i_min) / &
			(dt/dx*(max(ut(k,i),0._sp)*q_k_old(k,i)-min(ut(k,i-1),0._sp)*q_k_old(k,i)) + &
			dt/dz*(max(wt(k,i),0._sp)*q_k_old(k,i)-min(wt(k-1,i),0._sp)*q_k_old(k,i)) &
				+small)
		beta_ip_up=(psi_ip_max-q_k_old(k,i+1)) / &
			(dt/dx*(max(ut(k,i),0._sp)*q_k_old(k,i)-min(ut(k,i+1),0._sp)*q_k_old(k,i+2)) + &
			dt/dz*(max(wt(k,i),0._sp)*q_k_old(k,i)-min(wt(k+1,i),0._sp)*q_k_old(k+2,i)) &
				+small)
		beta_ip_down=(q_k_old(k,i+1)-psi_ip_min) / &
			(dt/dx*(max(ut(k,i+1),0._sp)*q_k_old(k,i+1)-min(ut(k,i),0._sp)*q_k_old(k,i+1)) + &
			dt/dz*(max(wt(k+1,i),0._sp)*q_k_old(k+1,i)-min(wt(k,i),0._sp)*q_k_old(k+1,i)) &
				+small)
				
 		ut_sav(k,i)=min(1._sp,beta_i_down,beta_ip_up)*max(ut(k,i),0._sp) + &
 				min(1._sp,beta_i_up,beta_ip_down)*min(ut(k,i),0._sp)

 		
		! z direction - note: should the last q in the max/min be q+1?
		psi_i_max=max(q_k(k-1,i),q_k(k,i),q_k(k+1,i), &
					q_k_old(k-1,i),q_k_old(k,i),q_k_old(k-1,i))
					
		psi_i_min=min(q_k(k-1,i),q_k(k,i),q_k(k+1,i), &
					q_k_old(k-1,i),q_k_old(k,i),q_k_old(k-1,i))
					
		psi_ip_max=max(q_k(k,i),q_k(k+1,i),q_k(k+2,i), &
					q_k_old(k,i),q_k_old(k+1,i),q_k_old(k,i))
					
		psi_ip_min=min(q_k(k,i),q_k(k+1,i),q_k(k+2,i), &
					q_k_old(k,i),q_k_old(k+1,i),q_k_old(k,i))




		! for 3-d add another term to the betas
		beta_i_up=(psi_i_max-q_k_old(k,i)) / &
			(dt/dx*(max(ut(k,i-1),0._sp)*q_k_old(k,i-1)-min(ut(k,i),0._sp)*q_k_old(k,i+1)) +&
			dt/dz*(max(wt(k-1,i),0._sp)*q_k_old(k-1,i)-min(wt(k,i),0._sp)*q_k_old(k+1,i)) &
				+small)
		beta_i_down=(q_k_old(k,i)-psi_i_min) / &
			(dt/dx*(max(ut(k,i),0._sp)*q_k_old(k,i)-min(ut(k,i-1),0._sp)*q_k_old(k,i)) + &
			dt/dz*(max(wt(k,i),0._sp)*q_k_old(k,i)-min(wt(k-1,i),0._sp)*q_k_old(k,i)) &
				+small)
		beta_ip_up=(psi_ip_max-q_k_old(k+1,i)) / &
			(dt/dx*(max(ut(k,i),0._sp)*q_k_old(k,i)-min(ut(k,i+1),0._sp)*q_k_old(k,i+2)) + &
			dt/dz*(max(wt(k,i),0._sp)*q_k_old(k,i)-min(wt(k+1,i),0._sp)*q_k_old(k+2,i)) &
				+small)
		beta_ip_down=(q_k_old(k+1,i)-psi_ip_min) / &
			(dt/dx*(max(ut(k,i+1),0._sp)*q_k_old(k,i+1)-min(ut(k,i),0._sp)*q_k_old(k,i+1)) +&
			dt/dz*(max(wt(k+1,i),0._sp)*q_k_old(k+1,i)-min(wt(k,i),0._sp)*q_k_old(k+1,i)) &
				+small)
				
 		wt_sav(k,i)=min(1._sp,beta_i_down,beta_ip_up)*max(wt(k,i),0._sp) + &
 				min(1._sp,beta_i_up,beta_ip_down)*min(wt(k,i),0._sp)
 			enddo
 		enddo
 		
 		ut=ut_sav
 		wt=wt_sav
 		endif
 		
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! upstream scheme using ut, wt                                                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do i=-o_halo+2,ip+o_halo-1
			do k=-o_halo+2,kp+o_halo-1
			  fip = 0._sp
			  ! i direction
			  fip = fip +( ( ut(k,i)+abs(ut(k,i)) )*q_k_old(k,i) + &
						 ( ut(k,i)-abs(ut(k,i)) )*q_k_old(k,i+1) )*dt/(2._sp*dx)
			  ! k direction
			  fip = fip +( ( wt(k,i)+abs(wt(k,i)) )*q_k_old(k,i) + &
						 ( wt(k,i)-abs(wt(k,i)) )*q_k_old(k+1,i) )*dt/(2._sp*dz)
	
			  fim = 0._sp
			  ! i direction
			  fim = fim +( ( ut(k,i-1)+abs(ut(k,i-1)) )*q_k_old(k,i-1) + &
						 ( ut(k,i-1)-abs(ut(k,i-1)) )*q_k_old(k,i) )*dt/(2._sp*dx)
			  ! k direction
			  fim = fim +( ( wt(k-1,i)+abs(wt(k-1,i)) )*q_k_old(k-1,i) + &
						 ( wt(k-1,i)-abs(wt(k-1,i)) )*q_k_old(k,i) )*dt/(2._sp*dz) 

			  q_kp1(k,i)=q_k_old(k,i)-(fip-fim)
			enddo
		enddo
		q_k_old(-o_halo+2:kp+o_halo-1,-o_halo+2:ip+o_halo-1)= &
		       q_kp1(-o_halo+2:kp+o_halo-1,-o_halo+2:ip+o_halo-1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	enddo
	q_k(1:kp,1:ip)=q_k_old(1:kp,1:ip)


	end subroutine mpdata
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end module advection_2d
