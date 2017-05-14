	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>variables for the thermal cloud model
    module variables
    use nrtype
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>variables and types for the thermal cloud model

    implicit none
    
		!>@brief
		!>main model prognostic variables
        type grid
            ! variables for grid
            integer(i4b) :: n_levels
            real(sp) :: dx,dz, dt
            real(sp), dimension(:,:,:), allocatable :: q, qold, precip
            real(sp), dimension(:,:), allocatable :: theta, p, rho, t, u, w
            real(sp), dimension(:), allocatable :: x,z,xn,zn
        end type grid

		!>@brief
		!>variables for sounding input
        type sounding
            ! variables for grid
            integer(i4b) :: n_levels
            real(sp), dimension(:,:), allocatable :: q
            real(sp), dimension(:), allocatable :: theta, p, z, rh
        end type sounding



        
		!>@brief
		!>variables for NetCDF file output
        type io
            ! variables for io
            integer(i4b) :: ncid, varid, x_dimid, y_dimid, z_dimid, &
                            dimids(2), a_dimid, xx_dimid, yy_dimid, &
                            zz_dimid, i_dimid, j_dimid, k_dimid, nq_dimid, nprec_dimid
            integer(i4b) :: icur=1
            logical :: new_file=.true.
        end type io


        ! declare a grid type
        type(grid) :: grid1
        ! declare a sounding type
        type(sounding) :: sounding1
        ! declare an io type
        type(io) :: io1

        ! constants
        integer(i4b), parameter :: nq = 11, nlevels_r=1000
        integer(i4b), parameter :: qv=1,qc=2,qr=3,nqc=4,nqr=5,qi=6,qs=7,qg=8, &
                                   nqi=9,nqs=10,nqg=11
        ! the type of q-variable. 0 vapour, 1 mass, 2 number conc.
        integer(i4b), dimension(nq) :: q_type=(/0,1,1,1,1,1,2,2,2,2,2/)
        ! the type of q-variable. 0 vapour, 1 mass, 2 number conc.
        logical, dimension(nq) :: q_init=(/.true.,.false.,.false.,.false., &
                                            .false.,.false.,.false.,.false., &
                                            .false.,.false.,.false./)
        logical :: micro_init=.true., adiabatic_prof=.false.
        logical :: monotone=.true.,microphysics_flag=.true.,theta_flag=.false., &
        			hm_flag=.true.
		integer(i4b) :: advection_scheme=0

        ! variables for model
        real(sp), dimension(nq,nlevels_r) :: q_read
        real(sp), dimension(nlevels_r) :: theta_read,rh_read, &
                  z_read
        real(sp) :: dx, dz,dt, runtime, psurf, theta_surf,tsurf, t_cbase, t_ctop, t_thresh, &
        			t_thresh2, w_peak, w_cb, theta_q_sat,t1old, p111, num_ice, mass_ice
        integer(i4b) :: ip,kp, n_levels_s, ord, o_halo,halo, updraft_type
        character (len=200) :: outputfile='output'
        
        
        ! variables for thermal properties
        real(sp) :: k, dsm_by_dz_z_eq_zc, b, del_gamma_mac, &
        			del_c_s, del_c_t, epsilon_therm
        logical :: therm_init=.true.
    end module variables




    module constants
        use nrtype
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>constants for the thermal cloud model

        implicit none
        real(sp), parameter :: ra=287.0_sp, cp=1005.0_sp, grav=9.81_sp, &
        						rv=461._sp, eps1=ra/rv, lv=2.5e6_sp, ttr=273.15_sp
    end module constants







