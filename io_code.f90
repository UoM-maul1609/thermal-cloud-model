	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>io routines
    module io_module
    use nrtype
    private
    public :: output_2d
    contains
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>output 1 time-step of model
	!>@param[in] time in seconds
	!>@param[in] nq number of q fields
	!>@param[in] nprec number of precipitation fields
	!>@param[in] ip number of horizontal levels
	!>@param[in] kp number of vertical levels
	!>@param[in] q_name: name of q-variables
	!>@param[in] q, precip, theta, pressure, x,xn,z,zn, temperature,u,w
	!>@param[inout] new_file
    subroutine output_2d(time,nq,nprec,ip,kp,q_name, &
                        q,precip,theta,p,x,xn,z,zn,t,u,w,tke,new_file)

    use nrtype
    use netcdf
    use variables, only : io1, outputfile

    implicit none
    real(sp), intent(in) :: time
    integer(i4b), intent(in) :: nq,nprec,kp,ip
    character(len=20), dimension(nq) :: q_name
    real(sp), dimension(kp,ip,nq), intent(in) :: q
    real(sp), dimension(kp,ip,nprec), intent(in) :: precip
    real(sp), dimension(kp,ip), intent(in) :: theta, p, t, u, w, tke
    real(sp), dimension(ip), intent(in) :: x,xn
    real(sp), dimension(kp), intent(in) :: z,zn
    logical, intent(inout) :: new_file


    ! output to netcdf file
    if(new_file) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! open / create the netcdf file                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call check( nf90_create(outputfile, NF90_CLOBBER, io1%ncid) )

        ! define dimensions (netcdf hands back a handle)
        call check( nf90_def_dim(io1%ncid, "times", NF90_UNLIMITED, io1%x_dimid) )
        call check( nf90_def_dim(io1%ncid, "nq", nq, io1%nq_dimid) )
        call check( nf90_def_dim(io1%ncid, "nprec", nprec, io1%nprec_dimid) )
        call check( nf90_def_dim(io1%ncid, "ip", ip, io1%i_dimid) )
        call check( nf90_def_dim(io1%ncid, "kp", kp, io1%k_dimid) )
        call check( nf90_def_dim(io1%ncid, "l_q_names", 20, io1%lq_dimid) )


        ! close the file, freeing up any internal netCDF resources
        ! associated with the file, and flush any buffers
        call check( nf90_close(io1%ncid) )


        ! now define some variables, units, etc
        call check( nf90_open(outputfile, NF90_WRITE, io1%ncid) )
        ! define mode
        call check( nf90_redef(io1%ncid) )


        ! define variable: q_names
        call check( nf90_def_var(io1%ncid, "q_names", NF90_CHAR, &
            (/io1%lq_dimid,io1%nq_dimid/), io1%varid) )


        ! define variable: time
        call check( nf90_def_var(io1%ncid, "time", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "time", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "seconds") )

        ! define variable: q
        call check( nf90_def_var(io1%ncid, "q", NF90_DOUBLE, &
            (/io1%nq_dimid, io1%k_dimid, io1%i_dimid, io1%x_dimid/), io1%varid) )
!        call check( nf90_def_var(io1%ncid, "q", NF90_DOUBLE, &
!            (/io1%k_dimid, io1%i_dimid,io1%nq_dimid,  io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "q", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "kg or number per kg") )

        ! define variable: precip
        call check( nf90_def_var(io1%ncid, "precip", NF90_DOUBLE, &
            (/io1%nprec_dimid, io1%k_dimid, io1%i_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "precip", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "mm hr-1") )



        ! define variable: theta
        call check( nf90_def_var(io1%ncid, "theta", NF90_DOUBLE, &
            (/io1%k_dimid, io1%i_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "theta", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "K") )


        ! define variable: p
        call check( nf90_def_var(io1%ncid, "p", NF90_DOUBLE, &
            (/io1%k_dimid, io1%i_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "p", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "Pa") )


        ! define variable: x
        call check( nf90_def_var(io1%ncid, "x", NF90_DOUBLE, &
            (/io1%i_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "x", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "m") )

        ! define variable: xn
        call check( nf90_def_var(io1%ncid, "xn", NF90_DOUBLE, &
            (/io1%i_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "xn", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "m") )


        ! define variable: z
        call check( nf90_def_var(io1%ncid, "z", NF90_DOUBLE, &
            (/io1%k_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "z", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "m") )

        ! define variable: zn
        call check( nf90_def_var(io1%ncid, "zn", NF90_DOUBLE, &
            (/io1%k_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "zn", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "m") )


        ! define variable: t
        call check( nf90_def_var(io1%ncid, "t", NF90_DOUBLE, &
            (/io1%k_dimid, io1%i_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "t", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "K") )

        ! define variable: u
        call check( nf90_def_var(io1%ncid, "u", NF90_DOUBLE, &
            (/io1%k_dimid, io1%i_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "u", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "ms-1") )

        ! define variable: w
        call check( nf90_def_var(io1%ncid, "w", NF90_DOUBLE, &
            (/io1%k_dimid, io1%i_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "w", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "ms-1") )

	!define variable tke
	call check( nf90_def_var(io1%ncid, "tke", NF90_DOUBLE, &
		(/io1%k_dimid, io1%i_dimid, io1%x_dimid/), io1%varid) )
	! get id to a_dimid
	call check( nf90_inq_varid(io1%ncid, "tke", io1%a_dimid) )
	! units
	call check (nf90_put_att(io1%ncid, io1%a_dimid, &
			"units", "m2 s-2") )

        call check( nf90_enddef(io1%ncid) )
        call check( nf90_close(io1%ncid) )

        new_file=.false.
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! write x,xn,z,zn data to file                                           !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call check( nf90_open(outputfile, NF90_WRITE, io1%ncid) )
		call check( nf90_inq_varid(io1%ncid, "q_names", io1%varid ) )
		call check( nf90_put_var(io1%ncid, io1%varid, q_name, start = (/1,1/)))

		call check( nf90_inq_varid(io1%ncid, "x", io1%varid ) )
		call check( nf90_put_var(io1%ncid, io1%varid, x, start = (/1/)))
        
		call check( nf90_inq_varid(io1%ncid, "xn", io1%varid ) )
		call check( nf90_put_var(io1%ncid, io1%varid, xn, start = (/1/)))
        
		call check( nf90_inq_varid(io1%ncid, "z", io1%varid ) )
		call check( nf90_put_var(io1%ncid, io1%varid, z, start = (/1/)))
        
		call check( nf90_inq_varid(io1%ncid, "zn", io1%varid ) )
		call check( nf90_put_var(io1%ncid, io1%varid, zn, start = (/1/)))
        
		call check( nf90_close(io1%ncid) )
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write data to file                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call check( nf90_open(outputfile, NF90_WRITE, io1%ncid) )
    ! write variable: time
    call check( nf90_inq_varid(io1%ncid, "time", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, time, &
                start = (/io1%icur/)))
! 

    ! write variable: q
    call check( nf90_inq_varid(io1%ncid, "q", io1%varid ) )
    ! had to reshape twice
    call check( nf90_put_var(io1%ncid, io1%varid, &
        reshape(reshape(q,[nq,ip,kp],order=[3,2,1]),[nq,kp,ip],order=[1,3,2] ), &
                start = (/1,1,1,io1%icur/)))
!    call check( nf90_put_var(io1%ncid, io1%varid, q, &
!                start = (/1,1,1,io1%icur/)))

    ! write variable: precip
    call check( nf90_inq_varid(io1%ncid, "precip", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        reshape(reshape(precip,[nprec,ip,kp],order=[3,2,1]),[nprec,kp,ip],order=[1,3,2] ), &
                start = (/1,1,1,io1%icur/)))
!    call check( nf90_put_var(io1%ncid, io1%varid, precip, &
!                start = (/1,1,1,io1%icur/)))

    ! write variable: theta
    call check( nf90_inq_varid(io1%ncid, "theta", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, theta, &
                start = (/1,1,io1%icur/)))

    ! write variable: p
    call check( nf90_inq_varid(io1%ncid, "p", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, p, &
                start = (/1,1,io1%icur/)))

    ! write variable: u
    call check( nf90_inq_varid(io1%ncid, "u", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, u, &
                start = (/1,1,io1%icur/)))

    ! write variable: w
    call check( nf90_inq_varid(io1%ncid, "w", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, w, &
                start = (/1,1,io1%icur/)))

    ! write variable: tke
    call check( nf90_inq_varid(io1%ncid, "tke", io1%varid) )
    call check( nf90_put_var(io1%ncid, io1%varid, tke, &
		start = (/1,1,io1%icur/)))

    ! write variable: t
    call check( nf90_inq_varid(io1%ncid, "t", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, t, &
                start = (/1,1,io1%icur/)))

    call check( nf90_close(io1%ncid) )
! 
! 
    io1%icur=io1%icur+1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine output_2d







    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! HELPER ROUTINE                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check(status)
    use netcdf
    use nrtype
    integer(I4B), intent ( in) :: status

    if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped"
    end if
    end subroutine check
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end module io_module
