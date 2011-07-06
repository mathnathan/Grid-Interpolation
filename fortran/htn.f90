module Interpolate
implicit none

contains

  ! Find the location of the largest number
  ! in the array that is less than value
  integer function find( value, array, array_size )
    implicit none

    integer, intent(in) :: array_size
    real*4, intent(in) :: value, array( array_size )
    real*4 :: tmp(array_size), minimum
    integer :: indx, i

    minimum = 1.e9
    tmp = abs( array - value )

    do i = 1, array_size 
      if( tmp(i) < minimum ) then
        minimum = tmp(i)
        indx = i
      endif
    enddo
      
    if( array( indx ) > value )  indx = indx -1
    find = indx

  end function find

  ! Find the location of the largest number
  ! in the array that is less than value
  integer function find_depth( value, array, array_size )
    implicit none

    integer, intent(in) :: array_size
    real*4, intent(in) :: value, array(array_size)
    real*4 :: tmp(array_size), minimum
    integer :: indx, i

    minimum = 1.e9
    tmp = abs( array - value )

    do i = 1, array_size 
      if( tmp(i) < minimum ) then
        minimum = tmp(i)
        indx = i
      endif
    enddo
      
    if( array( indx ) > value .or. indx == 1 ) then
      indx = indx + 1
    endif
    find_depth = indx

  end function find_depth

  ! This is the trilinear interpolation algorithm. Converts
  ! the HYCOM grid data to the NCOM grid data.
  subroutine hycom2ncom( ncom_lon, ncom_lat, ncom_z, ncom_r, &
  n, m, l, &
  hycom_lon, hycom_lat, hycom_z, hycom_r, &
  nh, mh, lh )

    real, parameter :: lrg = 1.e9
    integer, intent(in) :: nh, mh, lh, n, m, l

    ! loop variables, and indices
    integer :: i, j, k, x, y, z, nx, ny, nz

    ! Variable for the interpolation calculations
    real*4 :: dx, dy, dz, s1, s2, q1, q2, w1, w2

    ! Initialize NCOM data
    real*4, dimension(n), intent(in) :: ncom_lon
    real*4, dimension(m), intent(in) :: ncom_lat
    real*4, dimension(n,m,l), intent(in) :: ncom_z
    real*4, dimension(n,m,l), intent(out) :: ncom_r

    ! Initialize HYCOM data
    real*4, dimension(nh), intent(in) :: hycom_lon
    real*4, dimension(mh), intent(in) :: hycom_lat
    real*4, dimension(nh,mh,lh), intent(in) :: hycom_z, hycom_r

    ! Initialize nearest_hycom_point
    ! This stores the nearest south west
    ! hycom point of each ncom point
    real, allocatable, save :: nearest_hycom_point(:,:,:)

    ! 1st run flag
    logical, save :: first_run = .True.

! ----------------- TIMING JUNK --------------------
    integer :: total_start, total_stop, clock_start, clock_stop, clock_rate
    integer :: calls_to_find, calls_to_find_depth
    real :: time_in_find, time_in_find_depth, total_time
    character(len=30) :: fmt, fmt2, fmt3

    calls_to_find = 0
    calls_to_find_depth = 0
    time_in_find = 0
    time_in_find_depth = 0
    total_time = 0
! ------------- END TIMING JUNK -----------------

    print *, "ENTERING HYCOM2NCOM SUBROUTINE"

    call system_clock( COUNT=total_start )
    allocate(nearest_hycom_point(n,m,l))

    do i = 1,n
      do j = 1,m
        if( first_run ) then
          ! Find the indices of the nearest HYCOM point
          ! to the SW of the NCOM point
          call system_clock( COUNT=clock_start )
          calls_to_find = calls_to_find + 1
            nearest_hycom_point(i,j,1) = &
            find( ncom_lon(i), hycom_lon, nh )
          call system_clock( COUNT=clock_stop )
          time_in_find = time_in_find + (clock_stop-clock_start)

          call system_clock( COUNT=clock_start )
          calls_to_find = calls_to_find + 1
            nearest_hycom_point(i,j,2) = &
            find( ncom_lat(j), hycom_lat, mh )
          call system_clock( COUNT=clock_stop )
          time_in_find = time_in_find + (clock_stop-clock_start)
          if( i==n .and. j==m ) then
            first_run = .false.
          endif
        endif
        do k = 1,l

          x = nearest_hycom_point(i,j,1) 
          y = nearest_hycom_point(i,j,2) 
          call system_clock( COUNT=clock_start )
          calls_to_find_depth = calls_to_find_depth + 1
            z = find_depth( ncom_z(i,j,k), hycom_z(x,y,:), lh )
          call system_clock( COUNT=clock_stop )
          time_in_find_depth = time_in_find_depth + (clock_stop-clock_start)

          ! find the weights for parameterization
          dx = abs(hycom_lon(x)-ncom_lon(i))/abs(hycom_lon(x)-hycom_lon(x+1))

          ! Squish the hexahedrom into a plane
          if(  abs(hycom_r(x,y,z)) > lrg &
            .or. abs(hycom_r(x,y,z-1)) > lrg ) then
            ! just leave s1 equal to the previous value for now
          else
            dy = abs(hycom_lat(y)-ncom_lat(j))/ &
                 abs(hycom_lat(y)-hycom_lat(y+1))
            dz = abs(hycom_z(x,y,z)-ncom_z(i,j,k))/ &
                 abs(hycom_z(x,y,z)-hycom_z(x,y,z-1))
            s1 = hycom_r(x,y,z)*(1.-dz)+hycom_r(x,y,z-1)*dz
          endif

          if(  abs(hycom_r(x+1,y,z)) > lrg &
            .or. abs(hycom_r(x+1,y,z-1)) > lrg ) then
              ! just leave s2 equal to the previous value for now
          else
            dy = abs(hycom_lat(y)-ncom_lat(j))/ &
                 abs(hycom_lat(y)-hycom_lat(y+1))
            dz = abs(hycom_z(x+1,y,z)-ncom_z(i,j,k))/ &
            abs(hycom_z(x+1,y,z)-hycom_z(x+1,y,z-1))
            s2 = hycom_r(x+1,y,z)*(1.-dz)+hycom_r(x+1,y,z-1)*dz
          endif

          if(  abs(hycom_r(x,y+1,z)) > lrg &
            .or. abs(hycom_r(x,y+1,z-1)) > lrg ) then
              ! just leave q1 equal to the previous value for now
          else
            dy = abs(hycom_lat(y)-ncom_lat(j))/ &
                 abs(hycom_lat(y)-hycom_lat(y+1))
            dz = abs(hycom_z(x,y+1,z)-ncom_z(i,j,k))/ &
            abs(hycom_z(x,y+1,z)-hycom_z(x,y+1,z-1))
            q1 = hycom_r(x,y+1,z)*(1.-dz)+hycom_r(x,y+1,z-1)*dz
          endif

          if(  abs(hycom_r(x+1,y+1,z)) > lrg &
            .or. abs(hycom_r(x+1,y+1,z-1)) > lrg ) then
              ! just leave q2 equal to the previous value for now
          else
            dy = abs(hycom_lat(y)-ncom_lat(j))/ &
                 abs(hycom_lat(y)-hycom_lat(y+1))
            dz = abs(hycom_z(x+1,y+1,z)-ncom_z(i,j,k))/ &
            abs(hycom_z(x+1,y+1,z)-hycom_z(x+1,y+1,z-1))
            q2 = hycom_r(x+1,y+1,z)*(1.-dz)+hycom_r(x+1,y+1,z-1)*dz
          endif

          ! Squish the plane into a line
            w1 = s1*(1.-dx)+s2*dx
            w2 = q1*(1.-dx)+q2*dx

          ! Squish the line into a point
            ncom_r(i,j,k) = w1*(1.-dy)+w2*dy

        end do
      end do
    end do

    call system_clock( COUNT=total_stop )
    total_time = total_stop-total_start
    call system_clock( COUNT_RATE=clock_rate )
       
    open(unit=90,file='../data/timing/time.dat',form='formatted',status='old')
    fmt = "(A,F,A)"
    fmt2 = "(A,I,A)"
    fmt3 = "(A)"
    write(90,fmt) "Time in find =", time_in_find/clock_rate, " seconds"
    write(90,fmt2) "Number of calls to find =", calls_to_find, " times"
    write(90,fmt) "Average time per call to find =", &
    (time_in_find/clock_rate)/calls_to_find, " seconds"
    write(90,fmt3) ""
    write(90,fmt) "Time in find_depth =", time_in_find_depth/clock_rate, &
    " seconds"
    write(90,fmt2) "Number of calls to find_depth =", calls_to_find_depth, &
    " times"
    write(90,fmt) "Average time per call to find_depth =", &
    (time_in_find_depth/clock_rate)/calls_to_find_depth, " seconds"
    write(90,fmt3) ""
    write(90,fmt) "Total time =", total_time/clock_rate, " seconds"
    write(90,fmt) "Percentage of time in find =", &
    time_in_find/total_time, " seconds"
    write(90,fmt) "Percentage of time in find_depth =", &
    time_in_find_depth/total_time, " seconds"
    print *, "LEAVING HYCOM2NCOM SUBROUTINE"
  end subroutine hycom2ncom

end module Interpolate
          
