module Interpolate_2D
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

  ! This is the trilinear interpolation algorithm. Converts
  ! the HYCOM grid data to the NCOM grid data.
  subroutine hycom2ncom_2D( ncom_lon, ncom_lat, ncom_r, n, m, &
                          hycom_lon, hycom_lat, hycom_r, nh, mh )

    integer, intent(in) :: nh, mh, n, m

    ! loop variables, and indices
    integer :: i, j, x, y

    ! Variable for the interpolation calculations
    real*4 :: dx, dy, s1, s2

    ! Initialize NCOM data
    real*4, dimension(n), intent(in) :: ncom_lon
    real*4, dimension(m), intent(in) :: ncom_lat
    real*4, dimension(n,m), intent(out) :: ncom_r

    ! Initialize HYCOM data
    real*4, dimension(nh), intent(in) :: hycom_lon
    real*4, dimension(mh), intent(in) :: hycom_lat
    real*4, dimension(nh,mh), intent(in) :: hycom_r

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

    do i = 1,n
      do j = 1,m
          ! Find the indices of the nearest HYCOM point
          ! to the SW of the NCOM point
          call system_clock( COUNT=clock_start )
          calls_to_find = calls_to_find + 1
            x = find( ncom_lon(i), hycom_lon, nh )
          call system_clock( COUNT=clock_stop )
          time_in_find = time_in_find + (clock_stop-clock_start)

          call system_clock( COUNT=clock_start )
          calls_to_find = calls_to_find + 1
            y = find( ncom_lat(j), hycom_lat, mh )
          call system_clock( COUNT=clock_stop )
          time_in_find = time_in_find + (clock_stop-clock_start)

          ! find the weights for parameterization
          dx = abs(hycom_lon(x)-ncom_lon(i))/abs(hycom_lon(x)-hycom_lon(x+1))
          dy = abs(hycom_lat(y)-ncom_lat(j))/abs(hycom_lat(y)-hycom_lat(y+1))

          ! Squish the plane into a line
          s1 = hycom_r(x,y)*(1.-dx)+hycom_r(x+1,y)*dx
          s2 = hycom_r(x,y+1)*(1.-dx)+hycom_r(x+1,y+1)*dx

          ! Squish the line into a point
          ncom_r(i,j) = s1*(1.-dy)+s2*dy
      end do
    end do

    call system_clock( COUNT=total_stop )
    total_time = total_stop-total_start
    call system_clock( COUNT_RATE=clock_rate )
       
    open(unit=90,file='../data/timing/time_2d.dat',form='formatted',status='old')
    fmt = "(A,F,A)"
    fmt2 = "(A,I,A)"
    fmt3 = "(A)"
    write(90,fmt) "Time in find =", time_in_find/clock_rate, " seconds"
    write(90,fmt2) "Number of calls to find =", calls_to_find, " times"
    write(90,fmt) "Average time per call to find =", &
    (time_in_find/clock_rate)/calls_to_find, " seconds"
    write(90,fmt3) ""
    write(90,fmt) "Total time =", total_time/clock_rate, " seconds"
    write(90,fmt) "Percentage of time in find =", &
    time_in_find/total_time, " seconds"
    print *, "LEAVING HYCOM2NCOM SUBROUTINE"
  end subroutine hycom2ncom_2D

end module Interpolate_2D
          
