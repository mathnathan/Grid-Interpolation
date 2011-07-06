program test
implicit none

interface
  integer function find( value, array, array_size )
    integer, intent(in) :: array_size
    real, intent(in) :: value, array( array_size )
  end function find
end interface
  
  integer, parameter :: nh=2, mh=2, lh=2, n=1, m=1, l=1

  integer :: i, j, k, x, y, z

  real*4 :: dx, dy, dz, s1, s2, q1, q2, w1, w2

  real*4, dimension(n) :: ncom_lon
  real*4, dimension(m) :: ncom_lat
  real*4, dimension(n,m,l) :: ncom_z
  real*4, dimension(n,m,l) :: ncom_r

  real*4, dimension(nh) :: hycom_lon
  real*4, dimension(mh) :: hycom_lat
  real*4, dimension(nh,mh,lh) :: hycom_z, hycom_r


  ncom_lon(1) = 5
  ncom_lat(1) = 2
  ncom_z(1,1,1) = -1

  hycom_lon(1) = 0
  hycom_lon(2) = 6
  hycom_lat(1) = 0
  hycom_lat(2) = 6

  hycom_z(1,1,1) = 0
  hycom_z(1,2,1) = 0
  hycom_z(2,1,1) = 0
  hycom_z(2,2,1) = 0
  hycom_z(1,1,2) = -6
  hycom_z(1,2,2) = -6
  hycom_z(2,1,2) = -6
  hycom_z(2,2,2) = -6

  hycom_r(1,1,1) = 2
  hycom_r(1,2,1) = 1
  hycom_r(2,1,1) = 1
  hycom_r(2,2,1) = 6
  hycom_r(1,1,2) = 3
  hycom_r(1,2,2) = 2
  hycom_r(2,1,2) = 0
  hycom_r(2,2,2) = 4

  do i = 1,n
    do j = 1,m
      do k = 1,l

        x = find( ncom_lon(i), hycom_lon, nh )
        y = find( ncom_lat(j), hycom_lat, mh )
        z = find( ncom_z(i,j,k), hycom_z(x,y,:), lh )

        ! find the weights for parameterization
        dx = abs(hycom_lon(x)-ncom_lon(i))/abs(hycom_lon(x)-hycom_lon(x+1))

        ! Squish the hexahedrom into a plane
          if(  hycom_r(x,y,z) > 1000000000 &
          .or. hycom_r(x,y,z) < -1000000000 &
          .or. hycom_r(x,y,z-1) > 1000000000 &
          .or. hycom_r(x,y,z-1) < -1000000000 ) then
            ! just leave s1 equal to the previous value for now
          else
            dy = abs(hycom_lat(y)-ncom_lat(j))/ &
                 abs(hycom_lat(y)-hycom_lat(y+1))
            dz = abs(hycom_z(x,y,z)-ncom_z(i,j,k))/ &
                 abs(hycom_z(x,y,z)-hycom_z(x,y,z-1))
            s1 = hycom_r(x,y,z)*(1-dz)+hycom_r(x,y,z-1)*dz
          endif

          if(  hycom_r(x+1,y,z) > 1000000000 &
          .or. hycom_r(x+1,y,z) < -1000000000 &
          .or. hycom_r(x+1,y,z-1) > 1000000000 &
          .or. hycom_r(x+1,y,z-1) < -1000000000 ) then
            ! just leave s2 equal to the previous value for now
          else
            dy = abs(hycom_lat(y)-ncom_lat(j))/ &
                 abs(hycom_lat(y)-hycom_lat(y+1))
            dz = abs(hycom_z(x+1,y,z)-ncom_z(i,j,k))/ &
            abs(hycom_z(x+1,y,z)-hycom_z(x+1,y,z-1))
            s2 = hycom_r(x+1,y,z)*(1-dz)+hycom_r(x+1,y,z-1)*dz
          endif

          if(  hycom_r(x,y+1,z) > 1000000000 &
          .or. hycom_r(x,y+1,z) < -1000000000 &
          .or. hycom_r(x,y+1,z-1) > 1000000000 &
          .or. hycom_r(x,y+1,z-1) < -1000000000 ) then
            ! just leave q1 equal to the previous value for now
          else
            dy = abs(hycom_lat(y)-ncom_lat(j))/ &
                 abs(hycom_lat(y)-hycom_lat(y+1))
            dz = abs(hycom_z(x,y+1,z)-ncom_z(i,j,k))/ &
            abs(hycom_z(x,y+1,z)-hycom_z(x,y+1,z-1))
            q1 = hycom_r(x,y+1,z)*(1-dz)+hycom_r(x,y+1,z-1)*dz
          endif

          if(  hycom_r(x+1,y+1,z) > 1000000000 &
          .or. hycom_r(x+1,y+1,z) < -1000000000 &
          .or. hycom_r(x+1,y+1,z-1) > 1000000000 &
          .or. hycom_r(x+1,y+1,z-1) < -1000000000 ) then
            ! just leave q2 equal to the previous value for now
          else
            dy = abs(hycom_lat(y)-ncom_lat(j))/ &
                 abs(hycom_lat(y)-hycom_lat(y+1))
            dz = abs(hycom_z(x+1,y+1,z)-ncom_z(i,j,k))/ &
            abs(hycom_z(x+1,y+1,z)-hycom_z(x+1,y+1,z-1))
            q2 = hycom_r(x+1,y+1,z)*(1-dz)+hycom_r(x+1,y+1,z-1)*dz
          endif

        ! Squish the plane into a line
          w1 = s1*(1-dx)+s2*dx
          w2 = q1*(1-dx)+q2*dx

        ! Squish the line into a point
          ncom_r(i,j,k) = w1*(1-dy)+w2*dy

      end do
    end do
  end do

  ! If you do it out by hand it should 
  ! be 2277/972 or 2.342592593...
  print *, "ncom_r(1,1,1) =", ncom_r(1,1,1)

end program test

! Find the location of the largest number
! in the array that is less than value
integer function find( value,  array, array_size )
  implicit none

  integer, intent(in) :: array_size
  real, intent(in) :: value, array( array_size )
  integer, dimension(1) :: indx
  integer :: i 

  indx = minloc( abs(array-value) )
  if( array( indx(1) ) > value ) then
    if( value < 0 ) then
      indx(1) = indx(1) + 1
    else
      indx(1) = indx(1) - 1
    endif
  endif
  find = indx(1)

end function find

