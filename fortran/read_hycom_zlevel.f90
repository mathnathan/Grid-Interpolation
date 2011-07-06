program read_hycom_zlevel
use Interpolate
implicit none

! hycom vars
  integer :: nh,mh,lh,i,j,k
  real*4, allocatable :: xh(:), yh(:)
  real*4, allocatable :: zh(:,:,:), th(:,:,:)

! ncom vars
  integer, parameter :: n=540, m=420, l=80
  real*4 :: x(n), y(m), h(n,m)
  real*4 :: zw(n,m,l+1), z(n,m,l)
  real*4, dimension(n,m,l) :: t, SMt, dif

  open(unit=90,file='../data/f90data/hycom_zlevel_20090101.dat', form='unformatted', &
    status='old')
  read(unit=90) nh,mh,lh
  allocate(xh(nh))
  allocate(yh(mh))
  allocate(zh(nh,mh,lh))
  allocate(th(nh,mh,lh))
  read(unit=90) xh
  read(unit=90) yh
  read(unit=90) zh
  read(unit=90) th
  close(unit=90)

  open(unit=90,file='../data/f90data/ts20090101.dat', form='unformatted', &
    status='old')
  read(unit=90) SMt
  close(unit=90)
  
  do i=1,n
    x(i)=-92.25+(i-1)/120.
  enddo
  do j=1,m
    y(j)=25+(j-1)/120.
  enddo

  open(unit=90,file= '../data/f90data/sigsbee-topo2-3filt-blended.dat', &
    status='old',form='unformatted')
  read(unit=90) h
  close(unit=90)

  open(unit=90, file= '../data/f90data/zwt.dat', &
    status='old', form='unformatted')
  read(unit=90) zw
  close(unit=90)

  do k=1,l
    do j=1,m
      do i=1,n
        z(i,j,k)=0.5*(zw(i,j,k)+zw(i,j,k+1))
      enddo
    enddo
  enddo

  call hycom2ncom( x, y, z, t, n, m, l, &
                   xh, yh, zh, th, nh, mh, lh )
  
  dif = abs(SMt-t)

  open(unit=8, file='../data/f90data/SMinterp.dat', form='unformatted')
    write(8) SMt
  close(8)

  open(unit=8, file='../data/f90data/NCinterp.dat', form='unformatted')
    write(8) t
  close(8)

  open(unit=8, file='../data/f90data/DIFinterp.dat', form='unformatted')
    write(8) dif
  close(8)
  
  open(unit=8, file='../data/f90data/longitude.dat', form='unformatted')
    write(8) x
  close(8)

  open(unit=8, file='../data/f90data/latitude.dat', form='unformatted')
    write(8) y
  close(8)

  open(unit=8, file='../data/f90data/depth.dat', form='unformatted')
    write(8) z
  close(8)

end program read_hycom_zlevel
