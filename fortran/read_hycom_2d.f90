program read_hycom_2D
use Interpolate_2D
implicit none

! hycom vars
  integer :: nh,mh,lh,i,j,k
  real*4, allocatable :: xh(:), yh(:)
  real*4, allocatable :: zh(:,:,:), th(:,:,:)

! ncom vars
  integer, parameter :: n=540, m=420, l=80
  real*4 :: x(n), y(m), h(n,m) 
  real*4, dimension(n,m,l) :: SMt
  real*4, dimension(n,m) :: t, dif

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

  call hycom2ncom_2D( x, y, t, n, m, xh, yh, th(:,:,1), nh, mh )
  
  dif = abs(SMt(:,:,1)-t)

  open(unit=8, file='../data/f90data/SMinterp_2d.dat', form='unformatted')
    write(8) SMt(:,:,1)
  close(8)

  open(unit=8, file='../data/f90data/NCinterp_2d.dat', form='unformatted')
    write(8) t
  close(8)

  open(unit=8, file='../data/f90data/DIFinterp_2d.dat', form='unformatted')
    write(8) dif
  close(8)
  
  open(unit=8, file='../data/f90data/longitude_2d.dat', form='unformatted')
    write(8) x
  close(8)

  open(unit=8, file='../data/f90data/latitude_2d.dat', form='unformatted')
    write(8) y
  close(8)

end program read_hycom_2D
