  !> Test for ordinary differentiation
Program example_1
  use Mirana
  implicit none
  integer, parameter :: n = 1000
  double precision, parameter :: L = 2.0
  double precision :: h,x,y,e1,e2,e11,e22,e12,elp
  integer :: i,j
  double precision, dimension(:,:),allocatable :: p, px, py, pxx, pyy, pxy, pl
  allocate(p(n,n))
  allocate(px(n,n))
  allocate(py(n,n))
  allocate(pxx(n,n))
  allocate(pyy(n,n))
  allocate(pxy(n,n))
  allocate(pl(n,n))
  write(*,*) 'Hello, this is example_1, testing for ordinary differentiation'
  h = L / n
  do i = 1,n
     do j = 1,n
        x = h * i
        y = h * j
        p(i,j) = sin(2*x)*cos(3*y)
     end do
  end do

  call d1(p, h, px)
  call d2(p, h, py)
  call d12(p, h, pxy)
  call d11(p, h, pxx)
  call d22(p, h, pyy)
  call dlp(p, h, pl)

  e1 = 0
  e2 = 0
  e11 = 0
  e12 = 0
  e22 = 0
  elp = 0
  do i = 2,(n-1)
     do j = 2,(n-1)
        x = h * i
        y = h * j
        e1 = e1 + abs( px(i,j) - 2*cos(2*x)*cos(3*y) )
        e2 = e2 + abs( py(i,j) + 3*sin(2*x)*sin(3*y) )
        e11 = e11 + abs( pxx(i,j) + 4*sin(2*x)*cos(3*y) )
        e12 = e12 + abs( pxy(i,j) + 6*cos(2*x)*sin(3*y) )
        e22 = e22 + abs( pyy(i,j) + 9*sin(2*x)*cos(3*y) )
        elp = elp + abs( pl(i,j) + 4*sin(2*x)*cos(3*y) + 9*sin(2*x)*cos(3*y) )
     end do
  end do

  write(*,*) "==========================="
  write(*,*) "Error noms: "
  write(*,*) "D1: ", e1/(n-2)**2
  write(*,*) "D2: ", e2/(n-2)**2
  write(*,*) "D11: ", e11/(n-2)**2
  write(*,*) "D12: ", e12/(n-2)**2
  write(*,*) "D22: ", e22/(n-2)**2
  write(*,*) "Dlp: ", elp/(n-2)**2
  write(*,*) "==========================="

  deallocate( p ); deallocate( px ); deallocate( py )
  deallocate( pxx ); deallocate( pyy ); deallocate( pxy )  
end Program example_1

