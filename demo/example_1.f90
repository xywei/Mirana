  !> Test for ordinary differentiation
Program example_1
  use Mirana
  implicit none
  integer, parameter :: n1 = 1000, n2 = 500
  double precision, parameter :: L = 2.0
  double precision :: h1,h2,x,y,e1,e2,e1h1,e1h2,e2h1,e2h2,e11,e22,e12,elp
  integer :: i,j
  double precision, dimension(:,:),allocatable :: p, px, py, pxh1, pxh2, pyh1, pyh2, ppxx, ppyy, ppxy, pl
  allocate(p(n1,n2))
  allocate(px(n1,n2))
  allocate(pxh1(n1,n2))
  allocate(pxh2(n1,n2))
  allocate(py(n1,n2))
  allocate(pyh1(n1,n2))
  allocate(pyh2(n1,n2))
  allocate(ppxx(n1,n2))
  allocate(ppyy(n1,n2))
  allocate(ppxy(n1,n2))
  allocate(pl(n1,n2))
  write(*,*) 'Hello, this is example_1, testing for ordinary differentiation'
  h1 = L / n1
  h2 = L / n2
  do i = 1,n1
     do j = 1,n2
        x = h1 * i
        y = h2 * j
        p(i,j) = sin(2*x)*cos(3*y)
     end do
  end do

  call d1(p, h1, px)
  call d2(p, h2, py)
  call d1h1(p, h1, pxh1)
  call d1h2(p, h1, pxh2)
  call d2h2(p, h2, pyh2)
  call d2h1(p, h2, pyh1)
  call d12(p, h1, h2, ppxy)
  call d11(p, h1, ppxx)
  call d22(p, h2, ppyy)
  call dlp(p, h1, h2, pl)


  e1 = 0
  e2 = 0
  e1h1 = 0
  e1h2 = 0
  e2h1 = 0
  e2h2 = 0
  e11 = 0
  e12 = 0
  e22 = 0
  elp = 0
  do i = 1,n1
     do j = 1,n2
        x = h1 * i
        y = h2 * j
        if (i<n1) then
           e1h1 = e1h1 + abs( pxh1(i,j) - 2*cos(2*(x+h1*0.5))*cos(3*y) )
        end if
        if (j<n2) then
           e2h2 = e2h2 + abs( pyh2(i,j) + 3*sin(2*x)*sin(3*(y+h2*0.5)) )
        end if
        if (i>1 .and. i<n1) then
           e1 = e1 + abs( px(i,j) - 2*cos(2*x)*cos(3*y) )
           e11 = e11 + abs( ppxx(i,j) + 4*sin(2*x)*cos(3*y) )
        end if
        if (j>1 .and. j<n2) then
           e2 = e2 + abs( py(i,j) + 3*sin(2*x)*sin(3*y) )
           e22 = e22 + abs( ppyy(i,j) + 9*sin(2*x)*cos(3*y) )
        end if
        if (i>1 .and. j>1 .and. i<n1 .and. j<n2) then
           e1h2 = e1h2 + abs( pxh2(i,j) - 2*cos(2*(x))*cos(3*(y+0.5*h2)) )
           e2h1 = e2h1 + abs( pyh1(i,j) + 3*sin(2*x+h1)*sin(3*y) )
           e12 = e12 + abs( ppxy(i,j) + 6*cos(2*x)*sin(3*y) )
           elp = elp + abs( pl(i,j) + 4*sin(2*x)*cos(3*y) + 9*sin(2*x)*cos(3*y) )
        end if        
     end do
  end do

  write(*,*) "==========================="
  write(*,*) "Error noms: "
  write(*,*) "D1: ", e1/n1/n2
  write(*,*) "D2: ", e2/n1/n2
  write(*,*) "D1H1: ", e1h1/n1/n2
  write(*,*) "D1H2: ", e1h2/n1/n2
  write(*,*) "D2H1: ", e2h1/n1/n2
  write(*,*) "D1H1: ", e1h1/n1/n2
  write(*,*) "D11: ", e11/n1/n2
  write(*,*) "D12: ", e12/n1/n2
  write(*,*) "D22: ", e22/n1/n2
  write(*,*) "Dlp: ", elp/n1/n2
  write(*,*) "==========================="

  deallocate( p )
  deallocate( px )
  deallocate( py )
  deallocate( pyh2 )
  deallocate( pxh1 )
  deallocate( ppxx )
  deallocate( ppyy )
  deallocate( ppxy )  
end Program example_1

