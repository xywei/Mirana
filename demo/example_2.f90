  !> Test for mesh differentiation
Program example_2
  use Mirana
  implicit none
  integer, parameter :: n1 = 1000, n2 = 500
  double precision, parameter :: L = 2.0
  double precision :: h1,h2,ej,exi1,exi2,eet1,eet2,x,y
  integer :: i,j
  double precision, dimension(:,:),allocatable :: mesh1,mesh2,jacob,xi1,xi2,et1,et2
  allocate(mesh1(n1,n2))
  allocate(mesh2(n1,n2))
  allocate(jacob(n1,n2))
  allocate(xi1(n1,n2))
  allocate(xi2(n1,n2))
  allocate(et1(n1,n2))
  allocate(et2(n1,n2))
  write(*,*) 'Hello, this is example_2, testing for mesh differentiation.'
  h1 = L / n1
  h2 = L / n2
  do i = 1,n1
     do j = 1,n2
        x = h1 * i
        y = h2 * j
        mesh1(i,j) = (x-L/2.0)**3 + 1
        mesh2(i,j) = (y-L/2.0)**5 + 1
     end do
  end do

  call nmd(mesh1,mesh2,h1,h2,jacob,xi1,xi2,et1,et2)

  ej = 0
  exi1 = 0
  exi2 = 0
  eet1 = 0
  eet2 = 0
  do i = 1,n1
     do j = 1,n2
        x = h1 * i
        y = h2 * j
        exi1 = exi1 + 0.0
     end do
  end do

  write(*,*) "==========================="
  write(*,*) "Error noms: "
  write(*,*) "Jacobian: ", ej/n1/n2
  write(*,*) "xi_x: ", exi1/n1/n2
  write(*,*) "xi_y: ", exi2/n1/n2
  write(*,*) "eta_x: ", eet1/n1/n2
  write(*,*) "eta_y: ", eet2/n1/n2
  write(*,*) "==========================="

  deallocate( mesh1 )
  deallocate( mesh2 )
  deallocate( jacob )
  deallocate( xi1 )
  deallocate( xi2 )
  deallocate( et1 )
  deallocate( et2 )  
end Program example_2

