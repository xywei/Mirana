  !> Test for mesh differentiation
Program example_2
  use Mirana
  implicit none
  integer, parameter :: n1 = 30, n2 = 50
  integer :: st1 = -3, st2 = -1, ed1, ed2
  double precision, parameter :: L = 1.5 !> A square domain
  double precision :: h1,h2,ej,exi1,exi2,eet1,eet2,x,y
  integer :: i,j
  double precision, dimension(:,:),allocatable :: mesh1,mesh2,jacob,xi1,xi2,et1,et2
  ed1 = st1 + n1 - 1
  ed2 = st2 + n2 - 1
  allocate(mesh1(st1:ed1,st2:ed2))
  allocate(mesh2(st1:ed1,st2:ed2))
  allocate(jacob(st1:ed1,st2:ed2))
  allocate(xi1(st1:ed1,st2:ed2))
  allocate(xi2(st1:ed1,st2:ed2))
  allocate(et1(st1:ed1,st2:ed2))
  allocate(et2(st1:ed1,st2:ed2))
  write(*,*) 'Hello, this is example_2, testing for mesh differentiation.'
  h1 = L / n1
  h2 = L / n2
  do i = st1,ed1
     do j = st2,ed2
        x = h1 * i - L/2
        y = h2 * j - L/2
        !> - The x-mesh \f$\xi=tanh(5x)\f$.
        !> - The y-mesh \f$\eta=tanh(y)\f$.
        mesh1(i,j) = 5 * atanh(x)
        mesh2(i,j) = atanh(y)
     end do
  end do

  call nmd(mesh1,mesh2,h1,h2,st1,st2,jacob,xi1,xi2,et1,et2)
  
  ej = 0
  exi1 = 0
  exi2 = 0
  eet1 = 0
  eet2 = 0
  do i = st1+1,ed1-1
     do j = st2+1,ed2-1
        x = h1 * i - L/2
        y = h2 * j - L/2
        exi1 = exi1 + abs( xi1(i,j) -  0.2 * (cosh(mesh1(i,j)/5.0)**(-2)) )
        exi2 = exi2 + abs( xi2(i,j) )
        eet1 = eet1 + abs( et1(i,j) )
        eet2 = eet2 + abs( et2(i,j) - cosh(mesh2(i,j))**(-2) )
     end do
  end do

  write(*,*) "==========================="
  write(*,*) "Error noms: "
  !write(*,*) "Jacobian: ", ej/n1/n2
  write(*,*) "xi_x: ", exi1/n1/n2
  write(*,*) "xi_y: ", exi2/n1/n2
  write(*,*) "eta_x: ", eet1/n1/n2
  write(*,*) "eta_y: ", eet2/n1/n2
  write(*,*) "==========================="

  call save_mesh(mesh1,mesh2,st1,st2)

  deallocate( mesh1 )
  deallocate( mesh2 )
  deallocate( jacob )
  deallocate( xi1 )
  deallocate( xi2 )
  deallocate( et1 )
  deallocate( et2 )  
end Program example_2

