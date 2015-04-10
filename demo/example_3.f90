  !> Test for differentiation on the moving mesh
Program example_3
  use Mirana
  implicit none
  integer, parameter :: n1 = 1000, n2 = 1000
  integer :: st1 = -10, st2 = -1, ed1, ed2
  double precision, parameter :: L = 1.5 !> A square domain
  double precision :: h1,h2,e1,e2,elp,x,y
  integer :: i,j
  double precision, dimension(:,:),allocatable :: p,p1,p2,plp,mesh1,mesh2,jacob,xix,xiy,etx,ety,alpha,beta
  ed1 = st1 + n1 - 1
  ed2 = st2 + n2 - 1
  allocate(p(st1:ed1,st2:ed2))
  allocate(p1(st1:ed1,st2:ed2))
  allocate(p2(st1:ed1,st2:ed2))
  allocate(plp(st1:ed1,st2:ed2))
  allocate(mesh1(st1:ed1,st2:ed2))
  allocate(mesh2(st1:ed1,st2:ed2))
  allocate(jacob(st1:ed1,st2:ed2))
  allocate(xix(st1:ed1,st2:ed2))
  allocate(xiy(st1:ed1,st2:ed2))
  allocate(etx(st1:ed1,st2:ed2))
  allocate(ety(st1:ed1,st2:ed2))
  allocate(alpha(st1:ed1,st2:ed2))
  allocate(beta(st1:ed1,st2:ed2))
  write(*,*) 'Hello, this is example_3, testing for differentiation in computational domain.'
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
        p(i,j) = mesh1(i,j)**2 + 3.0*mesh1(i,j)*mesh2(i,j) + 2.0 * mesh2(i,j)**2
     end do
  end do

  call nmd(mesh1,mesh2,h1,h2,st1,st2,jacob,xix,xiy,etx,ety)
  call greek(mesh1,mesh2,xix,xiy,etx,ety,h1,h2,st1,st2,alpha,beta)

  call nmd1(p,xix,etx,h1,h2,st1,st2,p1)
  call nmd2(p,xiy,ety,h1,h2,st1,st2,p2)
  call nmdlp(p,xix,xiy,etx,ety,alpha,beta,h1,h2,st1,st2,plp) 

  e1 = 0
  e2 = 0
  elp = 0
  do i = st1+1,ed1-1
     do j = st2+1,ed2-1
        x = h1 * i - L/2
        y = h2 * j - L/2
        e1 = e1 + abs( p1(i,j) - ( 2*mesh1(i,j) + 3*mesh2(i,j) ) )
        e2 = e2 + abs( p2(i,j) - ( 3*mesh1(i,j) + 4*mesh2(i,j) ) )
        elp = elp + abs( plp(i,j) - 6 )
     end do
  end do

  write(*,*) "==========================="
  write(*,*) "Error noms: "
  write(*,*) "p_x: ", e1/n1/n2
  write(*,*) "p_y: ", e2/n1/n2
  write(*,*) "p_lp: ", elp/n1/n2
  write(*,*) "==========================="


  deallocate( mesh1 )
  deallocate( mesh2 )
  deallocate( jacob )
  deallocate( xix )
  deallocate( xiy )
  deallocate( etx )
  deallocate( ety )
  deallocate( p )
  deallocate( p1 )  
  deallocate( p2 )  
  deallocate( plp )  
  deallocate( alpha )  
  deallocate( beta )  
end Program example_3

