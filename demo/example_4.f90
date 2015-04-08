!> Test for calling C function
Program example_4
  !use itsol2
  implicit none
  integer :: n=10,nnz=10
  integer :: ierr,i
  double precision, allocatable, dimension(:) :: val,rhs,init
  double precision, allocatable, dimension(:) :: solu
  integer, allocatable, dimension(:) :: col_ind,row_ptr

  allocate( val(nnz) )
  allocate( rhs(n) )
  allocate( init(n) )
  allocate( solu(n) )
  allocate( col_ind(nnz) )
  allocate( row_ptr(n+1) )

  do i = 1,nnz
     val(i) = i**2
     col_ind(i) = i
  end do

  do i = 1,n
     rhs(i) = i**3
     row_ptr(i) = i
     init(i) = 0.0
  end do

  row_ptr(n+1) = n+1  

  write(*,*) "Matrix val:"
  write(*,*) val

  write(*,*) "Matrix col_ind:"
  write(*,*) col_ind

  write(*,*) "Matrix row_ptr:"
  write(*,*) row_ptr

  write(*,*) "Right hand side:"
  write(*,*) rhs

  write(*,*) "Initial guess:"
  write(*,*) init
  
  write(*,*) "Hi, this is example 4, testing for linear solvers."

  call arms_fgmres(n,val,col_ind,row_ptr,rhs,solu,init,ierr)

  write(*,*) "ierr = ",ierr,"."
  
end Program example_4
     
     
