! The main module
! Doxygen tips:
!   "!>" or "!<" starts a comment
!   "!!" or "!>" can be used to continue an one line comment into a multi-line comment.
! 

!> Provide the subroutines & functions for moving mesh computation.
!! @todo Add mesh derivatives
!! @todo add laplacian
Module Mirana
  implicit none
  
contains
  !> Handling exceptions
  !! @param n integer, exit code
  !! @param message string with length 128, error message.
  subroutine mirana_exception(message)
    implicit none
    character, intent(in) :: message(128)
    write(*,*) message
    stop 1
  end subroutine mirana_exception

  !> Allocate memory for new arrays
  !!
  subroutine mirana_memory(name, dim)
    implicit none
    double precision, allocatable, dimension(:,:), intent(inout) :: name
    integer, intent(in) :: dim
    if (allocated(name)) then
       call mirana_exception("Error in module Mirana allocating memory: repeated allocation!")
    end if
    
  end subroutine mirana_memory
    
  !> Differentiate with respect to the first dimension using central difference.
  !! @param phi real N-by-N matrix.
  !! @param h real number, step size.
  !! @return phi_x real N-by-N matrix, only points ranged in 2:(N-1) are updated.
  subroutine d1(phi, h, phi_x)
    implicit none
    integer :: n !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_x
    if (.not.allocated(phi) .or. .not.allocated(phi_x)) then
       call mirana_exception("Error in module Mirana: array not allocated!")
    end if
    n = size(phi,1)
    phi_x(2:(n-1),2:(n-1)) = 0.5 * ( phi(3:n,2:(n-1)) - phi(1:(n-2),2:(n-1)) ) / h
    return
  end subroutine d1

  !> Differentiate with respect to the second dimension using central difference.
  !! @param phi real N-by-N matrix
  !! @param h real number, step size.
  !! @return phi_y real N-by-N matrix, only points ranged in 2:(N-1) are updated.
  subroutine d2(phi, h, phi_y)
    implicit none
    integer :: n !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_y
    if (.not.allocated(phi) .or. .not.allocated(phi_y)) then
       call mirana_exception("Error in module Mirana: array not allocated!")
    end if
    n = size(phi,1)
    phi_y(2:(n-1),2:(n-1)) = 0.5 * ( phi(2:(n-1),3:n) - phi(2:(n-1),1:(n-2)) ) / h
    return
  end subroutine d2

  !> Compute the second order derivative along the first dimension using central difference.
  !! @param phi N-by-N matrix.
  !! @param h real number, step size.
  !! @param phi_xx real N-by-N matrix, only points ranged in 2:(N-1) are updated.
  subroutine d11(phi, h, phi_xx)
    implicit none
    integer :: n !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_xx
    if (.not.allocated(phi) .or. .not.allocated(phi_xx)) then
       call mirana_exception("Error in module Mirana: array not allocated!")
    end if
    n = size(phi,1)
    phi_xx(2:(n-1),2:(n-1)) = ( phi(3:n,2:(n-1)) - 2.0*phi(2:(n-1),2:(n-1)) + phi(1:(n-2),2:(n-1)) ) / (h**2)
    return
  end subroutine d11

  !> Comopute the mixed second order derivative using central difference.
  !! @param phi N-by-N matrix.
  !! @param h real number, step size.
  !! @param phi_xy real N-by-N matrix, only points ranged in 2:(N-1) are updated.
  subroutine d12(phi, h, phi_xy)
    implicit none
    integer :: n !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_xy
    if (.not.allocated(phi) .or. .not.allocated(phi_xy)) then
       call mirana_exception("Error in module Mirana: array not allocated!")
    end if
    n = size(phi,1)
    phi_xy(2:(n-1),2:(n-1)) = ( phi(1:(n-2),1:(n-2)) + phi(3:n,3:n) - phi(1:(n-2),3:n) - phi(3:n,1:(n-2)) ) / (4.0 * h**2)
  end subroutine d12

  !> Compute the second order derivative along the second dimension using central difference.
  !! @param phi N-by-N real matrix.
  !! @param h real number, step size.
  !! @param phi_yy real N-by-N matrix, only points ranged in 2:(N-1) are updated.
  subroutine d22(phi, h, phi_yy)
    implicit none
    integer :: n !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_yy
    if (.not.allocated(phi) .or. .not.allocated(phi_yy)) then
       call mirana_exception("Error in module Mirana: array not allocated!")
    end if
    n = size(phi,1)
    phi_yy(2:(n-1),2:(n-1)) = ( phi(2:(n-1),3:n) - 2.0*phi(2:(n-1),2:(n-1)) + phi(2:(n-1),1:(n-2)) ) / (h**2)
    return
  end subroutine d22

  !> Compute Laplacian using central difference.
  !! @param phi N-by-N matrix.
  !! @param h real number, step size.
  !! @param phi_lp real N-by-N matrix, only points ranged in 2:(N-1) are updated.
  subroutine dlp(phi, h, phi_lp)
    implicit none
    integer :: n !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_lp
    double precision, allocatable, dimension(:,:) :: pxx, pyy
    n = size(phi,1)
    allocate( pxx(n,n) ); allocate( pyy(n,n) )
    if (.not.allocated(phi) .or. .not.allocated(phi_lp)) then
       call mirana_exception("Error in module Mirana: array not allocated!")
    end if
    call d11(phi, h, pxx);
    call d22(phi, h, pyy);
    phi_lp(2:(n-1),2:(n-1)) = pxx(2:(n-1),2:(n-1)) + pyy(2:(n-1),2:(n-1))
    deallocate( pxx ); deallocate( pyy )
    return
  end subroutine dlp

  !> (Non-conservative Mesh Differentiation) Compute derivatives of mesh grid with respect to computational space and inverse using non-conservative form central difference.
  !! @param mx N-by-N real matrix, x-coordinate of mesh grid.
  !! @param my N-by-N real matrix, y-coordinate of mesh grid.
  !! @param hxi real number, step size in computational domain along \f$\xi\f$ direction.
  !! @param het real number, step size in computational domain along \f$\eta\f$ direction.
  !! @param jcb real N-by-N matrix, only points ranged in 2:(N-1) are updated.
  !! @param xi_x real N-by-N matrix, only points ranged in 2:(N-1) are updated.
  !! @param xi_y real N-by-N matrix, only points ranged in 2:(N-1) are updated.
  !! @param et_x real N-by-N matrix, only points ranged in 2:(N-1) are updated.
  !! @param et_y real N-by-N matrix, only points ranged in 2:(N-1) are updated.  
  subroutine nmd(mx,my,hxi,het,jcb,xi_x,xi_y,et_x,et_y)
    implicit none
    integer :: n !> dimension of arrays
    double precision, intent(in) :: het, hxi
    double precision, allocatable, dimension(:,:), intent(in) :: mx, my
    double precision, allocatable, dimension(:,:), intent(inout) :: jcb, xi_x, xi_y, et_x, et_y
    double precision, allocatable, dimension(:,:) :: x1, x2, y1, y2
    integer, allocatable, dimension(:) :: indx
    n = size(phi,1)
    allocate( indx(n-2) )
    indx = 2:(n-1)
    allocate( x1(n,n) ); allocate( x2(n,n) );
    allocate( y1(n,n) ); allocate( y2(n,n) );
    if (.not.allocated(mx) .or. .not.allocated(my) .or. .not.allocated(jcb) .or. &
         .not.allocated(xi_x) .or. .not.allocated(xi_y) .or. .not.allocated(et_x) .or. .not.allocated(et_y) ) then
       call mirana_exception("Error in module Mirana: array not allocated!")
    end if
    call d1(mx,hxi,x1)
    call d2(mx,het,x2)
    call d1(my,hxi,y1)
    call d2(my,het,y2)
    jcb(2:(n-1),2:(n-1)) = x1 * y2 - x2 * y1
    deallocate( x1 ); deallocate( x2 )
    deallocate( y1 ); deallocate( y2 )
    return
  end subroutine nmd1
  
end Module Mirana
