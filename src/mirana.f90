! Doxygen tips:
!   "!>" or "!<" starts a comment.
!   "!!" or "!>" can be used to continue an one line comment into a multi-line comment.
!
! Renaming facility:
! The names of any object in a module canbe renamed to resolve name clashes
! Syntax:
!      USE module_name, newName1 => Name1, &          
!                       newName2 => Name2, ...
!
! NOTE:
!      Once a new name newName1 is assigned, the object Name1 cannot be accessed using Name1 any longer ! 


!> Provide the subroutines & functions for moving mesh computation.
!! @todo Change to M-by-N matrices
Module Mirana
  implicit none

  public
    
contains
  !> Handling exceptions.
  !! @param n integer, exit code.
  !! @param message string with length 128, error message.
  subroutine mirana_exception(message)
    implicit none
    character, intent(in) :: message(128)
    write(*,*) message
    stop 1
  end subroutine mirana_exception

  !> Allocate memory for new arrays.
  !! @param name real allocatable 2-dim array (not allocated).
  !! @param dim integer, dimension.
  !! @return name dim1-by-dim2 real array with memory allocated.
  subroutine mirana_memory(name, dim1, dim2)
    implicit none
    double precision, allocatable, dimension(:,:), intent(inout) :: name
    integer, intent(in) :: dim1, dim2
    integer :: allocate_err
    if (allocated(name)) then
       call mirana_exception("Error in module Mirana allocating memory: repeated allocation!")
    end if
    allocate( name(dim1,dim2), stat=allocate_err )
    if (allocate_err > 0) then
       call mirana_exception("Error in module Mirana allocating memory: not enough memory!")
    end if
    return
  end subroutine mirana_memory
    
  !> Differentiate with respect to the first dimension using central difference.
  !! @param phi real N1-by-N2 matrix.
  !! @param h real number, step size.
  !! @return phi_x real N1-by-N2 matrix, only points ranged in [2:(N1-1)]X[1:N2] are updated.
  subroutine d1(phi, h, phi_x)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_x
    if (.not.allocated(phi) .or. .not.allocated(phi_x)) then
       call mirana_exception("Error in module Mirana.d1: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    phi_x(2:(n1-1),1:n2) = 0.5 * ( phi(3:n1,1:n2) - phi(1:(n1-2),1:n2) ) / h
    return
  end subroutine d1

  !> Differentiate with respect to the second dimension using central difference.
  !! @param phi real N1-by-N2 matrix
  !! @param h real number, step size.
  !! @return phi_y real N1-by-N2 matrix, only points ranged in [1:N1]X[2:(N2-1)] are updated.
  subroutine d2(phi, h, phi_y)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_y
    if (.not.allocated(phi) .or. .not.allocated(phi_y)) then
       call mirana_exception("Error in module Mirana.d2: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    phi_y(1:n1,2:(n2-1)) = 0.5 * ( phi(1:n1,3:n2) - phi(1:n1,1:(n2-2)) ) / h
    return
  end subroutine d2

  !> Compute the second order derivative along the first dimension using central difference.
  !! @param phi N1-by-N2 matrix.
  !! @param h real number, step size.
  !! @param phi_xx real N1-by-N2 matrix, only points ranged in [2:(N-1)]X[1:N] are updated.
  subroutine d11(phi, h, phi_xx)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_xx
    if (.not.allocated(phi) .or. .not.allocated(phi_xx)) then
       call mirana_exception("Error in module Mirana.d11: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    phi_xx(2:(n1-1),1:n2) = ( phi(3:n1,1:n2) - 2.0*phi(2:(n1-1),1:n2) + phi(1:(n1-2),1:n2) ) / (h**2)
    return
  end subroutine d11

  !> Comopute the mixed second order derivative using central difference.
  !! @param phi N1-by-N2 matrix.
  !! @param h1 real number, step size 1.
  !! @param h2 real number, step size 2.
  !! @param phi_xy real N1-by-N2 matrix, only points ranged in [2:(N1-1)]X[2:(N2-1)] are updated.
  subroutine d12(phi, h1, h2, phi_xy)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    double precision, intent(in) :: h1,h2
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_xy
    if (.not.allocated(phi) .or. .not.allocated(phi_xy)) then
       call mirana_exception("Error in module Mirana.d12: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    phi_xy(2:(n1-1),2:(n2-1)) = ( phi(1:(n1-2),1:(n2-2)) + phi(3:n1,3:n2) - &
                                phi(1:(n1-2),3:n2) - phi(3:n1,1:(n2-2)) ) / (4.0 * h1 * h2)
  end subroutine d12

  !> Compute the second order derivative along the second dimension using central difference.
  !! @param phi N1-by-N2 real matrix.
  !! @param h real number, step size.
  !! @param phi_yy real N1-by-N2 matrix, only points ranged in [1:N1]X[2:(N2-1)] are updated.
  subroutine d22(phi, h, phi_yy)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_yy
    if (.not.allocated(phi) .or. .not.allocated(phi_yy)) then
       call mirana_exception("Error in module Mirana.d22: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    phi_yy(1:n1,2:(n2-1)) = ( phi(1:n1,3:n2) - 2.0*phi(1:n1,2:(n2-1)) + phi(1:n1,1:(n2-2)) ) / (h**2)
    return
  end subroutine d22

  !> Compute Laplacian using central difference.
  !! @param phi N1-by-N2 matrix.
  !! @param h1 real number, step size 1.
  !! @param h2 real number, step size 2.
  !! @param phi_lp real N1-by-N2 matrix, only points ranged in [2:(N1-1)]X[2:(N2-1)] are updated.
  subroutine dlp(phi, h1, h2, phi_lp)
    implicit none
    integer :: n1, n2 !> dimension of arrays
    double precision, intent(in) :: h1,h2
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_lp
    double precision, allocatable, dimension(:,:) :: pxx, pyy
    n1 = size(phi,1)
    n2 = size(phi,2)
    allocate( pxx(n1,n2) ); allocate( pyy(n1,n2) )
    if (.not.allocated(phi) .or. .not.allocated(phi_lp)) then
       call mirana_exception("Error in module Mirana.dlp: array not allocated!")
    end if
    call d11(phi, h1, pxx);
    call d22(phi, h2, pyy);
    phi_lp(2:(n1-1),2:(n2-1)) = pxx(2:(n1-1),2:(n2-1)) + pyy(2:(n1-1),2:(n2-1))
    deallocate( pxx ); deallocate( pyy )
    return
  end subroutine dlp

  !> (Non-conservative Mesh Differentiation) Compute derivatives of mesh grid with respect to computational space and inverse using non-conservative form central difference.
  !! @param mx N1-by-N2 real matrix, x-coordinate of mesh grid (fully determined).
  !! @param my N1-by-N2 real matrix, y-coordinate of mesh grid (fully determined).
  !! @param hxi real number, step size in computational domain along \f$\xi\f$ direction.
  !! @param het real number, step size in computational domain along \f$\eta\f$ direction.
  !! @return jcb real N1-by-N2 matrix, the Jacobian \f$J\f$.
  !! @return xi_x real N1-by-N2 matrix, \f$\partial\xi/\partial x\f$.
  !! @return xi_y real N1-by-N2 matrix, \f$\partial\xi/\partial y\f$.
  !! @return et_x real N1-by-N2 matrix, \f$\partial\eta/\partial x\f$.
  !! @return et_y real N1-by-N2 matrix, \f$\partial\eta/\partial y\f$.
  subroutine nmd(mx,my,hxi,het,jcb,xi_x,xi_y,et_x,et_y)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    integer :: i
    double precision, intent(in) :: het, hxi
    double precision, allocatable, dimension(:,:), intent(in) :: mx, my
    double precision, allocatable, dimension(:,:), intent(inout) :: jcb, xi_x, xi_y, et_x, et_y
    double precision, allocatable, dimension(:,:) :: x1, x2, y1, y2
    n1 = size(mx,1)
    n2 = size(mx,2)
    allocate( x1(n1,n2) ); allocate( x2(n1,n2) );
    allocate( y1(n1,n2) ); allocate( y2(n1,n2) );
    if (.not.allocated(mx) .or. .not.allocated(my) .or. .not.allocated(jcb) .or. &
         .not.allocated(xi_x) .or. .not.allocated(xi_y) .or. .not.allocated(et_x) .or. .not.allocated(et_y) ) then
       call mirana_exception("Error in module Mirana.nmd: array not allocated!")
    end if
    call d1(mx,hxi,x1)
    call d2(mx,het,x2)
    call d1(my,hxi,y1)
    call d2(my,het,y2)
    call mbc_d(x1,x2,y1,y2)
    jcb = x1 * y2 - x2 * y1
    xi_x = y2 / jcb
    xi_y = - x2 / jcb
    et_x = - y1 / jcb
    et_y = x1 / jcb
    deallocate( x1 ); deallocate( x2 )
    deallocate( y1 ); deallocate( y2 )
    return
  end subroutine nmd

  !> (Mesh Boundary Condition for Derivatives) Set the boundary condition of the moving mesh.
  !! @note Be sure to call Mirana::d1 and Mirana::d2 on the mesh before calling this subroutine!
  !! @param x1, N1-by-N2 real array, resulted from call d1(mx,hxi,x1).
  !! @param x2, N1-by-N2 real array, resulted from call d2(mx,het,x2).
  !! @param y1, N1-by-N2 real array, resulted from call d1(my,hxi,y1).
  !! @param y2, N1-by-N2 real array, resulted from call d2(my,het,y2).
  subroutine mbc_d(x1,x2,y1,y2)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    double precision, allocatable, dimension(:,:), intent(inout) :: x1,x2,y1,y2
    n1 = size(x1,1)
    n2 = size(x1,2)
    if (.not.allocated(x1) .or. .not.allocated(x2) .or. &
         .not.allocated(y1) .or. .not.allocated(y2) ) then
       call mirana_exception("Error in module Mirana.mbc_d: array not allocated!")
    end if
    !> - Use uniform mesh at left/right boundary.
    x1( (/1,n1/), 1:n2 ) = x1( (/2,n1-1/), 1:n2 )
    y1( (/1,n1/), 1:n2 ) = 0.0
    !> - \f$\xi\f$-mesh is parallel to the wall at floor/ceiling.
    x2( 1:n1, (/1,n2/) ) = 0.0
    y2( 1:n1, (/1,n2/) ) = y2( 1:n1, (/2,n2-1/) )
    return
  end subroutine mbc_d
  
end Module Mirana
