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
!! @todo test greek()
!! @todo finish nmd1()
Module Mirana
  implicit none

  ! public  ! public members go here.

  ! private ! private members go here.
    
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

  !> Write mesh data into .dat file for plot
  !! @param mx real N1-by-N2 matrix, x-coodinates.
  !! @param my real N1-by-N2 matrix, x-coodinates.
  !! @note After generating mesh.dat, call the following command in Gnuplot:
  !! - set view 0,0
  !! - splot "mesh.dat" using 1:2:3 with lines
  subroutine save_mesh(mx,my)
    implicit none
    double precision, allocatable, dimension(:,:), intent(in) :: mx,my
    integer :: n1,n2,i,j
    n1 = size(mx,1)
    n2 = size(mx,2)
    open(unit=233,file="mesh.dat",status="replace")
    do i = 1,n1
       do j = 1,n2
          write(233,*) mx(i,j), my(i,j), 0.0
       end do
       write(233,*) " "
    end do
  end subroutine save_mesh
  
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

  !> Compute the derivative along the 1st direction at the 1st kind of half points using central difference.
  !! @param phi real N1-by-N2 matrix.
  !! @param h real number, step size.
  !! @return phi_1h, real N1-by-N2 matrix, the (i,j) element stores the derivative at (i+1/2,j), only
  !!                 [1:N1-1]X[1:N2] part is used.
  subroutine d1h1(phi,h,phi_1h)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_1h
    if (.not.allocated(phi) .or. .not.allocated(phi_1h)) then
       call mirana_exception("Error in module Mirana.d1h: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    phi_1h(1:(n1-1),1:n2) = ( phi(2:n1,1:n2) - phi(1:(n1-1),1:n2) ) / h
    return
  end subroutine d1h1

  !> Compute the derivative along the 1st direction at the 2nd kind of half points using central difference.
  !! @param phi real N1-by-N2 matrix.
  !! @param h real number, step size (h1).
  !! @return phi_1h, real N1-by-N2 matrix, the (i,j) element stores the derivative at (i,j+1/2), only
  !!                 [2:N1-1]X[1:N2-1] part is used.
  subroutine d1h2(phi,h,phi_1h)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_1h
    double precision, allocatable, dimension(:,:) :: phi_x
    if (.not.allocated(phi) .or. .not.allocated(phi_1h)) then
       call mirana_exception("Error in module Mirana.d1h: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    allocate( phi_x(n1,n2) )
    call d1(phi,h,phi_x)
    phi_1h(2:(n1-1),1:(n2-1)) = 0.5 * ( phi_x(2:(n1-1),1:(n2-1)) + phi_x(2:(n1-1),2:n2) )
    deallocate( phi_x )
    return
  end subroutine d1h2

  !> Compute the derivative along the 2nd direction at the 1st kind of half points using central difference.
  !! @param phi real N1-by-N2 matrix.
  !! @param h real number, step size (h2).
  !! @return phi_2h, real N1-by-N2 matrix, the (i,j) element stores the derivative at (i+1/2,j), only
  !!                 [1:N1-1]X[2:N2-1] part is used.
  subroutine d2h1(phi,h,phi_2h)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_2h
    double precision, allocatable, dimension(:,:) :: phi_y
    if (.not.allocated(phi) .or. .not.allocated(phi_2h)) then
       call mirana_exception("Error in module Mirana.d1h: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    allocate( phi_y(n1,n2) )
    call d2(phi,h,phi_y)
    phi_2h(1:(n1-1),2:(n2-1)) = 0.5 * ( phi_y(1:(n1-1),2:(n2-1)) + phi_y(2:n1,2:(n2-1)) )
    deallocate( phi_y )
    return
  end subroutine d2h1

  !> Compute the derivative along the 2nd direction at the 2nd kind of half points using central difference.
  !! @param phi real N1-by-N2 matrix.
  !! @param h real number, step size.
  !! @return phi_2h, real N1-by-N2 matrix, the (i,j) element stores the derivative at (i,j+1/2), only
  !!                 [1:N1]X[1:N2-1] part is used.
  subroutine d2h2(phi,h,phi_2h)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    double precision, intent(in) :: h
    double precision, allocatable, dimension(:,:), intent(in) :: phi
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_2h
    if (.not.allocated(phi) .or. .not.allocated(phi_2h)) then
       call mirana_exception("Error in module Mirana.d2h: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    phi_2h(1:n1,1:(n2-1)) = ( phi(1:n1,2:n2) - phi(1:n1,1:(n2-1)) ) / h
    return
  end subroutine d2h2  

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
  !! @note nmd() and greek() computes all fundamental coefficients in non-conservative form.
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

  !> (Non-conservative Mesh D1) Compute \f$\partial\phi/\partial x\f$ using mesh variables
  !! @param phi N1-by-N2 real matrix.
  !! @param xix N1-by-N2 real matrix, \f$\partial\xi\partial x\f$.
  !! @param etx N1-by-N2 real matrix, \f$\partial\eta\partial x\f$.
  !! @param h1 real number, step size in \f$\xi\f$ axis.
  !! @param h2 real number, step size in \f$\eta\f$ axis.
  !! @return phi_x N1-by-N2 real matrix, only [2:N1-1]X[2:N2-1] are used.
  subroutine nmd1(phi,xix,etx,h1,h2,phi_x)
    implicit none
    integer :: n1,n2
    double precision, intent(in) :: h1,h2
    double precision, allocatable, dimension(:,:), intent(in) :: phi,xix,etx
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_x
    double precision, allocatable, dimension(:,:) :: phi1,phi2
    if ( .not.allocated(phi) .or. .not.allocated(xix) .or. &
         .not.allocated(etx) .or. .not.allocated(phi_x) ) then
       call mirana_exception("Error in module Mirana.nmd1: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    allocate( phi1(n1,n2) )
    allocate( phi2(n1,n2) )
    call d1(phi,h1,phi1)
    call d2(phi,h2,phi2)
    phi_x = phi1 * xix + phi2 * etx
    deallocate( phi1 )
    deallocate( phi2 )
    return
  end subroutine nmd1

  !> (Non-conservative Mesh D2) Compute \f$\partial\phi/\partial y\f$ using mesh variables
  !! @param phi N1-by-N2 real matrix.
  !! @param xiy N1-by-N2 real matrix, \f$\partial\xi\partial y\f$.
  !! @param ety N1-by-N2 real matrix, \f$\partial\eta\partial y\f$.
  !! @param h1 real number, step size in \f$\xi\f$ axis.
  !! @param h2 real number, step size in \f$\eta\f$ axis.
  !! @return phi_y N1-by-N2 real matrix, only [2:N1-1]X[2:N2-1] are used.
  subroutine nmd2(phi,xiy,ety,h1,h2,phi_y)
    implicit none
    integer :: n1,n2
    double precision, intent(in) :: h1,h2
    double precision, allocatable, dimension(:,:), intent(in) :: phi,xiy,ety
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_y
    double precision, allocatable, dimension(:,:) :: phi1,phi2
    if ( .not.allocated(phi) .or. .not.allocated(xiy) .or. &
         .not.allocated(ety) .or. .not.allocated(phi_y) ) then
       call mirana_exception("Error in module Mirana.nmd2: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    allocate( phi1(n1,n2) )
    allocate( phi2(n1,n2) )
    call d1(phi,h1,phi1)
    call d2(phi,h2,phi2)
    phi_y = phi1 * xiy + phi2 * ety
    deallocate( phi1 )
    deallocate( phi2 )
    return
  end subroutine nmd2

  !> (Non-conservative Mesh DLP) Compute \f$\Delta\phi\f$ using mesh variables
  !! @param phi N1-by-N2 real matrix.
  !! @param xix N1-by-N2 real matrix, \f$\partial\xi\partial x\f$.
  !! @param xiy N1-by-N2 real matrix, \f$\partial\xi\partial y\f$.
  !! @param etx N1-by-N2 real matrix, \f$\partial\eta\partial x\f$.
  !! @param ety N1-by-N2 real matrix, \f$\partial\eta\partial y\f$.
  !! @param alp N1-by-N2 real matrix, obtained from greek().
  !! @param bet N1-by-N2 real matrix, obtained from greek().
  !! @param h1 real number, step size in \f$\xi\f$ axis.
  !! @param h2 real number, step size in \f$\eta\f$ axis.
  !! @return phi_lp N1-by-N2 real matrix, only [2:N1-1]X[2:N2-1] are used.
  subroutine nmdlp(phi,xix,xiy,etx,ety,alp,bet,h1,h2,phi_lp)
    implicit none
    integer :: n1,n2
    double precision, intent(in) :: h1,h2
    double precision, allocatable, dimension(:,:), intent(in) :: phi,alp,bet,xix,xiy,etx,ety
    double precision, allocatable, dimension(:,:), intent(inout) :: phi_lp
    double precision, allocatable, dimension(:,:) :: phi1,phi2,phi11,phi12,phi22
    if ( .not.allocated(phi) .or. .not.allocated(xix) .or. &
         .not.allocated(etx) .or. .not.allocated(xiy) .or. &
         .not.allocated(ety) .or. .not.allocated(alp) .or. &
         .not.allocated(bet) .or. .not.allocated(phi_lp) ) then
       call mirana_exception("Error in module Mirana.nmdlp: array not allocated!")
    end if
    n1 = size(phi,1)
    n2 = size(phi,2)
    allocate( phi1(n1,n2) )
    allocate( phi2(n1,n2) )
    allocate( phi11(n1,n2) )
    allocate( phi12(n1,n2) )
    allocate( phi22(n1,n2) )
    call d1(phi,h1,phi1)
    call d2(phi,h2,phi2)
    call d11(phi,h1,phi11)
    call d12(phi,h1,h2,phi12)
    call d22(phi,h2,phi22)
    !===============================
    phi_lp = (xix**2 + xiy**2) * phi11 + (etx**2 + ety**2) * phi22 + &
         2.0 * (xix*etx + xiy*ety) * phi12 + alp * phi1 + bet * phi2
    !===============================    
    deallocate( phi1 )
    deallocate( phi2 )
    deallocate( phi11 )
    deallocate( phi12 )
    deallocate( phi22 )
    return
  end subroutine nmdlp

  !> Compute the two coefficients for transformed laplacian \f$\alpha,\beta\f$.
  !! @param mx real N1-by-N2 matrix, mesh coordinates - x.
  !! @param my real N1-by-N2 matrix, mesh coordinates - y.
  !! @param xix real N1-by-N2 matrix, \f$\partial\xi/\partial x\f$.
  !! @param xiy real N1-by-N2 matrix, \f$\partial\xi/\partial y\f$.
  !! @param etx real N1-by-N2 matrix, \f$\partial\eta/\partial x\f$.
  !! @param ety real N1-by-N2 matrix, \f$\partial\eta/\partial y\f$.
  !! @param h1 real number, step size in the first direction.
  !! @param h2 real number, step size in the second direction.
  !! @return alp real N1-by-N2 matrix, the \f$\alpha\f$, only [2:N1-1]X[2:N2-1] is used.
  !! @return bet real N1-by-N2 matrix, the \f$\beta\f$, only [2:N1-1]X[2:N2-1] is used.
  subroutine greek(mx,my,xix,xiy,etx,ety,h1,h2,alp,bet)
    implicit none
    integer :: n1,n2 !> dimension of arrays
    integer, allocatable, dimension(:) :: kk,ll
    double precision, intent(in) :: h1,h2
    integer :: i
    double precision, allocatable, dimension(:,:), intent(in) :: mx,my,xix,xiy,etx,ety
    double precision, allocatable, dimension(:,:), intent(inout) :: alp,bet
    double precision, allocatable, dimension(:,:) :: jh1,jh2,&
         x1h1,x1h2,x2h1,x2h2,y1h1,y1h2,y2h1,y2h2,&
         xixh1,xixh2,xiyh1,xiyh2,etxh1,etxh2,etyh1,etyh2
    n1 = size(mx,1)
    n2 = size(mx,2)
    allocate( kk(n1-2) )
    allocate( ll(n2-2) )
    kk = (/ (i, i=2,(n1-1), 1) /)
    ll = (/ (i, i=2,(n2-1), 1) /)
    if (.not.allocated(mx) .or. .not.allocated(my) .or. .not.allocated(xix) .or. &
         .not.allocated(xiy) .or. .not.allocated(etx) .or. .not.allocated(ety) .or. &
         .not.allocated(alp) .or. .not.allocated(bet) ) then
       call mirana_exception("Error in module Mirana.greek: array not allocated!")
    end if
    allocate( jh1(n1,n2) )
    allocate( jh2(n1,n2) )
    allocate( x1h1(n1,n2) )
    allocate( x1h2(n1,n2) )
    allocate( x2h1(n1,n2) )
    allocate( x2h2(n1,n2) )
    allocate( y1h1(n1,n2) )
    allocate( y1h2(n1,n2) )
    allocate( y2h1(n1,n2) )
    allocate( y2h2(n1,n2) )
    !=====================================
    call d1h1(mx,h1,x1h1)
    call d1h2(mx,h1,x1h2)
    call d2h1(mx,h2,x2h1)
    call d2h2(mx,h2,x2h2)
    call d1h1(my,h1,y1h1)
    call d1h2(my,h1,y1h2)
    call d2h1(my,h2,y2h1)
    call d2h2(my,h2,y2h2)
    !=====================================
    jh1 = x1h1 * y2h1 - y1h1 * x2h1
    jh2 = x1h2 * y2h2 - y1h2 * x2h2
    !=====================================
    allocate( xixh1(n1,n2) )
    allocate( xiyh1(n1,n2) )
    allocate( xixh2(n1,n2) )
    allocate( xiyh2(n1,n2) )
    allocate( etxh1(n1,n2) )
    allocate( etyh1(n1,n2) )
    allocate( etxh2(n1,n2) )
    allocate( etyh2(n1,n2) )
    xixh1(1:n1-1,ll) = + y2h1(1:n1-1,ll) / jh1(1:n1-1,ll) ! the index must be restriced to avoid division by zero!
    xiyh1(1:n1-1,ll) = - x2h1(1:n1-1,ll) / jh1(1:n1-1,ll)
    xixh2(kk,1:n2-1) = + y2h2(kk,1:n2-1) / jh2(kk,1:n2-1) 
    xiyh2(kk,1:n2-1) = - x2h2(kk,1:n2-1) / jh2(kk,1:n2-1) 
    etxh1(1:n1-1,ll) = - y1h1(1:n1-1,ll) / jh1(1:n1-1,ll)
    etyh1(1:n1-1,ll) = + x1h1(1:n1-1,ll) / jh1(1:n1-1,ll)
    etxh2(kk,1:n2-1) = - y1h2(kk,1:n2-1) / jh2(kk,1:n2-1) 
    etyh2(kk,1:n2-1) = + x1h2(kk,1:n2-1) / jh2(kk,1:n2-1)
    !=====================================
    alp(kk,ll) = xix(kk,ll)*( xixh1(kk,ll) - xixh1(kk-1,ll) ) / h1 + xiy(kk,ll)*( xiyh1(kk,ll) - xiyh1(kk-1,ll) ) / h1 &
         + etx(kk,ll)*( xixh2(kk,ll) - xixh2(kk,ll-1) ) / h2 + ety(kk,ll)*( xiyh2(kk,ll) - xiyh2(kk,ll-1) ) / h2
    bet(kk,ll) = xix(kk,ll)*( etxh1(kk,ll) - etxh1(kk-1,ll) ) / h1 + xiy(kk,ll)*( etyh1(kk,ll) - etyh1(kk-1,ll) ) / h1 &
         + etx(kk,ll)*( etxh2(kk,ll) - etxh2(kk,ll-1) ) / h2 + ety(kk,ll)*( etyh2(kk,ll) - etyh2(kk,ll-1) ) / h2
    !=====================================
    deallocate( jh1 )
    deallocate( jh2 )
    deallocate( x1h1 )
    deallocate( x1h2 )
    deallocate( x2h1 )
    deallocate( x2h2 )
    deallocate( y1h1 )
    deallocate( y1h2 )
    deallocate( y2h1 )
    deallocate( y2h2 )
    deallocate( xixh1 )
    deallocate( xixh2 )
    deallocate( xiyh1 )
    deallocate( xiyh2 )
    deallocate( etxh1 )
    deallocate( etxh2 )
    deallocate( etyh1 )
    deallocate( etyh2 )
    return
  end subroutine greek
  

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
