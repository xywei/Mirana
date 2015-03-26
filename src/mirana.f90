! The main module
! Doxygen tips:
!   "!>" or "!<" starts a comment
!   "!!" or "!>" can be used to continue an one line comment into a multi-line comment.
! 

!> Provide the subroutines & functions for moving mesh computation.
!! @todo Add differenciation
Module Mirana
  implicit none

contains
  !> Differentiate with respect to \f$x\f$ using central difference.
  !! @param phi real N-by-N matrix, index range 1:N
  !! @return phi_x real (N-2)-by-(N-2) matrix, index range 2:(N-1)
  subroutine ddx
    implicit none

  end subroutine ddx

  !> Differentiate with respect to \f$y\f$ using central difference.
  !! @param phi N-by-N matrix
  !! @return phi_y (N-2)-by(N-2) matrix
  subroutine ddy
    implicit none

  end subroutine ddy

end Module Mirana
