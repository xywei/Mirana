  ! Test for differentiation
Program example_1
  impicit none

  real(8), dimension(:,:),allocatable :: Phi
  allocate(Phi, (1:10,1:10))


  deallocate( Phi )
  
end Program example_1

