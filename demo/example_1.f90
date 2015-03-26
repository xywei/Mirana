  ! Test for differentiation
Program example_1
  Implicit none

  real(8), dimension(:,:),allocatable :: phi, px, py
  allocate(phi (1:10,1:10))
  allocate(px(2:9,2:9),py(2:9,2:9))

  write(*,*) 'Hello, this is example_1.'

  deallocate( phi, px, py )
  
end Program example_1

