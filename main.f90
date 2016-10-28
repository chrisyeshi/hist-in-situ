! include 'mpif.h'
program histogram
  use topology_m
  use histogram_m

  use variables_m
  integer :: x, y, z
  do x = 1, 12
    do y = 1, 12
      do z = 1, 12
        u(x,y,z,1) = x
        u(x,y,z,2) = y
        u(x,y,z,3) = z
      enddo
    enddo
  end do

  call MPI_INIT( ierr )
  call initialize_histogram( 6 )
  call print_histogram_configs( 6 )
  call generate_and_output_histograms( 6 )
  call MPI_FINALIZE( ierr )
end program histogram