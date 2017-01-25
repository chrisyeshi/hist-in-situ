! include 'mpif.h'
program histogram
  use topology_m
  use histogram_m
  use profiler_m
  use variables_m
  use tracer_m
  use chemkin_m

  ! initialize some testing data
  integer :: x, y, z
  species_name(1) = "H1"
  species_name(2) = "HE2"
  species_name(3) = "HER3"
  species_name(4) = "HERE4"
  do x = 1, 12
    do y = 1, 12
      do z = 1, 12
        u(x,y,z,1) = x
        u(x,y,z,2) = y
        u(x,y,z,3) = z
        yspecies(x,y,z,1) = 100.0
        yspecies(x,y,z,2) = 200.0
        yspecies(x,y,z,3) = 300.0
        yspecies(x,y,z,4) = 400.0
      enddo
    enddo
  end do

  ! make data directory for output
  call execute_command_line('mkdir -p ../data/tracer-0.0000E+00/')

  ! initialize MPI
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( gcomm, myid )

  !---------------------------------------------------
  ! initialization block in solve_driver
  call initialize_histogram( 6 )
  call profile_clear()
  !---------------------------------------------------

  !---------------------------------------------------
  ! (optional) printing the histogram configuration
  call print_histogram_configs( 6 )
  !---------------------------------------------------

  !---------------------------------------------------
  ! within the simulation loop
  ! do while (i_time < i_time_end)
    call profile_start("simulation")
    ! integrate()
    call profile_start("histogram")
    call generate_and_output_histograms( 6 )
    call profile_pause("histogram")
    call profile_start("tracer")
    call write_sorted_tracer_savefile( 6 )
    call profile_pause("tracer")
    ! other tasks within the simulation loop
    call profile_pause("simulation")
    call profile_mpi_advance()
  ! enddo
  !---------------------------------------------------

  !---------------------------------------------------
  ! profiler output
  call profile_mpi_appendfile()

  call MPI_FINALIZE( ierr )
end program histogram
