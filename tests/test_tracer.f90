program test_tracer
  use topology_m
  use histogram_m
  use param_m
  use tracer_m
  integer i

  call mpi_init(ierr)
  call mpi_comm_rank(gcomm, myid)

  nx = 20
  ny = 20
  nz = 20
  nhx = 4
  nhy = 4
  nhz = 4
  fill = 4

  loc(1,1) = 8.5
  loc(1,2) = 8.5
  loc(1,3) = 8.5
  ssn(1) = 22
  age(1) = 22.0
  dob(1) = 22.0
  state(1) = 22

  loc(2,1) = 2.5
  loc(2,2) = 2.5
  loc(2,3) = 2.5
  ssn(2) = 1
  age(2) = 1.0
  dob(2) = 1.0
  state(2) = 1

  loc(3,1) = 12.5
  loc(3,2) = 8.5
  loc(3,3) = 8.5
  ssn(3) = 23
  age(3) = 23.0
  dob(3) = 23.0
  state(3) = 23

  loc(4,1) = 2.5
  loc(4,2) = 18.5
  loc(4,3) = 2.5
  ssn(4) = 13
  age(4) = 13.0
  dob(4) = 13.0
  state(4) = 13

  call sort_by_histograms_and_output_offsets(6)

  ! write(*,*) "ssn:", ssn(1:fill)
  ! write(*,*) "age:", age(1:fill)
  ! write(*,*) "dob:", dob(1:fill)
  ! write(*,*) "state:", state(1:fill)
  ! do i = 1,fill
  !   write(*,*) "loc:", loc(i,:)
  ! enddo

  if (ssn(1) .ne. 1) then
    stop 1
  end if
  if (ssn(2) .ne. 13) then
    stop 1
  end if
  if (ssn(3) .ne. 22) then
    stop 1
  end if
  if (ssn(4) .ne. 23) then
    stop 1
  end if

  call mpi_finalize(ierr)
end program test_tracer
