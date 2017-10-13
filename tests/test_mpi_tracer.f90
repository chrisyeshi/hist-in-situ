program test_tracer
  use topology_m
  use histogram_m
  use param_m
  use tracer_m
  use reference_m
  integer i
  integer, dimension(32) :: offsets, correct_offsets
  character id_ext*5, time_ext*10, dirname*200, filename*200

  call mpi_init(ierr)
  call mpi_comm_rank(gcomm, myid)
  call mpi_comm_rank(gcomm, yid)
  call mpi_comm_size(gcomm, ypes)

  nx = 10
  ny = 10
  nz = 10
  nhx = 2
  nhy = 2
  nhz = 2
  fill = 3

  loc(1,1) = 8.5
  loc(1,2) = 8.5
  loc(1,3) = 8.5
  ssn(1) = 8
  age(1) = 8.0
  dob(1) = 8.0
  state(1) = 8

  loc(2,1) = 2.5
  loc(2,2) = 2.5
  loc(2,3) = 2.5
  ssn(2) = 1
  age(2) = 1.0
  dob(2) = 1.0
  state(2) = 1

  loc(3,1) = 8.5
  loc(3,2) = 2.5
  loc(3,3) = 8.5
  ssn(3) = 6
  age(3) = 6.0
  dob(3) = 6.0
  state(3) = 6

  call execute_command_line('mkdir -p ../data/tracer-0.0000E+00/')
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
  if (ssn(2) .ne. 6) then
    stop 1
  end if
  if (ssn(3) .ne. 8) then
    stop 1
  end if

  correct_offsets = (/ 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 2, 0, 2, 1, &
                       3, 1, 4, 0, 4, 0, 4, 0, 4, 0, 4, 1, 5, 0, 5, 1 /)

  if (yid.eq.0) then
    write(time_ext, '(1pe10.4)') time * time_ref
    write(id_ext, '(I5.5)') xz_id
    dirname = '../data/'//'tracer-'//trim(time_ext)//'/'
    filename = trim(dirname)//'sortedtraceroffsets.'//trim(id_ext)
    open(unit=840, file=trim(filename), status='unknown', form='unformatted')
    read(840) offsets
    do i = 1,32
      if (correct_offsets(i) .ne. offsets(i)) stop 1
    enddo
  end if

  call mpi_finalize(ierr)
end program test_tracer
