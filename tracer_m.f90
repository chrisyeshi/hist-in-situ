module tracer_m
  implicit none
  private
  public sort_by_histograms_and_output_offsets, write_sorted_tracer_savefile, &
      fill, ssn, age, dob, state, loc, trace_save_fctr

  integer, parameter :: ARRAY_SIZE = 1000
  integer :: fill = 0
  integer(kind=8), dimension(ARRAY_SIZE) :: ssn
  real, dimension(ARRAY_SIZE) :: age
  real, dimension(ARRAY_SIZE) :: dob
  integer, dimension(ARRAY_SIZE) :: state
  real, dimension(ARRAY_SIZE, 3) :: loc
  real trace_save_fctr

  contains

  subroutine mpi_make_directory(dirname)
    implicit none
    character(len=*), intent(in) :: dirname
  end subroutine mpi_make_directory

  subroutine interpolate(field, field_tracer)
    use param_m, only: nx, ny, nz
      implicit none
      real(kind=8), intent(in) :: field(nx, ny, nz)
      real, intent(out) :: field_tracer(1:fill)
  end subroutine interpolate

  subroutine grid_coords(loc, xloc)
      implicit none
      real, intent(in) :: loc(3)
      real, intent(out) :: xloc(3)
  end subroutine grid_coords

  subroutine sort_by_histograms_and_output_offsets(io)
    use histogram_m
    use param_m, only: nx,ny,nz
    use topology_m, only: ypes, xz_id, ycomm, myid, yid, ierr
    use runtime_m, only: time
    use reference_m, only: time_ref
    implicit none
    include 'mpif.h'
    integer io
    integer :: ihists(3), ngrid_hists(3), nhist, ihist, i, j
    integer, dimension(:), allocatable :: histidx, desttraceridx
    integer, dimension(:), allocatable :: ntracersperhist, offsetperhist, &
        filledhists, displ, recvcounts
    integer, dimension(:), allocatable :: offsets, offsetsg
    integer(kind=8), dimension(:), allocatable :: newssn
    real, dimension(:), allocatable :: newage, newdob
    integer, dimension(:), allocatable :: newstate
    real, dimension(:,:), allocatable :: newloc
    character dirname*200, filename*200, id_ext*5, time_ext*10
    ! ssn, age, dob, state, loc
    ngrid_hists(1) = nx / nhx
    ngrid_hists(2) = ny / nhy
    ngrid_hists(3) = nz / nhz
    nhist = nhx * nhy * nhz
    allocate(ntracersperhist(nhist))
    allocate(offsetperhist(nhist))
    allocate(filledhists(nhist))
    allocate(offsets(2*nhist))
    allocate(displ(ypes))
    allocate(recvcounts(ypes))
    allocate(histidx(fill))
    allocate(desttraceridx(fill))
    allocate(newssn(fill))
    allocate(newage(fill))
    allocate(newdob(fill))
    allocate(newstate(fill))
    allocate(newloc(fill,3))
    allocate(offsetsg(2*nhist*ypes))
    ! get sampling region indexes from tracer locations
    ntracersperhist(:) = 0
    do i = 1,fill
      ihists(1) = int(mod((int(loc(i,1))-1),nx)/ngrid_hists(1))
      ihists(2) = int(mod((int(loc(i,2))-1),nx)/ngrid_hists(2))
      ihists(3) = int(mod((int(loc(i,3))-1),nx)/ngrid_hists(3))
      ihist = ihists(1) + ihists(2) * nhx + ihists(3) * nhy * nhz + 1
      histidx(i) = ihist;
      ntracersperhist(ihist) = ntracersperhist(ihist) + 1
    enddo
    ! write(io,*) "histidx:", histidx
    ! write(io,*) "ntracersperhist:", ntracersperhist
    offsetperhist(1) = 1
    do i = 2,nhist
      offsetperhist(i) = ntracersperhist(i - 1) + offsetperhist(i - 1)
    enddo
    ! write(io,*) "offsetperhist:", offsetperhist
    ! compute the new array index for each tracer
    filledhists(:) = 0
    do i = 1,fill
      desttraceridx(i) = offsetperhist(histidx(i)) + filledhists(histidx(i))
      filledhists(histidx(i)) = filledhists(histidx(i)) + 1
    enddo
    ! write(io,*) "desttraceridx:", desttraceridx
    ! move to the new array indexes
    do i = 1,fill
      newssn(desttraceridx(i)) = ssn(i)
      newage(desttraceridx(i)) = age(i)
      newdob(desttraceridx(i)) = dob(i)
      newstate(desttraceridx(i)) = state(i)
      newloc(desttraceridx(i),1) = loc(i,1)
      newloc(desttraceridx(i),2) = loc(i,2)
      newloc(desttraceridx(i),3) = loc(i,3)
    enddo
    ! copy back to the original array
    do i = 1,fill
      ssn(i) = newssn(i)
      age(i) = newage(i)
      dob(i) = newdob(i)
      state(i) = newstate(i)
      loc(i,1) = newloc(i,1)
      loc(i,2) = newloc(i,2)
      loc(i,3) = newloc(i,3)
    enddo
    ! write(io,*) offsetperhist
    ! write(io,*) ntracersperhist
    ! the local offsets array
    do i = 1,nhist
      offsets(2 * i - 1) = offsetperhist(i) - 1
      offsets(2 * i - 0) = ntracersperhist(i)
    enddo
    ! write(io,*) offsets
    displ(1) = 0
    do i = 2,ypes
      displ(i) = displ(i - 1) + 2 * nhist
    enddo
    ! write(io,*) 2*nhist*ypes
    ! write(io,*) displ
    recvcounts(:) = 2 * nhist
    ! gather all the offset arrays
    call mpi_gatherv(offsets, 2 * nhist, MPI_INTEGER, offsetsg, recvcounts, displ, MPI_INTEGER, 0, ycomm, ierr)
    ! increment the gathered offsets
    do i = 2,ypes
      do j = 1,nhist
        offsetsg((i - 1) * 2 * nhist + 2 * (j - 1) + 1) = &
            offsetsg((i - 1) * 2 * nhist + 2 * (j - 1) + 1) + &
            offsetsg((i - 1) * 2 * nhist - 1) + &
            offsetsg((i - 1) * 2 * nhist)
      enddo
    enddo
    ! output the offsetsg array
    if (yid.eq.0) then
    !   write(io,*) offsetsg
      write(time_ext, '(1pe10.4)') time * time_ref
      write(id_ext, '(I5.5)') xz_id
      dirname = '../data/'//'tracer-'//trim(time_ext)//'/'
      filename = trim(dirname)//'sortedtraceroffsets.'//trim(id_ext)
      open(unit=840, file=trim(filename), status='unknown', form='unformatted')
      write(840) offsetsg
      close(unit=840)
    end if
    ! deallocate
    deallocate(offsetsg)
    deallocate(newloc)
    deallocate(newstate)
    deallocate(newdob)
    deallocate(newage)
    deallocate(newssn)
    deallocate(desttraceridx)
    deallocate(histidx)
    deallocate(recvcounts)
    deallocate(displ)
    deallocate(offsets)
    deallocate(filledhists)
    deallocate(offsetperhist)
    deallocate(ntracersperhist)
  end subroutine sort_by_histograms_and_output_offsets

  subroutine write_sorted_tracer_savefile(io)
    ! ==============================================================================
    !
    ! write_sorted_tracer_savefile
    ! ---------------
    !
    ! Author:
    ! G. Borghesi -- SNL
    !
    ! History:
    ! 02-March-2016 -- Created
    !
    ! Description:
    ! Wrapper subroutine for outputting tracer data. Use the following data:
    !
    !   1. loc:  double dimensional array of size [ARRAY_SIZE][3]. Holds particle
    !            grid locations. Data type is double (8 bytes floating point)
    !   2. myid: MPI processor ID
    !   3. ssn:  one dimensional array of size [ARRAY_SIZE]. Holds particle ID
    !            number. Data type is long integer (8 bytes integer)
    !   3. state: one dimensional array of size [ARRAY_SIZE]. Holds particle state.
    !             Data type is 4 bytes integer
    !
    ! ==============================================================================

    ! ------------------------------------------------------------------------------
    use histogram_m
    use topology_m, only: myid
    use param_m, only: nx,ny,nz
    use reference_m, only: l_ref,t_ref,time_ref
    use runtime_m, only : time
    use variables_m, only: u,yspecies,temp
    use profiler_m
    ! ------------------------------------------------------------------------------

    ! ------------------------------------------------------------------------------
    implicit none
    integer io
    character time_ext*10, dirname*200
    ! ------------------------------------------------------------------------------

    integer i,ntracer,nfield
    integer NMAXFIELD
    parameter (NMAXFIELD = 20)

    real :: simtime
    real, dimension(fill,3) :: xloc
    real, dimension(fill,NMAXFIELD) :: sfield

    ! ==============================================================================

    ! Initialization
    nfield = 0

    ! make directory
    write(time_ext, '(1pe10.4)')time*time_ref
    dirname = '../data/'//'tracer-'//trim(time_ext)//'/'
    call mpi_make_directory(dirname)

    call profile_start("newtracer internal")
    ! Absolute tracer position
    do i = 1,fill
      ! write(io,*)loc(i,:)
      call grid_coords(loc(i,:),xloc(i,:))
      xloc(i,:) = xloc(i,:)*l_ref
    enddo

    ! Template code for Eulerian-to-Lagrangian interpolation
    nfield = nfield + 1
    if (nfield.le.NMAXFIELD) call interpolate(temp,sfield(:,nfield))
    sfield(:,nfield) = sfield(:,nfield)*t_ref ! Adimensional field to dimensional one

    ! Call to C function
    call tracerssortedoutput(time*time_ref,myid,nhx,nhy,nhz,nx,ny,nz,ssn,loc,ARRAY_SIZE,xloc,sfield,fill,nfield)
    call profile_pause("newtracer internal")

    return
  end subroutine write_sorted_tracer_savefile

  ! ===============================================================================

end module tracer_m
