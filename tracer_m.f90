module tracer_m
  implicit none
  private
  public write_sorted_tracer_savefile

  integer, parameter :: ARRAY_SIZE = 1000
  integer :: fill = 0
  integer(kind=8), dimension(ARRAY_SIZE) :: ssn
  real, dimension(ARRAY_SIZE, 3) :: loc

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