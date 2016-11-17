!===========================================================================
module profiler_m
!===========================================================================
! module for profiler

  implicit none
!---------------------------------------------------------------------------
! integers
  contains

  subroutine profile_clear()
    implicit none
    call profiler_clear()
  end subroutine profile_clear

  subroutine profile_advance()
    use runtime_m, only : time
    use reference_m, only : time_ref
    implicit none
    call profiler_advance(time * time_ref)
  end subroutine profile_advance

  subroutine profile_mpi_advance()
    use topology_m, only : myid, gcomm
    use runtime_m, only : time
    use reference_m, only : time_ref
    implicit none
    call profiler_mpi_advance(time * time_ref, myid, 0, gcomm)
  end subroutine profile_mpi_advance

  subroutine profile_reset(title)
    implicit none
    character(len=*), intent(in) :: title
    call profiler_reset(title//CHAR(0))
  end subroutine profile_reset

  subroutine profile_start(title)
    implicit none
    character(len=*), intent(in) :: title
    call profiler_start(title//CHAR(0))
  end subroutine profile_start

  subroutine profile_pause(title)
    implicit none
    character(len=*), intent(in) :: title
    call profiler_pause(title//CHAR(0))
  end subroutine profile_pause

  subroutine profile_appendfile()
    implicit none
    call profiler_appendfile()
  end subroutine profile_appendfile

  subroutine profile_mpi_appendfile()
    use topology_m, only : myid, gcomm
    implicit none
    call profiler_mpi_appendfile(myid, 0, gcomm)
  end subroutine profile_mpi_appendfile
!---------------------------------------------------------------------------
end module profiler_m
