module histogram_m
  implicit none
  include "mpif.h"
  private
  public initialize_histogram, print_histogram_configs, generate_and_output_histograms
  public nhx, nhy, nhz

  !----------------------------------------
  !variables for histograms (PDFs)
  integer nhx, nhy, nhz           ! number of histograms in each direction per domain decomposition
  ! a map from string to pointers
  integer, parameter :: limit = 10
  integer :: nHistConfigs
  integer, dimension(limit) :: ndim_hist
  integer, dimension(limit,3) :: nbin_hist
  character, dimension(limit,3) :: vars_hist*30
  character, dimension(limit,3) :: methods_hist*30
  real*8, dimension(limit,3) :: mins_hist
  real*8, dimension(limit,3) :: maxs_hist

  contains

  ! initialize_histogram
  subroutine initialize_histogram( io )
    use topology_m
    implicit none
    integer, intent(in) :: io
    character filename*100
    integer, parameter :: io_inp = 92
    character dim1*100, dim2*100, dim3*100
    integer :: stat, ndim

    filename = "../input/histogram.in"
    ! call inquire_about_input_file(filename, io)
    if(myid .eq. 0) then
      open(unit=io_inp, file=trim(filename), status='old')
      ! read the header 3 lines
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      ! read histogram dimensions
      read(io_inp, *) nhx, nhy, nhz
      ! read the middle comments
      ! TODO: automatically skip comments
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      read(io_inp, *)
      ! read histogram variables
      nHistConfigs = 0
      do
        ! read number of dimension for this histogram configuration
        read(io_inp, *, IOSTAT=stat) ndim
        if (stat /= 0) then
          exit
        endif
        nHistConfigs = nHistConfigs + 1
        ndim_hist(nHistConfigs) = ndim
        ! read nBins per dimension
        if (ndim == 1) then
          read(io_inp, *) nbin_hist(nHistConfigs, 1)
        elseif (ndim == 2) then
          read(io_inp, *) nbin_hist(nHistConfigs, 1), nbin_hist(nHistConfigs, 2)
        elseif (ndim == 3) then
          read(io_inp, *) nbin_hist(nHistConfigs, 1), nbin_hist(nHistConfigs, 2), nbin_hist(nHistConfigs, 3)
        endif
        ! read variable names and their ranges
        if (ndim_hist(nHistConfigs) > 0) then
          read(io_inp, *) vars_hist(nHistConfigs, 1), &
              methods_hist(nHistConfigs, 1), mins_hist(nHistConfigs, 1), &
              maxs_hist(nHistConfigs, 1)
        endif
        if (ndim_hist(nHistConfigs) > 1) then
          read(io_inp, *) vars_hist(nHistConfigs, 2), &
              methods_hist(nHistConfigs, 2), mins_hist(nHistConfigs, 2), &
              maxs_hist(nHistConfigs, 2)
        endif
        if (ndim_hist(nHistConfigs) > 2) then
          read(io_inp, *) vars_hist(nHistConfigs, 3), &
              methods_hist(nHistConfigs, 3), mins_hist(nHistConfigs, 3), &
              maxs_hist(nHistConfigs, 3)
        endif
      enddo
    endif

    ! broadcast the configurations to all processes
    call mpi_bcast(nhx, 1, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(nhy, 1, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(nhz, 1, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(nHistConfigs, 1, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(ndim_hist, limit, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(nbin_hist, limit * 3, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(vars_hist, limit * 3 * 30, MPI_CHARACTER, 0, gcomm, ierr)
    call mpi_bcast(methods_hist, limit * 3 * 30, MPI_CHARACTER, 0, gcomm, ierr)
    call mpi_bcast(mins_hist, limit * 3, MPI_REAL8, 0, gcomm, ierr)
    call mpi_bcast(maxs_hist, limit * 3, MPI_REAL8, 0, gcomm, ierr)

  end subroutine initialize_histogram

  ! print the histogram configurations
  subroutine print_histogram_configs( io )
    implicit none
    integer :: i
    integer, intent(in) :: io

    write(io, *) "nhx =", nhx, "nhy =", nhy, "nhz =", nhz
    do i=1,nHistConfigs
      write(io, *) "nDim=", ndim_hist(i)
      write(io, *) "var1:", vars_hist(i, 1), methods_hist(i, 1), mins_hist(i, 1), maxs_hist(i, 1)
      write(io, *) "var2:", vars_hist(i, 2), methods_hist(i, 2), mins_hist(i, 2), maxs_hist(i, 2)
      write(io, *) "var3:", vars_hist(i, 3), methods_hist(i, 3), mins_hist(i, 3), maxs_hist(i, 3)
      write(io, *) "nbinx=", nbin_hist(i, 1), "nbiny=", nbin_hist(i, 2), "nbinz=", nbin_hist(i, 3)
    enddo
  end subroutine print_histogram_configs

  ! generate all histograms according to the configuration
  subroutine generate_and_output_histograms( io )
    implicit none
    integer, intent(in) :: io
    integer :: i

    do i=1,nHistConfigs
      call generate_and_output_histogram(io, i)
    enddo
  end subroutine generate_and_output_histograms

  ! generate one histogram configuration
  subroutine generate_and_output_histogram( io, iConfig )
    use topology_m, only: myid
    use param_m, only: nx, ny, nz
    use reference_m, only: time_ref
    use runtime_m, only: time
    implicit none
    integer, intent(in) :: io, iConfig
    real*8, dimension(3) :: mins, maxs
    character, dimension(3) :: vars*30, methods*30
    character :: var*30
    integer :: nDim, iVar
    integer, dimension(3) :: nbin
    real*8, allocatable :: data(:,:,:,:)

    nDim = ndim_hist(iConfig)
    mins = mins_hist(iConfig,:)
    maxs = maxs_hist(iConfig,:)
    nbin = nbin_hist(iConfig,:)
    vars = vars_hist(iConfig,:)
    methods = methods_hist(iConfig,:)
    allocate(data(nx,ny,nz,nDim)); data=0.0
    do iVar=1, nDim
      var = vars(iVar)
      call mapStrToVar(io, var, data(:,:,:,iVar))
    enddo
    call generate_and_output_histogram_nd(io, nDim, data, methods, mins, maxs, nbin, iConfig)
    deallocate(data)
  end subroutine generate_and_output_histogram

  ! helper function to map variable names to actual variables
  subroutine mapStrToVar( io, str, ret )
    use variables_m
    implicit none
    integer, intent(in) :: io
    character, intent(in) :: str*30
    real*8, intent(out) :: ret(:,:,:)

    ! TODO: also map the variables in chem.asc
    if (str == 'xvel') then
      ret = u(:,:,:,1)
    elseif (str == 'yvel') then
      ret = u(:,:,:,2)
    elseif (str == 'zvel') then
      ret = u(:,:,:,3)
    elseif (str == 'volum') then
      ret = volum
    elseif (str == 'pressure') then
      ret = pressure
    elseif (str == 'temp') then
      ret = temp
    endif
  end subroutine mapStrToVar

  ! generate and output a histogram configuration of n-dimensional
  subroutine generate_and_output_histogram_nd( io, nDim, data, methods, mins, maxs, nbin, iConfig )
    use topology_m, only: myid, ycomm
    use param_m, only: nx, ny, nz
    use reference_m, only: time_ref
    use runtime_m, only: time
    implicit none
    integer, intent(in) :: io, nDim, nbin(:), iConfig
    character, dimension(3), intent(in) :: methods*30
    real*8, intent(in) :: data(:,:,:,:), mins(:), maxs(:)
    if (nDim == 1) then
      call c_generate_and_output_histogram_1d(data(:,:,:,1), &
                                              0.0, &       ! log base for scaling
                                              methods(1), &
                                              mins(1), &
                                              maxs(1), &
                                              time * time_ref, &
                                              nx, ny, nz, &
                                              nhx, nhy, nhz, &
                                              nbin(1), &
                                              myid, 0, ycomm, &
                                              iConfig)
    elseif (nDim == 2) then
      call c_generate_and_output_histogram_2d(data(:,:,:,1), data(:,:,:,2), &
                                              0.0, 0.0, &       ! log base for scaling
                                              methods(1), methods(2), &
                                              mins(1), mins(2), &
                                              maxs(1), maxs(2), &
                                              time * time_ref, &
                                              nx, ny, nz, &
                                              nhx, nhy, nhz, &
                                              nbin(1), nbin(2), &
                                              myid, 0, ycomm, &
                                              iConfig)
    elseif (nDim == 3) then
      call c_generate_and_output_histogram_3d(data(:,:,:,1), data(:,:,:,2), data(:,:,:,3), &
                                              0.0, 0.0, 0.0, &       ! log base for scaling
                                              methods(1), methods(2), methods(3), &
                                              mins(1), mins(2), mins(3), &
                                              maxs(1), maxs(2), maxs(3), &
                                              time * time_ref, &
                                              nx, ny, nz, &
                                              nhx, nhy, nhz, &
                                              nbin(1), nbin(2), nbin(3), &
                                              myid, 0, ycomm, &
                                              iConfig)
    endif

  end subroutine generate_and_output_histogram_nd

end module histogram_m
