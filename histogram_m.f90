module histogram_m
  implicit none
  include "mpif.h"
  private
  public initialize_histogram, print_histogram_configs, generate_and_output_histograms
  public nhx, nhy, nhz

  !----------------------------------------
  !variables for histograms (PDFs)
  integer nhx, nhy, nhz           ! number of histograms in each direction per domain decomposition
  integer, parameter :: TOTAL_NUMBER_OF_HISTS = 32768
  ! a map from string to pointers
  integer, parameter :: limit = 10
  integer :: nHistConfigs
  integer, dimension(limit) :: ndim_hist
  character, dimension(limit,3) :: nbin_hist*30
  character, dimension(limit,3) :: vars_hist*30
  character, dimension(limit,3) :: methods_hist*30
  real*8, dimension(limit,3) :: mins_hist
  real*8, dimension(limit,3) :: maxs_hist

  contains

  subroutine mpi_make_directory(dirname)
    use topology_m, only : myid, gcomm, ierr
    implicit none
    character(len=*), intent(in) :: dirname
    if(myid == 0) then
#ifdef SYSTEMCALLWONTWORK
      call makedirectory(trim(dirname)//char(0))
#else
      call execute_command('mkdir -p'//trim(dirname))
#endif
    end if
    ! All processes need to wait for the directory to be created
    call MPI_Barrier(gcomm, ierr)
  end subroutine mpi_make_directory

  ! initialize_histogram
  subroutine initialize_histogram( pdf_save_fctr, io )
    use topology_m
    use param_m
    use runtime_m
    use grid_m
    implicit none
    integer, intent(in) :: io
    real, intent(in) :: pdf_save_fctr
    character filename*100, dirname*100
    integer, parameter :: io_inp = 92
    character dim1*100, dim2*100, dim3*100
    integer :: stat, ndim, ihist, idim

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
        ! read variable names, number of bins, and their ranges
        if (ndim_hist(nHistConfigs) > 0) then
          read(io_inp, *) vars_hist(nHistConfigs, 1), &
              nbin_hist(nHistConfigs, 1), methods_hist(nHistConfigs, 1), &
              mins_hist(nHistConfigs, 1), maxs_hist(nHistConfigs, 1)
        endif
        if (ndim_hist(nHistConfigs) > 1) then
          read(io_inp, *) vars_hist(nHistConfigs, 2), &
              nbin_hist(nHistConfigs, 2), methods_hist(nHistConfigs, 2), &
              mins_hist(nHistConfigs, 2), maxs_hist(nHistConfigs, 2)
        endif
        if (ndim_hist(nHistConfigs) > 2) then
          read(io_inp, *) vars_hist(nHistConfigs, 3), &
              nbin_hist(nHistConfigs, 3), methods_hist(nHistConfigs, 3), &
              mins_hist(nHistConfigs, 3), maxs_hist(nHistConfigs, 3)
        endif
      enddo
      ! create the data_pdf directory
      dirname = "../data_pdf/"
#ifdef SYSTEMCALLWONTWORK
        call makedirectory(trim(dirname)//char(0))
#else
        call execute_command('mkdir -p'//trim(dirname))
#endif
      ! write the pdf.config file in the data_pdf directory
      filename = trim(dirname)//'pdf.config'
      open(unit=630, file=trim(filename))
      write(630,*) 'num grid points', nx_g, ny_g, nz_g
      write(630,*) 'num processors', xpes, ypes, zpes
      write(630,*) 'num timesteps', i_time_end
      write(630,*) 'save field frequency', i_time_save
      write(630,*) 'physical bounding box lower', xmin, ymin, zmin
      write(630,*) 'physical bounding box upper', xmax, ymax, zmax
      write(630,*) 'save tracer frequency', pdf_save_fctr
      write(630,*) 'num pdfs per domain', nhx, nhy, nhz
      write(630,*) 'histograms'
      do ihist=1,nHistConfigs
        write(630,*) 'dimension', ndim_hist(ihist)
        do idim=1,ndim_hist(ihist)
          write(630,*) vars_hist(ihist, idim), nbin_hist(ihist, idim), &
              methods_hist(ihist, idim), mins_hist(ihist, idim), &
              maxs_hist(ihist, idim)
        enddo
      enddo
      close(630)
    endif

    ! broadcast the configurations to all processes
    call mpi_bcast(nhx, 1, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(nhy, 1, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(nhz, 1, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(nHistConfigs, 1, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(ndim_hist, limit, MPI_INTEGER, 0, gcomm, ierr)
    call mpi_bcast(nbin_hist, limit * 3 * 30, MPI_CHARACTER, 0, gcomm, ierr)
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
    use reference_m, only: time_ref
    use runtime_m, only: time
    implicit none
    integer, intent(in) :: io
    integer :: i
    character time_ext*10

    ! create directory
    write(time_ext, '(1pe10.4)') time * time_ref
    call mpi_make_directory('../data_pdf/pdf-'//trim(time_ext)//'/')

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
    character, dimension(3) :: nbin*30
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
    use chemkin_m, only: n_species, species_name
    implicit none
    integer, intent(in) :: io
    character, intent(in) :: str*30
    real*8, intent(out) :: ret(:,:,:)
    integer ispecies

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
    do ispecies=1,n_species
      if (species_name(ispecies) == str) then
        ret = yspecies(:,:,:,ispecies)
      end if
    enddo
  end subroutine mapStrToVar

  ! generate and output a histogram configuration of n-dimensional
  subroutine generate_and_output_histogram_nd( io, nDim, data, methods, mins, maxs, nbin, iConfig )
    use topology_m, only: yid, xz_id, ypes, ycomm
    use param_m, only: nx, ny, nz
    use reference_m, only: time_ref
    use runtime_m, only: time
    implicit none
    integer, intent(in) :: io, nDim, iConfig
    character, dimension(3), intent(in) :: methods*30, nbin*30
    real*8, intent(in) :: data(:,:,:,:), mins(:), maxs(:)
    real*8, dimension(3) :: log_bases

    log_bases(1) = 0.0
    log_bases(2) = 0.0
    log_bases(3) = 0.0

    call c_generate_and_output_histogram( &
        data(:,:,:,:), log_bases(:), methods, nbin, &
        mins, maxs, &
        time * time_ref, &
        nDim, &
        nx, ny, nz, &
        nhx, nhy, nhz, &
        yid, 0, ypes, ycomm, xz_id, iConfig)

  end subroutine generate_and_output_histogram_nd

end module histogram_m
