module variables_m
  implicit none
  real*8, dimension(12,12,12,4,1) :: q
  real*8, dimension(12,12,12,4) :: yspecies
  real*8, dimension(12,12,12,3) :: u
  real*8, dimension(12,12,12) :: volum
  real*8, dimension(12,12,12) :: pressure
  real*8, dimension(12,12,12) :: temp
end module variables_m

module param_m
  implicit none
  integer :: nx = 12, ny = 12, nz = 12, n_spec = 4
end module param_m

module reference_m
  implicit none
  real*8 a_ref        
  real*8 l_ref        
  real*8 rho_ref      
  real*8 lambda_ref   
  real*8 mu_ref       
  real*8 t_ref        
  real*8 t_o          
  real*8 p_ref        
  real*8 time_ref
  real*8 cp_ref       
  real*8 univ_gascon  
  real*8 g_ref        
  real*8 pres_atm     
  real*8 re      
  real*8 re_real 
  real*8 mach_no 
  real*8 rr_ref  
  real*8 hr_ref  
end module reference_m

module runtime_m
  implicit none
  real*8 :: time
  integer :: i_time
end module runtime_m

module topology_m
  implicit none
  include 'mpif.h'
  integer :: myid = 0, yid = 0, ypes = 5, xz_id = 0
  integer :: gcomm = MPI_COMM_WORLD, ycomm = MPI_COMM_WORLD
  integer :: ierr

  contains
  subroutine makedirectory(path)
    implicit none
    character, dimension(200) :: path
  end subroutine makedirectory
  subroutine execute_command(command)
    implicit none
    character, dimension(200) :: command
  end subroutine execute_command

end module topology_m

module mixfrac_m
  implicit none
  real*8, dimension(4,4,4) :: mixfrac

  contains
  subroutine specToMixfr( massfr )
    use param_m, only : nx, ny, nz
    implicit none
    real*8, intent(in), dimension(nx,ny,nz,4) :: massfr
  end subroutine specToMixfr
  subroutine computeScalarDissipationRate(io, chiSDR)
    use param_m, only : nx, ny, nz
    implicit none
    integer, intent(in) :: io
    real*8, intent(out), dimension(nx,ny,nz) :: chiSDR
  end subroutine computeScalarDissipationRate
end module mixfrac_m

subroutine write_header(io,char_in)
  implicit none
  integer io
  character*1 char_in

  return
end subroutine write_header