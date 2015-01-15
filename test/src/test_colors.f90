program test_colors

  use colors_def
  use load_bc_tables
  use eval_bc_tables
  
  implicit none
    
  integer :: ierr, i
  type(bc_table) :: t
  real(sp), allocatable :: res(:)
  real(sp) :: logT, logg, Av, Rv

  t% filename = '/home/dotter/science/Spectra/Zp0d0ap0.UBVRIJHKsKp'
  call load_one_bc(t,ierr)
  write(*,*) ' ierr = ', ierr
  write(*,*) 'num_Av = ', t% num_Av
  write(*,*) 'num_filter = ', t% num_filter
  write(*,*) 'num_lines = ', t% num_lines
  write(*,*) 'num_g     = ', t% num_g
  write(*,*) 'num_T     = ', t% num_T

  Av   = 0.111
  Rv   = 2.01
  logT = 3.7777
  logg = 4.4444

  allocate(res(t% num_filter))
  call eval_one_bc(t,logT, logg, Av, Rv, res, ierr)

  write(*,*) res
  
  write(*,*)
  write(*,*) ' ierr = ', ierr
    
end program test_colors
  
