program test_colors

  use colors_def
  use load_bc_tables
  use eval_bc_tables
  
  implicit none
    
  integer :: ierr, i
  type(bc_table) :: t
  real(sp), allocatable :: res(:)
  real(sp) :: logT, logg, Av

  t% filename = '/home/dotter/science/Spectra/Zp0d0ap0.UBVRIJHKsKp.Rv31'
  call load_one_bc(t,ierr)
  write(*,*) ' ierr = ', ierr
  write(*,*) 'num_Av = ', t% num_Av
  write(*,*) 'num_filter = ', t% num_filter
  write(*,*) 'num_lines = ', t% num_lines
  write(*,*) 'BC(1,1,1) = ', t% BC(1,1,1)
  write(*,*) 'num_g     = ', t% num_g
  write(*,*) 'num_T     = ', t% num_T

  Av   = 0.111
  logT = 3.7777
  logg = 4.4444

  allocate(res(t% num_filter))
  call eval_one_bc(t,logT, logg, Av, res, ierr)

  write(*,*) res
  
  write(*,*)
  write(*,*) ' ierr = ', ierr
    
end program test_colors
  
