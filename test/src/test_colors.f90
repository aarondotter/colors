program test_colors

  use colors_def
  use load_bc_tables
  use eval_bc_tables
  
  implicit none
    
  integer :: ierr, i, n
  type(bc_table), allocatable :: t(:)
  real(sp), allocatable :: res(:), mags(:)
  real(sp) :: logT, logg, logL, Av, Rv

  n=9
  allocate(t(n))

  t(1)% filename = '/home/dotter/science/colors/data/feh-2.1.FSPS'
  t(2)% filename = '/home/dotter/science/colors/data/feh-1.8.FSPS'
  t(3)% filename = '/home/dotter/science/colors/data/feh-1.5.FSPS'
  t(4)% filename = '/home/dotter/science/colors/data/feh-1.2.FSPS'
  t(5)% filename = '/home/dotter/science/colors/data/feh-0.9.FSPS'
  t(6)% filename = '/home/dotter/science/colors/data/feh-0.6.FSPS'
  t(7)% filename = '/home/dotter/science/colors/data/feh-0.3.FSPS'
  t(8)% filename = '/home/dotter/science/colors/data/feh+0.0.FSPS'
  t(9)% filename = '/home/dotter/science/colors/data/feh+0.3.FSPS'

  do i=8,8
     call load_one_bc(t(i),ierr)
  enddo
  allocate(res(t(8)% num_filter),mags(t(8)% num_filter))
  
  do i=8,8
     write(*,*) ' i    = ', i
     write(*,*) '[Fe/H]= ', t(i)% FeH
     write(*,*) ' ierr = ', ierr
     write(*,*) 'num_Av = ', t(i)% num_Av
     write(*,*) 'num_filter = ', t(i)% num_filter
     write(*,*) 'num_lines = ', t(i)% num_lines
     write(*,*) 'num_g     = ', t(i)% num_g
     write(*,*) 'num_T     = ', t(i)% num_T

     if(ierr/=0) stop

     Av   = 0.0
     Rv   = 3.1
     logT = 3.15
     logg = 0.
     call eval_one_bc(t(i), logT, logg, Av, Rv, res, ierr)
     write(*,*) logT, logg, res(1:5)
  enddo

  Av=0.9
  Rv=3.09

  if(.true.)then
     open(1,file='iso.txt')
     open(2,file='iso.out')
     do i=1,2111
        read(1,*) logT, logg, logL
        call eval_one_bc(t(8),logT, logg, Av, Rv, res, ierr)
        if(ierr==0)then
           mags = SolBol - 2.5*logL - res
           write(2,'(99f14.7)') logT, logg, logL, mags(1:8)
        endif
     enddo
     close(1)
     close(2)
  endif

  deallocate(t)

  write(*,*)
  write(*,*) ' ierr = ', ierr
  
end program test_colors
  
