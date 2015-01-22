program test_colors

  use colors_def
  use load_bc_tables
  use eval_bc_tables
  
  implicit none

  integer :: ierr, i, n, m
  type(bc_table), allocatable :: t(:)
  real(sp), allocatable :: res(:), mags(:)
  real(sp) :: logT, logg, logL, Av, Rv
  character(len=256) :: output

  m=8
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

  do i=m,m
     call load_one_bc(t(i),ierr)
  enddo
  allocate(res(t(m)% num_filter),mags(t(m)% num_filter))
  
  do i=m,m
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
     call eval_one_bc(t(m), logT, logg, Av, Rv, res, ierr)
     write(*,*) logT, logg, res(1:5)
  enddo

  Av=0.1
  Rv=3.09

  if(.false.)then
     open(1,file='iso2.txt')
     open(2,file='iso2.out')
     do i=1,580
        read(1,*) logT, logg, logL
        call eval_one_bc(t(m),logT, logg, Av, Rv, res, ierr)
        if(ierr==0)then
           mags = SolBol - 2.5*logL - res
           write(2,'(99f14.7)') logT, logg, logL, &
                mags(35), mags(39), mags(41), mags(45), mags(48)
        endif
     enddo
     close(1)
     close(2)
  endif

  !output = 'junk'
  !call write_one_ascii(t(8),output,ierr)

  deallocate(t)

  write(*,*)
  write(*,*) ' ierr = ', ierr
  
  call interp_koester

contains

  subroutine interp_koester
    type(bc_table) :: k1, k2
    character(len=256) :: output
    integer, parameter :: num_T = 66, num_g=8
    real(sp) :: master_Teff(num_T), master_logg(num_g), logT, logg, Av, Rv
    integer :: i,j,k,r,a
    real(sp), allocatable :: mags(:)

    k1% filename = '../preprocessor/data/Koester.FSPS'
    k2% filename = 'Koester_interp.FSPS'

    print *, 'this  sucks!'

    call load_one_bc(k1,ierr)
    if(ierr/=0) return

    print *, 'still sucks!'

    k2% num_Av = k1% num_Av
    k2% num_Rv = k1% num_Rv
    k2% num_filter = k1% num_filter
    k2% num_lines = num_T*num_g
    k2% num_T = num_T
    k2% num_g = num_g
    k2% FeH = k1% FeH
    k2% alphaFe = k1% alphaFe
    k2% labels = k1% labels

    print *, 'ok, made it here'
    
    allocate(k2% Av(k2% num_Av), k2% Rv(k2% num_Rv))
    allocate(k2% logT(k2% num_T), k2% logg(k2% num_g))
    allocate(k2% grid(2,k2% num_lines))
    allocate(k2% index(k2% num_T, k2% num_g))
    allocate(k2% BC(k2% num_filter, k2% num_lines, k2% num_Av, k2% num_Rv))
    
    k2% Av = k1% Av
    k2% Rv = k1% Rv
    k2% index = 0

    master_logg(1) = 6.0
    do i=2,num_g
       master_logg(i) = master_logg(i-1) + 0.5
    enddo

    master_Teff(1) = 6.0e3
    !  6,000 to    13,000
    do i=2,29
       master_Teff(i) = master_Teff(i-1) + 2.5e2
    enddo

    ! 14,000 to    50,000
    do i=30,66
       master_Teff(i) = master_Teff(i-1) + 1.0e3
    enddo
    
    k2% logT = log10(master_Teff)
    k2% logg = master_logg

    do i=1,num_T
       do j=1,num_g
          k=(i-1)*num_g + j
          k2% grid(1,k) = k2% logT(i)
          k2% grid(2,k) = k2% logg(j)
       enddo
    enddo

    k2% BC = 0.

    allocate(mags(k2% num_filter))
    do r=1,k2% num_Rv
       do a=1,k2% num_Av
          do k=1,k2% num_lines
             logT = k2% grid(1,k)
             logg = k2% grid(2,k)
             Av   = k2% Av(a)
             Rv   = k2% Rv(r)
             call eval_one_bc(k1,logT,logg,Av,Rv,mags,ierr)
             if(ierr==0) k2% BC(:,k,a,r) = mags
          enddo
       enddo
    enddo

    write(*,*) master_logg
    write(*,*) master_Teff
    
    print *, 'made it here!'

    call write_one_ascii(k2,k2% filename,ierr)
    
  end subroutine interp_koester
end program test_colors
