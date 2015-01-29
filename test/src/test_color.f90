program test_color

  use color_def
  use color_lib
  
  implicit none

  integer :: ierr, i, j, n, m, iRv, iAv
  type(bc_table), allocatable :: t(:)
  type(bc_table) :: u
  real(sp), allocatable :: res(:), mags(:)
  real(sp) :: logT, logg, logL, Av, Rv
  character(len=256) :: output, table_list

  table_list = 'bc_table.list'
  call color_init(table_list,t,ierr)
  if(ierr/=0) then
     write(*,*) 'color_init failed with ierr = ', ierr
  endif

  n=size(t)
  m=6

  allocate(res(t(m)% num_filter),mags(t(m)% num_filter))
  
  do i=1,n
     write(*,*) ' i    = ', i
     write(*,*) '[Fe/H]= ', t(i)% FeH
     write(*,*) ' ierr = ', ierr
     write(*,*) 'num_Av = ', t(i)% num_Av
     write(*,*) 'num_filter = ', t(i)% num_filter
     write(*,*) 'num_lines = ', t(i)% num_lines
     write(*,*) 'num_g     = ', t(i)% num_g
     write(*,*) 'num_T     = ', t(i)% num_T

     iAv   = 1
     iRv   = 2
     logT = 3.9
     logg = 8.7
     call color_get(t(m), logT, logg, iAv, iRv, res, ierr)
     write(*,*) logT, logg, res(1:5)
  enddo


  Av=0.10
  Rv=3.1
  call color_create_fixed_Av_Rv(t(8),u,Av,Rv)
  u% filename = 'test.FSPS'
  call color_write_ascii(u,u% filename,ierr)
  

  if(.false.)then
     open(1,file='iso.txt')
     open(2,file='iso_cubic.out')

     do i=1,982
        read(1,*) logT, logg, logL
        call color_get(t(m),logT, logg, iAv, iRv, res, ierr)
        if(ierr==0)then
           mags = SolBol - 2.5*logL - res
           write(2,'(99f14.7)') logT, logg, logL, &
                mags(35), mags(39), mags(41), mags(45), mags(48)
        endif
     enddo
     close(1)
     close(2)
  endif

  if(.true.)then
     m=8
     open(1,file='iso8.txt')
     open(2,file='iso8.out')
     write(2,'(99a14)') 'logT', 'logg', 'logL', t(8)% labels
     do i=1,581
        read(1,*) logT, logg, logL
        call color_get(t(m), logT, logg, iAv, iRv, res, ierr)
        !res holds bolometric corrections
        mags = SolBol - 2.5*logL - res
        write(2,'(99f14.7)') logT, logg, logL, mags
     enddo
     close(1)
     close(2)

     open(1,file='iso85.txt')
     open(2,file='iso85.out')
     write(2,'(99a14)') 'logT', 'logg', 'logL', t(8)% labels
     do i=1,570
        read(1,*) logT, logg, logL
        call color_get(t(m), logT, logg, iAv, iRv, res, ierr)
        mags = SolBol - 2.5*logL - res
        write(2,'(99f14.7)') logT, logg, logL, mags
     enddo
     close(1)
     close(2)

     open(1,file='iso9.txt')
     open(2,file='iso9.out')
     write(2,'(99a14)') 'logT', 'logg', 'logL', t(8)% labels
     do i=1,559
        read(1,*) logT, logg, logL
        call color_get(t(m), logT, logg, iAv, iRv, res, ierr)
        mags = SolBol - 2.5*logL - res
        write(2,'(99f14.7)') logT, logg, logL, mags
     enddo
     close(1)
     close(2)
  endif


  deallocate(t)

  write(*,*)
  write(*,*) ' ierr = ', ierr
  
  if(.false.) call interp_koester
  
contains

  subroutine interp_koester
    type(bc_table) :: k1, k2
    character(len=256) :: output
    integer, parameter :: num_T = 66, num_g=4
    real(sp) :: master_Teff(num_T), master_logg(num_g), logT, logg
    integer :: i,j,k,r,a
    real(sp), allocatable :: mags(:)

    k1% filename = '../preprocessor/data/Koester.FSPS'
    k2% filename = 'Koester_interp.FSPS'

    call load_one_bc(k1,ierr)
    if(ierr/=0) return

    k2% num_Av = k1% num_Av
    k2% num_Rv = k1% num_Rv
    k2% num_filter = k1% num_filter
    k2% num_lines = num_T*num_g
    k2% num_T = num_T
    k2% num_g = num_g
    k2% FeH = k1% FeH
    k2% alphaFe = k1% alphaFe
    k2% labels = k1% labels
    
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
       master_logg(i) = master_logg(i-1) + 1.0
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
             call eval_one_bc(k1,k2% grid(1,k),k2% grid(2,k),a,r,mags,ierr)
             if(ierr==0) k2% BC(:,k,a,r)=mags
          enddo
       enddo
    enddo

    call write_one_ascii(k2,k2% filename,ierr)
    
  end subroutine interp_koester
end program test_color
