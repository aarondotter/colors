program table_merge

  implicit none

  integer, parameter :: sp = selected_real_kind(p=5)
  integer, parameter :: dp = selected_real_kind(p=15)

  character(len=256) :: outfile, infile
  integer :: ierr, i, j, k
  
  type bc_table
     real(dp), allocatable :: mags(:,:,:,:), Teff(:), logg(:), Av(:), Rv(:)
     real(dp) :: FeH
     character(len=12), allocatable :: header(:)
     integer :: num_Rv, num_Av, num_spectra, num_filters
  end type bc_table
  type(bc_table) :: bb, kur

  integer, parameter :: nT=106, ng=23
  real(dp) :: master_Teff(nT), master_logg(ng), master_mags(nT,ng)


!set the temperature and gravity scales for MIST bolometric correction tables
  master_Teff(1) = 1.0d3; master_Teff(2) = 1.5d3; master_Teff(3) = 2.0d3; master_Teff(4)=2.5d3
  master_Teff(5) = 2.8d3; master_Teff(6) = 3.0d3; master_Teff(7) = 3.2d3; master_Teff(8)=3.5d3
       
  !  3,750 to    13,000
  do i=9,46
     master_Teff(i) = master_Teff(i-1) + 2.5d2
  enddo

  ! 14,000 to    50,000
  do i=47,83
     master_Teff(i) = master_Teff(i-1) + 1.0d3
  enddo
  
  ! 60,000 to   200,000
  do i=84,98
     master_Teff(i) = master_Teff(i-1) + 1.0d4
  enddo

  !200,000 to 1,000,000
  do i=99,nT
     master_Teff(i) = master_Teff(i-1) + 1.0d5
  enddo
  
  do i=1,ng
     master_logg(i) = -1.0d0 + 0.5d0*real(i-1,kind=dp)
  enddo

  write(*,*) master_logg
  write(*,*) master_Teff


  infile='/home/dotter/science/colors/preprocessor/data/blackbody.UBVRIJHKsKp'
  call readBC(bb,infile)

  infile='/home/dotter/science/colors/preprocessor/data/Zp0d0ap0.UBVRIJHKsKp'
  call readBC(kur,infile)

  outfile='junk.bb'
  call writeBC(bb,outfile)

  outfile='junk.kur'
  call writeBC(kur,outfile)


  master_mags = -99.
  do j=1,nT
     !for blackbody only, we ignore logg
     do k=1,bb% num_spectra
        if( abs(master_Teff(j) - bb% Teff(k)) < 1d-1 ) master_mags(j,:) = bb% mags(1,1,1,k)
     enddo

     do i=1,ng
        !forr ATLAS, check logg
        do k=1,kur% num_spectra
           if( abs(master_Teff(j) - kur% Teff(k)) < 1d-1 .and. &
                abs(master_logg(i) - kur% logg(k)) < 1d-3 ) then 
              master_mags(j,i) = kur% mags(1,1,1,k)
              !write(*,*) kur% logg(k), kur% Teff(k)
           endif
        enddo
     enddo
  enddo

  open(1,file='master')
  do j=1,nT
     write(1,'(99f18.8)') master_Teff(j), master_mags(j,:)
  enddo
  close(1)

contains

  subroutine readBC(bc,filename)
    type(bc_table), intent(out) :: bc
    character(len=256) :: filename
    open(1,file=trim(filename),iostat=ierr)
    if(ierr/=0) then
       write(*,*) 'failed to open input file'
       return
    endif
    read(1,*)  !skip header
    read(1,'(2x,4i8)') bc% num_filters, bc% num_spectra, bc% num_Av, bc% num_Rv
    read(1,*)  !skip header
    allocate(bc% header(bc% num_filters+5))
    allocate(bc% Teff(bc% num_spectra),bc% logg(bc% num_spectra),bc% Rv(bc% num_Rv),bc% Av(bc% num_Av))
    allocate(bc% mags(bc% num_filters,bc% num_Av,bc% num_Rv,bc% num_spectra))
    do j=1,bc% num_Rv
       do k=1,bc% num_Av
          read(1,*) !skip header
          read(1,'(1x,a7,a5,3a6,99a12)') bc% header
          do i=1,bc% num_spectra
             read(1,'(f8.0,f5.1,3f6.2,99f12.6)') bc% Teff(i), bc% logg(i), bc% FeH, bc% Av(k), &
                  bc% Rv(j), bc% mags(:,k,j,i)
          enddo
          if(k<bc% num_Av.and.j<bc% num_Rv)then
             read(1,*) !2 blank lines
             read(1,*)
          endif
       enddo
    enddo
    close(1)
  end subroutine readBC

  subroutine writeBC(bc,filename)
    type(bc_table), intent(in) :: bc
    character(len=256), intent(in) :: filename
    open(1,file=trim(filename),iostat=ierr)
    if(ierr/=0) then
       write(*,*) 'failed to write output file'
       write(*,*) 'ierr = ', ierr
       return
    endif
    write(1,'(a1,1x,4a8)') '#', 'filters', 'spectra', 'num Av', 'num Rv'
    write(1,'(a1,1x,4i8)') '#', bc% num_filters, bc% num_spectra, bc% num_Av, bc% num_Rv
    write(1,'(a1)') '#'
    do j=1,bc% num_Rv
       do k=1,bc% num_Av
          write(1,'(a1,i7, i5, 3i6, 99(5x,i2,5x))') '#', (i,i=1,bc% num_filters+5)
          write(1,'(a1,a7, a5, 3a6 ,99a12)') '#', bc% header
          do i=1,bc% num_spectra
             write(1,'(f8.0,f5.1,3f6.2,99f12.6)') bc% Teff(i), bc% logg(i), bc% FeH, bc% Av(k), &
                  bc% Rv(j), bc% mags(:,k,j,i)
          enddo
          if(k<bc% num_Av.and.j<bc% num_Rv)then
             write(1,*)
             write(1,*)
          endif
       enddo
    enddo
    close(1)
    write(0,*) '   output to ', trim(outfile)
  end subroutine writeBC

end program table_merge
