program table_merge

  implicit none

  integer, parameter :: sp = selected_real_kind(p=5)
  integer, parameter :: dp = selected_real_kind(p=15)

  character(len=256) :: outfile, infile, wd_file, bb_file, metallicity
  integer :: ierr, i, j, k
  
  type bc_table
     real(sp), allocatable :: mags(:,:,:,:), Teff(:), logg(:), Av(:), Rv(:)
     real(sp) :: FeH
     character(len=12), allocatable :: header(:)
     integer :: num_Rv, num_Av, num_spectra, num_filters
  end type bc_table
  type(bc_table) :: blackbody, atlas, rauch, final, koester

  integer, parameter :: nT=106, ng=20
  real(sp) :: master_Teff(nT), master_logg(ng)

  if(command_argument_count()<2)then
     stop 'not enough command arguments. specify input ATLAS file and output file'
  endif

  call get_command_argument(1, infile)
  call get_command_argument(2,outfile)

  !base layer is blackbody spectra
  bb_file='/home/dotter/science/colors/preprocessor/data/blackbody.FSPS'
  call readBC(blackbody,bb_file)

  !add Rauch and Koester post-AGB/WD models
  wd_file='/home/dotter/science/colors/preprocessor/data/rauch_solar.FSPS'
  call readBC(rauch,wd_file)
  wd_file='/home/dotter/science/colors/preprocessor/data/Koester_interp.FSPS'
  call readBc(koester,wd_file)

  !final layer from ATLAS
  call readBC(atlas,infile)

  !make sure all the tables are compatible
  if(incompatible(blackbody,atlas).or.incompatible(blackbody,rauch) &
       .or.incompatible(blackbody,koester)) stop 'tables incompatible'

  !set up the final results
  call master_bc_init(final,blackbody)

  !merge the tables based on common Teff,logg points
  final% mags = -99.
  call merge_one(blackbody,final,.false.)
  call merge_one(rauch,final,.true.)
  call merge_one(koester,final,.true.)
  call merge_one(atlas,final,.true.)

  if(command_argument_count()==3)then
     call get_command_argument(3,metallicity)
     read(metallicity,*) final% FeH
  endif

  !write the final result
  call writeBC(final,outfile)

contains

  subroutine merge_one(b,m,check_logg)
    type(bc_table), intent(in) :: b
    type(bc_table), intent(inout) :: m
    logical, intent(in) :: check_logg
    integer :: i,j,k,l,lo,hi

    m% FeH = b% FeH

    do j=1,nT
       if(check_logg)then
          do i=1,ng
             !for ATLAS and Rauch, check logg
             do k=1,b% num_spectra
                if( abs(master_Teff(j) - b% Teff(k)) < 1e-1 .and. &
                     abs(master_logg(i) - b% logg(k)) < 1e-3 ) then 
                   m% mags(:,:,:,(j-1)*ng+i) = b% mags(:,:,:,k)
                endif
             enddo
          enddo
       else
          lo=(j-1)*ng+1
          hi=j*ng
          !for blackbody only, we ignore logg
          do k=1,b% num_spectra
             if( abs(master_Teff(j) - b% Teff(k)) < 1e-1 ) then
                do l=lo,hi
                   m% mags(:,:,:,l) = b% mags(:,:,:,k)
                enddo
             endif
          enddo
       endif
    enddo

  end subroutine merge_one

  subroutine master_bc_init(m,b)
    type(bc_table), intent(out) :: m
    type(bc_table), intent(in)  :: b
    integer :: i,j,k

    m% num_filters = b% num_filters
    m% num_Av = b% num_Av
    m% num_Rv = b% num_Rv
    m% FeH = b% FeH
    m% num_spectra = nT * ng

    allocate(m% header(m% num_filters+5))
    allocate(m% Teff(m% num_spectra),m% logg(m% num_spectra),m% Rv(m% num_Rv),m% Av(m% num_Av))
    allocate(m% mags(m% num_filters,m% num_Av,m% num_Rv,m% num_spectra))

    m% header = b% header
    m% Av = b% Av
    m% Rv = b% Rv

    call set_master_Teff_logg

    k=1
    do i=1,nT
       do j=1,ng
          m% Teff(k) = master_Teff(i)
          m% logg(k) = master_logg(j)
          k = k+1
       enddo
    enddo


  end subroutine master_bc_init

  logical function incompatible(a,b)
    type(bc_table), intent(in) :: a,b
    incompatible = (a% num_Av /= b% num_Av) .or. (a% num_Rv /= b% num_Rv) .or. (a% num_filters /= b% num_filters)
  end function incompatible

  subroutine set_master_Teff_logg
    !set the temperature and gravity scales for MIST bolometric correction tables
    master_Teff(1) = 1.0e3; master_Teff(2) = 1.5e3; master_Teff(3) = 2.0e3; master_Teff(4)=2.5e3
    master_Teff(5) = 2.8e3; master_Teff(6) = 3.0e3; master_Teff(7) = 3.2e3; master_Teff(8)=3.5e3
       
    !  3,750 to    13,000
    do i=9,46
       master_Teff(i) = master_Teff(i-1) + 2.5e2
    enddo

    ! 14,000 to    50,000
    do i=47,83
       master_Teff(i) = master_Teff(i-1) + 1.0e3
    enddo
  
    ! 60,000 to   200,000
    do i=84,98
       master_Teff(i) = master_Teff(i-1) + 1.0e4
    enddo

    !200,000 to 1,000,000
    do i=99,nT
       master_Teff(i) = master_Teff(i-1) + 1.0e5
    enddo

    !logg from -2 to +10 in 0.5 dex
    master_logg(1) = -2.0
    do i=2,15
       master_logg(i) = master_logg(i-1) + 0.5
    enddo
    do i=16,ng
       master_logg(i) = master_logg(i-1) + 1.0
    enddo

  end subroutine set_master_Teff_logg

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
