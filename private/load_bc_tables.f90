  module load_bc_tables

    use colors_def

    implicit none

    contains

    subroutine load_one_bc(t,ierr)
      type(bc_table), intent(inout) :: t
      integer, intent(out) :: ierr
      character(len=256) :: binfile
      logical :: binary_exists
      ierr=0
      !if .bin file exists, read it that way
      binfile = trim(t% filename) // '.bin'
      inquire(file=trim(binfile),exist=binary_exists)
      if(binary_exists)then
         call read_one_bin(t,ierr)
      else !read it the old-fashioned way and write a binary 
         call read_one_ascii(t,ierr)
         call write_one_bin(t,ierr)
      endif
    end subroutine load_one_bc

    subroutine read_one_bin(t,ierr)
      type(bc_table), intent(inout) :: t
      integer, intent(out) :: ierr
      character(len=256) :: binfile
      ierr=0
      binfile = trim(t% filename) // '.bin'
      open(1,file=trim(binfile),status='old',form='unformatted',iostat=ierr)
      read(1) t% num_filter, t% num_lines, t% num_Av, t% num_Rv, t% num_T, t% num_g
      allocate(t% grid(2,t% num_lines),t% Av(t% num_Av),t% Rv(t% num_Rv))
      allocate(t% BC(t% num_filter, t% num_lines, t% num_Av, t% num_Rv))
      allocate(t% labels(t% num_filter))
      allocate(t% logg(t% num_g), t% logT(t% num_T))
      allocate(t% index(t% num_T, t% num_g))
      read(1) t% labels
      read(1) t% index
      read(1) t% FeH, t% alphaFe
      read(1) t% logT, t% logg, t% Av, t% Rv
      read(1) t% grid
      read(1) t% BC
      close(1)
      write(*,*) ' read binary file: ', trim(binfile)
    end subroutine read_one_bin

    subroutine write_one_bin(t,ierr)
      type(bc_table), intent(inout) :: t
      integer, intent(out) :: ierr
      character(len=256) :: binfile
      ierr=0
      binfile = trim(t% filename) // '.bin'
      open(1,file=trim(binfile),form='unformatted',iostat=ierr)
      write(1) t% num_filter, t% num_lines, t% num_Av, t% num_Rv, t% num_T, t% num_g
      write(1) t% labels
      write(1) t% index
      write(1) t% FeH, t% alphaFe
      write(1) t% logT, t% logg, t% Av, t% Rv
      write(1) t% grid
      write(1) t% BC
      close(1)
      write(*,*) ' wrote binary file: ', trim(binfile)
    end subroutine write_one_bin

    subroutine read_one_ascii(t,ierr)
      type(bc_table), intent(inout) :: t
      integer, intent(out) :: ierr
      integer :: i, j, k, pass, ng=1, nT=1, r
      real(sp) :: Teff, logg, logT

      open(1,file=trim(t% filename),status='old',iostat=ierr)
      if(ierr/=0) return
      read(1,*) !skip the header
      read(1,'(2x,4i8)') t% num_filter, t% num_lines, t% num_Av, t% num_Rv
      read(1,*)

      allocate(t% grid(2,t% num_lines),t% Av(t% num_Av),t% Rv(t% num_Rv))
      allocate(t% BC(t% num_filter, t% num_lines, t% num_Av, t% num_Rv))
      allocate(t% labels(t% num_filter))

      do r=1,t% num_Rv
         do i=1,t% num_Av
            read(1,*) !skip the column numbers
            read(1,'(31x,199a12)') t% labels(1:t% num_filter)
            do j=1,t% num_lines
               read(1,'(f8.0,f5.1,3f6.2,99f12.6)') Teff, logg, t% FeH, t% Av(i), t% Rv(r), t% BC(:,j,i,r)
               t% grid(1,j) = log10(Teff)
               t% grid(2,j) = logg
            enddo
            !skip the blank lines after each table
            if(i < t% num_Av .and. r < t% num_Rv) then
               read(1,*)
               read(1,*)
            endif
         enddo
      enddo

      close(1)

      !these loops determine a set of unique logT and logg values
      do pass=1,2

         if(pass==2)then
            t% num_T = nT; t% num_g = ng
            allocate(t% logg(t% num_g), t% logT(t% num_T))
         endif
         logT = -99.; nT = 0
         logg = -99.; ng = 0

         do i = 1, t% num_lines

            if( logT < t% grid(1,i) )then
               nT = nT + 1
               logT = t% grid(1,i)
            endif
            
            if( logg < t% grid(2,i) )then
               ng = ng + 1
               logg = t% grid(2,i)
            endif

            if(pass==2)then
               t% logT(nT) = logT
               t% logg(ng) = logg
            endif

         enddo
      enddo


      ! these loops create an index that translates between
      ! unique logT,logg values and the original grid
      allocate(t% index(t% num_T, t% num_g))
      t% index = 0
      do i=1,t% num_lines
         do j=1,t% num_T
            do k=1,t% num_g
               if(t% logT(j)==t% grid(1,i) .and. t% logg(k)==t% grid(2,i))then
                  t% index(j,k) = i
               endif
            enddo
         enddo
      enddo

    end subroutine read_one_ascii

      
  end module load_bc_tables
