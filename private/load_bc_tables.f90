  module load_bc_tables

    use colors_def

    implicit none

    contains

    subroutine load_one_bc(t,ierr)
      type(bc_table), intent(inout) :: t
      integer, intent(out) :: ierr
      integer :: i, j, k, pass, ng=1, nT=1
      real(sp) :: Teff, logg, logT

      open(1,file=trim(t% filename),status='old',iostat=ierr)
      read(1,'(2x,3i4)') t% num_lines, t% num_Av, t% num_filter

      allocate(t% grid(2,t% num_lines),t% Av(t% num_Av))
      allocate(t% BC(t% num_lines,t% num_filter,t% num_Av))
      allocate(t% labels(t% num_filter+5))
      
      do i=1,t% num_Av
         read(1,*) !skip the column numbers
         read(1,'(1x, a7, a5, 3a6 ,99a12)') t% labels(1:t% num_filter+5)
         do j=1,t% num_lines
            read(1,'(f8.0,f5.1,3f6.2,99f12.6)') Teff, logg, t% FeH, t% Av(i), t% Rv, t% BC(j,:,i)
            t% grid(1,j) = log10(Teff)
            t% grid(2,j) = logg
         enddo
         !skip the blank lines after each table
         if(i<t% num_Av) then
            read(1,*)
            read(1,*)
         endif
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

    end subroutine load_one_bc

      
  end module load_bc_tables
