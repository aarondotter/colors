module load_new_tables

  use const_def, only: sp
  use interp_1d_def
  use interp_1d_lib_sg
  use interp_2d_lib_sg!, only: interp_mkbicub_sg, interp_evbicub_sg

  implicit none

  type BC
     real(sp), pointer :: data(:)
  end type BC
  
  type table
     character(len=256) :: filename
     character(len=12), allocatable :: labels(:)
     integer :: num_Av, num_Rv, num_filter, num_T, num_g, iling, ilinT
     real(sp) :: Rv, FeH, alphaFe
     real(sp), allocatable :: logT(:), logg(:), Av(:)
     type(BC), allocatable :: bcs(:,:) !Av, filters
  end type table
  
  contains

  subroutine load_one_bc(t,ierr)
      type(table), intent(inout) :: t
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
         if(ierr==0) call write_one_bin(t,ierr)
      endif
    end subroutine load_one_bc

    subroutine read_one_bin(t,ierr)
      type(table), intent(inout) :: t
      integer, intent(out) :: ierr
      character(len=256) :: binfile
      integer :: i, j, num_lines
      ierr=0
      binfile = trim(t% filename) // '.bin'
      open(1,file=trim(binfile),status='old',form='unformatted',iostat=ierr)
      read(1)  t% num_g, t% num_T, t% num_Av, t% num_filter
      allocate(t% Av(t% num_Av))
      allocate(t% bcs(t% num_filter,t% num_Av))
      allocate(t% labels(t% num_filter))
      allocate(t% logg(t% num_g), t% logT(t% num_T))
      read(1) t% labels
      read(1) t% FeH, t% alphaFe, t% Rv
      read(1) t% Av
      read(1) t% logg
      read(1) t% logT
      num_lines = t% num_g * t% num_T
      do i=1,t% num_filter
         do j=1,t% num_Av
            allocate(t% bcs(i,j)% data(4*num_lines))
            read(1) t% bcs(i,j)% data
         enddo
      enddo
      close(1)
    end subroutine read_one_bin

    subroutine write_one_bin(t,ierr)
      type(table), intent(in) :: t
      integer, intent(out) :: ierr
      character(len=256) :: binfile
      integer :: i, j
      ierr=0
      binfile = trim(t% filename) // '.bin'
      open(1,file=trim(binfile),action='write',form='unformatted',iostat=ierr)
      write(1)  t% num_g, t% num_T, t% num_Av, t% num_filter
      write(1) t% labels
      write(1) t% FeH, t% alphaFe, t% Rv
      write(1) t% Av
      write(1) t% logg
      write(1) t% logT
      do i=1,t% num_Filter
         do j=1,t% num_Av
            write(1) t% bcs(i,j)% data
         enddo
      enddo
      close(1)
    end subroutine write_one_bin

    subroutine read_one_ascii(t,ierr)
      type(table), intent(inout) :: t
      integer, intent(out) :: ierr
      integer :: i, j, k, pass, ng, nT, num_lines, ibcTmin, ibcTmax, ibcgmin, ibcgmax
      real(sp) :: Teff, logg, logT
      real(sp), allocatable :: grid(:,:), bcTmin(:), bcTmax(:), bcgmin(:), bcgmax(:), BCdata(:)

      write(*,*) ' opening ', trim(t% filename)
      open(1,file=trim(t% filename),status='old',action='read',iostat=ierr)
      if(ierr/=0) then
         write(*,*) ' problem opening ascii file, ierr = ', ierr
         return
      endif
      read(1,*) !skip the header
      read(1,'(2x,3i8)') t% num_filter, num_lines, t% num_Av
      read(1,*)
      
      allocate(t% Av(t% num_Av), grid(2+t% num_filter,num_lines))
      allocate(t% labels(t% num_filter), t% bcs(t% num_filter, t% num_Av))
      allocate(BCdata(t% num_filter))

      do i=1,t% num_Av
         read(1,*) !skip the column numbers
         read(1,'(31x,199a12)') t% labels(1:t% num_Filter)

         do k=1,t% num_filter
            allocate(t% bcs(k,i)% data(4*num_lines))
         enddo

         do j=1, num_lines
            read(1,'(f8.0,f5.1,3f6.2,199f12.6)') Teff, logg, t% FeH, t% Av(i), t% Rv, BCdata
            grid(1,j) = log10(Teff)
            grid(2,j) = logg
            do k=1,t% num_filter
               t% bcs(k,i)% data(4*(j-1)+1) = BCdata(k)
               !t% bcs(k,i)% data(j) = BCdata(k)
            enddo
         enddo
      enddo
      close(1)

      write(*,*) ' finished reading . . .'

      !these loops determine a set of unique logT and logg values
      nT=0; ng=0
      do pass=1,2
         
         if(pass==2)then
            if( ng*nT == num_lines) then !pass
               t% num_T = nT; t% num_g = ng
               allocate(t% logg(t% num_g), t% logT(t% num_T))
            else
               stop 'fail!'
            endif
         endif
         logT = -99.; nT = 0
         logg = -99.; ng = 0

         do i = 1, num_lines

            if( logT < grid(1,i) )then
               nT = nT + 1
               logT = grid(1,i)
            endif
            
            if( logg < grid(2,i) )then
               ng = ng + 1
               logg = grid(2,i)
            endif

            if(pass==2)then
               t% logT(nT) = logT
               t% logg(ng) = logg
            endif

         enddo
      enddo

      !process 2-d interpolation data once
      ! use "not a knot" bc's
      allocate(bcgmin(ng),bcgmax(ng),bcTmin(nT),bcTmax(nT))
      ibcTmin = 0; bcTmin(:) = 0d0
      ibcTmax = 0; bcTmax(:) = 0d0
      ibcgmin = 0; bcgmin(:) = 0d0
      ibcgmax = 0; bcgmax(:) = 0d0

      do i=1,t% num_filter       
         do j=1,t% num_Av
            call interp_mkbicub_sg(t% logg, t% num_g, t% logT, t% num_T, &
                 t% bcs(i,j)% data, t% num_g, ibcgmin, bcgmin, ibcgmax, bcgmax, &
                 ibcTmin, bcTmin, ibcTmax, bcTmax, t% iling, t% ilinT, ierr)
         enddo
      enddo

    end subroutine read_one_ascii

    function table_interp_filter_fixed_Av(t,logg,logT,iflt,iAv,ierr) result(res)
      type(table), intent(inout) :: t
      real(sp), intent(in) :: logT, logg
      integer, intent(in) :: iflt, iAv
      real(sp) :: res, result(6)
      integer, intent(out) :: ierr
      integer :: ict(6)=[1,0,0,0,0,0]     
      call interp_evbicub_sg(logg, logT, t% logg, t% num_g, &
           t% logT, t% num_T, t% iling, t% ilinT, t% bcs(iflt,iAv)% data, &
           t% num_g, ict, result, ierr)
      res=result(1)
    end function table_interp_filter_fixed_Av

    function table_interp_filter(t,logg,logT,iflt,Av,ierr) result(res)
      type(table), intent(inout) :: t
      real(sp), intent(in) :: logT, logg, Av
      integer, intent(in) :: iflt
      integer, intent(out) :: ierr
      real(sp) :: res, xnew(1), ynew(1)
      real(sp), allocatable :: yold(:)
      integer :: i
      integer, parameter :: nwork = max(mp_work_size,pm_work_size)
      real(sp), pointer :: work(:)=>NULL()
      character(len=256) :: str
      ierr = 0
      res = 0.
      if(t% num_Av==1)then
         res=table_interp_filter_fixed_Av(t,logg,logT,iflt,1,ierr)
         return
      endif
      allocate(yold(t% num_Av),work(nwork*t%num_Av))
      do i=1, t% num_Av
         yold(i) = table_interp_filter_fixed_Av(t,logg,logT,iflt,i,ierr)
      enddo
      xnew(1)=Av
      call interpolate_vector_sg( t% num_Av, t% Av, 1, xnew, yold, ynew, &
           interp_m3a_sg, nwork, work, str, ierr)  
      if(ierr/=0) then
         write(*,*) trim(str)
      else
         res=ynew(1)
      endif
      deallocate(work)
    end function table_interp_filter

    subroutine table_interp_filters(t,logg,logT,Av,res,ierr)
      type(table), intent(inout) :: t
      real(sp), intent(in) :: logT, logg, Av
      real(sp), intent(out) :: res(:)
      integer, intent(out) :: ierr
      integer :: i
      if(t% num_filter /= size(res)) then
         ierr=-1
         return
      endif
      do i=1,t% num_filter
         res(i)=table_interp_filter(t,logg,logT,i,Av,ierr)
      enddo
    end subroutine table_interp_filters

    subroutine table_interp_filters_fixed_Av(t,logg,logT,iAv,res,ierr)
      type(table), intent(inout) :: t
      real(sp), intent(in) :: logT, logg
      integer, intent(in) :: iAv
      real(sp) :: res(:)
      integer, intent(out) :: ierr
      integer :: i
      if(t% num_filter/=size(res)) then
         ierr=-1
         return
      endif
      do i=1,t% num_filter
         res(i)=table_interp_filter_fixed_Av(t,logg,logT,i,iAv,ierr)
      enddo
    end subroutine table_interp_filters_fixed_Av
    
end module load_new_tables


program test_new_tables
  use load_new_tables
  real(sp) :: logT, logg
  real(sp), allocatable :: res1(:), res2(:)
  type(table) :: x
  integer :: ierr

  x% filename = '/home/dotter/science/colors/data/fehp000_afep00.FSPS'

  call load_one_bc(x,ierr)

  allocate(res1(x% num_filter),res2(x% num_filter))

  logg=4.0
  logT=3.5
  do while(logT<5)
     call table_interp_filters_fixed_Av(x,logg,logT,1,res1,ierr)
     call table_interp_filters(x,logg,logT,1e-5,res2,ierr)
     res1=res1-res2
     write(*,*) logT, res1(1:5)
     logT = logT+0.1
  enddo

end program test_new_tables
