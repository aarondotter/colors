module color_lib

  use const_def, only: sp
  use color_def
  use load_bc_tables
  use eval_bc_tables

  implicit none

contains

  subroutine color_init(table_list,t,ierr)
    character(len=256), intent(in) :: table_list
    type(bc_table), allocatable, intent(inout) :: t(:)
    integer, intent(out) :: ierr
    integer :: i,n
    character(len=256) :: line
    open(99,file=trim(table_list),iostat=ierr)
    if(ierr/=0) return
    read(99,*,iostat=ierr) n
    allocate(t(n))
    if(ierr/=0) return
    i=1
    do while(i<=n)
       read(99,'(a)',iostat=ierr) line
       if(ierr/=0) exit
       if(line=='' .or. line(1:1)=='#' .or. line(1:1)=='!') cycle
       read(line,'(a)') t(i)% filename
       call load_one_bc(t(i),ierr)
       if(ierr/=0) exit
       i=i+1
    enddo
    close(99)
  end subroutine color_init

  subroutine color_get(t,logT,logg,iAv,iRv,res,ierr)
    type(BC_table), intent(in) :: t
    real(sp) :: logT, logg
    integer, intent(in) :: iAv, iRv
    real(sp), intent(out) :: res(:) !res(t% num_filter)
    integer, intent(out) :: ierr
    call eval_one_BC(t,logT,logg,iAv,iRv,res,ierr)
  end subroutine color_get

  subroutine color_create_fixed_Av_Rv(t,u,Av,Rv)
    type(bc_table), intent(in) :: t
    type(bc_table), intent(out) :: u
    real(sp), intent(in) :: Av, Rv
    call create_fixed_Av_Rv(t,u,Av,Rv)
  end subroutine color_create_fixed_Av_Rv

  subroutine color_write_ascii(t,output,ierr)
    type(bc_table), intent(in) :: t
    character(len=256) :: output
    integer, intent(out) :: ierr
    call write_one_ascii(t,output,ierr)
  end subroutine color_write_ascii

end module color_lib
