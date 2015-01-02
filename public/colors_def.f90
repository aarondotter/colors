module colors_def

  implicit none

  integer, parameter :: sp = selected_real_kind(p=5)
  integer, parameter :: dp = selected_real_kind(p=15)
  
  type bc_table
     character(len=256) :: filename
     character(len=12), allocatable :: labels(:)
     integer :: num_Av, num_filter, num_lines, num_T, num_g
     integer, allocatable :: index(:,:)
     real(sp), allocatable :: logT(:), logg(:), Av(:), grid(:,:), BC(:,:,:)
     real(sp) :: Rv, FeH, alphaFe
  end type bc_table
    
end module colors_def
