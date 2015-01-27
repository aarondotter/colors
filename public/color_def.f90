module color_def

  use const_def, only: sp

  implicit none

  real(sp), parameter :: SolBol = 4.75
  
  type bc_table
     character(len=256) :: filename
     character(len=12), allocatable :: labels(:)
     integer :: num_Av, num_Rv, num_filter, num_lines, num_T, num_g
     integer, allocatable :: index(:,:)
     real(sp), allocatable :: logT(:), logg(:), Av(:), Rv(:), grid(:,:), BC(:,:,:,:)
     real(sp) :: FeH, alphaFe
  end type bc_table
    
end module color_def
