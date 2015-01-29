module eval_bc_tables

  use const_def, only: sp
  use color_def
  
  implicit none
  
  integer, parameter :: interp_linear = 1, interp_cubic = 3
  
contains

  !create a table with one fixed value of Av and Rv
  subroutine create_fixed_Av_Rv(t,u,Av,Rv)
    type(bc_table), intent(in) :: t
    type(bc_table), intent(out) :: u
    real(sp), intent(in) :: Av, Rv
    real(sp), parameter :: eps = 0.001
    real(sp) :: q(4)
    real(sp), allocatable :: res(:)
    integer :: interp, i11, i, iT, ig, iAv=1, iRv=1

    !interp = interp_linear
    interp = interp_cubic

    !set table parameters
    u% filename = 'null and void'
    u% num_Av  = 1
    u% num_Rv  = 1
    u% num_T   = t% num_T
    u% num_g   = t% num_g
    u% num_filter = t% num_filter
    u% num_lines = t% num_lines
    u% FeH     = t% FeH
    u% alphaFe = t% alphaFe

    !allocates
    allocate(u% labels(u% num_filter), u% Av(u% num_Av), u% Rv(u% num_Rv))
    allocate(u% grid(2,u% num_lines))
    allocate(u% BC(u% num_filter, u% num_lines, u% num_Av, u% num_Rv))
    allocate(u% logT(u% num_T), u% logg(u% num_g), u% index(u% num_T, u% num_g))
    allocate(res(u% num_filter))

    !set static arrays
    u% labels  = t% labels
    u% index   = t% index
    u% logT    = t% logT
    u% logg    = t% logg
    u% grid    = t% grid
    u% Av(1)   = Av
    u% Rv(1)   = Rv

    do i=1,t% num_Av - 1
       if(Av>=t% Av(i) .and. Av < t% Av(i+1) )then
          iAv = i
          exit
       endif
    enddo

    do i=1,t% num_Rv - 1
       if(Rv>=t% Rv(i) .and. Rv < t% Rv(i+1) )then
          iRv = i
          exit
       endif
    enddo

    !do the interpolation; even if cubic is chosen, revert 
    !to linear if on the edge of the grid.

    do iT = 1, u% num_T
       do ig = 1, u% num_g

          i11 = t% index(iT,ig)

          !do binlinear interpolation for logT, logg and then for Av, Rv
          do i=1,t% num_filter
             !first steps are bilinear interpolation in logT, logg
             q(1) = t% BC(i,i11,iAv  ,iRv  )
             q(2) = t% BC(i,i11,iAv+1,iRv  )
             q(3) = t% BC(i,i11,iAv  ,iRv+1)
             q(4) = t% BC(i,i11,iAv+1,iRv+1)
             u% BC(i,i11,1,1) = bilinear(t% Av(iAv),t% Av(iAv+1),Av,t% Rv(iRv),t% Rv(iRv+1),Rv,q)
          enddo



       enddo
    enddo
    
    
  end subroutine create_fixed_Av_Rv

  
  !interpolate in logT,logg at fixed Av,Rv 
  subroutine eval_one_bc(t,logT,logg,iAv,iRv,res,ierr)
    type(bc_table), intent(in) :: t
    real(sp), intent(in) :: logT, logg
    integer, intent(in) :: iAv, iRv
    real(sp), intent(out) :: res(:) !res(t% num_filter)
    integer, intent(out) :: ierr
    integer :: i, ig=1, iT=1, interp, j, k
    integer :: i11, i12, i21, i22, i31, i41
    real(sp) :: q(4), r(4)

    if(t% num_filter /= size(res)) stop 'bad news!'
    ierr = 0
    res = 0.

    !interp = interp_linear
    interp = interp_cubic

    !check to make sure logT,logg is within the grid
    if( logT > t% logT(t% num_T) .or. logT < t% logT(1) .or. &
         logg > t% logg(t% num_g) .or. logg < t% logg(1) ) then  
       res=0.
       ierr = -1
       return
    endif

    !now locate the point within the grid
    do i=1, t% num_T - 1
       if(logT >= t% logT(i)  .and. logT < t% logT(i+1))then
          iT = i
          exit
       endif
    enddo

    do i=1, t% num_g - 1
       if(logg >= t% logg(i)  .and. logg < t% logg(i+1))then
          ig = i
          exit
       endif
    enddo

    !do the interpolation; even if cubic is chosen, revert 
    !to linear if on the edge of the grid.
    if(interp==interp_cubic .and. iT>1 .and. ig>1 .and. &
         iT < t% num_T -1 .and. ig < t% num_g -1 )then
       !first option is bicubic, do four interpolation in logT
       !followed by one in logg
       do i = 1, t% num_filter
          do j=1,4
             k=ig+j-2             

             q(1) = t% BC(i, t% index(iT-1,k), iAv, iRv)
             q(2) = t% BC(i, t% index(iT  ,k), iAv, iRv)
             q(3) = t% BC(i, t% index(iT+1,k), iAv, iRv)
             q(4) = t% BC(i, t% index(iT+2,k), iAv, iRv)
        
             r(j) = cubic(t% logT(iT-1:iT+2), logT, q)
          enddo

          res(i) = cubic(t% logg(ig-1:ig+2), logg, r)
       enddo

    else
       !for bilinear interpolation in logT and logg
       i11 = t% index(iT  ,ig  )
       i21 = t% index(iT+1,ig  )
       i12 = t% index(iT  ,ig+1)
       i22 = t% index(iT+1,ig+1)

       !do binlinear interpolation for logT, logg and then for Av, Rv
       do i=1,t% num_filter
          !first steps are bilinear interpolation in logT, logg
          q(1) = t% BC(i,i11,iAv,iRv)
          q(2) = t% BC(i,i21,iAv,iRv)
          q(3) = t% BC(i,i12,iAv,iRv)
          q(4) = t% BC(i,i22,iAv,iRv)
          res(i) = bilinear(t% logT(iT),t% logT(iT+1),logT,t% logg(ig),t% logg(ig+1),logg,q)
       enddo

    endif

  end subroutine eval_one_bc

  function linear(x1,x2,x,y1,y2) result(y)
    real(sp), intent(in) :: x1, x2, x, y1, y2
    real(sp) :: y, slope
    slope = (y2-y1)/(x2-x1)
    y = y1 + slope*(x-x1)
  end function linear

  function bilinear(x1,x2,x,y1,y2,y,q) result(p)
    ! q(1,1) => q(1)  q(2,1) => q(2) 
    ! q(1,2) => q(3)  q(2,2) => q(4)
    real(sp), intent(in) :: x1, x2, y1, y2, x, y, q(4)
    real(sp) :: dx, dy, x1dx, x2dx, y1dy, y2dy, r1, r2, p
    dx = x2 - x1
    dy = y2 - y1
    x2dx = (x2 - x)/dx; x1dx = (x - x1)/dx
    y2dy = (y2 - y)/dy; y1dy = (y - y1)/dy
    r1 = q(1)*x2dx + q(2)*x1dx
    r2 = q(3)*x2dx + q(4)*x1dx
     p =   r1*y2dy +   r2*y1dy
  end function bilinear

  function cubic(x,x0,y) result(y0)
    real(sp), intent(in) :: x(4), y(4), x0
    real(sp) :: a(3), dx, y0
    call interp_4pt_pm_sg(x,y,a)
    dx=x0-x(2)
    y0 = y(2) + dx*( a(1) + dx*(a(2)+dx*a(3)))
  end function cubic
      
  subroutine interp_4pt_pm_sg(x, y, a) !from MESA/interp_1d by Bill Paxton
    ! returns coefficients for monotonic cubic interpolation from x(2) to x(3)
    real(sp), intent(in)    :: x(4)    ! junction points, strictly monotonic
    real(sp), intent(in)    :: y(4)    ! data values at x's
    real(sp), intent(out)   :: a(3)    ! coefficients
    real(sp) :: h1, h2, h3, s1, s2, s3, p2, p3, as2, ss2, yp2, yp3
    ! for x(2) <= x <= x(3) and dx = x-x(2), 
    ! y(x) = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
    h1 = x(2)-x(1)
    h2 = x(3)-x(2)
    h3 = x(4)-x(3)
    s1 = (y(2)-y(1))/h1
    s2 = (y(3)-y(2))/h2
    s3 = (y(4)-y(3))/h3
    p2 = (s1*h2+s2*h1)/(h1+h2)
    p3 = (s2*h3+s3*h2)/(h2+h3)
    as2 = abs(s2)
    ss2 = sign(1.0, s2)
    yp2 = (sign(1.0, s1)+ss2)*min(abs(s1), as2, 0.5*abs(p2))
    yp3 = (ss2+sign(1.0, s3))*min(as2, abs(s3), 0.5*abs(p3))
    a(1) = yp2
    a(2) = (3*s2-2*yp2-yp3)/h2
    a(3) = (yp2+yp3-2*s2)/(h2*h2)
  end subroutine interp_4pt_pm_sg      

end module eval_bc_tables
