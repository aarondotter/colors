module eval_bc_tables

  use colors_def
    
  implicit none
    
contains
      
  subroutine eval_one_bc(t,logT,logg,Av,Rv,res,ierr)
    type(bc_table), intent(in) :: t
    real(sp), intent(in) :: logT, logg, Av, Rv
    real(sp), intent(out) :: res(:) !t% num_filter
    integer, intent(out) :: ierr
    integer :: i, i11, i12, i21, i22, ig=1, iT=1, iAv=1, iRv=1
    real(sp) :: q(4), r(4), mag1, mag2, a, b
   
    if(t% num_filter /= size(res)) stop 'bad news!'
    ierr = 0
    res = 0.

    if( logT > t% logT(t% num_T) .or. logT < t% logT(1) .or. &
         logg > t% logg(t% num_g) .or. logg < t% logg(1) .or. &
         Av > t% Av( t% num_Av) .or. Av < t% Av(1) ) then  
       ierr = -1
       return
    endif

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

    !for bilinear interpolation in logT and logg
    i11 = t% index(iT  ,ig  )
    i21 = t% index(iT+1,ig  )
    i12 = t% index(iT  ,ig+1)
    i22 = t% index(iT+1,ig+1)

    do i=1, t% num_Av - 1
       if(Av >= t% Av(i) .and. Av < t% Av(i+1))then
          iAv=i
          exit
       endif
    enddo

    do i=1, t% num_Rv - 1
       if(Rv >= t% Rv(i) .and. Rv < t% Av(i+1))then
          iRv=i
          exit
       endif
    enddo

    !do binlinear interpolation for logT, logg and then for Av, Rv
    do i=1,t% num_filter
       q=0.
       r=0.

       !first steps are bilinear interpolation in logT, logg
       q(1) = t% BC(i,i11,iAv  ,iRv  ); q(2) = t% BC(i,i21,iAv  ,iRv  )
       q(3) = t% BC(i,i12,iAv  ,iRv  ); q(4) = t% BC(i,i22,iAv  ,iRv  )
       r(1) = bilinear(t% logT(iT), t% logT(iT+1), logT, t% logg(ig), t% logg(ig+1), logg, q)

       q(1) = t% BC(i,i11,iAv+1,iRv  ); q(2) = t% BC(i,i21,iAv+1,iRv  )
       q(3) = t% BC(i,i12,iAv+1,iRv  ); q(4) = t% BC(i,i22,iAv+1,iRv  )
       r(2) = bilinear(t% logT(iT), t% logT(iT+1), logT, t% logg(ig), t% logg(ig+1), logg, q)

       q(1) = t% BC(i,i11,iAv  ,iRv+1); q(2) = t% BC(i,i21,iAv  ,iRv+1)
       q(3) = t% BC(i,i12,iAv  ,iRv+1); q(4) = t% BC(i,i22,iAv  ,iRv+1)
       r(3) = bilinear(t% logT(iT), t% logT(iT+1), logT, t% logg(ig), t% logg(ig+1), logg, q)

       q(1) = t% BC(i,i11,iAv+1,iRv+1); q(2) = t% BC(i,i21,iAv+1,iRv+1)
       q(3) = t% BC(i,i12,iAv+1,iRv+1); q(4) = t% BC(i,i22,iAv+1,iRv+1)
       r(4) = bilinear(t% logT(iT), t% logT(iT+1), logT, t% logg(ig), t% logg(ig+1), logg, q)

       !final step is bilinear in Av, Rv
       res(i) = bilinear(t% Av(iAv), t% Av(iAv+1), Av, t% Rv(iRv), t% Rv(iRv+1), Rv, r)
    enddo

  end subroutine eval_one_bc

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
      
end module eval_bc_tables
