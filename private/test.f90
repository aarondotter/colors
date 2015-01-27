program test

  implicit none

  real :: x(4), y(4), x0

  x = [1.,2.,3.,4.]
  y = [-3.,11.,67.,189.]

  x0=2.1
  do while(x0<3.)
     write(*,*) x0, linear(x(2), x(3), x0, y(2), y(3)), cubic(x,y,x0)
     x0=x0+0.1
  enddo

contains

  function cubic(x,y,x0) result(y0)
    real, intent(in) :: x(4), y(4), x0
    real :: a(3), dx, y0
    call interp_4pt_pm_sg(x,y,a)
    dx=x0-x(2)
    y0 = y(2) + dx*( a(1) + dx*(a(2)+dx*a(3)))
  end function cubic

  function linear(x1,x2,x,y1,y2) result(y)
    real, intent(in) :: x1, x2, x, y1, y2
    real :: y, slope
    slope = (y2-y1)/(x2-x1)
    y = y1 + slope*(x-x1)
  end function linear
    
  subroutine interp_4pt_pm_sg(x, y, a) !from MESA/interp_1d by Bill Paxton
    ! returns coefficients for monotonic cubic interpolation from x(2) to x(3)
    real, intent(in)    :: x(4)    ! junction points, strictly monotonic
    real, intent(in)    :: y(4)    ! data values at x's
    real, intent(out)   :: a(3)    ! coefficients
    real :: h1, h2, h3, s1, s2, s3, p2, p3, as2, ss2, yp2, yp3
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

end program test
