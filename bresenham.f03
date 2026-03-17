MODULE bresenham_mod

CONTAINS

SUBROUTINE bresenham( x1,y1,x2,y2 )
IMPLICIT NONE
INTEGER, INTENT(IN) :: x1,y1,x2,y2
INTEGER :: dx, dy, i, e
INTEGER :: incx, incy, inc1, inc2
INTEGER :: x,y

dx = x2 - x1
dy = y2 - y1
IF ( dx .LT. 0 ) dx = -dx
IF ( dy .LT. 0 ) dy = -dy
incx=1
IF ( x2 .LT. x1 ) incx = -1
incy = 1
IF ( y2 .LT. y1 ) incy = -1
x = x1
y = y1

IF ( dx .GT. dy ) THEN
  PRINT *, x,y
  e = 2*dy - dx
  inc1 = 2*(dy-dx)
  inc2 = 2*dy
  DO i=0,dx-1
    IF ( e .GE. 0 ) THEN
      y = y + incy
      e = e + inc1
    ELSE
      e = e + inc2
    END IF
    x = x + incx
    PRINT *, x,y
  END DO
ELSE
  PRINT *, x,y
  e = 2*dx - dy
  inc1 = 2*(dx-dy)
  inc2 = 2*dx
  DO i=0,dy-1
    IF ( e .GE. 0 ) THEN
      x = x + incx
      e = e + inc1
    ELSE
      e = e + inc2
    END IF
    y = y + incy
    PRINT *, x,y
  END DO
END IF

END SUBROUTINE bresenham

END MODULE bresenham_mod
