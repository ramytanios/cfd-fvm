module boundariesParametrized
implicit none 
! This module contains the parametrized equations of the 
! four boundaries of the 2D domain. 
! To be changed when necessary.
! ************************************************
! The parameter s lies in the interval [0 1].
contains 

function topBoundary(s)
	double precision x,y
	double precision, intent(in) :: s
	double precision topBoundary(2)
	x = -1+4*s
	y = 3+2*s
	topBoundary = [ x , y ]
end function

function leftBoundary(s)
	double precision x,y
	double precision, intent(in) :: s
	double precision leftBoundary(2)
	x = -1*s
	y = 3*s
	leftBoundary = [ x , y ]
end function

function bottomBoundary(s)
	double precision x,y
	double precision, intent(in) :: s
	double precision bottomBoundary(2)
	x = s
	y = 0
	bottomBoundary = [ x , y ]
end function

function rightBoundary(s)
	double precision x,y
	double precision, intent(in) :: s
	double precision rightBoundary(2)
	x = 2+s
	y = 5*s
	rightBoundary = [ x , y ]
end function


end module 
