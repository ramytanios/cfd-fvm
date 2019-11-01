! This subroutine generates a 2D grid based on the 
! Transfinite Interpolation method.
! ************************************************
! Developed by Ramy Tanios - American University of Beirut.
! ************************************************

! Description of arguments.
! -------------------------
! Inputs : m,n : The number of discrete points in xi and eta directions.
! Outputs : X,Y : 2D matrices containing the x and y coordinates of grid vertices.
! ************************************************

subroutine transfiniteInterpolation(m,n,X,Y)
use boundariesParametrized ! Parametrized equations of the boundaries. 
implicit none 

integer, intent(in) :: m,n  
double precision, intent(out) :: X(m,n),Y(m,n)
double precision dumb(2)
double precision xi(m),eta(n)
integer i,j

xi(1:m) = (/ (dble(i)/(m-1), i=0,m-1) /)
eta(1:n) = (/ (dble(i)/(n-1), i=0,n-1) /)

do i = 1,m
	do j = 1,n
		dumb =  (dble(1)-eta(j))*bottomBoundary(xi(i))+eta(j)*topBoundary(xi(i))+ & 
				(dble(1)-xi(i))*leftBoundary(eta(j))+xi(i)*rightBoundary(eta(j)) &
				-(xi(i)*eta(j)*topBoundary(dble(1))+xi(i)*(dble(1)-eta(j))*& 
				bottomBoundary(dble(1))+eta(j)*(dble(1)-xi(i))*topBoundary(dble(0))+ &
				(dble(1)-xi(i))*(dble(1)-eta(j))*bottomBoundary(dble(0))) ;
		X(i,j) = dumb(1)
		Y(i,j) = dumb(2)
	end do
end do

end subroutine transfiniteInterpolation
