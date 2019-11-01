! This subroutine solves the system Ax=b where A is a tridiagonal 
! matrix, using the thomas algorithm.
! ***************************************************************
! Developed by Ramy Tanios - American University of Beirut.
! ************************************************

! Description of arguments.
! -------------------------
! Inputs : a: The lower diagonal of A.
! 		   b: The main diagonal of A.
! 		   c: The upper diagonal of A.
! 		   d: The source vector of Ax=b.
! 		   n: number of equations.
! Outputs : x: solution of the system Ax=b.
! ************************************************
subroutine myTDMA(a,b,c,d,x,n)
implicit none

integer n
double precision, intent(in) :: a(n),b(n),c(n),d(n)
double precision, intent(out):: x(n)
double precision :: cPrime(n),dPrime(n)
integer i,k

cPrime(1) = c(1)/b(1)
dPrime(1) = d(1)/b(1)
do i = 2,n-1
	cPrime(i) = c(i) / ( b(i) - a(i) * cPrime(i-1) )
end do

do k=2,n
	dPrime(k) = ( d(k) - a(k)*dPrime(k-1) ) / ( b(k) - a(k)*cPrime(k-1) )
end do

x(n) = dPrime(n)

do i=n-1,1,-1
	x(i) = dPrime(i) - cPrime(i) * x(i+1)
end do


end subroutine 
