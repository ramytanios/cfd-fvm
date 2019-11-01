! The 2D diffusion convection equation.
! The finite volume method.
! **************************************************
program main 
use variablesDeclaration
implicit none 
integer i,j

m = 26
n = 26
t = 0.9 ! The desired time.
smartSchemeS0lver = .FALSE. 

call allocateArrays()

call transfiniteInterpolation(m,n,X,Y)

call geometricQuantities()
			
!call getAnalyticalSolution()

call surfaceVectors()

call getFields()

call getDiffusiveCoeff()

call getConvectiveCoeff()

call iterativeSolver()

open(unit=1,file="dataCrankN.txt")

do i = 1,m-1
	do j = 1,n-1
		write(1,*) x0fCellCentroid(i,j),achar(9),y0fCellCentroid(i,j),achar(9),phi(i,j)
	end do
end do

close(unit=1)
call system("python plot2D.py") 
!call system("python plot3D.py")

end program 
