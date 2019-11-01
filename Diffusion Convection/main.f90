! The 2D diffusion convection equation discretization.
! using the finite volume method.
! **************************************************
program main 
use variablesDeclaration
implicit none 
integer i,j

m = 30
n = 30

smartSchemeS0lver = .TRUE.

call allocateArrays()

call transfiniteInterpolation(m,n,X,Y)

call geometricQuantities()
			
call surfaceVectors()

call getFields()

call getDiffusiveCoeff()

call getConvectiveCoeff()

call iterativeSolver()

open(unit=1,file="test.txt")
do i = 1,m-1
	do j = 1,n-1
		write(1,*) x0fCellCentroid(i,j),",",y0fCellCentroid(i,j),",",phi(i,j)
	end do
end do
close(unit=1)

call system("python plot2D.py") 
!call system("python plot3D.py")

end program 
