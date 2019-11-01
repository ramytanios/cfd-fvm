program main 
use variablesDeclaration
use userInput
implicit none 
integer i,j

! Number of discrete points in the eta direction.
m = 5
! Number of discrete points in the xi direction. 
n = 5

call allocateMemoryForArrays()

call transfiniteInterpolation(m,n,X,Y)

call constructBoundaryConditions()

call geometricQuantities()
			
call surfaceVectors()

call getFields()


pCorrection = 100.0D0

!Initial pressure field guess. 
p = 20

! u and v initial guesses.
u = 2.0D0
v = 2.0D0

! Initial velocity gradient field guess.
x0fGradu = 1
y0fGradu = 1
x0fGradv = 1
y0fGradv = 1

do while (abs(maxval(pCorrection)).gt.tolerance)
	call solveNavierStokesf0rVelocityInXdiretion()
	call solveNavierStokesf0rVelocityInYdiretion()
	call rhieCh0wInterpolation()
	call solvePressureCorrectionEquation()
	call correctPressureAndVelocity()
	call computePressureGradient()
end do


call solvef0rPhi()


! Write data to a text file. 
open(unit=1,file="test.txt")
do i = 1,m-1
	do j = 1,n-1
		write(1,*) x0fCellCentroid(i,j),",",y0fCellCentroid(i,j),",",v(i,j)
	end do
end do
close(unit=1)

! Write the number of cells for python plotting. 
open(unit=2,file="mMinus0ne.txt")
write(2,*) m-1
close(unit=2)

call system("python plot2D.py") 

end program 
