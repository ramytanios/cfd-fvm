! This subroutine modify the coefficients to solve
! the transient problem. The method used is Backward Euler.
! ********************************************************
subroutine modifyCoefficientsForTransientSolving()
use variablesDeclaration
use parameters
use boundaryConditions
implicit none

integer i,j

do i = 1,m-1
	do j = 1,n-1
		aCt(i,j) = aCt(i,j) + densityField*cellVolume(i,j)/dt*2.0D0
		bCt(i,j) = bCt(i,j) + densityField*cellVolume(i,j)/dt*2.0D0 * phiPreviousTime(i,j)
	end do
end do

end subroutine 
