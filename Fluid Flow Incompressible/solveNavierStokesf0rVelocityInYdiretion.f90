! This subroutine solves the Y-Component of Navier-Stokes
! Equation on the 2D collocated grid. 
! *****************************************************
subroutine solveNavierStokesf0rVelocityInYdiretion()
use variablesDeclaration
use userInput
implicit none 

integer i,j

! Initial iteration.
vIteration = 0

! Initial u error.
vError = 100.0D0

do while (uError.ge.tolerance .AND. uIteration.le.maxIterations)
	
	call getDiffusiveCoeffForYmomentumEq() !

	call getConvectiveCoeffForYmomentumEq() !

	call getUpwindFlux0nAllFacesForYmomentumEq() !

	call getSmartFlux0nAllFacesForYmomentumEq() !

	 !The modification of the source term due to the deferred correction approach.
	do i = 1,m-1
		do j  = 1,n-1
			bCv(i,j) = bCv(i,j) + md0teve(i,j) + md0twvw(i,j) + md0tsvs(i,j) + md0tnvn(i,j) &
				- md0teveSmart(i,j) - md0twvwSmart(i,j) - md0tsvsSmart(i,j) - md0tnvnSmart(i,j)
		end do 
	end do
! -------------------------------------

	call getCrossDiffusionForYmomentumEq() !
	
	call addToSourcePressureAndShearContributionForYmomentumEq()

	call sweepingTDMAforYmomentumEq() !

	vError = abs(maxval(vNew-v))
	v = vNew

	call getGradientOfvForYmomentumEq() !
	
	vIteration = vIteration + 1
end do

end subroutine 
