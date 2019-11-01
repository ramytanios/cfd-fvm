! This subroutine solves the X-Component of Navier-Stokes
! Equation on the 2D collocated grid. 
! *****************************************************
subroutine solveNavierStokesf0rVelocityInXdiretion()
use variablesDeclaration
use userInput
implicit none 

integer i,j

! Initial iteration.
uIteration = 0

! Initial u error.
uError = 100.0D0

do while (uError.ge.tolerance .AND. uIteration.le.maxIterations)
	
	call getDiffusiveCoeffForXmomentumEq() !

	call getConvectiveCoeffForXmomentumEq() !

	call getUpwindFlux0nAllFacesForXmomentumEq() !

	call getSmartFlux0nAllFacesForXmomentumEq() !

	 !The modification of the source term due to the deferred correction approach.
	do i = 1,m-1
		do j  = 1,n-1
			bCu(i,j) = bCu(i,j) + md0teue(i,j) + md0twuw(i,j) + md0tsus(i,j) + md0tnun(i,j) &
				- md0teueSmart(i,j) - md0twuwSmart(i,j) - md0tsusSmart(i,j) - md0tnunSmart(i,j)
		end do 
	end do
! -------------------------------------

	call getCrossDiffusionForXmomentumEq() !
	
	call addToSourcePressureAndShearContributionForXmomentumEq()

	call sweepingTDMAforXmomentumEq() !

	uError = abs(maxval(uNew-u))
	u = uNew

	call getGradientOfuForXmomentumEq() !
	
	uIteration = uIteration + 1
end do

end subroutine 
