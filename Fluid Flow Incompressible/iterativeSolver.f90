! This subroutine solves for phi using a line by line TDMA.
! *****************************************************
subroutine solvef0rPhi()
use variablesDeclaration
use userInput
implicit none

integer i,j

! Phi field initial guess.
phi = dble(0) 

! Gradient field initial guess.
x0fGradientAtCell = dble(0) ! Gradient initial guess
y0fGradientAtCell = dble(0) ! Gradient initial guess

phiNew = phi

! Initial iteration.
iteration = 0

! Initial phi error.
phiError = 100.0D0

do while (phiError.ge.tolerance .AND. iteration.le.maxIterations)
	
	call getDiffusiveCoeff() !
	call getConvectiveCoeff() !
	call getUpwindFlux0nAllFaces() !
	call getSmartFlux0nAllFaces() !
	
	! The modification of the source term due to the deferred correction approach.
	do i = 1,m-1
		do j  = 1,n-1
			bC(i,j) = bC(i,j) + md0tePhie(i,j) + md0twPhiw(i,j) + md0tsPhis(i,j) + md0tnPhin(i,j) &
			- md0tePhieSmart(i,j) - md0twPhiwSmart(i,j) - md0tsPhisSmart(i,j) - md0tnPhinSmart(i,j)
		end do 
	end do
	
	call getSourcef0rRobin() !
	call getCrossDiffusion() !
	call sweepingTDMA()
	
	phiError = abs(maxval(phiNew-phi))
	phi = phiNew
	
	call getGradient() !
	
	iteration = iteration + 1
	
end do 
end subroutine 
