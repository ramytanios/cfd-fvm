! This subroutine solves for phi using a line by line TDMA.
! *****************************************************
subroutine iterativeSolver()
use boundaryConditions
use variablesDeclaration
use parameters
implicit none

integer i,j
double precision relaxationCoeff,time

time = 0.0D0
phiPreviousTime = phiInitial

do while (time.le.t)
	aCt = aC
	bCt = bC
	call modifyCoefficientsForTransientSolving()
! ------------------------------------------------------------------
	phi = 0.0D0 ! Phi field initial guess.
	x0fGradientAtCell = 0.0D0 ! Gradient initial guess.
	y0fGradientAtCell = 0.0D0 ! Gradient initial guess.
	phiNew = phi
	iteration = 0
	phiError = 100.0D0
	do while (phiError.ge.tolerance .AND. iteration.le.maxIterations)
		bCmodified = bCt
		
		if (smartSchemeS0lver) then
			call getUpwindFlux0nAllFaces()
			call getSmartFlux0nAllFaces()
			do i = 1,m-1
				do j  = 1,n-1
					bCmodified(i,j) = bCmodified(i,j) + md0tePhie(i,j) + md0twPhiw(i,j) + md0tsPhis(i,j) + md0tnPhin(i,j) &
					- md0tePhieSmart(i,j) - md0twPhiwSmart(i,j) - md0tsPhisSmart(i,j) - md0tnPhinSmart(i,j)
				end do 
			end do
		end if
		
		call getSourcef0rRobin() !
		
		call getCrossDiffusion() !
		
		do i = 1,m-1
			do j = 1,n-1 
			    ldVerSweep(j) = aS(i,j)
	            mdVerSweep(j) = aCt(i,j)/relaxationFact0r
	            udVerSweep(j) = aN(i,j)
	            relaxationCoeff = (1-relaxationFact0r)/relaxationFact0r * phiNew(i,j) * aC(i,j)
				if (i.eq.1) then	
					sourceVerSweep(j) = bCmodified(i,j)-aE(i,j)*phiNew(i+1,j) + relaxationCoeff
				else if (i.eq.m-1) then
					sourceVerSweep(j) = bCmodified(i,j)-aW(i,j)*phiNew(i-1,j) + relaxationCoeff
				else
					sourceVerSweep(j) = bCmodified(i,j)-aE(i,j)*phiNew(i+1,j)-aW(i,j)*phiNew(i-1,j) + relaxationCoeff
				end if
			end do
			call myTDMA(ldVerSweep,mdVerSweep,udVerSweep,sourceVerSweep,phiNew(i,:),n-1)
		end do
		
		do i = 1,m-1
			do j = 1,n-1
			    ldHorSweep(i) = aW(i,j)
	            mdHorSweep(i) = aCt(i,j)/relaxationFact0r
	            udHorSweep(i) = aE(i,j)
	            relaxationCoeff = (1-relaxationFact0r)/relaxationFact0r * phiNew(i,j) * aC(i,j)
				if (j.eq.1) then	
					sourceHorSweep(i) = bCmodified(i,j)-aN(i,j)*phiNew(i,j+1) + relaxationCoeff
				else if (j.eq.n-1) then
					sourceHorSweep(i) = bCmodified(i,j)-aS(i,j)*phiNew(i,j-1) + relaxationCoeff
				else
					sourceHorSweep(i) = bCmodified(i,j)-aN(i,j)*phiNew(i,j+1)-aS(i,j)*phiNew(i,j-1) + relaxationCoeff
				end if
			end do
			call myTDMA(ldHorSweep,mdHorSweep,udHorSweep,sourceHorSweep,phiNew(:,j),m-1)
		end do
		
		phiError = abs(maxval(phiNew-phi))
		phi = phiNew
		call getGradient() !
		iteration = iteration + 1
	end do 
! ------------------------------------------------------------------
	phi = 2.0D0*phi-phiPreviousTime
	phiPreviousTime = phi
	time = time + dt
end do


end subroutine 
