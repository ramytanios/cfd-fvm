! Line by Line TDMA solving in horizontal and vertical directions.
subroutine sweepingTDMAforPressureCorrection()
use variablesDeclaration
use userInput
implicit none

double precision relaxationCoeff 
integer i,j

do i = 1,m-1
	do j = 1,n-1 
	    ldVerSweep(j) = aS(i,j)
		mdVerSweep(j) = aC(i,j)/relaxationFact0r
		udVerSweep(j) = aN(i,j)
		relaxationCoeff = (1-relaxationFact0r)/relaxationFact0r * pCorrectionNew(i,j) * aC(i,j)
		if (i.eq.1) then	
			sourceVerSweep(j) = bC(i,j)-aE(i,j)*pCorrection(i+1,j) + relaxationCoeff
		else if (i.eq.m-1) then
			sourceVerSweep(j) = bC(i,j)-aW(i,j)*pCorrection(i-1,j) + relaxationCoeff
		else
			sourceVerSweep(j) = bC(i,j)-aE(i,j)*pCorrection(i+1,j)-aW(i,j)*pCorrection(i-1,j) + relaxationCoeff
		end if
	end do
	call myTDMA(ldVerSweep,mdVerSweep,udVerSweep,sourceVerSweep,pCorrectionNew(i,:),n-1)
end do

do i = 1,m-1
	do j = 1,n-1
	    ldHorSweep(i) = aW(i,j)
		mdHorSweep(i) = aC(i,j)/relaxationFact0r
		udHorSweep(i) = aE(i,j)
		relaxationCoeff = (1-relaxationFact0r)/relaxationFact0r * pCorrectionNew(i,j) * aC(i,j)
		if (j.eq.1) then	
			sourceHorSweep(i) = bC(i,j)-aN(i,j)*pCorrection(i,j+1) + relaxationCoeff
		else if (j.eq.n-1) then
			sourceHorSweep(i) = bC(i,j)-aS(i,j)*pCorrection(i,j-1) + relaxationCoeff
		else
			sourceHorSweep(i) = bC(i,j)-aN(i,j)*pCorrection(i,j+1)-aS(i,j)*pCorrection(i,j-1) + relaxationCoeff
		end if
	end do
	call myTDMA(ldHorSweep,mdHorSweep,udHorSweep,sourceHorSweep,pCorrectionNew(:,j),m-1)
end do

end subroutine 
