! Line by Line TDMA solving in horizontal and vertical directions.
subroutine sweepingTDMAforYmomentumEq
use variablesDeclaration
use userInput
implicit none

double precision relaxationCoeff 
integer i,j

do i = 1,m-1
	do j = 1,n-1 
	    ldVerSweep(j) = aSv(i,j)
		mdVerSweep(j) = aCv(i,j)/relaxationFact0r
		udVerSweep(j) = aNv(i,j)
		relaxationCoeff = (1-relaxationFact0r)/relaxationFact0r * vNew(i,j) * aCv(i,j)
		if (i.eq.1) then	
			sourceVerSweep(j) = bCv(i,j)-aEv(i,j)*v(i+1,j) + relaxationCoeff
		else if (i.eq.m-1) then
			sourceVerSweep(j) = bCv(i,j)-aWv(i,j)*v(i-1,j) + relaxationCoeff
		else
			sourceVerSweep(j) = bCv(i,j)-aEv(i,j)*v(i+1,j)-aWv(i,j)*v(i-1,j) + relaxationCoeff
		end if
	end do
	call myTDMA(ldVerSweep,mdVerSweep,udVerSweep,sourceVerSweep,vNew(i,:),n-1)
end do

do i = 1,m-1
	do j = 1,n-1
	    ldHorSweep(i) = aWv(i,j)
		mdHorSweep(i) = aCv(i,j)/relaxationFact0r
		udHorSweep(i) = aEv(i,j)
		relaxationCoeff = (1-relaxationFact0r)/relaxationFact0r * vNew(i,j) * aCv(i,j)
		if (j.eq.1) then	
			sourceHorSweep(i) = bCv(i,j)-aNv(i,j)*v(i,j+1) + relaxationCoeff
		else if (j.eq.n-1) then
			sourceHorSweep(i) = bCv(i,j)-aSv(i,j)*v(i,j-1) + relaxationCoeff
		else
			sourceHorSweep(i) = bCv(i,j)-aNv(i,j)*v(i,j+1)-aSv(i,j)*v(i,j-1) + relaxationCoeff
		end if
	end do
	call myTDMA(ldHorSweep,mdHorSweep,udHorSweep,sourceHorSweep,vNew(:,j),m-1)
end do

end subroutine 
