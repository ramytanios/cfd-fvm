! Line by Line TDMA solving in horizontal and vertical directions.
subroutine sweepingTDMAforXmomentumEq
use variablesDeclaration
use userInput
implicit none

double precision relaxationCoeff 
integer i,j

do i = 1,m-1
	do j = 1,n-1 
	    ldVerSweep(j) = aSu(i,j)
		mdVerSweep(j) = aCu(i,j)/relaxationFact0r
		udVerSweep(j) = aNu(i,j)
		relaxationCoeff = (1-relaxationFact0r)/relaxationFact0r * uNew(i,j) * aCu(i,j)
		if (i.eq.1) then	
			sourceVerSweep(j) = bCu(i,j)-aEu(i,j)*u(i+1,j) + relaxationCoeff
		else if (i.eq.m-1) then
			sourceVerSweep(j) = bCu(i,j)-aWu(i,j)*u(i-1,j) + relaxationCoeff
		else
			sourceVerSweep(j) = bCu(i,j)-aEu(i,j)*u(i+1,j)-aWu(i,j)*u(i-1,j) + relaxationCoeff
		end if
	end do
	call myTDMA(ldVerSweep,mdVerSweep,udVerSweep,sourceVerSweep,uNew(i,:),n-1)
end do

do i = 1,m-1
	do j = 1,n-1
	    ldHorSweep(i) = aWu(i,j)
		mdHorSweep(i) = aCu(i,j)/relaxationFact0r
		udHorSweep(i) = aEu(i,j)
		relaxationCoeff = (1-relaxationFact0r)/relaxationFact0r * uNew(i,j) * aCu(i,j)
		if (j.eq.1) then	
			sourceHorSweep(i) = bCu(i,j)-aNu(i,j)*u(i,j+1) + relaxationCoeff
		else if (j.eq.n-1) then
			sourceHorSweep(i) = bCu(i,j)-aSu(i,j)*u(i,j-1) + relaxationCoeff
		else
			sourceHorSweep(i) = bCu(i,j)-aNu(i,j)*u(i,j+1)-aSu(i,j)*u(i,j-1) + relaxationCoeff
		end if
	end do
	call myTDMA(ldHorSweep,mdHorSweep,udHorSweep,sourceHorSweep,uNew(:,j),m-1)
end do

end subroutine 
