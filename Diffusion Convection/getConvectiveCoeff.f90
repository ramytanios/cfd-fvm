! This subroutine generates the convective coefficients 
! for each control volume based on the upwind scheme.
! *****************************************************
subroutine getConvectiveCoeff() 
use variablesDeclaration
use boundaryConditions
use parameters
implicit none 

integer i,j
double precision md0te,md0tw,md0ts,md0tn
double precision convCoeff,convCoefff


do i = 1,m-1
	do j = 1,n-1
	
		md0te = md0tVerP0sDir(i+1,j)
		md0tw = -md0tVerP0sDir(i,j)
		md0tn = md0tH0rP0sDir(i,j+1)
		md0ts = -md0tH0rP0sDir(i,j)
		
		! North.
		if (j.eq.n-1) then
			if (inletCondition(1)) then
				convCoeff = md0tn*phiAtInlet0ne 
				bC(i,j) = bC(i,j) - convCoeff
			else if (outletCondition(1)) then
				convCoeff = dmax1(dble(0),md0tn)
				aC(i,j) = aC(i,j)+ convCoeff
			end if
		else if (j.ne.n-1) then
			convCoeff = -dmax1(dble(0),-md0tn)
			convCoefff = dmax1(dble(0),md0tn)
			aC(i,j) =  aC(i,j)+ convCoefff
			aN(i,j) = aN(i,j)+ convCoeff
		end if
		
		! West.
		if (i.eq.1) then
			if (inletCondition(2)) then
				convCoeff = md0tw*phiAtInletTw0 
				bC(i,j) = bC(i,j) -convCoeff
			else if (outletCondition(2)) then
				convCoeff = dmax1(dble(0),md0tw)
				aC(i,j) = aC(i,j) + convCoeff
			end if
		else if (i.ne.1) then
			convCoeff = -dmax1(dble(0),-md0tw)
			convCoefff = dmax1(dble(0),md0tw)
			aC(i,j) = aC(i,j) + convCoefff
			aW(i,j) = aW(i,j)+convCoeff
		end if
		
		! South.
		if (j.eq.1) then
			if (inletCondition(3)) then
				convCoeff = md0ts*phiAtInletThree 
				bC(i,j) = bC(i,j) - convCoeff
			else if (outletCondition(3)) then
				convCoeff = dmax1(dble(0),md0ts)
				aC(i,j) = aC(i,j) + convCoeff
			end if
		else if (j.ne.1) then
			convCoeff = -dmax1(dble(0),-md0ts)
			convCoefff = dmax1(dble(0),md0ts)
			aC(i,j) = aC(i,j) + convCoefff
			aS(i,j) = aS(i,j)+convCoeff
		end if
		
		! East.
		if (i.eq.m-1) then
			if (inletCondition(4)) then
				convCoeff = md0te*phiAtInletF0ur 
				bC(i,j) = bC(i,j) -convCoeff
			else if (outletCondition(4)) then
				convCoeff = dmax1(md0te,dble(0))
				aC(i,j) = aC(i,j)+convCoeff
			end	if
		else if (i.ne.m-1) then
			convCoeff = -dmax1(dble(0),-md0te)
			convCoefff = dmax1(dble(0),md0te)
			aC(i,j) = aC(i,j) + convCoefff
			aE(i,j) = aE(i,j)+convCoeff
		end if
		
	end do
end do


end subroutine 
