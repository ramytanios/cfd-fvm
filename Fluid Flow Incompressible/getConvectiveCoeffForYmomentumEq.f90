! This subroutine generates the convective coefficients 
! for eaCvh control volume baSved on the upwind scheme.
! *****************************************************
subroutine getConvectiveCoeffForYmomentumEq() 
use variablesDeclaration
use userInput
implicit none 

integer i,j
double precision md0te,md0tw,md0ts,md0tn
double precision convCoeff,convCoefff


do i = 1,m-1
	do j = 1,n-1
	
		md0te =  verFaCeDensity(i+1,j)*verFaCeArea(i+1,j) * & 
					(u(i+1,j)*x0fUnitSe(i+1,j)+v(i+1,j)*y0fUnitSe(i+1,j))
		md0tw = - verFaCeDensity(i,j)*verFaCeArea(i,j) * & 
					(u(i,j)*x0fUnitSe(i,j)+v(i,j)*y0fUnitSe(i,j))
		md0tn = horFaCeDensity(i,j+1)*horFaCeArea(i,j+1) * & 
					(u(i,j+1)*x0fUnitSn(i,j+1)+v(i,j+1)*y0fUnitSn(i,j+1))
		md0ts = -horFaCeDensity(i,j)*horFaCeArea(i,j) * & 
					(u(i,j)*x0fUnitSn(i,j)+v(i,j)*y0fUnitSn(i,j))
		
		! North.
		if (j.eq.n-1) then
			if (outletConditionNorth(i)) then
				convCoeff = dmax1(dble(0),md0tn)
				aCv(i,j) = aCv(i,j)+ convCoeff
			end if
		else if (j.ne.n-1) then
			convCoeff = -dmax1(dble(0),-md0tn)
			convCoefff = dmax1(dble(0),md0tn)
			aCv(i,j) =  aCv(i,j)+ convCoefff
			aNv(i,j) = aNv(i,j)+ convCoeff
		end if
		
		! West.
		if (i.eq.1) then
			if (inletConditionWest(j)) then
				convCoeff = md0tw*vAtInletWest(j)
				bCv(i,j) = bCv(i,j) - convCoeff 
			end if
		else if (i.ne.1) then
			convCoeff = -dmax1(dble(0),-md0tw)
			convCoefff = dmax1(dble(0),md0tw)
			aCv(i,j) = aCv(i,j) + convCoefff
			aWv(i,j) = aWv(i,j)+convCoeff
		end if
		
		! South.
		if (j.ne.1) then
			convCoeff = -dmax1(dble(0),-md0ts)
			convCoefff = dmax1(dble(0),md0ts)
			aCv(i,j) = aCv(i,j) + convCoefff
			aSv(i,j) = aSv(i,j)+convCoeff
		end if
		
		! EaSvt.
		if (i.ne.m-1) then
			convCoeff = -dmax1(dble(0),-md0te)
			convCoefff = dmax1(dble(0),md0te)
			aCv(i,j) = aCv(i,j) + convCoefff
			aEv(i,j) = aEv(i,j)+convCoeff
		end if
		
	end do
end do


end subroutine 
