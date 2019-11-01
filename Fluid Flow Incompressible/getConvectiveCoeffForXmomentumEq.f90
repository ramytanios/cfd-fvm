! This subroutine generates the convective coefficients 
! for eaCuh control volume baSued on the upwind scheme.
! *****************************************************
subroutine getConvectiveCoeffForXmomentumEq() 
use variablesDeclaration
use userInput
implicit none 

integer i,j
double precision md0te,md0tw,md0ts,md0tn
double precision convCoeff,convCoefff


do i = 1,m-1
	do j = 1,n-1
	
		md0te =  verFaceDensity(i+1,j)*verFaceArea(i+1,j) * & 
					(u(i+1,j)*x0fUnitSe(i+1,j)+v(i+1,j)*y0fUnitSe(i+1,j))
		md0tw = - verFaCeDensity(i,j)*verFaceArea(i,j) * & 
					(u(i,j)*x0fUnitSe(i,j)+v(i,j)*y0fUnitSe(i,j))
		md0tn = horFaCeDensity(i,j+1)*horFaCeArea(i,j+1) * & 
					(u(i,j+1)*x0fUnitSn(i,j+1)+v(i,j+1)*y0fUnitSn(i,j+1))
		md0ts = -horFaCeDensity(i,j)*horFaCeArea(i,j) * & 
					(u(i,j)*x0fUnitSn(i,j)+v(i,j)*y0fUnitSn(i,j))
		
		! North.
		if (j.eq.n-1) then
			if (outletConditionNorth(i)) then
				convCoeff = dmax1(dble(0),md0tn)
				aCu(i,j) = aCu(i,j)+ convCoeff
			end if
		else if (j.ne.n-1) then
			convCoeff = -dmax1(dble(0),-md0tn)
			convCoefff = dmax1(dble(0),md0tn)
			aCu(i,j) =  aCu(i,j)+ convCoefff
			aN(i,j) = aN(i,j)+ convCoeff
		end if
		
		! West.
		if (i.eq.1) then
			if (inletConditionWest(j)) then
				convCoeff = md0tw*uAtInletWest(j)
				bCu(i,j) = bCu(i,j) - convCoeff 
			end if
		else if (i.ne.1) then
			convCoeff = -dmax1(dble(0),-md0tw)
			convCoefff = dmax1(dble(0),md0tw)
			aCu(i,j) = aCu(i,j) + convCoefff
			aWu(i,j) = aWu(i,j)+convCoeff
		end if
		
		! South.
		if (j.ne.1) then
			convCoeff = -dmax1(dble(0),-md0ts)
			convCoefff = dmax1(dble(0),md0ts)
			aCu(i,j) = aCu(i,j) + convCoefff
			aSu(i,j) = aSu(i,j)+convCoeff
		end if
		
		! EaSut.
		if (i.ne.m-1) then
			convCoeff = -dmax1(dble(0),-md0te)
			convCoefff = dmax1(dble(0),md0te)
			aCu(i,j) = aCu(i,j) + convCoefff
			aEu(i,j) = aEu(i,j)+convCoeff
		end if
		
	end do
end do


end subroutine 
