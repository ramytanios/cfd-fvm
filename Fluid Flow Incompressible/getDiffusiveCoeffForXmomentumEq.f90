! This subroutine generates the diffusive coefficients 
! for eaCuh control volume, for the X-Component of the 
! Navier-Stokes equation.
! *****************************************************
subroutine getDiffusiveCoeffForXmomentumEq()
use variablesDeclaration
use userInput
implicit none

integer i,j
double precision dCb,aDumb,term1,term2,rDumb

do i = 1,m-1
	do j = 1,n-1
	
	if (i.eq.m-1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i+1,j))**2 + (y0fCellCentroid(i,j)-y0fVerFaceCenter(i+1,j))**2 )
		aDumb = -viscosity*normEe(i+1,j) / dCb
		aEu(i,j) = 0.0D0
		bCu(i,j) = -aDumb * vWall
		aCu(i,j) = - aDumb
	else if (i.ne.m-1) then
		aEu(i,j) = -viscosity*normEe(i+1,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i+1,j))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i+1,j))**2 )
		aCu(i,j) =  - aEu(i,j)
		bCu(i,j) = 0.0D0
	end if
! --------------------------
		if (i.eq.1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i,j))**2 + (y0fCellCentroid(i,j)-y0fVerFaceCenter(i,j))**2 )
		if (InletConditionWest(j)) then
			aDumb = -viscosity*normEe(i,j) / dCb
			aWu(i,j) = 0.0D0
			bCu(i,j) = bCu(i,j) -aDumb * uAtInletWest(j)
			aCu(i,j) = aCu(i,j) - aDumb
		else if (wallConditionWest(j)) then
			aDumb = -viscosity*normEe(i,j) / dCb
			aWu(i,j) = 0.0D0
			bCu(i,j) = bCu(i,j) -aDumb * vWall
			aCu(i,j) = aCu(i,j) - aDumb
		end if
	else if (i.ne.1) then
		aWu(i,j) = -viscosity*normEe(i,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i-1,j))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i-1,j))**2 )
		aCu(i,j) = aCu(i,j) - aWu(i,j)
	end if
! --------------------------
		if (j.eq.1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fHorFaceCenter(i,j))**2 + (y0fCellCentroid(i,j)-y0fHorFaceCenter(i,j))**2 )
			aDumb = -viscosity*normEn(i,j) / dCb
			aSu(i,j) = 0.0D0
			bCu(i,j) = bCu(i,j) -aDumb * vWall
			aCu(i,j) = aCu(i,j) - aDumb
	else if (j.ne.1) then
		aSu(i,j) = -viscosity*normEn(i,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i,j-1))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i,j-1))**2 )
		aCu(i,j) = aCu(i,j) - aSu(i,j)
	end if
! --------------------------
		if (j.eq.n-1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fHorFaceCenter(i,j+1))**2 + (y0fCellCentroid(i,j)-y0fHorFaceCenter(i,j+1))**2 )
		if (wallConditionNorth(i)) then
			aDumb = -viscosity*normEn(i,j+1) / dCb
			aNu(i,j) = 0.0D0
			bCu(i,j) = bCu(i,j) -aDumb * vWall
			aCu(i,j) = aCu(i,j) - aDumb
		else if (outletConditionNorth(i)) then
			aNu(i,j) = 0.0D0
			bCu(i,j) = bCu(i,j) + 0.0D0
	else if (j.ne.n-1) then
		aNu(i,j) = -viscosity*normEn(i,j+1) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i,j+1))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i,j+1))**2 )
		aCu(i,j) = aCu(i,j) - aNu(i,j)
	end if
end if
	bCu(i,j) = bCu(i,j) + sourceField(i,j) * cellVolume(i,j)
	
	end do
end do

end subroutine 
