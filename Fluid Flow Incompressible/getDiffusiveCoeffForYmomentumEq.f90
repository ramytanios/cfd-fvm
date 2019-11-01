! This subroutine generates the diffusive coefficients 
! for eaCvh control volume, for the Y-Component of the 
! Navier-Stokes equation.
! *****************************************************
subroutine getDiffusiveCoeffForYmomentumEq()
use variablesDeclaration
use userInput
implicit none

integer i,j
double precision dCb,aDumb,term1,term2,rDumb

do i = 1,m-1
	do j = 1,n-1
	
	if (i.eq.m-1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaCeCenter(i+1,j))**2 + (y0fCellCentroid(i,j)-y0fVerFaCeCenter(i+1,j))**2 )
		aDumb = -viscosity*normEe(i+1,j) / dCb
		aEv(i,j) = 0.0D0
		bCv(i,j) = -aDumb * vWall
		aCv(i,j) = - aDumb
	else if (i.ne.m-1) then
		aEv(i,j) = -viscosity*normEe(i+1,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i+1,j))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i+1,j))**2 )
		aCv(i,j) =  - aEv(i,j)
		bCv(i,j) = 0.0D0
	end if
! --------------------------
		if (i.eq.1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaCeCenter(i,j))**2 + (y0fCellCentroid(i,j)-y0fVerFaCeCenter(i,j))**2 )
		if (InletConditionWest(j)) then
			aDumb = -viscosity*normEe(i,j) / dCb
			aWv(i,j) = 0.0D0
			bCv(i,j) = bCv(i,j) -aDumb * vAtInletWest(j)
			aCv(i,j) = aCv(i,j) - aDumb
		else if (wallConditionWest(j)) then
			aDumb = -viscosity*normEe(i,j) / dCb
			aWv(i,j) = 0.0D0
			bCv(i,j) = bCv(i,j) -aDumb * vWall
			aCv(i,j) = aCv(i,j) - aDumb
		end if
	else if (i.ne.1) then
		aWv(i,j) = -viscosity*normEe(i,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i-1,j))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i-1,j))**2 )
		aCv(i,j) = aCv(i,j) - aWv(i,j)
	end if
! --------------------------
		if (j.eq.1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fHorFaCeCenter(i,j))**2 + (y0fCellCentroid(i,j)-y0fHorFaCeCenter(i,j))**2 )
			aDumb = -viscosity*normEn(i,j) / dCb
			aSv(i,j) = 0.0D0
			bCv(i,j) = bCv(i,j) -aDumb * vWall
			aCv(i,j) = aCv(i,j) - aDumb
	else if (j.ne.1) then
		aSv(i,j) = -viscosity*normEn(i,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i,j-1))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i,j-1))**2 )
		aCv(i,j) = aCv(i,j) - aSv(i,j)
	end if
! --------------------------
		if (j.eq.n-1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fHorFaCeCenter(i,j+1))**2 + (y0fCellCentroid(i,j)-y0fHorFaCeCenter(i,j+1))**2 )
		if (wallConditionNorth(i)) then
			aDumb = -viscosity*normEn(i,j+1) / dCb
			aNv(i,j) = 0.0D0
			bCv(i,j) = bCv(i,j) -aDumb * vWall
			aCv(i,j) = aCv(i,j) - aDumb
		else if (outletConditionNorth(i)) then
			aNv(i,j) = 0.0D0
			bCv(i,j) = bCv(i,j) + 0.0D0
	else if (j.ne.n-1) then
		aNv(i,j) = -viscosity*normEn(i,j+1) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i,j+1))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i,j+1))**2 )
		aCv(i,j) = aCv(i,j) - aNv(i,j)
	end if
end if
	bCv(i,j) = bCv(i,j) + sourceField(i,j) * cellVolume(i,j)
	
	end do
end do

end subroutine 
