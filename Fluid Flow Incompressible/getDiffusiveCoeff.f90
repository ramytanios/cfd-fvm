! This subroutine generates the diffusive coefficients 
! for each control volume.
! *****************************************************
subroutine getDiffusiveCoeff()
use variablesDeclaration
use userInput
implicit none 

integer i,j
double precision dCb,aDumb,term1,term2,rDumb

do i = 1,m-1
	do j = 1,n-1
	
	if (i.eq.m-1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i+1,j))**2 + (y0fCellCentroid(i,j)-y0fVerFaceCenter(i+1,j))**2 )
		if (dirichletEast(j)) then
			aDumb = -verInterfaceDiff(i+1,j)*normEe(i+1,j) / dCb
			aE(i,j) = 0.0D0
			bC(i,j) = -aDumb * dirichletValueEast(j)
			aC(i,j) = - aDumb
		else if (vonNeumannEast(j)) then
			aE(i,j) = 0.0D0
			bC(i,j) =  verInterfaceDiff(i+1,j)*verFaceArea(i+1,j) * &
									( x0fNeumannGradEast(j)*x0fUnitSe(i+1,j) + y0fNeumannGradEast(j)*y0fUnitSe(i+1,j) )
			aC(i,j) = 0.0D0
		else if (robinEast(j)) then
			aE(i,j) = 0.0D0
			term1 = hInfEast(j)*verFaceArea(i+1,j)
			term2 = verInterfaceDiff(i+1,j)*normEe(i+1,j) / dCb
			rDumb = term1*term2 / (term1+term2)
			aC(i,j) =  rDumb
			bC(i,j) =  rDumb*phiInfEast(j)
		end if
	else if (i.ne.m-1) then
		aE(i,j) = -verInterfaceDiff(i+1,j)*normEe(i+1,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i+1,j))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i+1,j))**2 )
		aC(i,j) =  - aE(i,j)
		bC(i,j) = 0.0D0
	end if
! --------------------------
		if (i.eq.1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i,j))**2 + (y0fCellCentroid(i,j)-y0fVerFaceCenter(i,j))**2 )
		if (dirichletWest(j)) then
			aDumb = -verInterfaceDiff(i,j)*normEe(i,j) / dCb
			aW(i,j) = 0.0D0
			bC(i,j) = bC(i,j) -aDumb * dirichletValueWest(j)
			aC(i,j) = aC(i,j) - aDumb
		else if (vonNeumannWest(j)) then
			aW(i,j) = 0.0D0
			bC(i,j) = bC(i,j) + verInterfaceDiff(i,j)*verFaceArea(i,j) * &
									( x0fNeumannGradWest(j)*(-x0fUnitSe(i,j)) + y0fNeumannGradWest(j)*(-y0fUnitSe(i,j)) )
	
		else if (robinWest(j)) then
			aW(i,j) = 0.0D0
			term1 = hInfWest(j)*verFaceArea(i,j)
			term2 = verInterfaceDiff(i,j)*normEe(i,j) / dCb
			rDumb = term1*term2 / (term1+term2)
			aC(i,j) = aC(i,j) + rDumb
			bC(i,j) =  bC(i,j) + rDumb*phiInfWest(j)
		end if
	else if (i.ne.1) then
		aW(i,j) = -verInterfaceDiff(i,j)*normEe(i,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i-1,j))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i-1,j))**2 )
		aC(i,j) = aC(i,j) - aW(i,j)
	end if
! --------------------------
	
		if (j.eq.1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fHorFaceCenter(i,j))**2 + (y0fCellCentroid(i,j)-y0fHorFaceCenter(i,j))**2 )
		if (dirichletSouth(i)) then
			aDumb = -horInterfaceDiff(i,j)*normEn(i,j) / dCb
			aS(i,j) = 0.0D0
			bC(i,j) = bC(i,j) -aDumb * dirichletValueSouth(i)
			aC(i,j) = aC(i,j) - aDumb
		else if (vonNeumannSouth(i)) then
			aS(i,j) = 0.0D0
			bC(i,j) = bC(i,j) + horInterfaceDiff(i,j)*horFaceArea(i,j) * &
									( x0fNeumannGradSouth(i)*(-x0fUnitSn(i,j)) + y0fNeumannGradSouth(i)*(-y0fUnitSn(i,j)) )
	
		else if (robinSouth(i)) then
			aS(i,j) = 0.0D0
			term1 = hInfSouth(i)*horFaceArea(i,j)
			term2 = horInterfaceDiff(i,j)*normEn(i,j) / dCb
			rDumb = term1*term2 / (term1+term2)
			aC(i,j) = aC(i,j) + rDumb
			bC(i,j) =  bC(i,j) + rDumb*phiInfSouth(i)
		end if
	else if (j.ne.1) then
		aS(i,j) = -horInterfaceDiff(i,j)*normEn(i,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i,j-1))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i,j-1))**2 )
		aC(i,j) = aC(i,j) - aS(i,j)
	end if
	
! --------------------------
		if (j.eq.n-1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fHorFaceCenter(i,j+1))**2 + (y0fCellCentroid(i,j)-y0fHorFaceCenter(i,j+1))**2 )
		if (dirichletNorth(i)) then
			aDumb = -horInterfaceDiff(i,j+1)*normEn(i,j+1) / dCb
			aN(i,j) = 0.0D0
			bC(i,j) = bC(i,j) -aDumb * dirichletValueNorth(i)
			aC(i,j) = aC(i,j) - aDumb
		else if (vonNeumannNorth(i)) then
			aN(i,j) = 0.0D0
			bC(i,j) = bC(i,j) + horInterfaceDiff(i,j+1)*horFaceArea(i,j+1) * &
									( x0fNeumannGradNorth(i)*(x0fUnitSn(i,j+1)) + y0fNeumannGradNorth(i)*(y0fUnitSn(i,j+1)) )
	
		else if (robinNorth(i)) then
			aN(i,j) = 0.0D0
			term1 = hInfNorth(i)*horFaceArea(i,j+1)
			term2 = horInterfaceDiff(i,j+1)*normEn(i,j+1) / dCb
			rDumb = term1*term2 / (term1+term2)
			aC(i,j) = aC(i,j) + rDumb
			bC(i,j) =  bC(i,j) + rDumb*phiInfNorth(i)
		end if
	else if (j.ne.n-1) then
		aN(i,j) = -horInterfaceDiff(i,j+1)*normEn(i,j+1) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i,j+1))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i,j+1))**2 )
		aC(i,j) = aC(i,j) - aN(i,j)
	end if
	
	bC(i,j) = bC(i,j) + sourceField(i,j) * cellVolume(i,j)
	
	end do
end do

end subroutine 
