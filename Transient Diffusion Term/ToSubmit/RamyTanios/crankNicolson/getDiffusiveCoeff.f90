! This subroutine generates the diffusive coefficients 
! for each control volume.
! *****************************************************
subroutine getDiffusiveCoeff()
use variablesDeclaration
use parameters
use boundaryConditions
implicit none 

integer i,j
double precision dCb,aDumb,term1,term2,rDumb

do i = 1,m-1
	do j = 1,n-1
	
	if (i.eq.m-1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i+1,j))**2 + (y0fCellCentroid(i,j)-y0fVerFaceCenter(i+1,j))**2 )
		if (dirichlet(4)) then
			aDumb = -verInterfaceDiff(i+1,j)*normEe(i+1,j) / dCb
			aE(i,j) = aE(i,j) + 0.0D0
			bC(i,j) = bC(i,j) -aDumb * dirichletValue(4)
			aC(i,j) = aC(i,j) - aDumb
		else if (vonNeumann(4)) then
			aE(i,j) = aE(i,j)+ 0.0D0
			bC(i,j) = bC(i,j) + verInterfaceDiff(i+1,j)*verFaceArea(i+1,j) * &
									( x0fNeumannGrad(4)*x0fUnitSe(i+1,j) + y0fNeumannGrad(4)*y0fUnitSe(i+1,j) )
		else if (robin(4)) then
			aE(i,j) = aE(i,j)+0.0D0
			term1 = hInf(4)*verFaceArea(i+1,j)
			term2 = verInterfaceDiff(i+1,j)*normEe(i+1,j) / dCb
			rDumb = term1*term2 / (term1+term2)
			aC(i,j) = aC(i,j) + rDumb
			bC(i,j) =  bC(i,j) + rDumb*phiInf(4)
		end if
	else if (i.ne.m-1) then
		aE(i,j) = aE(i,j)-verInterfaceDiff(i+1,j)*normEe(i+1,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i+1,j))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i+1,j))**2 )
		aC(i,j) = aC(i,j) - aE(i,j)
	end if
! --------------------------
		if (i.eq.1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i,j))**2 + (y0fCellCentroid(i,j)-y0fVerFaceCenter(i,j))**2 )
		if (dirichlet(2)) then
			aDumb = -verInterfaceDiff(i,j)*normEe(i,j) / dCb
			aW(i,j) = aW(i,j)+0.0D0
			bC(i,j) = bC(i,j) -aDumb * dirichletValue(2)
			aC(i,j) = aC(i,j) - aDumb
		else if (vonNeumann(2)) then
			aW(i,j) = aW(i,j)+0.0D0
			bC(i,j) = bC(i,j) + verInterfaceDiff(i,j)*verFaceArea(i,j) * &
									( x0fNeumannGrad(2)*(-x0fUnitSe(i,j)) + y0fNeumannGrad(2)*(-y0fUnitSe(i,j)) )
	
		else if (robin(2)) then
			aW(i,j) = aW(i,j)+0.0D0
			term1 = hInf(2)*verFaceArea(i,j)
			term2 = verInterfaceDiff(i,j)*normEe(i,j) / dCb
			rDumb = term1*term2 / (term1+term2)
			aC(i,j) = aC(i,j) + rDumb
			bC(i,j) =  bC(i,j) + rDumb*phiInf(2)
		end if
	else if (i.ne.1) then
		aW(i,j) = aW(i,j)-verInterfaceDiff(i,j)*normEe(i,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i-1,j))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i-1,j))**2 )
		aC(i,j) = aC(i,j) - aW(i,j)
	end if
! --------------------------
	
		if (j.eq.1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fHorFaceCenter(i,j))**2 + (y0fCellCentroid(i,j)-y0fHorFaceCenter(i,j))**2 )
		if (dirichlet(3)) then
			aDumb = -horInterfaceDiff(i,j)*normEn(i,j) / dCb
			aS(i,j) = aS(i,j)+0.0D0
			bC(i,j) = bC(i,j) -aDumb * dirichletValue(3)
			aC(i,j) = aC(i,j) - aDumb
		else if (vonNeumann(3)) then
			aS(i,j) = aS(i,j)+0.0D0
			bC(i,j) = bC(i,j) + horInterfaceDiff(i,j)*horFaceArea(i,j) * &
									( x0fNeumannGrad(3)*(-x0fUnitSn(i,j)) + y0fNeumannGrad(3)*(-y0fUnitSn(i,j)) )
	
		else if (robin(3)) then
			aS(i,j) = aS(i,j)+0.0D0
			term1 = hInf(3)*horFaceArea(i,j)
			term2 = horInterfaceDiff(i,j)*normEn(i,j) / dCb
			rDumb = term1*term2 / (term1+term2)
			aC(i,j) = aC(i,j) + rDumb
			bC(i,j) =  bC(i,j) + rDumb*phiInf(3)
		end if
	else if (j.ne.1) then
		aS(i,j) = aS(i,j)-horInterfaceDiff(i,j)*normEn(i,j) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i,j-1))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i,j-1))**2 )
		aC(i,j) = aC(i,j) - aS(i,j)
	end if
	
! --------------------------
		if (j.eq.n-1) then
		dCb = sqrt( (x0fCellCentroid(i,j)-x0fHorFaceCenter(i,j+1))**2 + (y0fCellCentroid(i,j)-y0fHorFaceCenter(i,j+1))**2 )
		if (dirichlet(1)) then
			aDumb = -horInterfaceDiff(i,j+1)*normEn(i,j+1) / dCb
			aN(i,j) = aN(i,j)+0.0D0
			bC(i,j) = bC(i,j) -aDumb * dirichletValue(1)
			aC(i,j) = aC(i,j) - aDumb
		else if (vonNeumann(1)) then
			aN(i,j) = aN(i,j)+0.0D0
			bC(i,j) = bC(i,j) + horInterfaceDiff(i,j+1)*horFaceArea(i,j+1) * &
									( x0fNeumannGrad(1)*(x0fUnitSn(i,j+1)) + y0fNeumannGrad(1)*(y0fUnitSn(i,j+1)) )
	
		else if (robin(1)) then
			aN(i,j) = aN(i,j)+0.0D0
			term1 = hInf(1)*horFaceArea(i,j+1)
			term2 = horInterfaceDiff(i,j+1)*normEn(i,j+1) / dCb
			rDumb = term1*term2 / (term1+term2)
			aC(i,j) = aC(i,j) + rDumb
			bC(i,j) =  bC(i,j) + rDumb*phiInf(1)
		end if
	else if (j.ne.n-1) then
		aN(i,j) = aN(i,j)-horInterfaceDiff(i,j+1)*normEn(i,j+1) / sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i,j+1))**2 & 
															+ (y0fCellCentroid(i,j)-y0fCellCentroid(i,j+1))**2 )
		aC(i,j) = aC(i,j) - aN(i,j)
	end if
	
	bC(i,j) = bC(i,j) + sourceField(i,j) * cellVolume(i,j)
	
	end do
end do

end subroutine 
