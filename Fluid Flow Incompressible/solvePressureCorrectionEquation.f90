! This subroutine solves the pressure correction equation.
! *******************************************************
subroutine solvePressureCorrectionEquation()
use variablesDeclaration
use userInput
implicit none 

integer i,j

double precision Df,bar_Duf,bar_Dvf,Sfx,Sfy,dxCF,dyCF,gC,gF

do i = 1,m-1
	do j = 1,n-1
	
	! North Coefficients
	if (j.eq.n-1) then 
		if (wallConditionNorth(i)) then
			aN(i,j) = 1
		else if (outletConditionNorth(i)) then
			aN(i,j) = 0
		end if
	else
		gF = NorthGeomFactor(i,j+1)
		gC = 1 - gF
		dxCF = NormEn(i,j+1) * x0fUnitEn(i,j+1)
		dyCF = NormEn(i,j+1) * y0fUnitEn(i,j+1)
		Sfx = horFaceArea(i,j+1) * x0fUnitSn(i,j+1)
		Sfy = horFaceArea(i,j+1) * y0fUnitSn(i,j+1)
		bar_Duf = gC * (cellVolume(i,j)/aCu(i,j) ) + gF * (cellVolume(i,j+1)/aCu(i,j+1)) 
		bar_Duf = gC * (cellVolume(i,j)/aCv(i,j) ) + gF * (cellVolume(i,j+1)/aCv(i,j+1)) 
		Df = (dxCF*bar_Duf*Sfx + dyCF*bar_Dvf*Sfy ) / (dxCF**2 + dyCF**2) 
		aN(i,j) = -densityField * Df
	end if
	
	! West Coefficients
	if (i.eq.1) then
		aW(i,j) = 1
	else
		gF = WestGeomFactor(i,j)
		gC = 1 - gF
		dxCF = -normEe(i,j) * x0fUnitEe(i,j)
		dyCF = -normEe(i,j) * y0fUnitEe(i,j)
		Sfx = -verFaceArea(i,j) * x0fUnitSe(i,j)
		Sfy = -verFaceArea(i,j) * y0fUnitSe(i,j)
		bar_Duf = gC * (cellVolume(i,j)/aCu(i,j) ) + gF * (cellVolume(i-1,j)/aCu(i-1,j)) 
		bar_Duf = gC * (cellVolume(i,j)/aCv(i,j) ) + gF * (cellVolume(i-1,j)/aCv(i-1,j)) 
		Df = (dxCF*bar_Duf*Sfx + dyCF*bar_Dvf*Sfy ) / (dxCF**2 + dyCF**2) 
		aW(i,j) = -densityField * Df
	end if
	
	! South Coefficients
	if (j.eq.1) then
		aS(i,j) = 1
	else
		gF = SouthGeomFactor(i,j)
		gC = 1 - gF
		dxCF = -NormEn(i,j) * x0fUnitEn(i,j)
		dyCF = -NormEn(i,j) * y0fUnitEn(i,j)
		Sfx = -horFaceArea(i,j) * x0fUnitSn(i,j)
		Sfy = -horFaceArea(i,j) * y0fUnitSn(i,j)
		bar_Duf = gC * (cellVolume(i,j)/aCu(i,j) ) + gF * (cellVolume(i,j-1)/aCu(i,j-1)) 
		bar_Duf = gC * (cellVolume(i,j)/aCv(i,j) ) + gF * (cellVolume(i,j-1)/aCv(i,j-1)) 
		Df = (dxCF*bar_Duf*Sfx + dyCF*bar_Dvf*Sfy ) / (dxCF**2 + dyCF**2) 
		aS(i,j) = -densityField * Df
	end if
	
	! East Coefficients
	if (i.eq.m-1) then
		aE(i,j) = 1
	else
		gF = EastGeomFactor(i+1,j)
		gC = 1 - gF
		dxCF = normEe(i+1,j) * x0fUnitEe(i+1,j)
		dyCF = normEe(i+1,j) * y0fUnitEe(i+1,j)
		Sfx = verFaceArea(i+1,j) * x0fUnitSe(i+1,j)
		Sfy = verFaceArea(i+1,j) * y0fUnitSe(i+1,j)
		bar_Duf = gC * (cellVolume(i,j)/aCu(i,j) ) + gF * (cellVolume(i+1,j)/aCu(i+1,j)) 
		bar_Duf = gC * (cellVolume(i,j)/aCv(i,j) ) + gF * (cellVolume(i+1,j)/aCv(i+1,j)) 
		Df = (dxCF*bar_Duf*Sfx + dyCF*bar_Dvf*Sfy ) / (dxCF**2 + dyCF**2) 
		aE(i,j) = -densityField * Df
	end if
	
	aC(i,j) = aN(i,j) - aW(i,j) - aS(i,j) - aE(i,j)
	bC(i,j) = - (md0tH0rP0sDir(i+1,j)+md0tH0rP0sDir(i,j)+md0tVerP0sDir(i,j+1)+md0tVerP0sDir(i,j)) 
	end do
end do

call sweepingTDMAforPressureCorrection()

pCorrection = pCorrectionNew
end subroutine 
