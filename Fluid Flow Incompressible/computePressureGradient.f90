! This subroutine calculates the pressure gradient 
! at the centroid of each control volume. 
! *****************************************************
subroutine computePressureGradient()
use variablesDeclaration
use userInput
implicit none 

integer i,j
double precision peSe,pwSw,pnSn,psSs,pe,pw,pn,ps

do i = 1,m-1
	do j = 1,n-1
	
	if (i.eq.1) then
		pe = p(i+1,j)*EastGeomFactor(i,j) + p(i,j)*(1-EastGeomFactor(i,j))
		pw = p(i,j)
	else if (i.eq.m-1) then
		pe = p(i,j)
		pw = p(i-1,j)*WestGeomFactor(i,j) + p(i,j)*(1-WestGeomFactor(i,j))
	else if (i.ne.1 .AND. i.ne.m-1) then
        pe = p(i-1,j)*EastGeomFactor(i,j) + p(i,j)*(1-EastGeomFactor(i,j))
        pw = p(i-1,j)*WestGeomFactor(i,j) + p(i,j)*(1-WestGeomFactor(i,j))
	end if
	
	if (j.eq.1) then
		pn = p(i,j+1)*NorthGeomFactor(i,j) + p(i,j)*(1-NorthGeomFactor(i,j))
		ps = p(i,j) 
	else if (j.eq.n-1) then
		if (WallConditionNorth(i)) then
		pn = p(i,j)
		else if (OutletConditionNorth(i)) then
		pn = pOutlet
		end if
		ps = p(i,j-1)*SouthGeomFactor(i,j) + p(i,j)*(1-SouthGeomFactor(i,j))
	else if (j.ne.1 .AND. j.ne.n-1) then
        ps = p(i,j-1)*SouthGeomFactor(i,j) + p(i,j)*(1-SouthGeomFactor(i,j))
        pn = p(i,j+1)*NorthGeomFactor(i,j) + p(i,j)*(1-NorthGeomFactor(i,j))
    end if

	x0fGradp(i,j) = pe*verFaceArea(i+1,j)*x0fUnitSe(i+1,j) + pw*verFaceArea(i,j)*(-x0fUnitSe(i,j)) + & 
					pn*horFaceArea(i,j+1)*x0fUnitSn(i,j+1) + ps*horFaceArea(i,j)*(-x0fUnitSn(i,j))
	y0fGradp(i,j) = pe*verFaceArea(i+1,j)*y0fUnitSe(i+1,j) + pw*verFaceArea(i,j)*(-y0fUnitSe(i,j)) + & 
					pn*horFaceArea(i,j+1)*y0fUnitSn(i,j+1) + ps*horFaceArea(i,j)*(-y0fUnitSn(i,j))

	x0fGradp(i,j) = x0fGradp(i,j)/cellVolume(i,j)
	y0fGradp(i,j) = y0fGradp(i,j)/cellVolume(i,j)
	end do
end do

end subroutine 
