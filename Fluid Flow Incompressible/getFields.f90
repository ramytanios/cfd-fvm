! This subroutine generates the density, velocity values 
! at the faces of the control volumes. It also generates the 
! source term and the diffusivity 
! at the centroid of each control volume and the mass
! flow rates for each control volume.

! --> The velocity field is assumed to be known at any point
! in the domain.
! ***********************************************************
subroutine getFields()
use variablesDeclaration
use userInput 
implicit none

integer i,j

do i = 1,m-1
	do j = 1,n-1
		sourceField(i,j) = sourceFunction(x0fCellCentroid(i,j),y0fCellCentroid(i,j))
		gammaField(i,j) = gammaFunction(x0fCellCentroid(i,j),y0fCellCentroid(i,j))
	end do
end do

! ------ Velocities and densities at faces. ------
do i = 1,m-1
	do j = 1,n
		if (j.eq.1) then
			horInterfaceDiff(i,j) = gammaField(i,j)
		else if (j.eq.n) then
			horInterfaceDiff(i,j) = gammaField(i,j-1)
		else 
			horInterfaceDiff(i,j) = 1/ ( (1-NorthGeomFactor(i,j-1))/gammaField(i,j) + NorthGeomFactor(i,j-1)/gammaField(i,j-1) )
		end if
		horFaceDensity(i,j) = densityField
		!md0tH0rP0sDir(i,j) = horFaceDensity(i,j)*horFaceArea(i,j) * & 
					!(x0fHorVelocity(i,j)*x0fUnitSn(i,j)+y0fHorVelocity(i,j)*y0fUnitSn(i,j))
		
	end do
end do

do i = 1,m
	do j = 1,n-1
		if (i.eq.1) then
			verInterfaceDiff(i,j) = gammaField(i,j)
		else if (i.eq.m) then
			verInterfaceDiff(i,j) = gammaField(i-1,j)
		else 
			verInterfaceDiff(i,j) = 1/ ( (1-EastGeomFactor(i-1,j))/gammaField(i,j) + EastGeomFactor(i-1,j)/gammaField(i-1,j) )
		end if
		verFaceDensity(i,j) = densityField
		!md0tVerP0sDir(i,j) = verFaceDensity(i,j)*verFaceArea(i,j) * & 
					!(x0fVerVelocity(i,j)*x0fUnitSe(i,j)+y0fVerVelocity(i,j)*y0fUnitSe(i,j))
	end do
end do



end subroutine 
