! This subroutine adds to the source term the contribution of the 
! pressure gradient and the shear stress in the X-Component of
! the Navier-Stokes Equation. 
! *****************************************************
subroutine addToSourcePressureAndShearContributionForXmomentumEq()
use variablesDeclaration
use userInput
implicit none

integer i,j
double precision peSe,pwSw,pnSn,psSs,pe,pw,pn,ps
double precision DuDx_w,DuDy_w,DuDx_e,DuDy_e,DuDx_s,DuDy_s,DuDx_n,DuDy_n

do i = 1,m-1
	do j = 1,n-1
	
		! Pressure contribution.
		call computePressureGradient()
		bCu(i,j) = bCu(i,j) - x0fGradp(i,j) * cellVolume(i,j)
		
		! Shear stress contribution.
		if (i.eq.1) then 
			DuDx_w = x0fGradu(i,j)
			DuDy_w = y0fGradu(i,j)
		else
			DuDx_w = x0fGradu(i-1,j)*WestGeomFactor(i,j) + x0fGradu(i,j)*(1-WestGeomFactor(i,j))
			DuDy_w = y0fGradu(i-1,j)*WestGeomFactor(i,j) + y0fGradu(i,j)*(1-WestGeomFactor(i,j))
		end if
		
		if (i.eq.m-1) then 
			DuDx_e = x0fGradu(i,j)
			DuDy_e = y0fGradu(i,j)
		else
			DuDx_e = x0fGradu(i-1,j)*EastGeomFactor(i,j) + x0fGradu(i,j)*(1-EastGeomFactor(i,j))
			DuDy_e = y0fGradu(i-1,j)*EastGeomFactor(i,j) + y0fGradu(i,j)*(1-EastGeomFactor(i,j))
		end if
		
		if (j.eq.1) then
			DuDx_s = x0fGradu(i,j)
			DuDy_s = y0fGradu(i,j)
		else
			DuDx_s = x0fGradu(i,j-1)*SouthGeomFactor(i,j) + x0fGradu(i,j)*(1-SouthGeomFactor(i,j))
			DuDy_s = y0fGradu(i,j-1)*SouthGeomFactor(i,j) + y0fGradu(i,j)*(1-SouthGeomFactor(i,j))
		end if
		
		if (j.eq.n-1) then
			DuDx_n = x0fGradu(i,j)
			DuDy_n = y0fGradu(i,j)
		else
			DuDx_n = x0fGradu(i,j+1)*NorthGeomFactor(i,j) + x0fGradu(i,j)*(1-NorthGeomFactor(i,j))
			DuDy_n = y0fGradu(i,j+1)*NorthGeomFactor(i,j) + y0fGradu(i,j)*(1-NorthGeomFactor(i,j))
		end if
		
		bCu(i,j) = bCu(i,j) + viscosity * ( &
						   DuDx_n * horFaceArea(i,j+1) * x0fUnitSn(i,j+1) + DuDy_n * horFaceArea(i,j+1) * y0fUnitSn(i,j+1) &
						  + DuDx_s * horFaceArea(i,j)   * (-x0fUnitSn(i,j))   + DuDy_s * horFaceArea(i,j)   * (-y0fUnitSn(i,j)) &
						  + DuDx_e * verFaceArea(i+1,j) * x0fUnitSe(i+1,j) + DuDy_e * verFaceArea(i+1,j) * y0fUnitSe(i+1,j) &
						  + DuDx_w * verFaceArea(i,j)   * (-x0fUnitSe(i,j))   + DuDy_w * verFaceArea(i,j)   * (-y0fUnitSe(i,j)) )
	end do
end do





end subroutine 
