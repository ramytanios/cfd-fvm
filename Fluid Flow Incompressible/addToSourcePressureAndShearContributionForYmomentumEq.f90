! This svbroutine adds to the source term the contribution of the 
! pressure gradient and the shear stress in the Y-Component of
! the Navier-Stokes Equation. 
! *****************************************************
subroutine addToSourcePressureAndShearContributionForYmomentumEq()
use variablesDeclaration
use userInput
implicit none

integer i,j
double precision peSe,pwSw,pnSn,psSs,pe,pw,pn,ps
double precision DvDx_w,DvDy_w,DvDx_e,DvDy_e,DvDx_s,DvDy_s,DvDx_n,DvDy_n

do i = 1,m-1
	do j = 1,n-1
	
		! Pressvre contribvtion.
		call computePressureGradient()
		bCv(i,j) = bCv(i,j) - x0fGradp(i,j) * cellVolume(i,j)
		
		! Shear stress contribvtion.
		if (i.eq.1) then 
			DvDx_w = x0fGradv(i,j)
			DvDy_w = y0fGradv(i,j)
		else
			DvDx_w = x0fGradv(i-1,j)*WestGeomFactor(i,j) + x0fGradv(i,j)*(1-WestGeomFactor(i,j))
			DvDy_w = y0fGradv(i-1,j)*WestGeomFactor(i,j) + y0fGradv(i,j)*(1-WestGeomFactor(i,j))
		end if
		
		if (i.eq.m-1) then 
			DvDx_e = x0fGradv(i,j)
			DvDy_e = y0fGradv(i,j)
		else
			DvDx_e = x0fGradv(i-1,j)*EastGeomFactor(i,j) + x0fGradv(i,j)*(1-EastGeomFactor(i,j))
			DvDy_e = y0fGradv(i-1,j)*EastGeomFactor(i,j) + y0fGradv(i,j)*(1-EastGeomFactor(i,j))
		end if
		
		if (j.eq.1) then
			DvDx_s = x0fGradv(i,j)
			DvDy_s = y0fGradv(i,j)
		else
			DvDx_s = x0fGradv(i,j-1)*SouthGeomFactor(i,j) + x0fGradv(i,j)*(1-SouthGeomFactor(i,j))
			DvDy_s = y0fGradv(i,j-1)*SouthGeomFactor(i,j) + y0fGradv(i,j)*(1-SouthGeomFactor(i,j))
		end if
		
		if (j.eq.n-1) then
			DvDx_n = x0fGradv(i,j)
			DvDy_n = y0fGradv(i,j)
		else
			DvDx_n = x0fGradv(i,j+1)*NorthGeomFactor(i,j) + x0fGradv(i,j)*(1-NorthGeomFactor(i,j))
			DvDy_n = y0fGradv(i,j+1)*NorthGeomFactor(i,j) + y0fGradv(i,j)*(1-NorthGeomFactor(i,j))
		end if
		
		bCv(i,j) = bCv(i,j) + viscosity * ( &
						   DvDx_n * horFaceArea(i,j+1) * x0funitSn(i,j+1) + DvDy_n * horFaceArea(i,j+1) * y0funitSn(i,j+1) &
						  + DvDx_s * horFaceArea(i,j)   * (-x0funitSn(i,j))   + DvDy_s * horFaceArea(i,j)   * (-y0funitSn(i,j)) &
						  + DvDx_e * verFaceArea(i+1,j) * x0funitSe(i+1,j) + DvDy_e * verFaceArea(i+1,j) * y0funitSe(i+1,j) &
						  + DvDx_w * verFaceArea(i,j)   * (-x0funitSe(i,j))   + DvDy_w * verFaceArea(i,j)   * (-y0funitSe(i,j)) )
	end do
end do





end subroutine 
