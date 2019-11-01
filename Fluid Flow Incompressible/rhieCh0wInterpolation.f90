! This subroutine generates the mass flow rates on the elements 
! faces using the Rhi-Chow interpolation method. 
! *****************************************************
subroutine rhieCh0wInterpolation()
use variablesDeclaration
use userInput
implicit none 

integer i,j
double precision uf,vf,bar_uf,bar_vf,Dfv,Dfu,dCF,eCFx,eCFy,gC,gF,bar_xGradp_f,bar_yGradp_f

do i = 1,m
	do j=1,n-1
		if (i.eq.1) then
			if (inletConditionWest(j)) then
				uf = uAtInletWest(j)
				vf = vAtInletWest(j)
			else
				uf = vWall 
				vf = vWall
			end if
		else if (i.eq.m) then
			uf = vWall
			vf = vWall
		else
			gF = EastGeomFactor(i,j)
			gC = 1 - gF
			bar_uf = gC * u(i-1,j) + gF * u(i,j)
			bar_vf = gC * v(i-1,j) + gF * v(i,j)
			Dfu =  gC * (cellVolume(i-1,j)/aCu(i-1,j)) + gF * (cellVolume(i,j)/aCu(i,j)) 
			Dfv =  gC * (cellVolume(i-1,j)/aCv(i-1,j)) + gF * (cellVolume(i,j)/aCv(i,j)) 
			dCF = sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i,j))**2 + (y0fCellCentroid(i-1,j)-y0fCellCentroid(i-1,j))**2 )
			bar_xGradp_f = gC * x0fGradp(i-1,j) + gF * x0fGradp(i,j)
			bar_yGradp_f = gC * y0fGradp(i-1,j) + gF * y0fGradp(i,j)
			eCFx = x0fUnitEe(i,j)
			eCFy = y0fUnitEe(i,j)
			
			uf = bar_uf - Dfu* ( (p(i,j)-p(i-1,j))/dCF - (bar_xGradp_f*eCFx + bar_yGradp_f*eCFy ) ) * eCFx 
			vf = bar_vf - Dfv* ( (p(i,j)-p(i-1,j))/dCF - (bar_xGradp_f*eCFx + bar_yGradp_f*eCFy ) ) * eCFy 	
			md0tH0rP0sDir(i,j) = densityField*verFaceArea(i,j)*(uf*x0fUnitSe(i,j)+vf*y0fUnitSe(i,j))
		end if
	end do
end do
	
	
do i = 1,m-1
	do j=1,n
		if (j.eq.1) then 
			uf = vWall
			vf = vWall
		else if (j.eq.n) then
			if (outletConditionNorth(i)) then
				uf = u(i,j-1)
				vf = v(i,j-1)
			else if (wallConditionNorth(i)) then
				uf = vWall
				vf = vWall 
			end if
		else
			gF = NorthGeomFactor(i,j)
			gC = 1 - gF
			bar_uf = gC * u(i,j-1) + gF * u(i,j)
			bar_vf = gC * v(i,j-1) + gF * v(i,j)
			Dfu = gC * (cellVolume(i,j-1)/aCu(i,j-1)) + gF * (cellVolume(i,j)/aCu(i,j)) 
			Dfv = gC * (cellVolume(i,j-1)/aCv(i-1,j)) + gF * (cellVolume(i,j)/aCv(i,j)) 
			dCF = sqrt( (x0fCellCentroid(i,j)-x0fCellCentroid(i,j))**2 + (y0fCellCentroid(i,j-1)-y0fCellCentroid(i,j-1))**2 )
			bar_xGradp_f = gC * x0fGradp(i,j-1) + gF * x0fGradp(i,j)
			bar_yGradp_f = gC * y0fGradp(i,j-1) + gF * y0fGradp(i,j)
			eCFx = x0fUnitEn(i,j)
			eCFy = y0fUnitEn(i,j)
	
			uf = bar_uf - Dfu* ( (p(i,j)-p(i,j-1))/dCF - (bar_xGradp_f*eCFx + bar_yGradp_f*eCFy ) ) * eCFx 
			vf = bar_vf - Dfv* ( (p(i,j)-p(i,j-1))/dCF - (bar_xGradp_f*eCFx + bar_yGradp_f*eCFy ) ) * eCFy 
			md0tVerP0sDir(i,j) = densityField*horFaceArea(i,j)*(uf*x0fUnitSn(i,j)+vf*y0fUnitSn(i,j))
		end if
	end do
end do

end subroutine 
