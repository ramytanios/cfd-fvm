! This subroutine caluclates the gradient of u at the 
! centroid of each control volume element using the 
! Gauss method.
! *************************************************
subroutine getGradientOfvForYmomentumEq()
use variablesDeclaration
use userInput
implicit none 

integer i,j
double precision ve,vw,vs,vn,term1,term2,non0rth,dCb

do i = 1,m-1
	do j = 1,n-1
		if (i.eq.1) then
			ve = v(i+1,j)*EastGeomFactor(i,j) + v(i,j)*(1-EastGeomFactor(i,j))
			if (InletConditionWest(j)) then
				vw = vAtInletWest(j)
			else if (wallConditionWest(j)) then
				vw = vWall
			end if
			
			else if (i.eq.m-1) then
				vw = v(i-1,j)*WestGeomFactor(i,j) + v(i,j)*(1-WestGeomFactor(i,j))
				ve = vWall
			else if (i.ne.1 .AND. i.ne.m-1) then
                ve = v(i-1,j)*EastGeomFactor(i,j) + v(i,j)*(1-EastGeomFactor(i,j))
                vw = v(i-1,j)*WestGeomFactor(i,j) + v(i,j)*(1-WestGeomFactor(i,j))
			end if

		if (j.eq.1) then
			vn = v(i,j+1)*NorthGeomFactor(i,j) + v(i,j)*(1-NorthGeomFactor(i,j))
			vs = vWall
		else if (j.eq.n-1) then
			vs = v(i,j-1)*SouthGeomFactor(i,j) + v(i,j)*(1-SouthGeomFactor(i,j)) 
			if (wallConditionNorth(i)) then
				vn = vWall
			else if (outletConditionNorth(i)) then
				vn = v(i,j)
			end if
		else if (j.ne.1 .AND. j.ne.n-1) then
            vs = v(i,j-1)*SouthGeomFactor(i,j) + v(i,j)*(1-SouthGeomFactor(i,j))
            vn = v(i,j+1)*NorthGeomFactor(i,j) + v(i,j)*(1-NorthGeomFactor(i,j))
        end if
            x0fGradv(i,j) = ve*verFaceArea(i+1,j)*x0fUnitSe(i+1,j) &
                + vw*verFaceArea(i,j)*(-x0fUnitSe(i,j)) + vn*horFaceArea(i,j+1)*x0fUnitSn(i,j+1)+ &
                vs*horFaceArea(i,j)*(-x0fUnitSn(i,j));
            y0fGradv(i,j) = ve*verFaceArea(i+1,j)*y0fUnitSe(i+1,j) &
                + vw*verFaceArea(i,j)*(-y0fUnitSe(i,j)) + vn*horFaceArea(i,j+1)*y0fUnitSn(i,j+1)+ &
                vs*horFaceArea(i,j)*(-y0fUnitSn(i,j));
            x0fGradv(i,j) =  x0fGradu(i,j)/cellVolume(i,j)
            y0fGradv(i,j) =  y0fGradu(i,j)/cellVolume(i,j)  
    end do
end do

end subroutine 
