! This subroutine caluclates the gradient of u at the 
! centroid of each control volume element using the 
! Gauss method.
! *************************************************
subroutine getGradientOfuForXmomentumEq()
use variablesDeclaration
use userInput
implicit none 

integer i,j
double precision ue,uw,us,un,term1,term2,non0rth,dCb

do i = 1,m-1
	do j = 1,n-1
		if (i.eq.1) then
			ue = u(i+1,j)*EastGeomFactor(i,j) + u(i,j)*(1-EastGeomFactor(i,j))
			if (InletConditionWest(j)) then
				uw = uAtInletWest(j)
			else if (wallConditionWest(j)) then
				uw = vWall
			end if
			
			else if (i.eq.m-1) then
				uw = u(i-1,j)*WestGeomFactor(i,j) + u(i,j)*(1-WestGeomFactor(i,j))
				ue = vWall
			else if (i.ne.1 .AND. i.ne.m-1) then
                ue = u(i-1,j)*EastGeomFactor(i,j) + u(i,j)*(1-EastGeomFactor(i,j))
                uw = u(i-1,j)*WestGeomFactor(i,j) + u(i,j)*(1-WestGeomFactor(i,j))
			end if

		if (j.eq.1) then
			un = u(i,j+1)*NorthGeomFactor(i,j) + u(i,j)*(1-NorthGeomFactor(i,j))
			us = vWall
		else if (j.eq.n-1) then
			us = u(i,j-1)*SouthGeomFactor(i,j) + u(i,j)*(1-SouthGeomFactor(i,j)) 
			if (wallConditionNorth(i)) then
				un = vWall
			else if (outletConditionNorth(i)) then
				un = u(i,j)
			end if
		else if (j.ne.1 .AND. j.ne.n-1) then
            us = u(i,j-1)*SouthGeomFactor(i,j) + u(i,j)*(1-SouthGeomFactor(i,j))
            un = u(i,j+1)*NorthGeomFactor(i,j) + u(i,j)*(1-NorthGeomFactor(i,j))
        end if
            x0fGradu(i,j) = ue*verFaceArea(i+1,j)*x0fUnitSe(i+1,j) &
                + uw*verFaceArea(i,j)*(-x0fUnitSe(i,j)) + un*horFaceArea(i,j+1)*x0fUnitSn(i,j+1)+ &
                us*horFaceArea(i,j)*(-x0fUnitSn(i,j));
            y0fGradu(i,j) = ue*verFaceArea(i+1,j)*y0fUnitSe(i+1,j) &
                + uw*verFaceArea(i,j)*(-y0fUnitSe(i,j)) + un*horFaceArea(i,j+1)*y0fUnitSn(i,j+1)+ &
                us*horFaceArea(i,j)*(-y0fUnitSn(i,j));
            x0fGradu(i,j) =  x0fGradu(i,j)/cellVolume(i,j)
            y0fGradu(i,j) =  y0fGradu(i,j)/cellVolume(i,j)  
    end do
end do

end subroutine 
