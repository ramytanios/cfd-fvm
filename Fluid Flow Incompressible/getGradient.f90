! This subroutine caluclates the gradient of phi at the 
! centroid of each control volume element using the 
! Gauss method.
! *************************************************
subroutine getGradient()
use variablesDeclaration
use userInput
implicit none 

integer i,j
double precision phie,phiw,phis,phin,term1,term2,non0rth,dCb

do i = 1,m-1
	do j = 1,n-1
		if (i.eq.1) then
			phie = phi(i+1,j)*EastGeomFactor(i,j) + phi(i,j)*(1-EastGeomFactor(i,j))
			if (dirichletWest(j)) then
				phiw = dirichletValueWest(2)
			else if (robinWest(j)) then
				term1 = hInfWest(j)*verFaceArea(i,j)
				dCb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i,j))**2 &
												+(y0fCellCentroid(i,j)-y0fVerFaceCenter(i,j))**2 )
				term2 = verInterfaceDiff(i,j) * normEe(i,j)/dCb
				non0rth = verInterfaceDiff(i,j) * ( x0fGradientAtCell(i,j)*normTe(i,j)*x0fUnitTe(i,j) &
												+ y0fGradientAtCell(i,j)*normTe(i,j)*y0fUnitTe(i,j) )
				phiw = ( term1*phiInfWest(j) + term2*phi(i,j) - non0rth ) / (term1 + term2)
			else
				phiw = phi(i,j)
			end if
			
			else if (i.eq.m-1) then
				phiw = phi(i-1,j)*WestGeomFactor(i,j) + phi(i,j)*(1-WestGeomFactor(i,j))
				if (dirichletEast(j)) then
					phie = dirichletValueEast(j)
				else if (robinEast(j)) then
					term1 = hInfEast(j)*verFaceArea(i,j)
					dCb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i+1,j))**2 &
												+(y0fCellCentroid(i,j)-y0fVerFaceCenter(i+1,j))**2 )
					term2 = verInterfaceDiff(i+1,j) * normEe(i+1,j) / dCb
					non0rth = verInterfaceDiff(i+1,j) * ( x0fGradientAtCell(i,j)*normTe(i+1,j)*x0fUnitTe(i+1,j) &
												+ y0fGradientAtCell(i,j)*normTe(i+1,j)*y0fUnitTe(i+1,j) )
					phie = ( term1*phiInfEast(j) + term2*phi(i,j) - non0rth ) / (term1 + term2)
                else
					phie = phi(i,j)
				end if
				
			else if (i.ne.1 .AND. i.ne.m-1) then
                phie = phi(i-1,j)*EastGeomFactor(i,j) + phi(i,j)*(1-EastGeomFactor(i,j))
                phiw = phi(i-1,j)*WestGeomFactor(i,j) + phi(i,j)*(1-WestGeomFactor(i,j))
			end if

		if (j.eq.1) then
			phin = phi(i,j+1)*NorthGeomFactor(i,j) + phi(i,j)*(1-NorthGeomFactor(i,j))
			if (dirichletSouth(i)) then
				phis = dirichletValueSouth(i)
			else if (robinSouth(i)) then
				term1 = hInfSouth(i)*horFaceArea(i,j)
				dCb = sqrt( (x0fCellCentroid(i,j)-x0fHorFaceCenter(i,j))**2 &
												+(y0fCellCentroid(i,j)-y0fHorFaceCenter(i,j))**2 )
				term2 = horInterfaceDiff(i,j) * normEn(i,j) / dCb
				non0rth = horInterfaceDiff(i+1,j) * ( x0fGradientAtCell(i,j)*normTn(i,j)*x0fUnitTn(i,j) &
												+ y0fGradientAtCell(i,j)*normTn(i,j)*y0fUnitTn(i,j) )
				phin = ( term1*phiInfSouth(i) + term2*phi(i,j) - non0rth ) / (term1 + term2)
			else
				phin = phi(i,j)
			end if
			
		else if (j.ne.1 .AND. j.ne.n-1) then
            phis = phi(i,j-1)*SouthGeomFactor(i,j) + phi(i,j)*(1-SouthGeomFactor(i,j))
            phin = phi(i,j+1)*NorthGeomFactor(i,j) + phi(i,j)*(1-NorthGeomFactor(i,j))
        end if
            x0fGradientAtCell(i,j) = phie*verFaceArea(i+1,j)*x0fUnitSe(i+1,j) &
                + phiw*verFaceArea(i,j)*(-x0fUnitSe(i,j)) + phin*horFaceArea(i,j+1)*x0fUnitSn(i,j+1)+ &
                phis*horFaceArea(i,j)*(-x0fUnitSn(i,j));
            y0fGradientAtCell(i,j) = phie*verFaceArea(i+1,j)*y0fUnitSe(i+1,j) &
                + phiw*verFaceArea(i,j)*(-y0fUnitSe(i,j)) + phin*horFaceArea(i,j+1)*y0fUnitSn(i,j+1)+ &
                phis*horFaceArea(i,j)*(-y0fUnitSn(i,j));
            x0fGradientAtCell(i,j) =  x0fGradientAtCell(i,j)/cellVolume(i,j)
            y0fGradientAtCell(i,j) =  y0fGradientAtCell(i,j)/cellVolume(i,j)  
    end do
end do

end subroutine 
