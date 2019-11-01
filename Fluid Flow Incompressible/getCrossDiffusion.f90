! This subroutine adds the cross diffusion terms to
! the source term.
! *************************************************
subroutine getCrossDiffusion()
use variablesDeclaration
use userInput
implicit none 

integer i,j
double precision x0fgradE,x0fgradW,x0fgradS,x0fgradN
double precision y0fgradE,y0fgradW,y0fgradS,y0fgradN
double precision crossDiffTermWest,crossDiffTermSouth,crossDiffTermEast,crossDiffTermNorth,crossDiffTerm

do i = 1,m-1
	do j = 1,n-1
	
	if (i.eq.1) then
		x0fgradE = x0fGradientAtCell(i+1,j)*EastGeomFactor(i,j) + x0fGradientAtCell(i,j)*(1-EastGeomFactor(i,j))
		y0fgradE = y0fGradientAtCell(i+1,j)*EastGeomFactor(i,j) + y0fGradientAtCell(i,j)*(1-EastGeomFactor(i,j))
		if (vonNeumannWest(j) .OR. robinWest(j)) then
			x0fgradW = 0.0D0
			y0fgradW = 0.0D0
		else
			x0fgradW = x0fGradientAtCell(i,j)
			y0fgradW = y0fGradientAtCell(i,j)
		end if
	else if (i.eq.m-1) then
		x0fgradW = x0fGradientAtCell(i-1,j)*WestGeomFactor(i,j) + x0fGradientAtCell(i,j)*(1-WestGeomFactor(i,j))
		y0fgradW = y0fGradientAtCell(i-1,j)*WestGeomFactor(i,j) + y0fGradientAtCell(i,j)*(1-WestGeomFactor(i,j))
		if (vonNeumannEast(j) .OR. robinEast(j)) then
			x0fgradE = 0.0D0
			y0fgradE = 0.0D0
		else
			x0fgradE = x0fGradientAtCell(i,j)
			y0fgradE = y0fGradientAtCell(i,j)
		end if
	else if (i.ne.1 .AND. i.ne.m-1) then
		x0fgradE = x0fGradientAtCell(i+1,j)*EastGeomFactor(i,j) + x0fGradientAtCell(i,j)*(1-EastGeomFactor(i,j))
		y0fgradE = y0fGradientAtCell(i+1,j)*EastGeomFactor(i,j) + y0fGradientAtCell(i,j)*(1-EastGeomFactor(i,j))
		x0fgradW = x0fGradientAtCell(i-1,j)*WestGeomFactor(i,j) + x0fGradientAtCell(i,j)*(1-WestGeomFactor(i,j))
		y0fgradW = y0fGradientAtCell(i-1,j)*WestGeomFactor(i,j) + y0fGradientAtCell(i,j)*(1-WestGeomFactor(i,j))
	end if
	
	if (j.eq.1) then
		x0fgradN = x0fGradientAtCell(i,j+1)*NorthGeomFactor(i,j) + x0fGradientAtCell(i,j)*(1-NorthGeomFactor(i,j))
		y0fgradN = y0fGradientAtCell(i,j+1)*NorthGeomFactor(i,j) + y0fGradientAtCell(i,j)*(1-NorthGeomFactor(i,j))
		if (vonNeumannSouth(i) .OR. robinSouth(i)) then
			x0fgradS = 0.0D0
			y0fgradS = 0.0D0
		else
			x0fgradS = x0fGradientAtCell(i,j)
			y0fgradS = y0fGradientAtCell(i,j)
		end if
	else if (j.eq.n-1) then
		x0fgradS = x0fGradientAtCell(i,j-1)*NorthGeomFactor(i,j) + x0fGradientAtCell(i,j)*(1-NorthGeomFactor(i,j))
		y0fgradS = y0fGradientAtCell(i,j-1)*NorthGeomFactor(i,j) + y0fGradientAtCell(i,j)*(1-NorthGeomFactor(i,j))
		if (vonNeumannNorth(i) .OR. robinNorth(i)) then
			x0fgradN = 0.0D0
			y0fgradN = 0.0D0
		else
			x0fgradN = x0fGradientAtCell(i,j)
			y0fgradN = y0fGradientAtCell(i,j)
		end if
	else if (j.ne.1 .AND. j.ne.n-1) then
		x0fgradN = x0fGradientAtCell(i,j+1)*NorthGeomFactor(i,j) + x0fGradientAtCell(i,j)*(1-NorthGeomFactor(i,j))
		y0fgradN = y0fGradientAtCell(i,j+1)*NorthGeomFactor(i,j) + y0fGradientAtCell(i,j)*(1-NorthGeomFactor(i,j))
		x0fgradS = x0fGradientAtCell(i,j-1)*NorthGeomFactor(i,j) + x0fGradientAtCell(i,j)*(1-NorthGeomFactor(i,j))
		y0fgradS = y0fGradientAtCell(i,j-1)*NorthGeomFactor(i,j) + y0fGradientAtCell(i,j)*(1-NorthGeomFactor(i,j))
	end if
	
	crossDiffTermNorth = horInterfaceDiff(i,j+1) * normTn(i,j+1) * (x0fgradN*x0fUnitTn(i,j+1)+y0fgradN*y0fUnitTn(i,j+1))
	crossDiffTermWest = verInterfaceDiff(i,j) * normTe(i,j) * (x0fgradW*(-x0fUnitTe(i,j))+y0fgradW*(-y0fUnitTe(i,j)))
	crossDiffTermSouth = horInterfaceDiff(i,j) * normTn(i,j) * (x0fgradS*(-x0fUnitTn(i,j))+y0fgradN*(-y0fUnitTn(i,j)))
	crossDiffTermEast = verInterfaceDiff(i,j) * normTe(i+1,j) * (x0fgradE*x0fUnitTe(i+1,j)+y0fgradE*y0fUnitTe(i+1,j))

	crossDiffTerm = crossDiffTermNorth+crossDiffTermWest+crossDiffTermSouth+crossDiffTermEast;
		
	bC(i,j) = bC(i,j) + crossDiffTerm
		
	end do
end do



end subroutine 
