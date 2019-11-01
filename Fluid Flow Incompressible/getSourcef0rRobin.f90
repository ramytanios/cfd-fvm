! This subroutine modifies the source term for 
! robin boundary conditions
!*********************************************
subroutine getSourcef0rRobin()
use variablesDeclaration
use userInput
implicit none

integer i,j
double precision term1,term2,term3

do i = 1,m-1
	do j = 1,n-1
	
		if (i.eq.1) then
			if (robinWest(j)) then
				term1 = hInfWest(j)*verFaceArea(i,j)
				term2 = verInterfaceDiff(i,j)*normEe(i,j) / sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i,j))**2 & 
													+ (y0fCellCentroid(i,j)-y0fVerFaceCenter(i,j))**2 )
				term3 = verInterfaceDiff(i,j)*normTe(i,j)*(x0fGradientAtCell(i,j)*(-x0fUnitTe(i,j)) &
														+ y0fGradientAtCell(i,j)*(-y0fUnitTe(i,j)) )
                bC(i,j)= bC(i,j) + term1*term3/(term1+term2)
            end if
            
         else if (i.eq.m-1) then
			if (robinEast(j)) then
				term1 = hInfEast(j)*verFaceArea(i+1,j)
				term2 = verInterfaceDiff(i+1,j)*normEe(i+1,j) / sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i+1,j))**2 & 
													+ (y0fCellCentroid(i,j)-y0fVerFaceCenter(i+1,j))**2 )
				term3 = verInterfaceDiff(i+1,j)*normTe(i+1,j)*(x0fGradientAtCell(i,j)*(x0fUnitTe(i+1,j)) &
														+ y0fGradientAtCell(i,j)*(y0fUnitTe(i+1,j)) )
                bC(i,j)= bC(i,j) + term1*term3/(term1+term2)
			end if
		
		else if (j.eq.1) then
			if (robinSouth(i)) then
				term1 = hInfSouth(i)*horFaceArea(i,j)
				term2 = horInterfaceDiff(i,j)*normEn(i,j) / sqrt( (x0fCellCentroid(i,j)-x0fHorFaceCenter(i,j))**2 & 
													+ (y0fCellCentroid(i,j)-y0fHorFaceCenter(i,j))**2 )
				term3 = horInterfaceDiff(i,j)*normTn(i,j)*(x0fGradientAtCell(i,j)*(-x0fUnitTn(i,j)) &
														+ y0fGradientAtCell(i,j)*(-y0fUnitTn(i,j)) )
                bC(i,j)= bC(i,j) + term1*term3/(term1+term2)
			end if
		
		else if (j.eq.n-1) then
			if (robinNorth(i)) then
				term1 = hInfNorth(i)*horFaceArea(i,j+1)
				term2 = horInterfaceDiff(i,j+1)*normEn(i,j+1) / sqrt( (x0fCellCentroid(i,j)-x0fHorFaceCenter(i,j+1))**2 & 
													+ (y0fCellCentroid(i,j)-y0fHorFaceCenter(i,j+1))**2 )
				term3 = horInterfaceDiff(i,j+1)*normTn(i,j+1)*(x0fGradientAtCell(i,j)*(x0fUnitTn(i,j+1)) &
														+ y0fGradientAtCell(i,j)*(y0fUnitTn(i,j+1)) )
                bC(i,j)= bC(i,j) + term1*term3/(term1+term2)
			end if

		end if
	end do
end do
end subroutine 
