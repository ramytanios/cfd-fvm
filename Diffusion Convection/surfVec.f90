! This subroutine calculates the surface vectors (S,T,E)
! in the East and North directions for each
! face of all control volumes.
! *******************************************************************

subroutine surfaceVectors()
use variablesDeclaration
implicit none

integer i,j
double precision dumb 

! North direction vectors.
do i = 1,m-1
	do j = 1,n
	
	x0fUnitSn(i,j) = ( -Y(i+1,j) + Y(i,j) ) / horFaceArea(i,j)
	y0fUnitSn(i,j) = ( X(i+1,j) - X(i,j) ) / horFaceArea(i,j)
	
	if (j.ne.n) then
		dumb=sqrt( (x0fCellCentroid(i,j)-x0fHorFaceCenter(i,j))**2 + (y0fCellCentroid(i,j)-y0fHorFaceCenter(i,j))**2 )
		x0fUnitEn(i,j) = ( x0fCellCentroid(i,j) - x0fHorFaceCenter(i,j) ) / dumb
		y0fUnitEn(i,j) = ( y0fCellCentroid(i,j) - y0fHorFaceCenter(i,j) ) / dumb
	else
		dumb=sqrt( (x0fCellCentroid(i,j-1)-x0fHorFaceCenter(i,j))**2 + (y0fCellCentroid(i,j-1)-y0fHorFaceCenter(i,j))**2 )
		x0fUnitEn(i,j) = (x0fHorFaceCenter(i,j) - x0fCellCentroid(i,j-1)) / dumb
		y0fUnitEn(i,j) = (y0fHorFaceCenter(i,j) - y0fCellCentroid(i,j-1)) / dumb
	end if
	
	normEn(i,j) = horFaceArea(i,j) / (x0fUnitSn(i,j)*x0fUnitEn(i,j)+y0fUnitSn(i,j)*y0fUnitEn(i,j))
	normTn(i,j) = sqrt (normEn(i,j)**2 - horFaceArea(i,j)**2)
	
	if (normTn(i,j).eq.0) then
		x0fUnitTn(i,j) = 0
		y0fUnitTn(i,j) = 0
	else
		x0fUnitTn(i,j) = ( x0fUnitSn(i,j)*horFaceArea(i,j) - x0fUnitEn(i,j)*normEn(i,j) ) / normTn(i,j)
		y0fUnitTn(i,j) = ( y0fUnitSn(i,j)*horFaceArea(i,j) - y0fUnitEn(i,j)*normEn(i,j) ) / normTn(i,j)
	end if
	
	end do
end do

! East direction vectors.
do i = 1,m
	do j = 1,n-1
	
	x0fUnitSe(i,j) = ( Y(i,j+1) - Y(i,j) ) / verFaceArea(i,j)
	y0fUnitSe(i,j) = ( -X(i,j+1) + X(i,j) ) / verFaceArea(i,j)
	
	if (i.ne.m) then
		dumb = sqrt( (x0fCellCentroid(i,j)-x0fVerFaceCenter(i,j))**2 + (y0fCellCentroid(i,j)-y0fVerFaceCenter(i,j))**2 )
		x0fUnitEe(i,j) = ( x0fCellCentroid(i,j) - x0fVerFaceCenter(i,j) ) / dumb
		y0fUnitEe(i,j) = ( y0fCellCentroid(i,j) - y0fVerFaceCenter(i,j) ) / dumb
	else
		dumb = sqrt( (-x0fCellCentroid(i-1,j)+x0fVerFaceCenter(i,j))**2 + (-y0fCellCentroid(i-1,j)+y0fVerFaceCenter(i,j))**2 )
		x0fUnitEe(i,j) = ( -x0fCellCentroid(i-1,j) + x0fVerFaceCenter(i,j) ) / dumb
		y0fUnitEe(i,j) = ( -y0fCellCentroid(i-1,j) + y0fVerFaceCenter(i,j) ) / dumb
	end if
	
	normEe(i,j) = verFaceArea(i,j) / (x0fUnitSe(i,j)*x0fUnitEe(i,j)+y0fUnitSe(i,j)*y0fUnitEe(i,j))
	normTe(i,j) = sqrt (normEe(i,j)**2 - verFaceArea(i,j)**2)
	
	if (normTe(i,j).eq.0) then
		x0fUnitTe(i,j) = 0
		y0fUnitTe(i,j) = 0
	else
		x0fUnitTe(i,j) = ( x0fUnitSe(i,j)*verFaceArea(i,j) - x0fUnitEe(i,j)*normEe(i,j) ) / normTe(i,j)
		y0fUnitTe(i,j) = ( y0fUnitSe(i,j)*verFaceArea(i,j) - y0fUnitEe(i,j)*normEe(i,j) ) / normTe(i,j)
	end if
	
	end do
end do

end subroutine 
