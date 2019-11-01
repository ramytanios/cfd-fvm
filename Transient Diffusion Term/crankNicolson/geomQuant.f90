! This subroutine calculates the geometric quantities of the 2D mesh:
! *******************************************************************
! Cell volume.
! Cell centroids coordinates.
! Horizontal and vertical faces centers coordinates.
! Horizonta and vertical faces areas.
! Geometric factors.
! *******************************************************************

subroutine geometricQuantities()
use variablesDeclaration
implicit none 

integer i,j

do i = 1,m-1
	do j = 1,n-1
		! The volume of a cell is approximated by the area of a parallelogram.
		cellVolume(i,j) = abs( (X(i+1,j)-X(i,j) )*( Y(i,j+1)-Y(i,j) ) &
                - ( X(i,j+1)-X(i,j) )*( Y(i+1,j)-Y(i,j) ) )
        ! CellCentroid coordinates.
        x0fCellCentroid(i,j) = dble(1)/4 * ( X(i,j)+X(i,j+1)+X(i+1,j)+X(i+1,j+1) )
        y0fCellCentroid(i,j) = dble(1)/4 * ( Y(i,j)+Y(i,j+1)+Y(i+1,j)+Y(i+1,j+1) )
	end do
end do

! Faces Areas and Faces Centers.
do i = 1,m
	do j = 1,n-1 
        x0fVerFaceCenter(i,j) =  dble(1)/2 * ( X(i,j) + X(i,j+1) ) 
        y0fVerFaceCenter(i,j) =  dble(1)/2 * ( Y(i,j) + Y(i,j+1) ) 

        verFaceArea(i,j) = sqrt( (X(i,j+1)-X(i,j))**2 + (Y(i,j+1)-Y(i,j))**2 )
	end do
end do

do i = 1,m-1
	do j = 1,n
	    x0fHorFaceCenter(i,j) =  dble(1)/2 * ( X(i,j) + X(i+1,j) ) 
        y0fHorFaceCenter(i,j) =  dble(1)/2 * ( Y(i,j) + Y(i+1,j) )
        
        horFaceArea(i,j) = sqrt( (X(i+1,j)-X(i,j))**2 + (Y(i+1,j)-Y(i,j))**2 )
    end do
end do

! Geometric factors calculations.
do i = 1,m-1
	do j = 1,n-1
		if (i.eq.m-1) then
			EastGeomFactor(i,j) = 0.0D0
		else
			EastGeomFactor(i,j) = cellVolume(i,j) / ( cellVolume(i,j) + cellVolume(i+1,j) ) 
		end if
		
		if (i.eq.1) then
			WestGeomFactor(i,j) = 0.0D0
		else
			WestGeomFactor(i,j) = cellVolume(i,j) / ( cellVolume(i,j) + cellVolume(i-1,j) )
		end if
		
		if (j.eq.n-1) then
			NorthGeomFactor(i,j) = 0.0D0
		else
			NorthGeomFactor(i,j) = cellVolume(i,j) / ( cellVolume(i,j+1) + cellVolume(i,j) )
		end if
		
		if (j.eq.1) then
			SouthGeomFactor(i,j) = 0.0D0
		else
			SouthGeomFactor(i,j) = cellVolume(i,j-1) / ( cellVolume(i,j) + cellVolume(i,j-1) )
		end if
	end do
end do

end subroutine 
