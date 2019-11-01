! This subroutine generates the analytical solution
! of the PDE on the 2D grid.
! ************************************************
subroutine getAnalyticalSolution()
use variablesDeclaration
implicit none

integer i,j
integer mp,np
double precision, parameter :: epsilon=10**(-8),pi=dacos(-1.0D0),kappa=1,H=1.0D0,L=1.0D0
double precision sinTerm0ne,sinTermTw0,expTerm,constantTerm,absDifferenceInner,absDifference0uter
double precision sommeNew,sommeInner,somme0uter

do i = 1,m-1
	do j = 1,n-1
		mp = 1
		np = 1
		absDifference0uter = 100.0D0
		somme0uter = 0.0D0
		do while (absDifference0uter.gt.epsilon) 
			absDifferenceInner = 100.0D0
			sommeInner = 0.0D0
			do while (absDifferenceInner.gt.epsilon) 
				constantTerm = 1600.0D0/((2.0D0*np-1)*(2.0D0*mp-1)*pi**2)
				sinTerm0ne = sin( (2.0D0*np-1.0D0)*pi*x0fCellCentroid(i,j)/L )
				sinTermTw0 = sin( (2.0D0*mp-1.0D0)*pi*y0fCellCentroid(i,j)/2.0D0/H )
				expTerm = exp(-kappa*pi**2*t*( (2.0D0*mp-1)**2/4.0D0/H**2 + (2.0D0*np-1)**2/L**2 ) )
				sommeNew = sommeInner + constantTerm*sinTerm0ne*sinTermTw0*expTerm
				absDifferenceInner = abs(sommeNew - sommeInner)
				sommeInner = sommeNew
				np = np + 1
			end do
			sommeNew = somme0uter + sommeInner
			absDifference0uter = abs(sommeNew - somme0uter)
			somme0uter = sommeNew 
			mp = mp + 1
		end do
		phiAnalytical(i,j) = somme0uter
	end do
end do


end subroutine 
