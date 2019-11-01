! This subroutine calculates sum(md0tfPhifUpwind) in the source term
! of the deferred correction approach.

! Note that in this source term, the boundary faces fluxes are 
! considered to be 0.
! *****************************************************
subroutine getUpwindFlux0nAllFaces()
use variablesDeclaration
implicit none

integer i,j
double precision md0te,md0tw,md0ts,md0tn

do i = 1,m-1
		do j = 1,n-1
			md0te = md0tVerP0sDir(i+1,j)
			md0tw = -md0tVerP0sDir(i,j)
			md0tn = md0tH0rP0sDir(i,j+1)
			md0ts = -md0tH0rP0sDir(i,j)
				
			if (i.eq.m-1) then
				md0tePhie(i,j) = 0.0D0
			else
				md0tePhie(i,j) = dmax1(md0te,dble(0))*phi(i,j) - &
												dmax1(dble(0),-md0te)*phi(i+1,j)
			end if
				
			if (i.eq.1) then
				md0twPhiw(i,j) = 0.0D0
			end if
			if (i.ne.m-1) then
					md0twPhiw(i+1,j) = -md0tePhie(i,j) 
			end if
!			--------------------------------------
			if (j.eq.n-1) then
				md0tnPhin(i,j) = 0.0D0
			else
				md0tnPhin(i,j) = dmax1(md0tn,dble(0))*phi(i,j) - &
												dmax1(dble(0),-md0tn)*phi(i,j+1)
			end if
				
			if (j.eq.1) then
				md0tsPhis(i,j) = 0.0D0
			end if
			if (j.ne.n-1) then
					md0tsPhis(i,j+1) = -md0tnPhin(i,j)
			end if
		end do
	end do

end subroutine 
