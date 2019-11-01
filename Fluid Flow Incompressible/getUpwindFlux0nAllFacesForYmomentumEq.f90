! This subroutine calculates sum(md0tfufUpwind) in the source term
! of the deferred correction approach.

! Note that in this source term, the boundary faces fluxes are 
! considered to be 0.
! *****************************************************
subroutine getUpwindFlux0nAllFacesForYmomentumEq()
use variablesDeclaration
implicit none

integer i,j
double precision md0te,md0tw,md0ts,md0tn

do i = 1,m-1
		do j = 1,n-1
		md0te =  verFaceDensity(i+1,j)*verFaceArea(i+1,j) * & 
					(u(i+1,j)*x0fUnitSe(i+1,j)+v(i+1,j)*y0fUnitSe(i+1,j))
		md0tw = - verFaceDensity(i,j)*verFaceArea(i,j) * & 
					(u(i,j)*x0fUnitSe(i,j)+v(i,j)*y0fUnitSe(i,j))
		md0tn = horFaceDensity(i,j+1)*horFaceArea(i,j+1) * & 
					(u(i,j+1)*x0fUnitSn(i,j+1)+v(i,j+1)*y0fUnitSn(i,j+1))
		md0ts = -horFaceDensity(i,j)*horFaceArea(i,j) * & 
					(u(i,j)*x0fUnitSn(i,j)+v(i,j)*y0fUnitSn(i,j))
				
			if (i.eq.m-1) then
				md0teue(i,j) = 0.0D0
			else
				md0teue(i,j) = dmax1(md0te,dble(0))*u(i,j) - &
												dmax1(dble(0),-md0te)*u(i+1,j)
			end if
				
			if (i.eq.1) then
				md0twuw(i,j) = 0.0D0
			end if 
			if (i.ne.m-1) then
					md0twuw(i+1,j) = -md0teue(i,j) 
			end if
!			--------------------------------------
			if (j.eq.n-1) then
				md0tnun(i,j) = 0.0D0
			else
				md0tnun(i,j) = dmax1(md0tn,dble(0))*u(i,j) - &
												dmax1(dble(0),-md0tn)*u(i,j+1)
			end if
				
			if (j.eq.1) then
				md0tsus(i,j) = 0.0D0
			end if
			if (j.ne.n-1) then
					md0tsus(i,j+1) = -md0tnun(i,j)
			end if
		end do
	end do

end subroutine 
