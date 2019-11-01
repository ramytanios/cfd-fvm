! This subroutine calculates sum(md0tfufSmart) in the source term
! of the deferred correction approach.
! *****************************************************
subroutine getSmartFlux0nAllFacesForXmomentumEq()
use variablesDeclaration
use userInput
implicit none 

integer i,j 
double precision uC,uU,uD,tildauC,tildaue,tildaun,tildaus,tildauw
double precision ueSmart,unSmart,uwSmart,usSmart
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
	
		uC = u(i,j)
! 		If the flow is entering from the West boundary.	
		if (inletConditionWest(j)) then
			if (i.eq.m-1) then
				tildauC = 0.0D0
				uU = 0.0D0
			else if (i.eq.1) then
				uU = uAtInletWest(j)
				uD = u(i+1,j)
				tildauC = (uC-uU)/(uD-uU + epsilon)
			else
				uU = u(i-1,j)
				uD = u(i+1,j)
				tildauC = (uC-uU)/(uD-uU + epsilon)
			end if
			
			if (tildauC.ge.0 .AND. tildauC.lt.dble(1)/6) then
				tildaue = dble(3)*tildauC
			else if (tildauC.ge.dble(1)/6 .AND. tildauC.lt.dble(7)/10) then
				tildaue = dble(3)/4*tildauC + dble(3)/8
			else if (tildauC.ge.dble(7)/10 .AND. tildauC.lt.1) then
				tildaue = dble(1)/3*tildauC + dble(2)/3
			else
				tildaue = tildauC
			end if 
			ueSmart = tildaue * (uD-uU)  + uU
			md0teueSmart(i,j) = md0te*ueSmart
			if (i.eq.1) then
				md0twuwSmart(i,j) = 0.0D0!
			end if
			if (i.ne.m-1) then
				md0twuwSmart(i+1,j) = -md0teueSmart(i,j)
			end if
		end if

	end do
end do

end subroutine 
