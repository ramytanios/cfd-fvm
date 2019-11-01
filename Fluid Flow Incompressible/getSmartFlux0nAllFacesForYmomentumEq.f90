! This subroutine calculates sum(md0tfufSmart) in the source term
! of the deferred correction approach.
! *****************************************************
subroutine getSmartFlux0nAllFacesForYmomentumEq()
use variablesDeclaration
use userInput
implicit none 

integer i,j 
double precision vC,vU,vD,tildavC,tildave,tildavn,tildavs,tildavw
double precision veSmart,vnSmart,vwSmart,vsSmart
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
	
		vC = v(i,j)
! 		If the flow is entering from the West boundary.	
		if (inletConditionWest(j)) then
			if (i.eq.m-1) then
				tildavC = 0.0D0
				vU = 0.0D0
			else if (i.eq.1) then
				vU = vAtInletWest(j)
				vD = v(i+1,j)
				tildavC = (vC-vU)/(vD-vU + epsilon)
			else
				vU = v(i-1,j)
				vD = v(i+1,j)
				tildavC = (vC-vU)/(vD-vU + epsilon)
			end if
			
			if (tildavC.ge.0 .AND. tildavC.lt.dble(1)/6) then
				tildave = dble(3)*tildavC
			else if (tildavC.ge.dble(1)/6 .AND. tildavC.lt.dble(7)/10) then
				tildave = dble(3)/4*tildavC + dble(3)/8
			else if (tildavC.ge.dble(7)/10 .AND. tildavC.lt.1) then
				tildave = dble(1)/3*tildavC + dble(2)/3
			else
				tildave = tildavC
			end if 
			veSmart = tildave * (vD-vU)  + vU
			md0teveSmart(i,j) = md0te*veSmart
			if (i.eq.1) then
				md0twvwSmart(i,j) = 0.0D0!
			end if
			if (i.ne.m-1) then
				md0twvwSmart(i+1,j) = -md0teveSmart(i,j)
			end if
		end if

	end do
end do

end subroutine 
