! This subroutine calculates sum(md0tfPhifSmart) in the source term
! of the deferred correction approach.
! *****************************************************
subroutine getSmartFlux0nAllFaces()
use variablesDeclaration
use userInput
implicit none 

integer i,j 
double precision phiC,phiU,phiD,tildaPhiC,tildaPhie,tildaPhin,tildaPhis,tildaPhiw
double precision phieSmart,phinSmart,phiwSmart,phisSmart
double precision md0te,md0tw,md0ts,md0tn
do i = 1,m-1
	do j = 1,n-1
	
		md0te = md0tVerP0sDir(i+1,j)
        md0tw = -md0tVerP0sDir(i,j)
        md0tn = md0tH0rP0sDir(i,j+1)
        md0ts = -md0tH0rP0sDir(i,j)
	
		phiC = phi(i,j)
! 		If the flow is entering from the West boundary.	
		if (inletConditionWest(j)) then
			if (i.eq.m-1) then
				tildaPhiC = 0.0D0
				phiU = 0.0D0
			else if (i.eq.1) then
				phiU = phiAtInletWest(j)
				phiD = phi(i+1,j)
				tildaPhiC = (phiC-phiU)/(phiD-phiU + epsilon)
			else
				phiU = phi(i-1,j)
				phiD = phi(i+1,j)
				tildaPhiC = (phiC-phiU)/(phiD-phiU + epsilon)
			end if
			
			if (tildaPhiC.ge.0 .AND. tildaPhiC.lt.dble(1)/6) then
				tildaPhie = dble(3)*tildaPhiC
			else if (tildaPhiC.ge.dble(1)/6 .AND. tildaPhiC.lt.dble(7)/10) then
				tildaPhie = dble(3)/4*tildaPhiC + dble(3)/8
			else if (tildaPhiC.ge.dble(7)/10 .AND. tildaPhiC.lt.1) then
				tildaPhie = dble(1)/3*tildaPhiC + dble(2)/3
			else
				tildaPhie = tildaPhiC
			end if 
			phieSmart = tildaPhie * (phiD-phiU)  + phiU
			md0tePhieSmart(i,j) = md0te*phieSmart
			if (i.eq.1) then
				md0twPhiwSmart(i,j) = 0.0D0!
			end if
			if (i.ne.m-1) then
				md0twPhiwSmart(i+1,j) = -md0tePhieSmart(i,j)
			end if
		
		end if
		
! 		If the flow is entering from the East boundary.		
		if (inletConditionEast(j)) then
			if (i.eq.1) then
				tildaPhiC = 0.0D0
				phiU = 0.0D0
			else if (i.eq.m-1) then
				phiU = phiAtInletEast(j)
				phiD = phi(i-1,j)
				tildaPhiC = (phiC-phiU)/(phiD-phiU + epsilon)
			else
				phiU = phi(i+1,j)
				phiD = phi(i-1,j)
				tildaPhiC = (phiC-phiU)/(phiD-phiU + epsilon)
			end if
			if (tildaPhiC.ge.0 .AND. tildaPhiC.lt.dble(1)/6) then
				tildaPhiw = dble(3)*tildaPhiC
			else if (tildaPhiC.ge.dble(1)/6 .AND. tildaPhiC.lt.dble(7)/10) then
				tildaPhiw = dble(3)/4*tildaPhiC + dble(3)/8
			else if (tildaPhiC.ge.dble(7)/10 .AND. tildaPhiC.lt.1) then
				tildaPhiw = dble(1)/3*tildaPhiC + dble(2)/3
			else
				tildaPhiw = tildaPhiC
			end if 
			phiwSmart = tildaPhiw * (phiD-phiU) + phiU
			md0twPhiwSmart(i,j) = md0tw*phiwSmart
			if (i.eq.m-1) then
				md0tePhieSmart(i,j) = 0.0D0
			end if
			if (i.ne.1) then
				md0tePhieSmart(i-1,j) = -md0twPhiwSmart(i,j)
			end if
		end if
		
!  		----------------------------------------------------------------------

! 		If the flow is entering from the South boundary.
		if (inletConditionSouth(i)) then
			if (j.eq.n-1) then
				tildaPhiC = 0.0D0
				phiU = 0.0D0
			else if (j.eq.1) then
				phiU = phiAtInletSouth(i)
				phiD = phi(i,j+1)
				tildaPhiC = (phiC-phiU)/(phiD-phiU + epsilon)
			else
				phiU=phi(i,j-1)
				phiD=phi(i,j+1)
				tildaPhiC = (phiC-phiU)/(phiD-phiU + epsilon)
			end if
			
			if (tildaPhiC.ge.0 .AND. tildaPhiC.lt.dble(1)/6) then
				tildaPhin = dble(3)*tildaPhiC
			else if (tildaPhiC.ge.dble(1)/6 .AND. tildaPhiC.lt.dble(7)/10) then
				tildaPhin = dble(3)/4*tildaPhiC + dble(3)/8
			else if (tildaPhiC.ge.dble(7)/10 .AND. tildaPhiC.lt.1) then
				tildaPhin = dble(1)/3*tildaPhiC + dble(2)/3
			else
				tildaPhin = tildaPhiC
			end if 
			phinSmart = tildaPhin * (phiD-phiU) + phiU
			md0tnPhinSmart(i,j) = md0tn*phinSmart
			if (j.eq.1) then
				md0tsPhisSmart(i,j) = 0.0D0 !
			end if
			if (j.ne.n-1) then
				md0tsPhisSmart(i,j+1) = -md0tnPhinSmart(i,j)
			end if
		end if
		
		
! 		If the flow is entering from the North boundary.		
		if (inletConditionNorth(i)) then
			if (j.eq.1) then
				tildaPhiC = 0.0D0
				phiU = 0.0D0
			else if (j.eq.n-1) then
				phiU = phiAtInletNorth(i)
				phiD = phi(i,j-1)
				tildaPhiC = (phiC-phiU)/(phiD-phiU + epsilon)
			else
				phiU=phi(i,j+1)
				phiD=phi(i,j-1)
				tildaPhiC = (phiC-phiU)/(phiD-phiU + epsilon)
			end if
			
			if (tildaPhiC.ge.0 .AND. tildaPhiC.lt.dble(1)/6) then
				tildaPhis = dble(3)*tildaPhiC
			else if (tildaPhiC.ge.dble(1)/6 .AND. tildaPhiC.lt.dble(7)/10) then
				tildaPhis = dble(3)/4*tildaPhiC + dble(3)/8
			else if (tildaPhiC.ge.dble(7)/10 .AND. tildaPhiC.lt.1) then
				tildaPhis = dble(1)/3*tildaPhiC + dble(2)/3
			else
				tildaPhis = tildaPhiC
			end if 
			phisSmart = tildaPhis * (phiD-phiU) + phiU
			md0tsPhisSmart(i,j) = md0ts*phisSmart
			if (j.eq.n-1) then
				md0tnPhinSmart(i,j) = 0.0D0 !
			end if
			if (j.ne.1) then
				md0tnPhinSmart(i,j-1) = -md0tsPhisSmart(i,j)
			end if	
		end if		
	end do
end do

end subroutine 
