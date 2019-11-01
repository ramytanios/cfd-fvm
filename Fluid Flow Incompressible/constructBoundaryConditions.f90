! This subroutine constructs the boundary conditions.
! *****************************************************
subroutine constructBoundaryConditions()
use variablesDeclaration
implicit none

double precision LdirichletNorth,L0fcvNorth,LrobinWest,L0fcvWest
integer numCvDirichlet,numCvRobin

! <<<<<<<<<<<<<<<< BOUNDARY CONDITIONS >>>>>>>>>>>>>>>>>>>
!

vWall = 0.0D0


!  -------------------------- North boundary conditions.
!

! <<<<<<<<<<Specifications>>>>>>>>>>
LdirichletNorth = 0.8
L0fcvNorth = 1.0D0 / (m-1)
numCvDirichlet = int(LdirichletNorth/L0fcvNorth)

!         - Inlet -
phiAtInletNorth = 0.0D0
inletConditionNorth=.FALSE.
vAtInletNorth = 0
uAtInletNorth = 0
!         - Outlet -
outletConditionNorth(1:numCvDirichlet) =.FALSE.
outletConditionNorth(numCvDirichlet+1:m-1) = .TRUE.
pOutlet = 1.0D-5
!         - Wall -
wallConditionNorth(1:numCvDirichlet)=.TRUE.
wallConditionNorth(numCvDirichlet+1:m-1)=.FALSE.

!         - Dirichlet -
dirichletNorth(1:numCvDirichlet)=.TRUE.
dirichletNorth(numCvDirichlet+1:m-1) = .FALSE.
dirichletValueNorth(1:numCvDirichlet) = 300

!         - Von Neumann -
vonNeumannNorth = .FALSE.
x0fNeumannGradNorth = 0
y0fNeumannGradNorth = 0

!         - Robin -
robinNorth = .FALSE.
phiInfNorth = 300
hInfNorth = 15

!  -------------------------- West boundary conditions. 
!

! <<<<<<<<<<Specifications>>>>>>>>>>
LrobinWest = 0.8
L0fcvWest = 1.0D0 / (n-1)
numCvRobin = int(LrobinWest/L0fcvWest)

!         - Inlet -
vAtInletWest = 0
uAtInletWest = 5
phiAtInletWest(1:n-1-numCvRobin) = 500.0D0 
inletConditionWest(1:n-1-numCvRobin)=.TRUE.
inletConditionWest(n-1-numCvRobin+1:n-1) =.FALSE.
!         - Outlet -
outletConditionWest =.FALSE.

!         - Wall -
wallConditionWest(1:n-1-numCvRobin)=.FALSE.
wallConditionWest(n-1-numCvRobin+1:n-1)=.TRUE.
!         - Dirichlet -
dirichletWest=.FALSE.
dirichletValueWest = 0

!         - Von Neumann -
vonNeumannWest = .FALSE.
x0fNeumannGradWest = 0
y0fNeumannGradWest = 0

!         - Robin -
robinWest(1:n-1-numCvRobin) = .FALSE.
robinWest(n-1-numCvRobin+1:n-1) = .TRUE.
phiInfWest(n-1-numCvRobin+1:n-1) = 300
hInfWest(n-1-numCvRobin+1:n-1) = 15

!  -------------------------- South boundary conditions.
!

! <<<<<<<<<<Specifications>>>>>>>>>>


!         - Inlet -
phiAtInletSouth= 0.0D0
inletConditionSouth=.FALSE.
vAtInletSouth = 0
uAtInletSouth = 0
!         - Outlet -
outletConditionSouth =.FALSE.

!         - wall -
wallConditionSouth=.TRUE.

!         - Dirichlet -
dirichletSouth=.TRUE.
dirichletValueSouth = 400

!         - Von Neumann -
vonNeumannSouth = .FALSE.
x0fNeumannGradSouth = 0
y0fNeumannGradSouth = 0

!         - Robin -
robinSouth = .FALSE.
phiInfSouth = 0
hInfSouth = 0

!  -------------------------- East boundary conditions.
!

! <<<<<<<<<<Specifications>>>>>>>>>>


!         - Inlet -
phiAtInletEast= 0.0D0
inletConditionEast=.FALSE.
vAtInletEast = 0
uAtInletEast = 0
!         - Outlet -
outletConditionEast =.FALSE.

!         - Wall -
wallConditionEast=.TRUE.

!         - Dirichlet -
dirichletEast=.FALSE.
dirichletValueEast = 0

!         - Von Neumann -
vonNeumannEast = .FALSE.
x0fNeumannGradEast = 0
y0fNeumannGradEast = 0

!         - Robin -
robinEast = .FALSE.
phiInfEast = 0
hInfEast = 00


end subroutine 
