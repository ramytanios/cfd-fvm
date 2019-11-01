! This module contains the declarations of all variables.
! ******************************************************
module variablesDeclaration
implicit none

! The number of discretization in eta and xi directions. 
integer m,n

! The coordinates of the grid points. 
double precision, allocatable :: X(:,:),Y(:,:)

! The volume of the cells.
double precision, allocatable :: cellVolume(:,:)

! The coordinates of the control volumes. 
double precision, allocatable :: x0fCellCentroid(:,:),y0fCellCentroid(:,:)

! The coordinates of the centers of horizontal and vertical faces.
double precision, allocatable :: x0fHorFaceCenter(:,:),y0fHorFaceCenter(:,:)
double precision, allocatable :: x0fVerFaceCenter(:,:),y0fVerFaceCenter(:,:)

! The geometric factors.
double precision, allocatable :: EastGeomFactor(:,:),WestGeomFactor(:,:)
double precision, allocatable :: NorthGeomFactor(:,:),SouthGeomFactor(:,:)

! The diffusivity at the interfaces.
double precision, allocatable :: horInterfaceDiff(:,:),verInterfaceDiff(:,:)

! The gamma field.
double precision, allocatable :: gammaField(:,:)

! Horizontal and vertical faces areas.
double precision, allocatable :: horFaceArea(:,:),verFaceArea(:,:)

! The norm of the E and T vectors in east and north directions.
double precision, allocatable :: normEe(:,:),normTe(:,:),normEn(:,:),normTn(:,:)

! The coordinates of the unit Surface vector in east and north directions. 
double precision, allocatable :: x0fUnitSe(:,:),y0fUnitSe(:,:)
double precision, allocatable :: x0fUnitSn(:,:),y0fUnitSn(:,:)

! The coordinates of the unit E vector in east and north directions. 
double precision, allocatable :: x0fUnitEe(:,:),y0fUnitEe(:,:)
double precision, allocatable :: x0fUnitEn(:,:),y0fUnitEn(:,:)

! The coordinates of the unit T vector in east and north directions. 
double precision, allocatable :: x0fUnitTe(:,:),y0fUnitTe(:,:)
double precision, allocatable :: x0fUnitTn(:,:),y0fUnitTn(:,:)

! The source term distribution of the grid. 
double precision, allocatable :: sourceField(:,:)

! The density at the horizontal and vertical faces.
double precision, allocatable :: horFaceDensity(:,:),verFaceDensity(:,:)

! The discrete equation coefficients. 
double precision, allocatable :: aE(:,:),aW(:,:),aS(:,:),aN(:,:),aC(:,:),bC(:,:)

! The number of iterations in the TDMA iterative solver. 
integer iteration,uIteration,vIteration

! The phi field at the centroids. 
double precision, allocatable :: phi(:,:),phiNew(:,:)

! Diagonals for the TDMA subroutine. 
double precision, allocatable :: ldVerSweep(:),mdVerSweep(:),udVerSweep(:),sourceVerSweep(:)
double precision, allocatable :: ldHorSweep(:),mdHorSweep(:),udHorSweep(:),sourceHorSweep(:)

! The mass fluxes of phi values on the faces of each control volume. 
double precision, allocatable :: md0tePhie(:,:),md0tnPhin(:,:),md0tsPhis(:,:),md0twPhiw(:,:)
double precision, allocatable :: md0tePhieSmart(:,:),md0tnPhinSmart(:,:),md0tsPhisSmart(:,:),md0twPhiwSmart(:,:)

! The gradient vector of phi at the centroids.
double precision, allocatable :: x0fGradientAtCell(:,:),y0fGradientAtCell(:,:)

! The error on phi between two iterations. 
double precision ::  phiError

! The velocity field at the centroids. 
double precision, allocatable :: x0fVelocityField(:,:),y0fVelocityField(:,:) 

! The velocity values at the horizontal and vertical faces velocities.
double precision, allocatable :: x0fHorVelocity(:,:),y0fHorVelocity(:,:)
double precision, allocatable :: x0fVerVelocity(:,:),y0fVerVelocity(:,:) 

! The mass flow rate values at the horizontal and vertical faces in the positive x and y.
double precision, allocatable :: md0tH0rP0sDir(:,:),md0tVerP0sDir(:,:)

! The value of phi at inlets and outlets.
double precision, allocatable :: phiAtInletNorth(:), phiAtInletWest(:)
double precision, allocatable :: phiAtInletSouth(:), phiAtInletEast(:)

! Outlet, inlet and wall conditions. 
logical,allocatable :: outletConditionNorth(:)
logical,allocatable :: outletConditionEast(:)
logical,allocatable :: outletConditionSouth(:)
logical,allocatable :: outletConditionWest(:)
logical,allocatable :: inletConditionNorth(:)
logical,allocatable :: inletConditionEast(:)
logical,allocatable :: inletConditionSouth(:)
logical,allocatable :: inletConditionWest(:)
logical,allocatable :: wallConditionNorth(:)
logical,allocatable :: wallConditionEast(:)
logical,allocatable :: wallConditionSouth(:)
logical,allocatable :: wallConditionWest(:)

! Dirichlet boundary conditions.
logical,allocatable :: dirichletNorth(:)
logical,allocatable :: dirichletEast(:)
logical,allocatable :: dirichletSouth(:)
logical,allocatable :: dirichletWest(:)
double precision,allocatable :: dirichletValueNorth(:)
double precision,allocatable :: dirichletValueEast(:)
double precision,allocatable :: dirichletValueSouth(:)
double precision,allocatable :: dirichletValueWest(:)
! Von Neumann boundary conditions.
logical,allocatable :: vonNeumannNorth(:)
logical,allocatable :: vonNeumannEast(:)
logical,allocatable :: vonNeumannSouth(:)
logical,allocatable :: vonNeumannWest(:)
double precision,allocatable :: x0fNeumannGradNorth(:)
double precision,allocatable :: x0fNeumannGradEast(:)
double precision,allocatable :: x0fNeumannGradSouth(:) 
double precision,allocatable :: x0fNeumannGradWest(:) 
double precision,allocatable :: y0fNeumannGradNorth(:) 
double precision,allocatable :: y0fNeumannGradEast(:)
double precision,allocatable :: y0fNeumannGradSouth(:)
double precision,allocatable :: y0fNeumannGradWest(:) 

! Robin boundary conditions. 
logical,allocatable :: robinNorth(:)
logical,allocatable :: robinEast(:)
logical,allocatable :: robinSouth(:)
logical,allocatable :: robinWest(:)
double precision,allocatable :: phiInfNorth(:)
double precision,allocatable :: phiInfEast(:)
double precision,allocatable :: phiInfSouth(:)
double precision,allocatable :: phiInfWest(:)
double precision,allocatable :: hInfNorth(:)
double precision,allocatable :: hInfEast(:)
double precision,allocatable :: hInfSouth(:)
double precision,allocatable :: hInfWest(:)

! Inlet and Outlet velocities and pressures conditions.
double precision,allocatable::uAtInletEast(:),vAtInletEast(:),uAtInletWest(:),vAtInletWest(:)
double precision,allocatable::uAtInletNorth(:),vAtInletNorth(:),uAtInletSouth(:),vAtInletSouth(:)
double precision :: pOutlet

! Initial and New Pressure Fields over the grid.
double precision,allocatable :: initialp(:,:),p(:,:) 

! The velocity field at the centroids. 
double precision, allocatable :: u(:,:),uNew(:,:),v(:,:),vNew(:,:)

! The gradient vector of u at the centroids.
double precision, allocatable :: x0fGradu(:,:),y0fGradu(:,:)

! The gradient vector of v at the centroids.
double precision, allocatable :: x0fGradv(:,:),y0fGradv(:,:)

! The gradient vector of p at the centroids.
double precision, allocatable :: x0fGradp(:,:),y0fGradp(:,:)

! The error on u and v between two iterations. 
double precision ::  uError, vError

! The mass fluxes of u on the faces of each control volume. 
double precision, allocatable :: md0teue(:,:),md0tnun(:,:),md0tsus(:,:),md0twuw(:,:)
double precision, allocatable :: md0teueSmart(:,:),md0tnunSmart(:,:),md0tsusSmart(:,:),md0twuwSmart(:,:)

! The mass fluxes of v on the faces of each control volume. 
double precision, allocatable :: md0teve(:,:),md0tnvn(:,:),md0tsvs(:,:),md0twvw(:,:)
double precision, allocatable :: md0teveSmart(:,:),md0tnvnSmart(:,:),md0tsvsSmart(:,:),md0twvwSmart(:,:)

! The wall velocity.
double precision vWall

! The pressure correction values matrix 
double precision, allocatable :: pCorrection(:,:),pCorrectionNew(:,:)

! The velocity correction values matrix 
double precision, allocatable :: uCorrection(:,:),uCorrectionNew(:,:)
double precision, allocatable :: vCorrection(:,:),vCorrectionNew(:,:)

! Initial velocities guesses.
double precision, allocatable :: initialu(:,:),initialv(:,:)

! The discrete X-Momentum equation coefficients. 
double precision, allocatable :: aEu(:,:),aWu(:,:),aSu(:,:),aNu(:,:),aCu(:,:),bCu(:,:)

! The discrete Y-Momentum equation coefficients. 
double precision, allocatable :: aEv(:,:),aWv(:,:),aSv(:,:),aNv(:,:),aCv(:,:),bCv(:,:)

end module 
