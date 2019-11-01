! This module contains the declarations of all variables.
! ******************************************************
module variablesDeclaration
implicit none

integer m,n
double precision, allocatable :: X(:,:),Y(:,:)
double precision, allocatable :: cellVolume(:,:)
double precision, allocatable :: x0fCellCentroid(:,:),y0fCellCentroid(:,:)
double precision, allocatable :: x0fHorFaceCenter(:,:),y0fHorFaceCenter(:,:)
double precision, allocatable :: x0fVerFaceCenter(:,:),y0fVerFaceCenter(:,:)
double precision, allocatable :: EastGeomFactor(:,:),WestGeomFactor(:,:)
double precision, allocatable :: NorthGeomFactor(:,:),SouthGeomFactor(:,:)
double precision, allocatable :: horInterfaceDiff(:,:),verInterfaceDiff(:,:),gammaField(:,:)
double precision, allocatable :: horFaceArea(:,:),verFaceArea(:,:)
double precision, allocatable :: normEe(:,:),normTe(:,:),normEn(:,:),normTn(:,:)
double precision, allocatable :: x0fUnitSe(:,:),y0fUnitSe(:,:)
double precision, allocatable :: x0fUnitEe(:,:),y0fUnitEe(:,:)
double precision, allocatable :: x0fUnitTe(:,:),y0fUnitTe(:,:)
double precision, allocatable :: x0fUnitSn(:,:),y0fUnitSn(:,:)
double precision, allocatable :: x0fUnitEn(:,:),y0fUnitEn(:,:)
double precision, allocatable :: x0fUnitTn(:,:),y0fUnitTn(:,:)


double precision, allocatable :: sourceField(:,:)
double precision, allocatable :: horFaceDensity(:,:),verFaceDensity(:,:)
double precision, allocatable :: x0fVelocityField(:,:),y0fVelocityField(:,:) ! Velocity at the centroids.
double precision, allocatable :: x0fHorVelocity(:,:),y0fHorVelocity(:,:) ! Velocity at horizontal faces.
double precision, allocatable :: x0fVerVelocity(:,:),y0fVerVelocity(:,:) ! Velocity at vertical faces.
double precision, allocatable :: md0tH0rP0sDir(:,:),md0tVerP0sDir(:,:)


double precision, allocatable :: aE(:,:),aW(:,:),aS(:,:),aN(:,:),aC(:,:),bC(:,:)
integer iteration
double precision, allocatable :: phi(:,:),phiNew(:,:),bCmodified(:,:)
double precision, allocatable :: ldVerSweep(:),mdVerSweep(:),udVerSweep(:),sourceVerSweep(:)
double precision, allocatable :: ldHorSweep(:),mdHorSweep(:),udHorSweep(:),sourceHorSweep(:)

logical smartSchemeS0lver
double precision, allocatable :: md0tePhie(:,:),md0tnPhin(:,:),md0tsPhis(:,:),md0twPhiw(:,:)
double precision, allocatable :: md0tePhieSmart(:,:),md0tnPhinSmart(:,:),md0tsPhisSmart(:,:),md0twPhiwSmart(:,:)

double precision, allocatable :: x0fGradientAtCell(:,:),y0fGradientAtCell(:,:)

double precision ::  phiError


end module 
