! This subroutine allocates memory for all arrays.
! ***********************************************
subroutine allocateArrays()
use variablesDeclaration
implicit none

allocate(X(m,n))
allocate(Y(m,n))
allocate(cellVolume(m-1,n-1))
allocate(x0fCellCentroid(m-1,n-1))
allocate(y0fCellCentroid(m-1,n-1))
allocate(x0fHorFaceCenter(m-1,n))
allocate(y0fHorFaceCenter(m-1,n))
allocate(x0fVerFaceCenter(m,n-1))
allocate(y0fVerFaceCenter(m,n-1))
allocate(EastGeomFactor(m-1,n-1))
allocate(WestGeomFactor(m-1,n-1))
allocate(NorthGeomFactor(m-1,n-1))
allocate(SouthGeomFactor(m-1,n-1))
allocate(horFaceArea(m-1,n))
allocate(verFaceArea(m,n-1))
allocate(x0fUnitSe(m,n-1))
allocate(y0fUnitSe(m,n-1))
allocate(x0fUnitEe(m,n-1))
allocate(y0fUnitEe(m,n-1))
allocate(x0fUnitTe(m,n-1))
allocate(y0fUnitTe(m,n-1))
allocate(x0fUnitSn(m-1,n))
allocate(y0fUnitSn(m-1,n))
allocate(x0fUnitEn(m-1,n))
allocate(y0fUnitEn(m-1,n))
allocate(x0fUnitTn(m-1,n))
allocate(y0fUnitTn(m-1,n))
allocate(normEe(m,n-1))
allocate(normTe(m,n-1))
allocate(normEn(m-1,n))
allocate(normTn(m-1,n))
allocate(horFaceDensity(m-1,n))
allocate(verFaceDensity(m,n-1))
allocate(x0fHorVelocity(m-1,n))
allocate(y0fHorVelocity(m-1,n))
allocate(x0fVerVelocity(m,n-1))
allocate(y0fVerVelocity(m,n-1))
allocate(sourceField(m-1,n-1))
allocate(md0tH0rP0sDir(m-1,n),md0tVerP0sDir(m,n-1))
allocate(aE(m-1,n-1),aW(m-1,n-1),aS(m-1,n-1),aN(m-1,n-1),bC(m-1,n-1),aC(m-1,n-1),bCt(m-1,n-1),aCt(m-1,n-1))
allocate(phi(m-1,n-1),phiNew(m-1,n-1),bCmodified(m-1,n-1), phiPreviousTime(m-1,n-1))
allocate(ldVerSweep(n-1),mdVerSweep(n-1),udVerSweep(n-1),sourceVerSweep(n-1))
allocate(ldHorSweep(m-1),mdHorSweep(m-1),udHorSweep(m-1),sourceHorSweep(m-1))
allocate(md0tePhie(m-1,n-1),md0tnPhin(m-1,n-1),md0twPhiw(m-1,n-1),md0tsPhis(m-1,n-1))
allocate(md0tePhieSmart(m-1,n-1),md0tnPhinSmart(m-1,n-1),md0twPhiwSmart(m-1,n-1),md0tsPhisSmart(m-1,n-1))
allocate(x0fGradientAtCell(m-1,n-1),y0fGradientAtCell(m-1,n-1))
allocate(horInterfaceDiff(m-1,n),verInterfaceDiff(m,n-1),gammaField(m-1,n-1))
allocate(phiAnalytical(m-1,n-1))
end subroutine 
