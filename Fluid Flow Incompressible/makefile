objects =  userInput.o getFields.o boundEq.o tfi.o variables.o geomQuant.o surfVec.o getConvectiveCoeff.o allocateMemoryForArrays.o myTDMA.o main.o solvef0rPhi.o getUpwindFlux0nAllFaces.o getSmartFlux0nAllFaces.o getGradient.o getCrossDiffusion.o getSourcef0rRobin.o getDiffusiveCoeff.o sweepingTDMA.o constructBoundaryConditions.o solveNavierStokesf0rVelocityInXdiretion.o getSmartFlux0nAllFacesForXmomentumEq.o getGradientOfuForXmomentumEq.o getCrossDiffusionForXmomentumEq.o  getDiffusiveCoeffForXmomentumEq.o sweepingTDMAForXmomentumEq.o getConvectiveCoeffForXmomentumEq.o   computePressureGradient.o addToSourcePressureAndShearContributionForXmomentumEq.o getUpwindFlux0nAllFacesForXmomentumEq.o solveNavierStokesf0rVelocityInYdiretion.o getSmartFlux0nAllFacesForYmomentumEq.o getGradientOfvForYmomentumEq.o getCrossDiffusionForYmomentumEq.o  getDiffusiveCoeffForYmomentumEq.o sweepingTDMAForYmomentumEq.o getConvectiveCoeffForYmomentumEq.o addToSourcePressureAndShearContributionForYmomentumEq.o getUpwindFlux0nAllFacesForYmomentumEq.o rhieCh0wInterpolation.o solvePressureCorrectionEquation.o sweepingTDMAForPressureCorrection.o correctPressureAndVelocity.o 
a.out: $(objects)
	gfortran $(objects)

myTDMA.o : myTDMA.f90
	gfortran -c myTDMA.f90
boundEq.mod: boundEq.o boundEq.f90
	gfortran -c boundEq.f90
boundEq.o: boundEq.f90
	gfortran -c boundEq.f90
userInput.mod : userInput.o userInput.f90
	gfortran -c userInput.f90
userInput.o : variables.mod userInput.f90
	gfortran -c userInput.f90
variables.mod: variables.o variables.f90
	gfortran -c variables.f90
variables.o: variables.f90
	gfortran -c variables.f90
tfi.o: boundEq.mod tfi.f90
	gfortran -c tfi.f90
allocateMemoryForArrays.o: variables.mod allocateMemoryForArrays.f90
	gfortran -c allocateMemoryForArrays.f90
geomQuant.o: variables.mod geomQuant.f90
	gfortran -c geomQuant.f90
surfVec.o: variables.mod surfVec.f90
	gfortran -c surfVec.f90
getFields.o: userInput.mod variables.mod getFields.f90 
	gfortran -c getFields.f90
getConvectiveCoeff.o : userInput.mod variables.mod getConvectiveCoeff.f90
	gfortran -c getConvectiveCoeff.f90
solvef0rPhi.o : userInput.mod variables.mod solvef0rPhi.f90
	gfortran -c solvef0rPhi.f90
getUpwindFlux0nAllFaces.o : variables.mod getUpwindFlux0nAllFaces.f90
	gfortran -c getUpwindFlux0nAllFaces.f90
getSmartFlux0nAllFaces.o : userInput.mod variables.mod getSmartFlux0nAllFaces.f90
	gfortran -c getSmartFlux0nAllFaces.f90
getGradient.o : variables.mod userInput.mod getGradient.f90
	gfortran -c getGradient.f90
getCrossDiffusion.o : variables.mod userInput.mod getCrossDiffusion.f90
	gfortran -c getCrossDiffusion.f90
getSourcef0rRobin.o : variables.mod userInput.mod getSourcef0rRobin.f90
	gfortran -c getSourcef0rRobin.f90
getDiffusiveCoeff.o : variables.mod userInput.mod getDiffusiveCoeff.f90
	gfortran -c getDiffusiveCoeff.f90
sweepingTDMA.o : variables.mod userInput.mod sweepingTDMA.f90
	gfortran -c sweepingTDMA.f90
constructBoundaryConditions.o : variables.mod constructBoundaryConditions.f90
	gfortran -c constructBoundaryConditions.f90
solveNavierStokesf0rVelocityInXdiretion.o: variables.mod userInput.mod solveNavierStokesf0rVelocityInXdiretion.f90
	gfortran -c solveNavierStokesf0rVelocityInXdiretion.f90
getSmartFlux0nAllFacesForXmomentumEq.o : userInput.mod variables.mod getSmartFlux0nAllFacesForXmomentumEq.f90
	gfortran -c getSmartFlux0nAllFacesForXmomentumEq.f90
getGradientOfuForXmomentumEq.o : variables.mod userInput.mod getGradientOfuForXmomentumEq.f90
	gfortran -c getGradientOfuForXmomentumEq.f90
getCrossDiffusionForXmomentumEq.o : variables.mod userInput.mod getCrossDiffusionForXmomentumEq.f90
	gfortran -c getCrossDiffusionForXmomentumEq.f90
getDiffusiveCoeffForXmomentumEq.o : variables.mod userInput.mod getDiffusiveCoeffForXmomentumEq.f90
	gfortran -c getDiffusiveCoeffForXmomentumEq.f90
sweepingTDMAForXmomentumEq.o : variables.mod userInput.mod sweepingTDMAForXmomentumEq.f90
	gfortran -c sweepingTDMAForXmomentumEq.f90
getConvectiveCoeffForXmomentumEq.o : userInput.mod variables.mod getConvectiveCoeffForXmomentumEq.f90
	gfortran -c getConvectiveCoeffForXmomentumEq.f90
computePressureGradient.o : userInput.mod variables.mod computePressureGradient.f90
	gfortran -c computePressureGradient.f90
addToSourcePressureAndShearContributionForXmomentumEq.o : userInput.mod variables.mod addToSourcePressureAndShearContributionForXmomentumEq.f90
	gfortran -c addToSourcePressureAndShearContributionForXmomentumEq.f90
getUpwindFlux0nAllFacesForXmomentumEq.o : userInput.mod variables.mod getUpwindFlux0nAllFacesForXmomentumEq.f90
	gfortran -c getUpwindFlux0nAllFacesForXmomentumEq.f90
solveNavierStokesf0rVelocityInYdiretion.o: variables.mod userInput.mod solveNavierStokesf0rVelocityInYdiretion.f90
	gfortran -c solveNavierStokesf0rVelocityInYdiretion.f90
getSmartFlux0nAllFacesForYmomentumEq.o : userInput.mod variables.mod getSmartFlux0nAllFacesForYmomentumEq.f90
	gfortran -c getSmartFlux0nAllFacesForYmomentumEq.f90
getGradientOfvForYmomentumEq.o : variables.mod userInput.mod getGradientOfvForYmomentumEq.f90
	gfortran -c getGradientOfvForYmomentumEq.f90
getCrossDiffusionForYmomentumEq.o : variables.mod userInput.mod getCrossDiffusionForYmomentumEq.f90
	gfortran -c getCrossDiffusionForYmomentumEq.f90
getDiffusiveCoeffForYmomentumEq.o : variables.mod userInput.mod getDiffusiveCoeffForYmomentumEq.f90
	gfortran -c getDiffusiveCoeffForYmomentumEq.f90
sweepingTDMAForYmomentumEq.o : variables.mod userInput.mod sweepingTDMAForYmomentumEq.f90
	gfortran -c sweepingTDMAForYmomentumEq.f90
getConvectiveCoeffForYmomentumEq.o : userInput.mod variables.mod getConvectiveCoeffForYmomentumEq.f90
	gfortran -c getConvectiveCoeffForYmomentumEq.f90
addToSourcePressureAndShearContributionForYmomentumEq.o : userInput.mod variables.mod addToSourcePressureAndShearContributionForYmomentumEq.f90
	gfortran -c addToSourcePressureAndShearContributionForYmomentumEq.f90
getUpwindFlux0nAllFacesForYmomentumEq.o : userInput.mod variables.mod getUpwindFlux0nAllFacesForYmomentumEq.f90
	gfortran -c getUpwindFlux0nAllFacesForYmomentumEq.f90
rhieCh0wInterpolation.o: userInput.mod variables.mod rhieCh0wInterpolation.f90
	gfortran -c rhieCh0wInterpolation.f90
solvePressureCorrectionEquation.o : userInput.mod variables.mod solvePressureCorrectionEquation.f90
	gfortran -c solvePressureCorrectionEquation.f90
sweepingTDMAForPressureCorrection.o : variables.mod userInput.mod sweepingTDMAForPressureCorrection.f90
	gfortran -c sweepingTDMAForPressureCorrection.f90
correctPressureAndVelocity.o : variables.mod userInput.mod correctPressureAndVelocity.f90
	gfortran -c correctPressureAndVelocity.f90
main.o: variables.mod userInput.mod tfi.o geomQuant.o surfVec.o 
	gfortran -c main.f90
clean:
	rm *.o *.mod a.out