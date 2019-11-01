objects = fields.o getFields.o boundEq.o tfi.o variables.o geomQuant.o surfVec.o getConvectiveCoeff.o allocate.o myTDMA.o main.o iterativeSolver.o getUpwindFlux0nAllFaces.o getSmartFlux0nAllFaces.o boundaryConditions.o parameters.o getGradient.o getCrossDiffusion.o getSourcef0rRobin.o getDiffusiveCoeff.o getAnalyticalSolution.o modifyCoefficientsForTransientSolving.o
a.out: $(objects)
	gfortran $(objects)

myTDMA.o : myTDMA.f90
	gfortran -c myTDMA.f90
boundEq.mod: boundEq.o boundEq.f90
	gfortran -c boundEq.f90
boundEq.o: boundEq.f90
	gfortran -c boundEq.f90
parameters.mod : parameters.o parameters.f90
	gfortran -c parameters.f90
parameters.o : parameters.f90
	gfortran -c parameters.f90
variables.mod: variables.o variables.f90
	gfortran -c variables.f90
variables.o: variables.f90
	gfortran -c variables.f90
boundaryConditions.mod: boundaryConditions.o boundaryConditions.f90
	gfortran -c boundaryConditions.f90
boundaryConditions.o : boundaryConditions.f90
	gfortran -c boundaryConditions.f90
tfi.o: boundEq.mod tfi.f90
	gfortran -c tfi.f90
allocate.o: variables.mod allocate.f90
	gfortran -c allocate.f90
geomQuant.o: variables.mod geomQuant.f90
	gfortran -c geomQuant.f90
surfVec.o: variables.mod surfVec.f90
	gfortran -c surfVec.f90
fields.mod:  fields.o fields.f90
	gfortran -c fields.f90
fields.o: parameters.mod variables.mod fields.f90
	gfortran -c fields.f90
getFields.o: fields.mod variables.mod getFields.f90 
	gfortran -c getFields.f90
getConvectiveCoeff.o : parameters.mod boundaryConditions.mod variables.mod getConvectiveCoeff.f90
	gfortran -c getConvectiveCoeff.f90
iterativeSolver.o : parameters.mod boundaryConditions.mod variables.mod iterativeSolver.f90
	gfortran -c iterativeSolver.f90
getUpwindFlux0nAllFaces.o : variables.mod getUpwindFlux0nAllFaces.f90
	gfortran -c getUpwindFlux0nAllFaces.f90
getSmartFlux0nAllFaces.o : parameters.mod boundaryConditions.mod variables.mod getSmartFlux0nAllFaces.f90
	gfortran -c getSmartFlux0nAllFaces.f90
getGradient.o : variables.mod parameters.mod boundaryConditions.mod getGradient.f90
	gfortran -c getGradient.f90
getCrossDiffusion.o : variables.mod parameters.mod boundaryConditions.mod getCrossDiffusion.f90
	gfortran -c getCrossDiffusion.f90
getSourcef0rRobin.o : variables.mod parameters.mod boundaryConditions.mod getSourcef0rRobin.f90
	gfortran -c getSourcef0rRobin.f90
getDiffusiveCoeff.o : variables.mod parameters.mod boundaryConditions.mod getDiffusiveCoeff.f90
	gfortran -c getDiffusiveCoeff.f90
getAnalyticalSolution.o : variables.mod getAnalyticalSolution.f90
	gfortran -c getAnalyticalSolution.f90
modifyCoefficientsForTransientSolving.o : variables.mod parameters.mod boundaryConditions.mod modifyCoefficientsForTransientSolving.f90
	gfortran -c modifyCoefficientsForTransientSolving.f90
main.o: parameters.mod variables.mod tfi.o geomQuant.o surfVec.o 
	gfortran -c main.f90
clean:
	rm *.o *.mod a.out