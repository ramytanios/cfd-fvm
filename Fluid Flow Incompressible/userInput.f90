! This module contains all the user inputs:
! Boundary Conditions, parameters and fields.
! ****************************************
module userInput
use variablesDeclaration
implicit none

!
! <<<<<<<<<<<<<<<< PARAMETERS >>>>>>>>>>>>>>>>>>>
!
double precision, parameter :: pi=dacos(-1.0D0)

! The density of the fluid for incompressible flows. 
double precision, parameter :: densityField = 0.6

! The viscosity of the fluid. 
double precision, parameter :: viscosity = 5.0D-5

! The specific heat of the fluid. 
double precision, parameter :: Cp = 1.03

! The thermal conductivity of the fluid.
double precision, parameter :: kappa = 0.036

! Tolerance for solving phi on the grid.
double precision, parameter :: tolerance = 10**(-4.0D0)

! The maximum number of iterations for solving phi on the grid. 
integer,parameter :: maxIterations = 1000

! Small number to be added to denominators to avoid division by 0.
double precision, parameter :: epsilon = 10**(-7.0D0)

! The relaxation factor.
double precision, parameter :: relaxationFact0r = 1
! ------------------------------------------------
! ------------------------------------------------
! ------------------------------------------------
! ------------------------------------------------
! ------------------------------------------------
!

! ------------------------------------------------
! ------------------------------------------------
! ------------------------------------------------
! ------------------------------------------------
! ------------------------------------------------
!
! <<<<<<<<<<<<<<<< FIELDS >>>>>>>>>>>>>>>>>>>
!
contains 
function x0fVelocityFunction(x,y) ! Scalar Field.
double precision x0fVelocityFunction,x,y
x0fVelocityFunction = 0!cos(pi/4.0D0)
end function

function y0fVelocityFunction(x,y) ! Scalar Field.
double precision  y0fVelocityFunction,x,y
y0fVelocityFunction = 0!sin(pi/4.0D0)
end function

! The source term field.
function sourceFunction(x,y) ! Scalar Field.
double precision sourceFunction,x,y
sourceFunction = 0
end function

! The diffusivity field gamma.
function gammaFunction(x,y) ! Scalar Field.
double precision gammaFunction,x,y
gammaFunction = kappa / densityField / Cp
end function
! ------------------------------------------------




end module 
