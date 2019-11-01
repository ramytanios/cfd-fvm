! This module contains all the equations of all fields:
! Velocity Field ,Source Field, Diffusivity Field.
! The density field is assumed to be constant.
! ******************************************************************

module fields
use variablesDeclaration
use parameters
implicit none

contains 

! The velocity field is a vector field.
function x0fVelocityFunction(x,y) ! Scalar Field.
double precision x0fVelocityFunction,x,y
x0fVelocityFunction = cos(pi/4.0D0)
end function

function y0fVelocityFunction(x,y) ! Scalar Field.
double precision  y0fVelocityFunction,x,y
y0fVelocityFunction = sin(pi/4.0D0)
end function

function sourceFunction(x,y) ! Scalar Field.
double precision sourceFunction,x,y
sourceFunction = (2*x-0.2*y)
end function

function gammaFunction(x,y) ! Scalar Field.
double precision gammaFunction,x,y
gammaFunction = x**2+exp(0.1*y)
end function

end module
