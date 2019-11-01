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
x0fVelocityFunction = 0.0D0
end function

function y0fVelocityFunction(x,y) ! Scalar Field.
double precision  y0fVelocityFunction,x,y
y0fVelocityFunction = 0.0D0
end function

function sourceFunction(x,y) ! Scalar Field.
double precision sourceFunction,x,y
sourceFunction = 0.0D0
end function

function gammaFunction(x,y) ! Scalar Field.
double precision gammaFunction,x,y
gammaFunction = 1.4D-3
end function

end module
