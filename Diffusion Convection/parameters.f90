! This module contains all parameters
! ***********************************
module parameters
implicit none


double precision, parameter :: pi=dacos(-1.0D0)
double precision, parameter :: densityField = 1.0D0 ! The density field is constant (incompressible flow).
double precision, parameter :: tolerance = 10**(-4.0D0),relaxationFact0r = 1
integer,parameter :: maxIterations = 3000
double precision, parameter :: epsilon = 10**(-50.0D0)


end module 
