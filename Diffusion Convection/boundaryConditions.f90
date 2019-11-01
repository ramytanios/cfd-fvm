! This module contains the boundary conditions.
! The order in each array is North-West-South-East.
! *************************************************
module boundaryConditions
implicit none

double precision, parameter :: phiAtInlet0ne = 0.0D0, phiAtInletTw0 = 50.0D0 ! Convection.

double precision, parameter :: phiAtInletThree= 0.0D0,phiAtInletF0ur= 0.0D0 ! Convection.

logical,parameter,dimension(4) :: outletCondition=[.FALSE., .FALSE., .FALSE., .TRUE.] ! Convection.

logical,parameter,dimension(4) :: inletCondition=[.FALSE., .TRUE., .FALSE., .FALSE. ] ! Convection.

logical,parameter,dimension(4) :: wallCondition=[.TRUE., .FALSE., .TRUE., .FALSE.] ! Convection.

logical,parameter,dimension(4) :: dirichlet=[.TRUE., .TRUE., .TRUE., .TRUE. ] ! Diffusion.
double precision,parameter,dimension(4) :: dirichletValue = [400.0D0,400.0D0,320.0D0,400.0D0]

logical,parameter,dimension(4) :: vonNeumann=[.FALSE., .FALSE., .FALSE., .FALSE.] ! Diffusion.
double precision,parameter,dimension(4) :: x0fNeumannGrad = [0.0D0,0.0D0,0.0D0,0.0D0]
double precision,parameter,dimension(4) :: y0fNeumannGrad = [0.0D0,0.0D0,0.0D0,0.0D0]

logical,parameter,dimension(4) :: robin(4)=[.FALSE., .FALSE., .FALSE., .FALSE.] ! Diffusion.
double precision,parameter :: phiInf(4)=[300.0D0,200.0D0,300.0D0,300.0D0] ! Diffusion.
double precision,parameter :: hInf(4)=[15.0D0,15.0D0,15.0D0,15.0D0] ! Diffusion.


end module 
