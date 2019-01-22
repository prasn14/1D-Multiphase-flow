Module Input
Use Grids
implicit none
public
! Defining the Grid details
double precision, parameter :: dx = L/xcells
double precision, parameter :: h_ref = L ! The reference height 
! Time steps
double precision , parameter :: time = 0.5 ! Time till which the solution needs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CORRECT BETWEEN VARIABLE TIME STEPPING AND A CONSTANT TIME STEPPING METHOD BY
! PROPER CHOICE OF EITHER DT OR CFL
!double precision, parameter :: dt = 0.15e-3/1.0 ! arbitrarily chosen intially
! to be computed (units of sec)
! Total number of timesteps
!integer, parameter :: timesteps = time/dt
double precision, parameter :: cfl = 0.5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer, parameter :: test_case = 1
! 1 - Two fluid shock tube problem 
! 2 - Moving contact discontinuity problem (same mass averaged velocity across
! the jump
! 3 - Stationay shock case
! 4 - Moving contact discontinuity (same volumetric velocity across the jump -
! Incompressible flow assumption
! 5 - Flow through a vertical duct with gas injection - Make sure to switch the gravity on 
! 6 - Sedimentation test case 
integer, parameter :: Flux_method = 3
! 1 - For LxF, 2 - FUP, 3 - HLL , 4 - HLLC
integer, parameter :: flux_order_switch = 1
! 0 - For getting the required first order flux 
! 1 - To have the FCT switch on
! 2 - To retain the second order Lax Wendroff Flux (obtained using two step L-W)
! Boundary method
integer, parameter :: boundary_variable = 0
! 0 - for Ghost cell  1 - for giving the fluxes at the boundaries 
! Switch for the boundary condition
integer, parameter :: slip_type = 1
! 0 - No slip, 1 - drift flux model
double precision, parameter :: flow_velocity_gas_kick = 0.0
! velocity of the liquid flowing initially without the gas being present
double precision, parameter :: gas_kick_time = 0.0
! Time when the gas starts kicking in  (for test case 5)
double precision, parameter :: shut_in_time = 30
! Time when the channel is shut in  (for test case 5)
double precision, parameter :: pi = 4.0*atan(1.0)
double precision, parameter :: pipe_dia = 0.1 ! units in m 
double precision, parameter :: Area = pi*(pipe_dia)**2 ! units in m^2
! Switch for the inclusion of the wall friction term
integer, parameter :: wall_friction_switch = 0
! 0 - without friction 1 - with wall friction
! Liquid properties 
double precision, parameter :: p_ref = 1e5 ! units in Pa (for both liquid and
! gas its the same as we have assumed mechanical equilibrium)
! AMBIENT ATMOSPHERIC CONDITION AT THE TOP OF THE WELL
double precision, parameter :: rho_liquid_ref = 850 ! units kg/m^3
double precision, parameter :: a_liquid = 1500
! Gas phase properties
double precision, parameter :: rho_gas_ref = 1.1111 ! units in kg/m^3
double precision, parameter :: gamma_gas = 1.0! Polytropic gas constant 
double precision, parameter :: a_gas = sqrt(gamma_gas*p_ref/rho_gas_ref)

! Acceleration due to gravity 
double precision, parameter :: g = 9.81 ! units m/s^2
! Inclination of the pipe 
double precision, parameter :: theta =  0.0 !pi/2.0  !for a horizontal duct
! Frequency interms of timesteps with which the data has to be written
integer, parameter :: frequency = 1
! Selection on the number of initial conditions 
integer, parameter :: initial_list = 0
! PARAMETERS TO BE INCORPORATED WITH THE PHASE CONDENSATION OF GAS TO LIQUID
double precision, parameter :: M_oil = 158 ! Molecular weight of the liquid 
double precision, parameter :: M_m = 16 ! Molecular weight of the gas
!component (methane is the molecular weight)
double precision, parameter :: p_crit = 2e7 !(critical pressure for the K
!factor)
double precision, parameter :: c_l_l = 1e-4
integer, parameter :: phase_transfer_switch = 1
! 0 for flow without phase transition    1 - for flow with condensation of the
! vapor into liquid phase






















end module Input

