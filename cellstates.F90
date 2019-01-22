module cellstates
use input
use Grids
implicit none
public :: U_vector_initial
public :: phase_transfer_properties
public :: pressure_from_U
public :: alpha_from_U
public :: cfl_computation
public :: wavespeed_calculator
public :: maximum_from_array
public :: minimum_from_array
public :: drift_flux_law
public :: drift_flux_law_uv
public :: appx_Roe_avg_waves
public :: flux_calculator
public :: flux_interface_calculator_first_order
public :: flux_interface_calculator_second_order
public :: omega_flux_cal
public :: source_vector
public :: boundary_condition
public :: flux_boundary_condition
public :: soft_boundary
public :: hard_boundary
public :: U_to_variables
public :: writeout
contains
!!!!!!!!!!!!!!! Direct initial distribution of conservative variables !!!!!!!!!!!
subroutine U_vector_initial(U_vector,x)
implicit none
double precision, dimension(1:3,0:xcells+1) :: U_vector
double precision, dimension(0:xcells+1) :: x
double precision, dimension(0:xcells+1) :: pressure
double precision, dimension(0:xcells+1) :: alpha
double precision :: alpha_left,alpha_right,u_M
double precision :: u_V,C_l,vs_l,C_r,vs_r,v_minus_u
double precision :: u_l,u_r,v_l,v_r
double precision :: alpha_constant,K,lambda
double precision :: exit_pressure
integer :: j
selectcase(test_case)
case(1) ! Two fluid shock tube problem    (Test case obtained from reference for
! Franck's code)
! Giving the left state variables
U_vector(1,0:xcells/2) = 453.197
U_vector(2,0:xcells/2) = 3.19503885*2
U_vector(3,0:xcells/2) = 453.197*24.8074
! Giving the right state variables 
U_vector(1,xcells/2+1:xcells+1) = 454.915
U_vector(2,xcells/2+1:xcells+1) = 4.913082 
U_vector(3,xcells/2+1:xcells+1) = 454.915*1.7461
case(2) ! Moving contact discontinuity problem with the same mass averaged
!velocity
! Here we assume that the pressure is constant initially with just the jump in
! alpha for both the left and the right states
alpha_left = 0.3
alpha_right = 0.8
u_M = 1.0
! Giving the left state variables
U_vector(1,0:xcells/2) = alpha_left*rho_gas_ref + (1-alpha_left)*rho_liquid_ref
U_vector(2,0:xcells/2) = alpha_left*rho_gas_ref
U_vector(3,0:xcells/2) = (alpha_left*rho_gas_ref + &
&(1-alpha_left)*rho_liquid_ref)*u_M
! Giving the right state variables 
U_vector(1,xcells/2+1:xcells+1) = alpha_right*rho_gas_ref + (1-alpha_right)*rho_liquid_ref
U_vector(2,xcells/2+1:xcells+1) = alpha_right*rho_gas_ref
U_vector(3,xcells/2+1:xcells+1) = (alpha_right*rho_gas_ref +&
&(1-alpha_right)*rho_liquid_ref)*u_M
case(3) ! Case of a stationary shock or moving shock 
! Giving the left state variables
U_vector(1,0:xcells/2) = 457.35
U_vector(2,0:xcells/2) = 0.045735
U_vector(3,0:xcells/2) = 457.35*(10+0)
! Giving the right state variables 
U_vector(1,xcells/2+1:xcells+1) = 624.07
U_vector(2,xcells/2+1:xcells+1) = 0.062407
U_vector(3,xcells/2+1:xcells+1) = 624.07*(7.3285+0)
case(4) ! Case for a jump in alpha with zero volumetric mean velocity(so that it
!can bie compared with scalar advection for alpha) - ASSUMPTION OF INCOMPRESSIBLE
!FLOW
u_V = 0.5
alpha_left = 1.0
alpha_right = 0.2! Defining the left state variables
U_vector(1,0:xcells/2) = alpha_left*rho_gas_ref + (1.0-alpha_left)*rho_liquid_ref
U_vector(2,0:xcells/2) = alpha_left*rho_gas_ref 
do j = 0,xcells/2
   call drift_flux_law_uv(v_minus_u,C_l,vs_l,u_V,U_vector(:,j))
   v_l = C_l*u_V + vs_l
   if(alpha_left == 1.0)then 
   u_l = 0.0
   elseif(alpha_left == 0.0)then
   u_l = vs_l + u_V*C_l
   else
   u_l = ((1.0-C_l*alpha_left)*u_V - alpha_left*vs_l)/(1.0 - alpha_left) 
   endif
!  print*,'vl-ul and v_minus_u is',v_l-u_l,v_minus_u
end do
! Giving the right state variables 
U_vector(1,xcells/2+1:xcells+1) = alpha_right*rho_gas_ref + (1.0-alpha_right)*rho_liquid_ref  
U_vector(2,xcells/2+1:xcells+1) = alpha_right*rho_gas_ref 
do j = xcells/2+1,xcells+1
   call drift_flux_law_uv(v_minus_u,C_r,vs_r,u_V,U_vector(:,j))
   v_r = C_r*u_V + vs_r
   if(alpha_right == 1.0)then
    u_r = 0.0
   elseif(alpha_right == 0.0)then
    u_r = C_r*u_V + vs_r
    else
    u_r = ((1-C_r*alpha_right)*u_V - alpha_right*vs_r)/(1.0 - alpha_right)
    endif
!  print*,'vr-ur and v-u is',v_r-u_r,v_minus_u
end do
case(5)
! Gas flow injection case - Initially water column is at rest 
!
exit_pressure = p_ref*10
U_vector(1,xcells+1) = rho_liquid_ref + (exit_pressure-p_ref)/(a_liquid**2)
do j = xcells,0,-1
   U_vector(1,j) = U_vector(1,j+1)*exp(g*(x(j+1)-x(j))/a_liquid**2)
   U_vector(2,j) = 0.0
   U_vector(3,j) = 0.0
end do
case(6)
alpha(:) = 0.5

do j = 0,xcells+1
   U_vector(1,j) = (1-alpha(j))*rho_liquid_ref*1.001 + alpha(j)*rho_gas_ref*20
   U_vector(2,j) = alpha(j)*rho_gas_ref*20
   U_vector(3,j) = 0.0 
end do

endselect
return
end subroutine U_vector_initial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Sub routine to compute phase tranfer properties!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine phase_transfer_properties(U_vector,alpha,p,rho_l,rho_g,w_m_l)
implicit none
double precision, dimension(1:3) :: U_vector
! Mole fraction of the two components (represented with subscripts o and m)
double precision :: z_o ! Mole fraction of the oil component
double precision :: z_m ! Mole fraction of the methane (gaseous component)
double precision :: N ! Total number of moles present in the mixture 
double precision :: y_o ! Mole fraction of oil in VAPOUR phase
double precision :: y_m ! Mole fraction of methane in the VAPOUR phase
double precision :: x_o ! Mole fraction of oil in the LIQUID phase
double precision :: x_m ! Mole fraction of methane in the LIQUID phase
double precision :: K_o ! Empirical derived K factors for the oil ( = y_o/x_o)
double precision :: K_m ! K factor for methane component
double precision :: n_g ! Mole fraction of the vapour phase in the entire
                        ! mixture
double precision :: p   ! pressure obtained for the mixture
double precision :: w_o_l,w_o_v,w_m_l,w_m_v ! The mass fraction of both the
! Components in their respective phases
double precision :: A_v,A_l ! Volume occupied by the respective phases
double precision :: alpha   ! Volume fraction of vapour phase
double precision :: rho_l,rho_g ! both the vapour and the liquid phase densities
double precision :: p_min,p_max,eps1,eps2,eps3,error,eps
double precision :: dng_dp,drho_ml_dp,drho_ol_dp,drho_mv_dp
double precision :: rho_ml,rho_ol,rho_mv
double precision :: dxm_dp,dxo_dp
double precision :: f_p, f_dash_p,p_buffer
integer :: iter_num,j
iter_num = 0.0
error = 1.0
! Since the area is the constant
N = U_vector(2)*Area/M_m + (U_vector(1)-U_vector(2))*Area/M_oil
z_o = (U_vector(1) - U_vector(2))*Area/M_oil *(1.0/N) 
z_m = U_vector(2)*Area/M_m*(1.0/N)

! FOR THE SIDE KICK MODEL WITH ONLY CONDENSTATION OF METHANE AND NON VOLATILE
! OIL
y_o = 0.0
y_m = 1.0
K_o = 0.0   ! as y_o = 0.0

eps = 1e-6
error = 1.0
p = p_ref/4 ! Initial data for the Newton Raphson method
! GETTING A NEWTON RAPSHON METHOD TO COMPUTE THE PRESSURE
do while(error>=eps)
   K_m = p_crit/p
   ! Computing n_g and its derivative
   if((z_m*K_m)<=1.0)then
     n_g = 0.0                   ! Region of single phase liquid
     dng_dp = 0.0

    x_m = z_m 
    dxm_dp = 0.0
    x_o = z_o
    dxo_dp = 0.0
   else
     n_g = (K_m*z_m-1.0)/(K_m-1.0)
     dng_dp = p_crit*(z_m-1.0)/(p_crit-p)**2  ! only valid for K_o = 0.0

     x_m = 1.0/K_m
     dxm_dp = 1.0/p_crit
     x_o = 1.0-1.0/K_m
     dxo_dp = -1.0/p_crit
   endif
! Methane Vapour
   rho_mv = rho_gas_ref*p/p_ref
   drho_mv_dp = rho_gas_ref/p_ref 
! Methane liquid
   rho_ml = rho_liquid_ref*(1.0 + (p-p_ref)*c_l_l/p_ref)
   drho_ml_dp =  c_l_l/p_ref*rho_liquid_ref
! Oil liquid 
   rho_ol = rho_liquid_ref + (p-p_ref)/a_liquid**2
   drho_ol_dp =  1.0/a_liquid**2


! function for pressure based on the volume filling condition (for non volatile
! liquid
!print*,"Pressuse before N-R",p
f_p = n_g*N*M_m/rho_mv + (1.0-n_g)*N*(M_m*x_m/rho_ml + M_oil*x_o/rho_ol) - Area
f_dash_p = dng_dp*N*M_m/rho_mv - n_g*N*M_m/rho_mv**2*drho_mv_dp - &
         & dng_dp*N*(M_m*x_m/rho_ml + M_oil*x_o/rho_ol) + (1-n_g)*N*&
         & (M_m/rho_ml*dxm_dp + M_oil/rho_ol*dxo_dp) - (1.0-n_g)*N*&
         & (M_m*x_m*drho_ml_dp/rho_ml**2 + M_oil*x_o*drho_ol_dp/rho_ol**2)     
p_buffer = p - f_p/f_dash_p
error = f_p
if(p<0.0)then
   print*," Getting negative pressure stop the iteration"
   p = p_buffer
endif
p = p_buffer
enddo
!print*,"Iterated pressure is ",p
! COMPUTING ALL THE DETAILS AFTER THE FINAL PRESSURE
   K_m = p_crit/p
 ! Computing n_g
   if((z_o*K_o+z_m*K_m)<=1.0)then
     n_g = 0.0
   else
     n_g = -((K_o-1.0)*z_o + (K_m-1.0)*z_m)/((K_o-1.0)*(K_m-1.0))
   endif
   ! Computing the x ratio
   if(p/p_crit<z_m)then
     x_m = 1.0/K_m
     x_o = 1.0-1.0/K_m
   elseif(p/p_crit>=z_m)then   ! Region of single phase liquid
     x_m = z_m 
     x_o = z_o
   endif  

rho_g = rho_gas_ref*p/p_ref
rho_l = (M_oil*x_o + M_m*x_m)/(M_oil*x_o/(rho_liquid_ref + &
&(p-p_ref)/a_liquid**2) + M_m*x_m/(rho_liquid_ref + c_l_l*(p-p_ref)/p_ref))

! Computation of volume occupied by the vapour phase
A_v = n_g*N*(M_m*y_m)/(p/p_ref*rho_gas_ref)

! Computatio of volume occupied by the liquid phase
A_l = (1.0-n_g)*N*(M_oil*x_o/(rho_liquid_ref + &
&(p-p_ref)/a_liquid**2) + M_m*x_m/(rho_liquid_ref + c_l_l*(p-p_ref)/p_ref)) 

alpha = A_v/(A_l+A_v)

! Getting the mass fractions of both the components in the respective phases#
! LIQUID PHASE
w_o_l = M_oil*x_o/(M_oil*x_o + M_m*x_m)
w_m_l = 1.0 - w_o_l  

! VAPOUR PHASE
w_o_v = 0.0 
w_m_v = 1.0-w_o_v 


return
end subroutine phase_transfer_properties
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Sub routine to calculate the pressure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pressure_from_U(p,rho_l,rho_g,w_m_l,U)
implicit none 
double precision, dimension(1:3) :: U
double precision :: p,rho_l,rho_g
double precision :: z,alpha,w_m_l
   selectcase(phase_transfer_switch)
   case(0)
   z = p_ref - a_liquid**2*(rho_liquid_ref) + U(2)*a_gas**2 + &
       &a_liquid**2*(U(1)-U(2)) 
   p = z/2.0 + 1.0/2.0*sqrt(z**2 - 4*U(2)*a_gas**2*(p_ref - &
          &a_liquid**2*rho_liquid_ref)) 
   w_m_l = 0.0 
   rho_g = rho_gas_ref*p/p_ref
   rho_l = rho_liquid_ref + (p-p_ref)/a_liquid**2 
   case(1)
   call phase_transfer_properties(U,alpha,p,rho_l,rho_g,w_m_l)
   end select
return 
end subroutine pressure_from_U
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Sub routine to calculate alpha from p,U!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alpha_from_U(alpha,U)
implicit none 
double precision, dimension(1:3) :: U
double precision :: p,alpha,eps,w_m_l,rho_l,rho_g
call pressure_from_U(p,rho_l,rho_g,w_m_l,U)
   selectcase(phase_transfer_switch)
   case(0)
   alpha = 1 - (U(1)-U(2))*a_liquid**2/(p - p_ref + rho_liquid_ref*a_liquid**2)
   case(1)
   call phase_transfer_properties(U,alpha,p,rho_l,rho_g,w_m_l)
   end select
if(alpha > 1.0)then
print*,'Alpha is greater than 1'
print*,'The pressure,U1,U2 is',p,U(1),U(2)
endif
return
end subroutine alpha_from_U
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! Sub routine to compute the cfl based on the cell data !!!!!!!!!
subroutine cfl_computation(U_vector,cfl_wave)
implicit none
double precision, dimension(1:3,0:xcells+1) :: U_vector
double precision, dimension(0:xcells+1) :: lambda_1 ! u-a wave
double precision, dimension(0:xcells+1) :: lambda_2 ! u wave
double precision, dimension(0:xcells+1) :: lambda_3 ! u+a wave 
double precision, dimension(0:xcells+1) :: max_lambda
double precision :: cfl_wave,lambda
integer :: j
do j = 0,xcells+1
   call wavespeed_calculator(U_vector(:,j),lambda_1(j),lambda_2(j),lambda_3(j))
   max_lambda(j) = max(dabs(lambda_1(j)),dabs(lambda_3(j)))
end do   
call maximum_from_array(max_lambda,lambda)
cfl_wave = lambda
return
end subroutine cfl_computation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Subroutine to compute the maximum in an array!!!!!!!!!!!!!!!!
subroutine maximum_from_array(l,l_max)
implicit none
double precision,dimension(0:xcells+1) :: l
double precision :: l_buffer, l_max
integer :: j
l_buffer = l(0)
do j = 1,xcells+1 
   if(l(j)>=l_buffer)then
      l_buffer = l(j)
   endif
end do
l_max = l_buffer
return
end subroutine maximum_from_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine minimum_from_array(l,l_max)
implicit none
double precision,dimension(0:xcells+1) :: l
double precision :: l_buffer, l_max
integer :: j
l_buffer = l(0)
do j = 1,xcells+1
   if(l(j)<=l_buffer)then
      l_buffer = l(j)
   endif
end do
l_max = l_buffer
return
end subroutine minimum_from_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Sub routine to compute the three wave speeds !!!!!!!!!!!!!!!!!!
subroutine wavespeed_calculator(U,l1,l2,l3)
implicit none
double precision, dimension(1:3) :: U
double precision :: l1,l2,l3
double precision :: alpha,c
call alpha_from_U(alpha,U)
if(alpha == 1.0)then
c = a_gas
elseif(alpha == 0.0)then
c = a_liquid
elseif(alpha .ne. 0.0 .or. alpha .ne. 1.0)then
c = sqrt(1.0/U(1)*((U(1)-U(2))*U(2)*a_liquid**2*a_gas**2)/&
&((1-alpha)**2*U(2)*a_gas**2 + alpha**2*(U(1)-U(2))*a_liquid**2))
endif
l1 = U(3)/U(1) - c
l2 = U(3)/U(1)
l3 = U(3)/U(1) + c
return
end subroutine wavespeed_calculator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Sub routine to calculate the drift flux(v-u)!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine drift_flux_law(v_minus_u,C,v_s,U_vector)
implicit none
double precision, dimension(1:3) :: U_vector
double precision :: alpha, C, v_s, v_minus_u
double precision :: u,v,u_V,p,rho_g,w_m_l,rho_l
select case(slip_type)
case(0)! The case of no slip
v_minus_u = 0.0
case(1) ! The case of slip
call alpha_from_U(alpha,U_vector)
! Computation of the coefficient C and v_s based on experimentation observations
if(alpha >=0.and.alpha <= 0.015)then
   C = 1.0
   v_s = 0.0
elseif(alpha>0.015.and.alpha<=0.1)then
   C = 1 + 4.1176*(alpha-0.015)
   v_s = 5.8823*(alpha-0.015)
elseif(alpha>0.1.and.alpha <= 0.4)then
   C = 1.35
   v_s = 0.5 
elseif(alpha>0.4.and.alpha<=1.0)then
   C = 1 + 0.58333*(1-alpha)
   v_s = 0.8333*(1-alpha)
endif
if(alpha == 1.0)then
        v_minus_u = U_vector(3)/(U_vector(1))
        elseif(alpha == 0.0)then
        v_minus_u = -U_vector(3)/U_vector(1) 
        else
!        selectcase(phase_transfer_switch)
!       case(0)
!        v_minus_u = v_s/(1.0-C*alpha) + (C-1)/(1-C*alpha)* (U_vector(3)*(1.0-C*alpha) - &
!        & U_vector(2)*v_s)/(U_vector(1)*(1.0-C*alpha) - U_vector(2)*(1.0-C))
!        case(1)
        call pressure_from_U(p,rho_l,rho_g,w_m_l,U_vector)
        v_minus_u = v_s/(1.0-C*alpha) + (C-1)/(1-C*alpha)* (U_vector(3)*(1.0-C*alpha) - &
        & alpha*rho_g*v_s)/(U_vector(1)*(1.0-C*alpha) - alpha*rho_g*(1.0-C))
!        endselect
endif 

end select

return
end subroutine drift_flux_law
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Sub routine to calculate the drift flux(v-u)!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine drift_flux_law_uv(v_minus_u,C,v_s,u_V,U_vector)
implicit none
double precision, dimension(1:3) :: U_vector
double precision :: alpha, C, v_s, v_minus_u
double precision :: u,v,u_V,p,rho_g,w_m_l,rho_l
select case(slip_type)
case(0)! The case of no slip
v_minus_u = 0.0
case(1) ! The case of slip
call alpha_from_U(alpha,U_vector)
! Computation of the coefficient C and v_s based on experimentation observations
if(alpha >=0.and.alpha <= 0.015)then
   C = 1.0
   v_s = 0.0
elseif(alpha>0.015.and.alpha<=0.1)then
   C = 1 + 4.1176*(alpha-0.015)
   v_s = 5.8823*(alpha-0.015)
elseif(alpha>0.1.and.alpha <= 0.4)then
   C = 1.35
   v_s = 0.5 
elseif(alpha>0.4.and.alpha<=1.0)then
   C = 1 + 0.58333*(1-alpha)
   v_s = 0.8333*(1-alpha)
endif
if(alpha == 1.0)then
        v = v_s
        u = 0.0
        else
        v = C*u_V + v_s
        u = ((1-C*alpha)*u_V - alpha*v_s)/(1-alpha)
endif 
v_minus_u = v-u
call pressure_from_U(p,rho_l,rho_g,w_m_l,U_vector)
U_vector(3) = U_vector(1)*(u) + rho_g*alpha*(v-u) 
end select
return
end subroutine drift_flux_law_uv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Subroutine to calculate the flux with the cell centered values!!!!!!!!!!!!!!
subroutine flux_calculator(F,U)
implicit none 
double precision, dimension(1:3) :: U
double precision, dimension(1:3) :: F
double precision :: pressure,G,I,v_minus_u,C,v_s
double precision :: alpha,w_m_l,rho_g,rho_l
call drift_flux_law(v_minus_u,C,v_s,U)
call pressure_from_U(pressure,rho_l,rho_g,w_m_l,U)
call alpha_from_U(alpha,U)

I = alpha*rho_g*(1.0-(alpha*rho_g)/U(1))
G = I*(v_minus_u)**2 

! Getting the flux variables
F(1) = U(3)
F(2) = U(2)*U(3)/U(1) + I*v_minus_u*(1.0-w_m_l)
F(3) = U(3)**2/U(1) + pressure + G


return
endsubroutine flux_calculator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Subroutine to compute the approximated Roe averages!!!!!!!!!!!!!!!!!!!
subroutine appx_Roe_avg_waves(U_left,U_right,SL,SM,SR)
implicit none
double precision, dimension(1:3) :: U_left
double precision, dimension(1:3) :: U_right
double precision :: rhom_l,alpha_l,rhol_l,rhog_l,um_l,u_l,v_l,p_l
double precision :: rhom_r,alpha_r,rhol_r,rhog_r,um_r,u_r,v_r,p_r
double precision :: SL, SM, SR
double precision :: u_ave,c_ave
double precision :: c_l,c_r,w_m_l
call U_to_variables(U_left,rhom_l,alpha_l,rhol_l,rhog_l,um_l,u_l,v_l,p_l,w_m_l)
call U_to_variables(U_right,rhom_r,alpha_r,rhol_r,rhog_r,um_r,u_r,v_r,p_r,w_m_l)

if(alpha_l==1.0)then
   c_l = a_gas
elseif(alpha_l == 0.0)then
   c_l = a_liquid
else
c_l = sqrt(1.0/U_left(1)*((U_left(1)-U_left(2))*U_left(2)*a_liquid**2*a_gas**2)/&
&((1-alpha_l)**2*U_left(2)*a_gas**2 + alpha_l**2*(U_left(1)-U_left(2))*a_liquid**2))
endif

if(alpha_r == 1.0)then
  c_r = a_gas
elseif(alpha_r == 0.0)then
  c_r = a_liquid
else
c_r = sqrt(1.0/U_right(1)*((U_right(1)-U_right(2))*U_right(2)*a_liquid**2*a_gas**2)/&
&((1-alpha_r)**2*U_right(2)*a_gas**2 + alpha_r**2*(U_right(1)-U_right(2))*a_liquid**2))
endif

u_ave = (sqrt(rhom_l)*um_l + sqrt(rhom_r)*um_r)/(sqrt(rhom_l)+sqrt(rhom_r))
c_ave = (sqrt(rhom_l)*c_l**2 + sqrt(rhom_r)*c_r**2)/(sqrt(rhom_l)+sqrt(rhom_r))&
&+0.5*sqrt(rhom_l*rhom_r)/(sqrt(rhom_l)+sqrt(rhom_r))**2*(um_r-um_l)**2 

c_ave = sqrt(c_ave)

!SL = min(u_ave-c_ave,um_l-c_l)
!SR = max(u_ave+c_ave,um_r+c_r)

SL = min(um_l-c_l,um_r-c_r)
SR = max(um_l+c_l,um_r+c_r)

!SL = um_l-c_l
!SR = um_r+c_r

SM = (p_r-p_l + rhom_l*um_l*(SL-um_l)-rhom_r*um_r*(SR-um_r))/&
    &(rhom_l*(SL-um_l)-rhom_r*(SR-um_r)) 

return
end subroutine appx_Roe_avg_waves
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Sub routine to calculate Fluxes at the interfaces !!!!!!!!!!!!!!!!!!!!!!!
subroutine flux_interface_calculator_first_order(U_left,U_right,F,dt)
implicit none 
double precision, dimension(1:3) :: U_left
double precision, dimension(1:3) :: U_right
double precision, dimension(1:3) :: U_star
double precision, dimension(1:3) :: F
double precision, dimension(1:3) :: F_L
double precision, dimension(1:3) :: F_W
double precision, dimension(1:3) :: F_left
double precision, dimension(1:3) :: F_right
double precision :: l_w_1,l_w_2,l_w_3
double precision :: r_w_1,r_w_2,r_w_3
double precision :: q_nu,ql_nu,qr_nu
double precision :: SL,SR,SM
double precision :: r,dt
! Calculation of the fluxes with the centered values 
call flux_calculator(F_left,U_left)
call flux_calculator(F_right,U_right)
select case(Flux_method)
case(1)
! For the Lax friedrich's method
   F(:) = 1.0/2.0*(F_left(:) + F_right(:)) - &
          &1.0/2.0*(dx/dt)*(U_right(:)-U_left(:))
case(2) 
! For a simple first order upwinding
call wavespeed_calculator(U_left,l_w_1,l_w_2,l_w_3)
call wavespeed_calculator(U_right,r_w_1,r_w_2,r_w_3)
ql_nu = max(dabs(l_w_1),dabs(l_w_2),dabs(l_w_3))
qr_nu = max(dabs(r_w_1),dabs(r_w_2),dabs(r_w_3))
q_nu = max(ql_nu,qr_nu)*dt/(dx)
   F(:) = 1.0/2.0*(F_left(:) + F_right(:)) - &
          &1.0/2.0*(dx/dt*q_nu)*(U_right(:)-U_left(:))
case(3) 
! A HLL scheme without resolving the contact
! Calculation of the average wave speed from Approximated Roe average
call appx_Roe_avg_waves(U_left,U_right,SL,SM,SR)
if(SL>=0)then 
  F(:) = F_left(:)
elseif(SR<=0)then
  F(:) = F_right(:)
elseif(SL<0.and.SR>0)then 
  F(:) = (SR*F_left(:) - SL*F_right(:) + SL*SR*(U_right(:)-U_left(:)))/&
         &(SR-SL)
end if
case(4)
! HLLC Scheme from TORO's work to preserve the contact discontinuity 
call wavespeed_calculator(U_left,l_w_1,l_w_2,l_w_3)
call wavespeed_calculator(U_right,r_w_1,r_w_2,r_w_3)
call appx_Roe_avg_waves(U_left,U_right,SL,SM,SR) 
! Using the standard model to capture the wave span from HLL(approximated Roe
! averages gives overshoot in alpha
!SL = min(l_w_1,r_w_1)
!SR = max(r_w_3,l_w_3)
if(SL>=0)then 
  F(:) = F_left(:)
elseif(SR<=0)then
  F(:) = F_right(:)
elseif(SL<0.and.SR>0)then 
  if(SM>=0)then 
  U_star(1) = U_left(1)*(SL-U_left(3)/U_left(1))/(SL-SM)    
  U_star(2) = U_star(1)*U_left(2)/U_left(1)
  U_star(3) = U_star(1)*SM
  ! Getting the fluxes 
  F(:) = F_left(:) + SL*(U_star(:)-U_left(:)) 
  elseif(SM<0)then
  U_star(1) = U_right(1)*(SR-U_right(3)/U_right(1))/(SR-SM)    
  U_star(2) = U_star(1)*U_right(2)/U_right(1)
  U_star(3) = U_star(1)*SM
  ! Getting the fluxes 
  F(:) = F_right(:) + SR*(U_star(:)-U_right(:)) 
  endif
endif
end select
return
end subroutine flux_interface_calculator_first_order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine flux_interface_calculator_second_order(U_left,U_right,F,dt)
implicit none
double precision, dimension(1:3) :: U_left
double precision, dimension(1:3) :: U_right
double precision, dimension(1:3) :: U_star
double precision, dimension(1:3) :: F
double precision, dimension(1:3) :: F_left
double precision, dimension(1:3) :: F_right
double precision :: dt
! Calculation of the fluxes with the centered values
call flux_calculator(F_left,U_left)
call flux_calculator(F_right,U_right)
! Two step Lax Wendroff method
! First step computation of a Uj+1/2
  U_star(:) = 1.0/2.0*(U_left(:)+U_right(:)) - &
&1.0/2.0*(dt/dx)*(F_right(:)-F_left(:))
  call flux_calculator(F,U_star)
endsubroutine flux_interface_calculator_second_order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Computation of the coefficient to determine the flux limiting!!!!!!!!!!!!
subroutine omega_flux_cal(U2_left,U1_left,U1_right,U2_right,omega_val)
implicit none
double precision, dimension(1:3) :: U2_left
double precision, dimension(1:3) :: U1_left
double precision, dimension(1:3) :: U1_right
double precision, dimension(1:3) :: U2_right
double precision, dimension(1:3) :: delta_UR
double precision, dimension(1:3) :: delta_UL
double precision, dimension(1:3) :: delta_UM
double precision, dimension(1:3) :: omega_val
double precision :: p2_left,p1_left,p1_right,p2_right
double precision :: alpha2_left,alpha1_left,alpha1_right,alpha2_right
double precision :: Y2_left,Y1_left,Y1_right,Y2_right
double precision :: um2_left,um1_left,um1_right,um2_right
double precision :: delta_pl,delta_pr,delta_pm
double precision :: delta_alphal,delta_alphar,delta_alpham
double precision :: delta_Yl,delta_Yr,delta_Ym
double precision :: delta_uml,delta_umr,delta_umm
double precision :: omega_scalar_p,omega_scalar_Y,omega_scalar
double precision :: omega_scalar_alpha
double precision :: q_plus,w_m_l,rho_l,rho_g
integer :: j
q_plus = 0.33 !(Taking from the scalar case of having it to be for FUP method) 

! GETTING THE PRESSRE AND OBTAINING ITS GRADIENTS
call pressure_from_U(p2_left,rho_l,rho_g,w_m_l,U2_left)
call pressure_from_U(p1_left,rho_l,rho_g,w_m_l,U1_left)
call pressure_from_U(p1_right,rho_l,rho_g,w_m_l,U1_right)
call pressure_from_U(p2_right,rho_l,rho_g,w_m_l,U2_right)
! GETTING THE VOLUME FRACTION AND ITS GRADIENTS 
call alpha_from_U(alpha2_left,U2_left)
call alpha_from_U(alpha1_left,U1_left)
call alpha_from_U(alpha1_right,U1_right)
call alpha_from_U(alpha2_right,U2_right)

delta_alphal = alpha1_left-alpha2_left
delta_alpham = alpha1_right-alpha1_left
delta_alphar = alpha2_right-alpha1_right

   if(delta_alphal*delta_alpham >=0 .and. delta_alphar*delta_alpham>=0)then
   omega_scalar_alpha = &
&  min(1.0,(2.0/q_plus*abs(delta_alphar/delta_alpham)),(2.0/q_plus*abs(delta_alphal/delta_alpham)))
   else
   omega_scalar_alpha = 0.0
   endif

delta_pl = p1_left-p2_left
delta_pm = p1_right-p1_left
delta_pr = p2_right-p1_right
   if(delta_pl*delta_pm >=0 .and. delta_pr*delta_pm>=0)then
   omega_scalar_p = &
&  min(1.0,(2.0/q_plus*abs(delta_pr/delta_pm)),(2.0/q_plus*abs(delta_pl/delta_pm)))
   else
   omega_scalar_p = 0.0
   endif

! GETTING THE MASS FRACTION AND OBTAINING ITS GRADIENTS
Y2_left = U2_left(2)/U2_left(1)
Y1_left = U1_left(2)/U1_left(1)
Y1_right = U1_right(2)/U1_right(1)
Y2_right = U2_right(2)/U2_right(1)

delta_Yl = Y1_left - Y2_left
delta_Ym = Y1_right - Y1_left
delta_Yr = Y2_right -Y1_right
   if(delta_Yl*delta_Ym >=0 .and. delta_Yr*delta_Ym>=0)then
   omega_scalar_Y = &
&  min(1.0,(2.0/q_plus*abs(delta_Yr/delta_Ym)),(2.0/q_plus*abs(delta_Yl/delta_Ym)))
   else
   omega_scalar_Y = 0.0
   endif
   
! LIMITING THE CONSERVATIVE VARIABLES 
delta_UL = U1_left-U2_left
delta_UR = U2_right-U1_right
delta_UM = U1_right-U1_left
do j = 1,3
   if(delta_UL(j)*delta_UM(j) >=0 .and. delta_UR(j)*delta_UM(j)>=0)then
   omega_val(j) = &
&  min(1.0,(2.0/q_plus*abs(delta_UR(j)/delta_UM(j))),(2.0/q_plus*abs(delta_UL(j)/delta_UM(j))))
   else
   omega_val(j) = 0.0
   endif
end do
omega_scalar = &
& min(omega_val(1),omega_val(2),omega_val(3))
!omega_scalar = min(omega_val(1),omega_scalar_alpha,omega_scalar_p)
!omega_scalar = 0.0
omega_val(1) = omega_scalar_alpha
omega_val(2) = omega_scalar_alpha
omega_val(3) = omega_scalar_alpha
return
end subroutine omega_flux_cal
!!!!!!!!!!!!!!! Computation of the source vector!!!!!!!!!!!!!!!!!!!!!!
subroutine source_vector(U_vector,S_vector)
implicit none
double precision, dimension(1:3) :: U_vector
double precision, dimension(1:3) :: S_vector
double precision :: pf ! the friction factor for the wall
double precision :: f ,rhom,um 
! Getting the details for the wall friction 
f = 0.05
rhom = U_vector(1)   
um = U_vector(3)/rhom
selectcase(wall_friction_switch)
case(0)
pf = 0
case(1)
pf = f*rhom*um*dabs(um)/sqrt(Area/pi)
end select

! For the time being no source in the propduction of gas and liquid from the
! formation is considered 
S_vector(1) = 0.0
S_vector(2) = 0.0
! For the initial computation no pipe friction is considered 
S_vector(3) = -U_vector(1)*g*cos(theta) - pf
return
end subroutine source_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine boundary_condition(U_vector,t)
implicit none
double precision, dimension(1:3,0:xcells+1) :: U_vector
double precision ::  t
double precision :: exit_pressure, BHP,rho_g_exit,rho_l_exit
double precision :: rho_l_inlet,rho_g_inlet,alpha_injection
double precision :: alpha_exit,w_m_l
! Defining the boundary condition based on the test cases
selectcase(test_case)
case(1)
     call soft_boundary(U_vector(:,1),U_vector(:,0))
     call soft_boundary(U_vector(:,xcells),U_vector(:,xcells+1))
case(2)
     call soft_boundary(U_vector(:,1),U_vector(:,0))
     call soft_boundary(U_vector(:,xcells),U_vector(:,xcells+1))
case(3)
     call soft_boundary(U_vector(:,1),U_vector(:,0))
     call soft_boundary(U_vector(:,xcells),U_vector(:,xcells+1))
case(4)
     call soft_boundary(U_vector(:,1),U_vector(:,0))
     call soft_boundary(U_vector(:,xcells),U_vector(:,xcells+1))
case(5) ! This is the case where we inject water through the system at a
! particular velocity (un pressurized liquid )
if(t <=gas_kick_time)then
     U_vector(1,0) = U_vector(1,0)
     U_vector(2,0) = 0.0
     U_vector(3,0) = U_vector(1,0)*flow_velocity_gas_kick
     
     U_vector(1,xcells+1) = U_vector(1,xcells+1)
     U_vector(2,xcells+1) = 0.0
     U_vector(3,xcells+1) = U_vector(1,xcells+1)*flow_velocity_gas_kick
elseif(t > gas_kick_time.and.t<=shut_in_time)then ! Injection of gas
call pressure_from_U(BHP,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,0)) 
     alpha_injection = 0.4
!     rho_l_inlet = rho_liquid_ref + (BHP-p_ref)/(a_liquid**2)
!     rho_g_inlet = BHP/(a_gas**2) 
     print*,BHP,w_m_l
     U_vector(1,0) = (1-alpha_injection)*rho_l_inlet +&
                     &alpha_injection*rho_g_inlet
     U_vector(2,0) = alpha_injection*rho_g_inlet +&
                   & (1.0-alpha_injection)*w_m_l*rho_l_inlet
     U_vector(3,0) = 1815.77
     
! Conditions at the exit having constant ambient pressure
!     call hard_boundary(U_vector(:,xcells),U_vector(:,xcells+1))
exit_pressure = p_ref*10
rho_g_exit = exit_pressure/(a_gas**2)
rho_l_exit = rho_liquid_ref +(exit_pressure-p_ref)/(a_liquid**2)  
call alpha_from_U(alpha_exit,U_vector(:,xcells)) 
     U_vector(1,xcells+1) = (1-alpha_exit)*rho_l_exit + alpha_exit*rho_g_exit
     U_vector(2,xcells+1) = alpha_exit*rho_g_exit
     U_vector(3,xcells+1) = U_vector(3,xcells)  
elseif(t>shut_in_time)then
     call hard_boundary(U_vector(:,1),U_vector(:,0))
     call hard_boundary(U_vector(:,xcells),U_vector(:,xcells+1))
end if
!print*," I am inside this loop",t

case(6) ! In this case after initial distribution of 
     call hard_boundary(U_vector(:,1),U_vector(:,0))
     call hard_boundary(U_vector(:,xcells),U_vector(:,xcells+1))
end select
return
end subroutine boundary_condition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine flux_boundary_condition(F_1,F_2,U_vector,x,t)
implicit none
double precision, dimension(1:3,0:xcells+1) :: U_vector
double precision, dimension(1:3,0:xcells) :: F_1
double precision, dimension(1:3,0:xcells) :: F_2
double precision, dimension(0:xcells) :: x
double precision :: u_inlet,v_inlet
double precision :: alpha_exit,alpha_injection,BHP
double precision :: t, p1,p2,w_m_l
double precision :: rho_l_inlet,rho_g_inlet
select case(test_case)
case(1)
! Giving the soft boundary conditions
  call flux_calculator(F_1(:,0),U_vector(:,1))
  call flux_calculator(F_1(:,xcells),U_vector(:,xcells))
case(2)
! Giving the soft boundary conditions
  call flux_calculator(F_1(:,0),U_vector(:,1))
  call flux_calculator(F_1(:,xcells),U_vector(:,xcells))
case(3)
! Giving the soft boundary conditions
  call flux_calculator(F_1(:,0),U_vector(:,1))
  call flux_calculator(F_1(:,xcells),U_vector(:,xcells))
case(4)
! Giving the soft boundary conditions
  call flux_calculator(F_1(:,0),U_vector(:,1))
  call flux_calculator(F_1(:,xcells),U_vector(:,xcells))
case(5) ! Case of having an initial water column and then injecting the gas
! STANDING WATER COLUMN
if(t>shut_in_time)then
! At the bottom
  F_1(1,0) = 0.0
  F_1(2,0) = 0.0 
  call pressure_from_U(p1,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,1))
  call pressure_from_U(p2,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,2))
! Interpolating the pressure to get the as the momentum is zero
  F_1(3,0) = (3.0*p1 - p2)/2.0
! Getting the conditions at the top of the surface
  F_1(1,xcells) = 0.0
  F_1(2,xcells) = 0.0
  call pressure_from_U(p1,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,xcells))
  call pressure_from_U(p2,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,xcells-1))
  F_1(3,xcells) = (3.0*p1-p2)/2.0
! INJECTION OF GAS-LIQUID MIXTURE FROM THE BOTTOM
elseif(t > gas_kick_time.and.t<=shut_in_time)then ! Injection of gas
!call pressure_from_U(p1,U_vector(:,1)) 
!call pressure_from_U(p2,U_vector(:,2)) 
  call pressure_from_U(p2,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,2))
  call pressure_from_U(p1,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,1))
!  call pressure_from_U(p2,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,2))
alpha_injection = 0.4
BHP = (3.0*p1 - p2)/2.0
     rho_l_inlet = rho_liquid_ref + (BHP-p_ref)/(a_liquid**2)
     rho_g_inlet = BHP/(a_gas**2) 
     F_1(2,0) = 99.09941 
     F_1(1,0) = 1815.77 
     v_inlet = F_1(2,0)/(alpha_injection*rho_g_inlet)
     u_inlet = (F_1(1,0)-F_1(2,0))/((1.0-alpha_injection)*rho_l_inlet)
     F_1(3,0) =  rho_l_inlet*(1-alpha_injection)*u_inlet**2 + &
                &rho_g_inlet*alpha_injection*v_inlet**2 + BHP
     F_2(:,0) = F_1(:,0)
! Exit boundary condition - Maintaining the exit pressure of 10 Bar
! and is obtained from the ghost cell values itself
    endif
case(6) ! Sedimentation case
! At the bottom
  F_1(1,0) = 0.0
  F_1(2,0) = 0.0 
  call pressure_from_U(p1,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,1))
  call pressure_from_U(p2,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,2))
! Interpolating the pressure to get the as the momentum is zero
  F_1(3,0) = (3.0*p1 - p2)/2.0
! Getting the conditions at the top of the surface
  F_1(1,xcells) = 0.0
  F_1(2,xcells) = 0.0
  call pressure_from_U(p1,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,xcells))
  call pressure_from_U(p2,rho_l_inlet,rho_g_inlet,w_m_l,U_vector(:,xcells-1))
  F_1(3,xcells) = (3.0*p1-p2)/2.0
endselect
F_2(:,0) = F_1(:,0)
F_2(:,xcells) = F_1(:,xcells)
return
endsubroutine flux_boundary_condition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Ghost cell Soft boundary condition
subroutine soft_boundary(U_inner,U_ghost)
implicit none 
double precision, dimension(1:3) :: U_inner
double precision, dimension(1:3) :: U_ghost
U_ghost(:) = U_inner(:)
return
end subroutine soft_boundary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ghost cell Hard boundary condition
subroutine hard_boundary(U_inner,U_ghost)
implicit none 
double precision, dimension(1:3) :: U_inner
double precision, dimension(1:3) :: U_ghost
U_ghost(1:2) = U_inner(1:2)
U_ghost(3) = -U_inner(3)
return
end subroutine hard_boundary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Computing the variables from the U vector!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine U_to_variables(U_vector,rho_m,alpha,rho_l,rho_g,u_m,u,v,p,w_m_l)
implicit none 
double precision, dimension(1:3) :: U_vector
double precision :: rho_m,alpha,rho_l,rho_g,u_m,u,v,p
double precision :: v_minus_u,C,v_s,w_m_l
rho_m = U_vector(1)
u_m = U_vector(3)/U_vector(1)
call drift_flux_law(v_minus_u,C,v_s,U_vector)
call pressure_from_U(p,rho_l,rho_g,w_m_l,U_vector)
call alpha_from_U(alpha,U_vector)
!rho_g = p*rho_gas_ref/p_ref
!if(w_m_l == 0.0)then
!rho_l = (p-p_ref)/(a_liquid**2) + rho_liquid_ref
!else
!rho_l = (U_vector(2) - alpha*rho_g)/(1-alpha)*w_m_l
!endif

u = (U_vector(3) - alpha*rho_g*v_minus_u)/(U_vector(1))
v = u + v_minus_u
return
end subroutine U_to_variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to writeout the data on the edges of the cell rather than the cell centroid
subroutine writeout(x,U_vector,iter)
implicit none
double precision,dimension(0:xcells+1) :: x
double precision,dimension(1:3,0:xcells+1):: U_vector
double precision :: rho_m,rho_l,rho_g,alpha,u_m,u,v,p,w_m_l
character*20 :: filname,fil2
integer :: nn,iter,j,imax
nn = iter/frequency
imax = xcells
filname ="out"
write(filname(4:8),"(i5.5)"),nn
filname(9:12) =".dat"
! Calculating the cell centered value of the quatities (Writing in binary
! format for efficient storage
open(unit = 1,file = filname)
write(1,*)'Variables = &
&x,rho_m,alpha,rho_l,rho_g,u_m,u,v,p,W_l,W_g,W_mix,W_mom,Y_g,Y_l'
write(1,*)'Zone T =','"',nn,'"','I = ',imax,',','F = point'
do j = 1,xcells
call U_to_variables(U_vector(:,j),rho_m,alpha,rho_l,rho_g,u_m,u,v,p,w_m_l) 
 write(1,*)x(j),rho_m,alpha,rho_l,rho_g,u_m,u,v,p,rho_l*(1-alpha),rho_g*alpha,&
           &rho_m,rho_m*u_m,alpha*rho_g/rho_m,1.0-alpha*rho_g/rho_m
end do
close(1)
end subroutine writeout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writeout_binary(x,U_vector,iter)
implicit none 
double precision,dimension(0:xcells+1) :: x
double precision,dimension(1:3,0:xcells+1):: U_vector
double precision,dimension(1:xcells) :: rho_m
double precision,dimension(1:xcells) :: rho_l
double precision,dimension(1:xcells) :: rho_g
double precision,dimension(1:xcells) :: alpha
double precision,dimension(1:xcells) :: u_m
double precision,dimension(1:xcells) :: u_v
double precision,dimension(1:xcells) :: u
double precision,dimension(1:xcells) :: v
double precision,dimension(1:xcells) :: p
double precision,dimension(1:xcells) :: W_l
double precision,dimension(1:xcells) :: W_g
double precision,dimension(1:xcells) :: W_mix
double precision,dimension(1:xcells) :: W_mom
double precision,dimension(1:xcells) :: Y_g
double precision,dimension(1:xcells) :: Y_l
double precision,dimension(1:xcells) :: w_m_l
integer :: i,j,imax,jmax,kmax
integer :: debug,ier,itot
integer :: tecini,tecdat,teczne,tecnod,tecfil,tecend
integer :: visdouble,disdouble
character*1 :: nulchar
character*12 :: filname,filname1
integer :: nn,iter
nn = iter/frequency
filname ="out"
write(filname(4:8),"(i5.5)"),nn
filname(9:12) =".plt"


nulchar = char(0)
debug = 0
visdouble = 0
disdouble = 1
imax = xcells
jmax = 1
kmax = 1

! Extracting the date from the conserved variables
do j = 1,xcells
call U_to_variables(U_vector(:,j),rho_m(j),alpha(j),rho_l(j),rho_g(j),u_m(j),u(j),v(j),p(j),w_m_l(j))
W_l(j) =  rho_l(j)*(1.0-alpha(j))*(1.0-w_m_l(j))
W_g(j) =  rho_g(j)*alpha(j)+(1.0-alpha(j))*w_m_l(j)*rho_l(j)
W_mix(j) = rho_m(j)
W_mom(j) = rho_m(j)*u_m(j)
Y_g(j) = alpha(j)*rho_g(j)/rho_m(j)
Y_l(j) = 1.0-Y_g(j)
u_v(j) = u(j)*(1.0-alpha(j)) + v(j)*alpha(j)
end do


! Open the file and write the Tecplot datafile header information
ier = &
&tecini('Variables'//nulchar,&
&'x,<greek>r</greek><sub>m</sub>,<greek>a</greek>,<greek>r</greek><sub>l</sub>,&
&<greek>r</greek><sub>g</sub>,u<sub>M</sub>,u<sub>V</sub>,u,v,p,W<sub>l</sub>,&
&W<sub>g</sub>,W<sub>mix</sub>,W<sub>mom</sub>,Y<sub>g</sub>,Y<sub>l</sub>,<greek>w</greek><sub>ml</sub>'//nulchar,&
&filname//nulchar,'.'//nulchar,debug,visdouble)

!
! Write the zone header information
!
ier = teczne('Variables'//nulchar,imax,jmax,kmax,'BLOCK'//nulchar,nulchar)

! Writing out the field data 
itot = imax*jmax*kmax
ier = tecdat(itot,x(1:xcells),disdouble)
ier = tecdat(itot,rho_m,disdouble)
ier = tecdat(itot,alpha,disdouble)
ier = tecdat(itot,rho_l,disdouble)
ier = tecdat(itot,rho_g,disdouble)
ier = tecdat(itot,u_m,disdouble)
ier = tecdat(itot,u_v,disdouble)
ier = tecdat(itot,u,disdouble)
ier = tecdat(itot,v,disdouble)
ier = tecdat(itot,p*0.00014,disdouble)
ier = tecdat(itot,W_l,disdouble)
ier = tecdat(itot,W_g,disdouble)
ier = tecdat(itot,W_mix,disdouble)
ier = tecdat(itot,W_mom,disdouble)
ier = tecdat(itot,Y_g,disdouble)
ier = tecdat(itot,Y_l,disdouble)
ier = tecdat(itot,w_m_l,disdouble)

! to close the file 
ier = tecend()

! Writing the similar file to be written in Matlab

!filname1 ="mat"
!write(filname1(4:8),"(i5.5)"),nn
!filname1(9:12) =".plt"
!open(1,file='filname1',form = 'unformatted', access = 'direct',RECL = 32)
!do j = 1,xcells
!write(1,REC =j)x(j),rho_m(j),alpha(j),rho_l(j),rho_g(j),u_M(j)
!    &u_V(j),u(j),v(j),p(j),W_l(j),W_g(j),W_g(j),W_mix(j),W_mom(j),&
!    &Y_g(j),Y_l(j)
!end do


return
end subroutine writeout_binary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module cellstates
