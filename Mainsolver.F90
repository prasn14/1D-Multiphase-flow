Program Mainsolver
use cellstates
use input
implicit none
! Grid Generation --- x dimension
double precision, dimension(0:xcells+1) :: x
double precision, dimension(1:3,0:xcells+1) :: U_vector ! contains the
double precision, dimension(1:3,0:xcells+1) :: U_star ! contains the
! conservative variables
double precision, dimension(1:3,0:xcells+1) :: S_vector ! contains the
double precision, dimension(1:3,0:xcells) :: Flux_vector_first_order
double precision, dimension(1:3,0:xcells) :: Flux_vector_second_order
double precision, dimension(1:3,0:xcells) :: Flux_vector
double precision, dimension(1:3,0:xcells) :: omega_flux
double precision, dimension(1:3) :: F_dummy
double precision :: current_time,max_wavespeed,rhom_max !(For getting the
! for the forward Euler scheme
double precision :: bhp1,bhp2,rho_l,rho_g
double precision :: w_ml
double precision :: dt
! double precision :: cfl SWITCH THIS ONE WHEN RUNNING FOR CONSTANT TIME STEPS
integer::i,j,iter,n
real ::  time1,time2
call cpu_time(time1)
! Coordinates of Cell centres
do i = 0,xcells+1
   x(i) = (i-0.5)*dx
end do  
current_time = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
iter = 0
! Initial conditions for the specified problem can be given in different
! variables
call U_vector_initial(U_vector,x)
print*,"Writing out the initial data"
                call writeout_binary(x,U_vector,iter)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!print*,'The total timesteps for this iteration is',timesteps
!do n = 1,timesteps
open(unit = 2,file = "BHP.dat")
write(2,*)'Variables = t,p'
call pressure_from_U(bhp1,rho_l,rho_g,w_ml,U_vector(:,0))
call pressure_from_U(bhp2,rho_l,rho_g,w_ml,U_vector(:,1))
write(2,*)current_time,(bhp1+bhp2)*0.00014/2.0
do while(current_time<=time)
        iter = iter + 1
        call cfl_computation(U_vector,max_wavespeed)
        call maximum_from_array(U_vector(1,:),rhom_max)
!        print*,'The iter number is',iter
!        print*,iter,max_wavespeed*dt/dx
! GETTING THE VARIABLE TIME STEP FROM THE CFL NUMBER 
         dt = min(cfl*dx/dabs(max_wavespeed),1.5/(rhom_max*g*cos(theta)))
         dt = cfl*dx/dabs(max_wavespeed)
        current_time = current_time + dt
        print*,iter,dt,current_time,max_wavespeed
! CALLING THE  BOUNDARY CONDITION BASED ON THE PROBLEM
      call boundary_condition(U_vector,current_time)

!Computation of the source terms in the cell centers
        do j = 0,xcells+1
                call source_vector(U_vector(:,j),S_vector(:,j))
        end do

! Computation of the Fluxes at the interfaces given the left and right states
        !(Getting the first order flux either with first order upwinding or Lax
        !Freidrichs flux method) 
!$OMP DO 
        do j = 1,xcells+1
                call flux_interface_calculator_first_order(U_vector(:,j-1),U_vector(:,j),Flux_vector_first_order(:,j-1),dt) 
        end do         

!$OMP END DO 

!$OMP  DO 
        ! Getting the flux based on two step Lax Wendroff method)         
        do j = 1,xcells+1
                call flux_interface_calculator_second_order(U_vector(:,j-1),U_vector(:,j),Flux_vector_second_order(:,j-1),dt) 
        end do         

!$OMP END DO 

! IF PRESCRIBING THE FLUXES AT THE BOUNDARIES THEN THE VARIABLES OF THE FLUX
! FUNTION NEEDS TO BE GIVEN AND THIS IS DONE FOR TEST CASES 5 and 6
selectcase(boundary_variable)
case(0)
! We retain the values of the ghost cell to compute the fluxes
case(1)
call flux_boundary_condition(Flux_vector_first_order,Flux_vector_second_order,U_vector,x,current_time)
end select
! Computation of the the factor to limit the flux based on the gradients
selectcase(flux_order_switch)
case(0)
omega_flux(1:3,:) = 0.0
case(1)
! Getting the first order solution to a star state and then contructing the
! gradients of that solution
        do j = 1,xcells
                U_star(:,j) = U_vector(:,j) - dt/dx*(Flux_vector_first_order(:,j) - Flux_vector_first_order(:,j-1)) &
                                & + dt*S_vector(:,j) 
        end do        
! Getting the boundary value for the star state
     call boundary_condition(U_star,current_time)
     ! Computation of the coefficients to be computed   
     do j = 1,xcells-1
               call omega_flux_cal(U_vector(:,j-1),U_vector(:,j),U_vector(:,j+1),U_vector(:,j+2),omega_flux(:,j))
!               call omega_flux_cal(U_star(:,j-1),U_star(:,j),U_star(:,j+1),U_star(:,j+2),omega_flux(:,j))
     end do
     omega_flux(:,0) = omega_flux(:,1)
     omega_flux(:,xcells) = omega_flux(:,xcells-1) 
case(2)
omega_flux(1:3,:) = 1.0
end select
! Getting the limited flux 
do i = 1,3
   do  j = 0,xcells
   Flux_vector(i,j) = Flux_vector_first_order(i,j) + &
                      &omega_flux(i,j)*(Flux_vector_second_order(i,j) - Flux_vector_first_order(i,j))
   end do
select case(boundary_variable)
case(0)
case(1)
Flux_vector(:,0) = Flux_vector_first_order(:,0)
Flux_vector(:,xcells) = Flux_vector_first_order(:,xcells)
end select

end do  

! Updating the solution in time (for the 1st order time stepping)
!$OMP DO 
        do j = 1,xcells
                U_vector(:,j) = U_vector(:,j) - dt/dx*(Flux_vector(:,j) - Flux_vector(:,j-1)) &
                                & + dt*S_vector(:,j) 
        end do        
!$OMP END DO 
        if(mod(iter,frequency) == 0)then
                call writeout_binary(x,U_vector,iter)
!                call writeout(x,U_vector,iter)
        endif
          

call pressure_from_U(bhp1,rho_l,rho_g,w_ml,U_vector(:,0))
call pressure_from_U(bhp2,rho_l,rho_g,w_ml,U_vector(:,1))
write(2,*)current_time,(bhp1+bhp2)*0.00014/2.0


end do
close(2)
                call writeout_binary(x,U_vector,iter)

call cpu_time(time2)
print*,"The computational time for this test case in sec is ",time2-time1
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program Mainsolver
