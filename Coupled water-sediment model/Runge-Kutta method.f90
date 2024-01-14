module Runge_Kutta_calculation
    
    use initialisation
    use variables
    use subcalculations
    
    implicit none
    
    contains
        
        !-------------------------------------------------------------------------
        !This subroutine contains the Runge-Kutta method
        !Partial derivatives are found at four different positions
        subroutine main_calculation(n,n_steps,dx,dt,q,dh_dt_1,ramp_up_type,t,t_ramp,q_ramp_end,q_ramp_0,&
                                    bed_friction_type,c_f,c,c_b,n_mannings,h,g,c_b_input,c_input,n_mannings_input,h_input,h_lim,&
                                    u,rho_w,dq_dt_1,tau_bx,west_bc,east_bc,h_1,h_2,h_3,q_1,q_2,q_3,h_east,h_west,q_east,q_west, &
                                    tau_bx_1,tau_bx_2,tau_bx_3,u_1,u_2,u_3,z_b,i, drive_flow)
        
            implicit none
        
            integer :: n,n_steps,ramp_up_type,bed_friction_type,west_bc,east_bc,i, drive_flow
            double precision,dimension(-1:n_steps + 2) :: q,c_f,c,c_b,n_mannings,h,u,z_b
            double precision,dimension(-1:n_steps + 2) :: dh_dt_1,dq_dt_1,dq_dt_2,dq_dt_3,dq_dt_4,dh_dt_2,dh_dt_3,dh_dt_4, u_1
            double precision,dimension(-1:n_steps + 2) :: tau_bx,tau_bx_1,tau_bx_2,tau_bx_3,rho_w
            double precision :: dx,dt,t,t_ramp,q_ramp_end,q_ramp_0,g,c_b_input,c_input,n_mannings_input,h_input,h_lim
            double precision :: h_1,h_2,h_3,q_1,q_2,q_3,h_east,h_west,q_east,q_west,u_2,u_3
            
            
            !Updates the ramped-up surface stress/inflows
            call ramp_up(n, n_steps,ramp_up_type,hyperconcentration_type, &
                            t,t_ramp,q_ramp_end,q_ramp_0,q_input, rho_w_clear, q_fixed, rho_w_multiplier, &
                            q, h, rho_w)
            
            !Drives the flow if applicable
            if (drive_flow == 1) then
                    
                q(1) = q_fixed 
            
            end if
            
            !The first part of the Runge-Kutta calculation
            do n = 1, n_steps
                
                !Calculates the friction coefficient
                call friction_coefficient(bed_friction_type,n,n_steps,c_f,c,c_b,n_mannings,h,g,c_b_input,&
                                        c_input,n_mannings_input,h_input,h_lim,d_median, von_karman_constant, z_0s)
                
                !Calculates the difference in head in each cell
                dh_dt_1(n) = -(q(n + 1)- q(n-1))/(2*dx)
                
                !Calculates the velocity in each cell
                u_current = q(n)/h(n)
                
                !Calculates the friction coefficient
                tau_bx(n) = rho_w(n) * abs(u_current(n)) * u_current(n) * c_f(n)
                
                !Calculates the derivative of flowrate with respect to time
                dq_dt_1(n) = -(q(n + 1)*q(n + 1)/h(n + 1) + g*h(n+1)*h(n+1)/2.0d0 - q(n - 1)*q(n - 1)/h(n - 1) - g*h(n - 1)*h(n + 1)/2.0d0)/2.0d0/dx & 
                        - g * h(n)*(z_b(n + 1) - z_b(n - 1))/2.0d0/dx - tau_bx(n)/rho_w(n)
                
            end do
            
            !Handles the western boundary being a uniform flow
            if (west_bc == 5) then
                
                dh_dt_1(0) = 0.0d0
                dq_dt_1(0) = dq_dt_1(1)
                
            !Handles the western boundary being the transmissive condition
            else if (west_bc == 2) then
                
                dh_dt_1(0) = dh_dt_1(1)
                dh_dt_1(-1) = dh_dt_1(0)
                dq_dt_1(0) = dq_dt_1(1)
                dq_dt_1(-1) = dq_dt_1(0)
                
                else
            
                end if
                
            !Handles the eastern boundary being a uniform flow
            if (east_bc == 5) then
            
                dh_dt_1(n_steps + 1) = 0.0d0
                dq_dt_1(n_steps + 1) = dq_dt_1(n_steps)
            
            !Handles the eastern boundary being the transmissive condition
            else if (east_bc == 2) then
                
                dh_dt_1(n_steps + 1) = dh_dt_1(n_steps)
                dh_dt_1(n_steps + 2) = dh_dt_1(n_steps + 1)
                dq_dt_1(n_steps + 1) = dq_dt_1(n_steps)
                dq_dt_1(n_steps + 2) = dq_dt_1(n_steps + 1)
                
            else
                
            end if
            
            !The second iteration over all spatial steps in the Runge-Kutta method    
            do n = 1, n_steps
                
                h_1 = h(n) + dt * dh_dt_1(n)/2.0d0
                h_east = h(n + 1) + dt * dh_dt_1(n + 1)/2.0d0
                h_west = h(n - 1) + dt * dh_dt_1(n - 1)/2.0d0
                
                q_1 = q(n) + dt * dq_dt_1(n)/2.0d0
                q_east = q(n + 1) + dt * dq_dt_1(n + 1)/2.0d0
                q_west = q(n - 1) + dt * dq_dt_1(n - 1)/2.0d0
                
                dh_dt_2(n) = -(q_east-q_west)/(dx * 2.0d0)
                
                call friction_coefficient_1D(bed_friction_type,n,n_steps,c_f_1,c,c_b,n_mannings,h,g,c_b_input,c_input,n_mannings_input,h_lim,h_f)
                
                !The flow velocity is calculated
                u_1(n) = q_1/h_1
                
                !The bed friction coefficient is calculated
                tau_bx_1(n) = rho_w(n) * c_f_1 * u_1(n) * abs(u_1(n))
                
                !Calculates the derivative of flowrate with respect to time
                dq_dt_2(n) = -(q_east**2/h_east - q_west**2/h_west + g*(h_east**2.0d0 - h_west**2.0d0)/2.0d0)/(2.0d0 * dx) & 
                        - g * h_1*(z_b(n + 1) - z_b(n - 1))/(2.0d0 * dx) - tau_bx_1(n)/rho_w(n)
                
            end do
            
            !Handles the western boundary being a uniform flow
            if (west_bc == 5) then
                
                dh_dt_2(0) = 0.0d0
                dq_dt_2(0) = dq_dt_2(1)
              
            !Handles the western boundary being the transmissive condition
            else if (west_bc == 2) then
                
                dh_dt_2(0) = dh_dt_2(1)
                dh_dt_2(-1) = dh_dt_2(0)
                dq_dt_2(0) = dq_dt_2(1)
                dq_dt_2(-1) = dq_dt_2(0)
                
            else
            
                end if
            
            !Handles the eastern boundary being a uniform flow
            if (east_bc == 5) then
                
                dh_dt_2(n_steps + 1) = 0.0d0
                dq_dt_2(n_steps + 1) = dq_dt_2(n_steps)
        
            !Handles the eastern boundary being the transmissive condition
            else if (east_bc == 2) then
                
                dh_dt_2(n_steps + 1) = dh_dt_2(n_steps)
                dh_dt_2(n_steps + 2) = dh_dt_2(n_steps + 1)
                dq_dt_2(n_steps + 1) = dq_dt_2(n_steps)
                dq_dt_2(n_steps + 2) = dq_dt_2(n_steps + 1)
                
            else
                
                end if
            
            !The third iteration over all spatial steps in the Runge-Kutta method    
            do n = 1, n_steps
                
                h_2 = h(n) + dt * dh_dt_2(n)/2.0d0
                h_east = h(n + 1) + dt * dh_dt_2(n + 1)/2.0d0
                h_west = h(n - 1) + dt * dh_dt_2(n - 1)/2.0d0
                
                q_2 = q(n) + dt * dq_dt_2(n)/2.0d0
                q_east = q(n + 1) + dt * dq_dt_2(n + 1)/2.0d0
                q_west = q(n - 1) + dt * dq_dt_2(n - 1)/2.0d0
                
                dh_dt_3(n) = -(q_east-q_west)/(2.0d0 * dx)
                
                call friction_coefficient_1D(bed_friction_type,n,n_steps,c_f_1,c,c_b,n_mannings,h,g,c_b_input,c_input,n_mannings_input,h_lim,h_f)
                !The flow velocity is calculated
                u_2 = q_2/h_2
                
                !The bed friction coefficient is calculated
                tau_bx_2(n) = rho_w(n) * c_f_1 * u_2 * abs(u_2)
                
                !Calculates the derivative of flowrate with respect to time
                dq_dt_3(n) = -(q_east**2/h_east - q_west**2/h_west + g/2.0d0*(h_east**2.0d0 - h_west**2.0d0))/(2.0d0 * dx) & 
                        - g * h_2*(z_b(n + 1) - z_b(n - 1))/(2.0d0 * dx) - tau_bx_2(n)/rho_w(n)
                
            end do
            
            !Handles the western boundary being a uniform flow
            if (west_bc == 5) then
                
                dh_dt_3(0) = 0.0d0
                dq_dt_3(0) = dq_dt_3(1)
                
            !Handles the western boundary being the transmissive condition
            else if (west_bc == 2) then
                
                dh_dt_3(0) = dh_dt_3(1)
                dh_dt_3(-1) = dh_dt_3(0)
                dq_dt_3(0) = dq_dt_3(1)
                dq_dt_3(-1) = dq_dt_3(0)
                
            else
            
                end if
            
            !Handles the eastern boundary being a uniform flow
            if (east_bc == 5) then
                
                dh_dt_3(n_steps + 1) = 0.0d0
                dq_dt_3(n_steps + 1) = dq_dt_3(n_steps)
               
            !Handles the eastern boundary being the transmissive condition
            else if (east_bc == 2) then
                
                dh_dt_3(n_steps + 1) = dh_dt_3(n_steps)
                dh_dt_3(n_steps + 2) = dh_dt_3(n_steps + 1)
                dq_dt_3(n_steps + 1) = dq_dt_3(n_steps)
                dq_dt_3(n_steps + 2) = dq_dt_3(n_steps + 1)
                    
            else
                
                end if
            
            !The fourth iteration over all spatial steps in the Runge-Kutta method    
            do n = 1, n_steps
                
                h_3 = h(n) + dt * dh_dt_3(n)
                h_east = h(n + 1) + dt * dh_dt_3(n + 1)
                h_west = h(n - 1) + dt * dh_dt_3(n - 1)
                
                q_3 = q(n) + dt * dq_dt_3(n)
                q_east = q(n + 1) + dt * dq_dt_3(n + 1)
                q_west = q(n - 1) + dt * dq_dt_3(n - 1)
                
                dh_dt_4(n) = -(q_east-q_west)/(2.0d0 * dx)
                
                call friction_coefficient_1D(bed_friction_type,n,n_steps,c_f_1,c,c_b,n_mannings,h,g,c_b_input,c_input,n_mannings_input,h_lim,h_f)
                
                !The flow velocity is calculated
                u_3 = q_3/h_3
                
                !The bed friction coefficient is calculated
                tau_bx_3(n) = rho_w(n) * c_f_1 * u_3 * abs(u_3)
                
                !Calculates the derivative of flowrate with respect to time
                dq_dt_4(n) = -(q_east**2/h_east - q_west**2/h_west + g/2.0d0*(h_east**2.0d0 - h_west**2.0d0))/(2.0d0 * dx) & 
                        - g * h_3*(z_b(n + 1) - z_b(n - 1))/(2.0d0 * dx) - tau_bx_3(n)/rho_w(n)
                
                !Now the actual Runge-Kutta method is deployed, using these four coefficients
                h_new(n) = h(n) + dt/6.0d0*(dh_dt_1(n) + 2.0d0 * dh_dt_2(n) + 2.0d0 * dh_dt_3(n) + dh_dt_4(n))
                q_new(n) = q(n) + dt/6.0d0*(dq_dt_1(n) + 2.0d0 * dq_dt_2(n) + 2.0d0 * dq_dt_3(n) + dq_dt_4(n))
                
                !If the western boundary is a uniform flow
                if (west_bc == 1) then
            
                    h_new(0) = h_new(1)
                    h_new(-1) = h_new(2)
                    q_new(0) = - q_new(1)
                    q_new(-1) = - q_new(2)
                
                !If the upstream boundary is a transmissive condition:
                else if (west_bc == 2) then
            
                    h_new(0) = 2.0d0 * h_new(1) - h_new(2)
                    h_new(-1) = 2.0d0 * h_new(0) - h_new(1)
                    q_new(0) = q_input
                    q_new(-1) = q_new(0)
                
                !If the upstream boundary has a prescribed depth and tranamissive flow
                else if (west_bc == 3) then
                
                    h_new(1) = h(1)
                    h_new(0) = h(0)
                    h_new(-1) = h_new(-1)
                    q_new(0) = 2.0d0 * q_new(1) - q_new(2)
                    q_new(-1) = 2.0d0 * q_new(0) - q_new(1)
            
                else
                    print*,'Invalid upstream boundary condition!'
            
                    call exit()
            
                end if
            
                !Handles the eastern (downstream) boundary condition using one of two
                !If the upstream boundary is a transmissive condition:
                if (east_bc == 2) then
                
                    h_new(n_steps + 1) = 2.0d0 * h_new(n_steps) - h_new(n_steps - 1)
                    h_new(n_steps + 2) = 2.0d0 * h_new(n_steps + 1) - h_new(n_steps)
                    q_new(n_steps + 1) = q_new(n_steps)
                    q_new(n_steps + 2) = q_new(n_steps + 1)
            
                !If the upstream boundary has uniform flow
                else if (east_bc == 5) then
                
                    h_new(n_steps + 1) = h(n_steps + 1)
                    h_new(n_steps + 2) = h(n_steps + 2)
                    q_new(n_steps + 1) = q_new(n_steps)
                    q_new(n_steps + 2) = q_new(n_steps + 1)
            
                else
                    print*,'Invalid upstream boundary condition!'
            
                    call exit()
            
                end if
                
            end do

            !Updates the cells to their new values
            do n = -1,n_steps + 2
                
                h(n) = h_new(n)
                q(n) = q_new(n)
                u(n) = q(n)/h(n)
              
            end do
            
        end subroutine main_calculation
                                    
end module Runge_kutta_calculation