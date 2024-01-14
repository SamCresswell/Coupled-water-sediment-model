module adams_bashford
    
    use initialisation
    use variables
    use subcalculations
    
    implicit none
    
    
    contains
        
        !-------------------------------------------------------------------------
        !This subroutine contains the Adams-Bashford method
        !Partial derivatives are found at two different positions
    
        subroutine adams_bashford_method(q,h,h_new,q_new,c,c_b,n_mannings,dh_dt_a,dq_dt_a,dh_dt_b,&
                                        dq_dt_b,c_f_1,h_1,g,c_b_input,c_input,n_mannings_input,h_lim,&
                                        n, n_steps,z_b,dt, d_median, von_karman_constant, z_0s, h_input,&
                                        q_fixed,drive_flow)
        
        
            double precision,dimension(0:n_steps) :: c,c_b,n_mannings
            double precision,dimension(0: n_steps + 1) :: q,h,h_new,q_new,dh_dt_a
            double precision,dimension(1: n_steps) :: dq_dt_a
            double precision,dimension(0:n_steps+1) :: z_b,z_0s, d_median
            double precision :: c_f_1,h_1,g,c_b_input,c_input,n_mannings_input,h_lim,dh_dt_b,dq_dt_b,dt,von_karman_constant, h_input, q_fixed
            integer :: n, n_steps, drive_flow
            
            q(1) = 2.0d0 * q(2) - q(3)
            h(1) = 2.0d0 * h(2) - h(3)
            
            !Drives the flow if applicable
            if (drive_flow == 1) then
                    
                q(1) = q_fixed 
            
            end if
            
            !This loop performs the calculation
            do n = 2, n_steps - 1
                
                !Calculates the change in the flow depth
                dh_dt_b = - (q(n + 1) - q(n - 1))/(2.0d0 * dx)
                
                !Calls the friction coefficient
                !call friction_coefficient_1D(bed_friction_type,n,n_steps,c_f_1,c,c_b,n_mannings,h,g,c_b_input,c_input,n_mannings_input,h_lim,h_f)
                call friction_coefficient(bed_friction_type,n,n_steps,c_f,c,c_b,n_mannings,h,g,c_b_input,&
                                        c_input,n_mannings_input,h_input,h_lim,d_median, von_karman_constant, z_0s)
    
                !Calculates the change in the specific flowrate
                dq_dt_b =    - (q(n + 1) * q(n + 1)/h(n + 1)- q(n - 1)*q(n-1)/h(n - 1) + 0.5d0 * g *(h(n + 1)*h(n + 1) - h(n - 1)*h(n - 1))&
                                +g*h(i)*(z_b(n + 1)- z_b(n - 1)))/(2.0d0 * dx) - c_f_1*q(n)**2.0d0/h(n)**2.0d0
                
                !Calculates the updated depth of flow
                h_new(n) = h(n) + dt * (1.50d0 * dh_dt_b - 0.50d0 * dh_dt_a(n))
                
                !print*,n,h(n),h_new(n)

                !Calculates the updated flowrates
                q_new(n) = q(n) + dt * (1.50d0 * dq_dt_b - 0.50d0 * dq_dt_a(n))
                
                !Sets the old depth of flow to the new one
                dh_dt_a(n) = dh_dt_b
            
                !Sets the old flowrate to the new one
                dq_dt_a(n) = dq_dt_b
                
            end do
            
            !Sets the boundary conditions for the depth of flow
            h_new(1) = h(1)
            h_new(n_steps) = h(n_steps)

            !Sets the boundary conditions for the flowrate                  
            q_new(1) = q_new(2)
            q_new(n_steps) = q_new(n_steps - 1)
            
            !Updates the flow depth, specific flowrate and velocity to their new values
            do n = 1,n_steps
                
                h(n) = h_new(n)
                q(n) = q_new(n)
                u(n) = q(n)/h(n)
                
                !print*,n,h(n),q(n)
                
            end do
            
        end subroutine adams_bashford_method    
    
end module adams_bashford