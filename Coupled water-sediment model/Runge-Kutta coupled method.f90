module Runge_Kutta_coupled_calculation
    
    use initialisation
    use variables
    use subcalculations
    
    implicit none
    
    contains
        
        !-------------------------------------------------------------------------
        !This subroutine contains the Runge-Kutta method
        !Partial derivatives are found at four different positions
        subroutine runge_kutta_coupled_calc(n,n_steps, west_bc,east_bc,&
                                            ramp_up_type,bed_friction_type, drive_flow, d_median_equation, use_entrainment_equation, use_deposition_equation, &
                                            upwinding_type, is_bed_fixed, bedload_method, hyperconcentration_type, update_grain_size, end_program_at_mobilisation, &
                                            h,u,z_b,conc, free_surface_elevation, x, S_b, bedrock_elevation, bed_slope_angle, u_1, a_1,a_2,a_3, &
                                            da1_dt_1,da1_dt_2,da1_dt_3, da1_dt_4, da2_dt_1, da2_dt_2,da2_dt_3,da2_dt_4, da3_dt_1,da3_dt_2,da3_dt_3,da3_dt_4, &
                                            dzb_dt_1,dzb_dt_2,dzb_dt_3, dzb_dt_4, a_1_new,a_2_new,a_3_new,z_b_new, rho, s, rho_w, rho_0, bed_porosity, nu, &
                                            tau_bx,tau_bx_1,tau_bx_2,tau_bx_3, c_f, c, c_b, n_mannings, z_0s, d_median, w_s_0, w_s_h, m_coefficient, d_dimensionless, &
                                            dimensionless_critical_shields_parameter, shields_parameter, modified_critical_shields_parameter, &
                                            Fr, Re_p, R_n, q_b, dq_b_dt, u_b, q_b_0, D, En, En_1, En_2, En_3, En_4, D_1, D_2, D_3, D_4, minmod_ab,minmod_a,minmod_b, &
                                            z_b_initial, d_median_initial, d_median_new, erosion_depth, deposition_depth, &
                                            a_1_1,a_1_east,a_1_west,a_2_1,a_2_east,a_2_west,a_3_1,a_3_east,a_3_west,z_b_east,z_b_west, h_1,h_east,h_west, &
                                            g, rho_s, von_karman_constant, rho_w_clear, angle_of_repose, dx, dt, &
                                            h_input, q_input, coeff_1_input,coeff_2_input, nu_input, c_b_input,c_input,n_mannings_input, &
                                            D_input,En_input, d_median_input, hyperconcentration_conc, h_lim, coeff_1,coeff_2, coeff_3, t, uniform_flow_time, &
                                            total_sediment_mass, bed_sediment_area, critical_shields_parameter,alpha_d,alpha_en, u_shear, m_d, &
                                            t_ramp,q_ramp_end,q_ramp_0, q_fixed, rho_w_multiplier, simulation_batch_name)
        
            implicit none
            
            integer :: n,n_steps !Mesh related integers
            integer :: west_bc,east_bc !Boundary values
            integer :: ramp_up_type,bed_friction_type, drive_flow, d_median_equation, use_entrainment_equation, use_deposition_equation !Program settings
            integer :: upwinding_type, is_bed_fixed, bedload_method, hyperconcentration_type, update_grain_size, end_program_at_mobilisation !More program settings

            double precision,dimension(-1:n_steps + 2) :: h,u,z_b,conc, free_surface_elevation !Key flow variables
            double precision,dimension(-1:n_steps + 2) :: x, S_b, bedrock_elevation, bed_slope_angle !Geometry-related variables
            double precision,dimension(-1:n_steps + 2) :: u_1, a_1,a_2,a_3 !Program parameters to solve for
            double precision,dimension(-1:n_steps + 2) :: da1_dt_1,da1_dt_2,da1_dt_3, da1_dt_4 !Fluxes in a1
            double precision,dimension(-1:n_steps + 2) :: da2_dt_1, da2_dt_2,da2_dt_3,da2_dt_4 !Fluxes in a2
            double precision,dimension(-1:n_steps + 2) :: da3_dt_1,da3_dt_2,da3_dt_3,da3_dt_4 !Fluxes in a3
            double precision,dimension(-1:n_steps + 2) :: dzb_dt_1,dzb_dt_2,dzb_dt_3, dzb_dt_4 !Fluxes in the bed elevation
            double precision,dimension(-1:n_steps + 2) :: a_1_new,a_2_new,a_3_new,z_b_new !Updated parameters
            double precision,dimension(-1:n_steps + 2) :: rho, s, rho_w, rho_0, bed_porosity, nu !Density related variables
            double precision,dimension(-1:n_steps + 2) :: tau_bx,tau_bx_1,tau_bx_2,tau_bx_3, c_f, c, c_b, n_mannings, z_0s !Friction coefficient related variables
            double precision,dimension(-1:n_steps + 2) :: d_median, w_s_0, w_s_h, m_coefficient, d_dimensionless !Particle related variables
            double precision,dimension(-1:n_steps + 2) :: dimensionless_critical_shields_parameter, shields_parameter, modified_critical_shields_parameter !Shields parameter variables
            double precision,dimension(-1:n_steps + 2) :: Fr, Re_p, R_n !Dimensionless parameters
            double precision,dimension(-1:n_steps + 2) :: q_b, dq_b_dt, u_b, q_b_0 !Bed load flux and bedload velocity
            double precision,dimension(-1:n_steps + 2) :: D, En, En_1, En_2, En_3, En_4, D_1, D_2, D_3, D_4 !Deposition and entrainment fluxes
            double precision,dimension(-1:n_steps + 2) :: minmod_ab,minmod_a,minmod_b !Minmod slope-limiting values
            double precision,dimension(-1:n_steps + 2) :: z_b_initial, d_median_initial, d_median_new, erosion_depth, deposition_depth !Variables related to changes in geometry
            
            double precision :: a_1_1,a_1_east,a_1_west,a_2_1,a_2_east,a_2_west,a_3_1,a_3_east,a_3_west,z_b_east,z_b_west !Intermediate Runge-Kutta values
            double precision :: h_1,h_east,h_west !Intermediate Runge-Kutta values
            double precision :: g, rho_s, von_karman_constant, rho_w_clear, angle_of_repose !Material or physics-related constants
            double precision :: dx, dt !Mesh-related constants
            double precision :: h_input, q_input, coeff_1_input,coeff_2_input, nu_input, c_b_input,c_input,n_mannings_input !Input values
            double precision :: D_input,En_input, d_median_input, hyperconcentration_conc, h_lim !More input values
            double precision :: coeff_1,coeff_2, coeff_3 !Model coefficients
            double precision :: t, uniform_flow_time !Time related variables
            double precision :: total_sediment_mass, bed_sediment_area !Verification values
            double precision :: critical_shields_parameter,alpha_d,alpha_en, u_shear,m_d !Particle related variables
            double precision :: t_ramp,q_ramp_end,q_ramp_0, q_fixed, rho_w_multiplier !Ramping-related variables
            
            character(len = 5) :: simulation_batch_name

            !Updates the ramped-up surface stress/inflows
            !If the bed is fixed, ignores the ramp up, allowing uniform flow to be established
            if (is_bed_fixed == 1) then
                
                q_fixed = q_input
                
                continue
                
            else
            
                 call ramp_up(n, n_steps,ramp_up_type,hyperconcentration_type, &
                            t,t_ramp,q_ramp_end,q_ramp_0,q_input, rho_w_clear, q_fixed, rho_w_multiplier, &
                            q, h, rho_w)
        
            end if
             
            !Drives the flow if applicable by fixing the upstream flowrate
            if (drive_flow == 1) then
                    
                a_2(0) = q_fixed * (rho_w(0) * (1.0d0 - a_3(0)/(rho_s * (h(0)))) + rho_s * a_3(0)/(rho_s * (h(0))))
                a_2(-1) = 2.0d0 * a_2(0) - a_2(1)
            
            end if

            
            !Calculates the bed slope and bed angle
            do n = 0, n_steps + 1
                
                !Calculates the bed slope
                S_b(n) = (z_b(n + 1) - z_b(n - 1))/(2.0d0 * dx)
                
                !Calculates the bed slope angle
                bed_slope_angle(n) = atan(S_b(n))
                
            end do
  
            !If the hyperconcentration is ramped alongside a sinusoidal flood, ramps it in proportion to the flowrate
            if (hyperconcentration_type == 2) then
                
                do n = -1, n_steps + 2
                    
                    call ramp_concentration(q, n, q_ramp_0, rho_w, rho_w_clear, rho_w_multiplier, pi, q_ramp_end)
            
                end do
                
            end if
                    
            !The initial part of the water-sediment Runge-Kutta calculation
            do n = -1, n_steps + 2
                
                !Calculates the velocity
                u_1(n) = a_2(n)/(a * a_1(n)) !!XX Flag for deletion?
                conc(n) = a_3(n)/(rho_s * h(n))
                 
            end do

            !---------------------------------------------
            !The first step of the Runge-Kutta calculation
            do n = 1, n_steps
                
                !Calculates the sediment-water mixture concentration
                rho(n) = rho_w(n) * (1.0d0 - conc(n)) + rho_s * conc(n)
                
                !Calculates the friction coefficient
                call friction_coefficient(bed_friction_type,n,n_steps,c_f,c,c_b,n_mannings,h,g,c_b_input,&
                                        c_input,n_mannings_input,h_input,h_lim,d_median, von_karman_constant, z_0s)

                !Calculates the friction coefficient
                tau_bx(n) = rho(n) * abs(u_1(n)) * u_1(n) * c_f(n)
                
                !Calculates the particle settling velocity for use initial_calculations calculating deposition
                call particle_transport_calcs(n, d_median_equation, use_deposition_equation, is_bed_fixed,&
                                            s,rho_s, g, d_median_input, rho_w_clear, &
                                            x, tau_bx, d_median, bed_slope_angle, &
                                            d_dimensionless, R_n, u_1, h,D,En, conc,bed_porosity, &
                                            Re_p, m_coefficient, dimensionless_critical_shields_parameter,rho_w, &
                                            Fr, shields_parameter,q, nu, rho, w_s_0, w_s_h, modified_critical_shields_parameter, &
                                            alpha_en, alpha_d, D_input, En_input, &
                                            u_shear, von_karman_constant, m_d, nu_input, angle_of_repose,&
                                            hyperconcentration_type, hyperconcentration_conc)
                
                !Applies boundary conditions to the dimensionless critical Shields parameter
                dimensionless_critical_shields_parameter(0) = dimensionless_critical_shields_parameter(1)
                dimensionless_critical_shields_parameter(-1) = dimensionless_critical_shields_parameter(1)
                dimensionless_critical_shields_parameter(n_steps + 1) = dimensionless_critical_shields_parameter(n_steps)
                dimensionless_critical_shields_parameter(n_steps + 2) = dimensionless_critical_shields_parameter(n_steps + 1)
                
            end do
            
            !Updates the bedload parameters
            do n = -1, n_steps + 2
             
                u_1(n) = a_2(n)/(a * a_1(n))
                call bedload_intermediate_calc(bedload_method, q_b, n, n_steps, u_b, coeff_1, coeff_2, s, rho_s,&
                                                rho_w,rho, conc, tau_bx, u_1, c_f, shields_parameter, g, d_median,&
                                                dimensionless_critical_shields_parameter,dq_b_dt, q_b_0, dt, coeff_3, &
                                                is_bed_fixed, end_program_at_mobilisation,simulation_batch_name, &
                                                hyperconcentration_conc)
                
            end do
            
            q_b(0) = q_b(1)
            q_b(-1) = q_b(0)
            u_b(0) = u_b(1)
            u_b(-1) = u_b(0)
            q_b(n_steps + 1)= q_b(n_steps)
            q_b(n_steps + 2)= q_b(n_steps + 1)
            u_b(n_steps + 1)= u_b(n_steps)
            u_b(n_steps + 2)= u_b(n_steps + 1)
            dq_b_dt(0) = dq_b_dt(1)
            dq_b_dt(-1) = dq_b_dt(0)
            dq_b_dt(n_steps + 1) = dq_b_dt(n_steps)
            dq_b_dt(n_steps + 2) = dq_b_dt(n_steps + 1)
            
            do n = 1, n_steps
                
                !Calculates the difference in bed elevation in each cell
                !Modifies the equation used based on the upwinding type:
                call bed_elevation_change(upwinding_type, n, n_steps, is_bed_fixed,&
                                        q_b,dzb_dt_1,D,En,bed_porosity, bedrock_elevation, &
                                        minmod_a, minmod_b, minmod_ab, dx)
                
                !Ensures no entrainment occurs if a fixed bed is in place
                if (bedrock_elevation_equation .ne. 0 .and. dzb_dt_1(n) == 0.0d0) then
                
                    En(n) = 0.0d0
                    
                end if
                
                !Saves the values of the fluxes to a variable so they can be averaged later on
                En_1(n) = En(n)
                D_1(n) = D(n)
                
                !Calculates the derivative of a_1
                da1_dt_1(n) = -(a_2(n + 1) - a_2(n-1))/(2.0d0 * dx) - rho_0(n) * dzb_dt_1(n) - rho_s*(q_b(n + 1) - q_b(n - 1))/(2.0d0 * dx)
                
                !Calculates the derivative of a_2
                da2_dt_1(n) =   -((b * a_2(n + 1)**2.0d0)/(a**2.0d0 * a_1(n + 1))+(g*a_1(n + 1) * h(n + 1))/(2.0d0) &
                                -(b * a_2(n - 1)**2.0d0)/(a**2.0d0 * a_1(n - 1))-(g*a_1(n - 1) * h(n - 1))/(2.0d0))/(2.0d0 * dx) &
                                - (rho(n) * g * h(n))*(z_b(n + 1) - z_b(n - 1))/(2.0d0 * dx) - tau_bx(n) &
                                - rho_s * dq_b_dt(n) - rho_s * (u_b(n + 1) * q_b(n + 1) - u_b(n - 1) * q_b(n - 1))/(2.0d0 * dx)
                
                !Calculates the derivative of a_3
                da3_dt_1(n) = -((a_3(n + 1) * a_2(n + 1))/(a_1(n + 1))-(a_3(n - 1) * a_2(n - 1))/(a_1(n - 1)))/(2.0d0 * dx) - rho_s * (D(n) - En(n))
                
            end do
        
            !Handles the western boundary being a uniform flow
            if (west_bc == 5) then

                da1_dt_1(0) = 0.0d0
                da2_dt_1(0) = da2_dt_1(1)
                da3_dt_1(0) = 0.0d0
                dzb_dt_1(0) = dzb_dt_1(1)
            
            !Handles the western boundary being transmissive 
            else if (west_bc == 2) then
                
                da1_dt_1(0) = da1_dt_1(1)
                da2_dt_1(0) = da2_dt_1(1)
                da3_dt_1(0) = da3_dt_1(1)
                dzb_dt_1(0) = dzb_dt_1(1)
                
                da1_dt_1(-1) = da1_dt_1(0)
                da2_dt_1(-1) = da2_dt_1(0)
                da3_dt_1(-1) = da3_dt_1(0)
                dzb_dt_1(-1) = dzb_dt_1(0)
                
            else if (west_bc == 3) then
                
                da1_dt_1(0) = da1_dt_1(1)
                da2_dt_1(0) = da2_dt_1(1)
                da3_dt_1(0) = da3_dt_1(1)
                dzb_dt_1(0) = dzb_dt_1(1)
                
                da1_dt_1(-1) = da1_dt_1(0)
                da2_dt_1(-1) = da2_dt_1(0)
                da3_dt_1(-1) = da3_dt_1(0)
                dzb_dt_1(-1) = dzb_dt_1(0)
                
            end if
                
            !Handles the eastern boundary being a uniform flow
            if (east_bc == 5) then
            
                da1_dt_1(n_steps + 1) = 0.0d0
                da2_dt_1(n_steps + 1) = da2_dt_1(n_steps)
                da3_dt_1(n_steps + 1) = 0.0d0
                dzb_dt_1(n_steps + 1) = dzb_dt_1(n_steps)
   
            !If the upstream boundary is a transmissive condition:
            else if (east_bc == 2) then
                
                da1_dt_1(n_steps + 1) = da1_dt_1(n_steps)
                da2_dt_1(n_steps + 1) = da2_dt_1(n_steps)
                da3_dt_1(n_steps + 1) = da3_dt_1(n_steps)
                dzb_dt_1(n_steps + 1) = dzb_dt_1(n_steps)
            
            else if (east_bc == 3) then
                
                da1_dt_1(n_steps + 1) = da1_dt_1(n_steps)
                da2_dt_1(n_steps + 1) = da2_dt_1(n_steps)
                da3_dt_1(n_steps + 1) = da3_dt_1(n_steps)
                dzb_dt_1(n_steps + 1) = dzb_dt_1(n_steps)
                
                da1_dt_1(n_steps + 2) = da1_dt_1(n_steps + 1)
                da2_dt_1(n_steps + 2) = da2_dt_1(n_steps + 1)
                da3_dt_1(n_steps + 2) = da3_dt_1(n_steps + 1)
                dzb_dt_1(n_steps + 2) = dzb_dt_1(n_steps + 1)
                
            else
                
                
                return
                
            end if
         
            !Updates the parameters for the full range of steps including the end steps
            do n = -1, n_steps + 2
                
                a_1_1 = a_1(n) + dt * da1_dt_1(n)/2.0d0
                a_2_1 = a_2(n) + dt * da2_dt_1(n)/2.0d0 
                u_1(n) = a_2_1/(a * a_1_1)

                call bedload_intermediate_calc(bedload_method, q_b, n, n_steps, u_b, coeff_1, coeff_2, s, rho_s,&
                                                rho_w,rho, conc, tau_bx, u_1, c_f, shields_parameter, g, d_median,&
                                                dimensionless_critical_shields_parameter,dq_b_dt, q_b_0, dt, coeff_3, &
                                                is_bed_fixed, end_program_at_mobilisation,simulation_batch_name,&
                                                hyperconcentration_conc)
            
            end do
            q_b(0) = q_b(1)
            q_b(-1) = q_b(0)
            u_b(0) = u_b(1)
            u_b(-1) = u_b(0)
            q_b(n_steps + 1)= q_b(n_steps)
            q_b(n_steps + 2)= q_b(n_steps + 1)
            u_b(n_steps + 1)= u_b(n_steps)
            u_b(n_steps + 2)= u_b(n_steps + 1)
            dq_b_dt(0) = dq_b_dt(1)
            dq_b_dt(-1) = dq_b_dt(0)
            dq_b_dt(n_steps + 1) = dq_b_dt(n_steps)
            dq_b_dt(n_steps + 2) = dq_b_dt(n_steps + 1)
            
            !---------------------------------------------
            !The second step of the Runge-Kutta method
            do n = 1, n_steps
                
                !print*,'a_3_1',n,da3_dt_1(n)
                
                a_1_1 = a_1(n) + dt * da1_dt_1(n)/2.0d0
                a_1_east = a_1(n + 1) + dt * da1_dt_1(n + 1)/2.0d0
                a_1_west = a_1(n - 1) + dt * da1_dt_1(n - 1)/2.0d0

                a_2_1 = a_2(n) + dt * da2_dt_1(n)/2.0d0
                a_2_east = a_2(n + 1) + dt * da2_dt_1(n + 1)/2.0d0
                a_2_west = a_2(n - 1) + dt * da2_dt_1(n - 1)/2.0d0
                
                a_3_1 = a_3(n) + dt * da3_dt_1(n)/2.0d0
                a_3_east = a_3(n + 1) + dt * da3_dt_1(n + 1)/2.0d0
                a_3_west = a_3(n - 1) + dt * da3_dt_1(n - 1)/2.0d0
                
                z_b_east = z_b(n + 1) + dt * dzb_dt_1(n + 1)/2.0d0
                z_b_west = z_b(n - 1) + dt * dzb_dt_1(n - 1)/2.0d0
                
                h_1 = (a_1_1 - (1.0d0 - rho_w(n)/rho_s) * a_3_1)/rho_w(n)
                h_east = (a_1_east - (1.0d0 - rho_w(n)/rho_s) * a_3_east)/rho_w(n)
                h_west = (a_1_west - (1.0d0 - rho_w(n)/rho_s) * a_3_west)/rho_w(n)
                
                !Calculates the sediment concentration from the problem coefficients
                conc(n) = a_3_1/(rho_s * h_1)
                
                !Calculates the sediment-water mixture concentration
                rho(n) = rho_w(n) * (1.0d0 - conc(n)) + rho_s * conc(n)
                
                !Calculates the friction coefficient
                call friction_coefficient(bed_friction_type,n,n_steps,c_f,c,c_b,n_mannings,h,g,c_b_input,&
                                        c_input,n_mannings_input,h_input,h_lim,d_median, von_karman_constant, z_0s)
                
                !Recalculates velocity
                u_1(n) = a_2_1/(a * a_1_1)
                
                !Calculates the bed shear stress
                tau_bx_1(n) = rho(n) * abs(u_1(n))*u_1(n)*c_f(n)
                
                !Calls the particle transport equations
                call particle_transport_calcs(n, d_median_equation, use_deposition_equation, is_bed_fixed,&
                                            s,rho_s, g, d_median_input, rho_w_clear, &
                                            x, tau_bx_1, d_median, bed_slope_angle, &
                                            d_dimensionless, R_n, u_1, h,D,En, conc,bed_porosity, &
                                            Re_p, m_coefficient, dimensionless_critical_shields_parameter,rho_w, &
                                            Fr, shields_parameter,q, nu, rho, w_s_0, w_s_h, modified_critical_shields_parameter, &
                                            alpha_en, alpha_d, D_input, En_input, &
                                            u_shear, von_karman_constant, m_d, nu_input, angle_of_repose,&
                                            hyperconcentration_type, hyperconcentration_conc)
                
                !Calculates the difference in bed elevation in each cell
                call bed_elevation_change(upwinding_type, n, n_steps, is_bed_fixed,&
                                        q_b,dzb_dt_2,D,En,bed_porosity, bedrock_elevation, &
                                        minmod_a, minmod_b, minmod_ab, dx)     
                
                !Ensures no entrainment occurs if a fixed bed is in place
                if (bedrock_elevation_equation .ne. 0 .and. dzb_dt_2(n) == 0.0d0) then
                
                    En(n) = 0.0d0
                    
                end if
                
                !Saves the values of the fluxes to a variable so they can be averaged later on
                En_2(n) = En(n)
                D_2(n) = D(n)
                dimensionless_critical_shields_parameter(0) = dimensionless_critical_shields_parameter(1)
                dimensionless_critical_shields_parameter(-1) = dimensionless_critical_shields_parameter(1)
                dimensionless_critical_shields_parameter(n_steps + 1) = dimensionless_critical_shields_parameter(n_steps)
                dimensionless_critical_shields_parameter(n_steps + 2) = dimensionless_critical_shields_parameter(n_steps + 1)
          
                !Calculates the derivative of a_1
                da1_dt_2(n) = -(a_2_east - a_2_west)/(2.0d0 * dx) - rho_0(n) * dzb_dt_2(n) - rho_s*(q_b(n + 1) - q_b(n - 1))/(2.0d0 * dx)
                
                !Calculates the derivative of a_2
                da2_dt_2(n) =   -((b * a_2_east**2.0d0)/(a**2.0d0 * a_1_east)+(g*a_1_east * h_east)/(2.0d0) &
                                -(b * a_2_west**2.0d0)/(a**2.0d0 * a_1_west)-(g*a_1_west * h_west)/(2.0d0))/(2.0d0 * dx) &
                                - (rho(n) * g * h_1)*(z_b_east - z_b_west)/(2.0d0 * dx) - tau_bx_1(n) &
                              - rho_s * dq_b_dt(n) - rho_s * (u_b(n + 1) * q_b(n + 1) - u_b(n - 1) * q_b(n - 1))/(2.0d0 * dx)
                
                !Calculates the derivative of a_3
                da3_dt_2(n) = -((a_3_east * a_2_east)/(a_1_east)-(a_3_west * a_2_west)/(a_1_west))/(2.0d0 * dx) - rho_s * (D(n) - En(n))

            end do
            
            !Handles the western boundary being a uniform flow
            if (west_bc == 5) then

                da1_dt_2(0) = 0.0d0
                da2_dt_2(0) = da2_dt_2(1)
                da3_dt_2(0) = 0.0d0
                dzb_dt_2(0) = dzb_dt_2(1)
                
                da1_dt_2(-1) = da1_dt_2(0)
                da2_dt_2(-1) = da2_dt_2(0)
                da3_dt_2(-1) = da3_dt_2(0)
                dzb_dt_2(-1) = dzb_dt_2(0)
                
            !If the upstream boundary is a transmissive condition:
            else if (west_bc == 2) then
                
                da1_dt_2(0) = da1_dt_2(1)
                da2_dt_2(0) = da2_dt_2(1)
                da3_dt_2(0) = da3_dt_2(1)
                dzb_dt_2(0) = dzb_dt_2(1)
                
                da1_dt_2(-1) = da1_dt_2(0)
                da2_dt_2(-1) = da2_dt_2(0)
                da3_dt_2(-1) = da3_dt_2(0)
                dzb_dt_2(-1) = dzb_dt_2(0)
                
            else if (west_bc == 3) then
                
                da1_dt_2(0) = da1_dt_2(1)
                da2_dt_2(0) = da2_dt_2(1)
                da3_dt_2(0) = da3_dt_2(1)
                dzb_dt_2(0) = dzb_dt_2(1)
                
                da1_dt_2(-1) = da1_dt_2(0)
                da2_dt_2(-1) = da2_dt_2(0)
                da3_dt_2(-1) = da3_dt_2(0)
                dzb_dt_2(-1) = dzb_dt_2(0)
            
            end if
                
            !Handles the eastern boundary being a uniform flow
            if (east_bc == 5) then
            
                da1_dt_2(n_steps + 1) = 0.0d0
                da2_dt_2(n_steps + 1) = da2_dt_2(n_steps)
                da3_dt_2(n_steps + 1) = 0.0d0
                dzb_dt_3(n_steps+1) = dzb_dt_3(n_steps)
            
            !If the upstream boundary is a transmissive condition:
            else if (east_bc == 2) then
                
                da1_dt_2(n_steps + 1) = da1_dt_2(n_steps)
                da2_dt_2(n_steps + 1) = da2_dt_2(n_steps)
                da3_dt_2(n_steps + 1) = da3_dt_2(n_steps)
                dzb_dt_2(n_steps + 1) = dzb_dt_2(n_steps)
                
                da1_dt_2(n_steps + 2) = da1_dt_2(n_steps + 1)
                da2_dt_2(n_steps + 2) = da2_dt_2(n_steps + 1)
                da3_dt_2(n_steps + 2) = da3_dt_2(n_steps + 1)
                dzb_dt_2(n_steps + 2) = dzb_dt_2(n_steps + 1)
                
            else if (east_bc == 3) then
                
                da1_dt_2(n_steps + 1) = da1_dt_2(n_steps)
                da2_dt_2(n_steps + 1) = da2_dt_2(n_steps)
                da3_dt_2(n_steps + 1) = da3_dt_2(n_steps)
                dzb_dt_2(n_steps + 1) = dzb_dt_2(n_steps)
                
                da1_dt_2(n_steps + 2) = da1_dt_2(n_steps + 1)
                da2_dt_2(n_steps + 2) = da2_dt_2(n_steps + 1)
                da3_dt_2(n_steps + 2) = da3_dt_2(n_steps + 1)
                dzb_dt_2(n_steps + 2) = dzb_dt_2(n_steps + 1)
                
                else
                
                end if
                
            !Updates the parameters for the full range of steps including the end steps
            do n = -1, n_steps + 2
                
                a_1_1 = a_1(n) + dt * da1_dt_2(n)/2.0d0
                a_2_1 = a_2(n) + dt * da2_dt_2(n)/2.0d0
                
                u_1(n) = a_2_1/(a * a_1_1)
                
                call bedload_intermediate_calc(bedload_method, q_b, n, n_steps, u_b, coeff_1, coeff_2, s, rho_s,&
                                                rho_w,rho, conc, tau_bx, u_1, c_f, shields_parameter, g, d_median,&
                                                dimensionless_critical_shields_parameter,dq_b_dt, q_b_0, dt, coeff_3, &
                                                is_bed_fixed, end_program_at_mobilisation,simulation_batch_name, &
                                                hyperconcentration_conc)
                
            end do    
            
            !Applies boundary conditions to the bedload fluxes
            q_b(0) = q_b(1)
            q_b(-1) = q_b(0)
            u_b(0) = u_b(1)
            u_b(-1) = u_b(0)
            q_b(n_steps + 1)= q_b(n_steps)
            q_b(n_steps + 2)= q_b(n_steps + 1)
            u_b(n_steps + 1)= u_b(n_steps)
            u_b(n_steps + 2)= u_b(n_steps + 1)
            dq_b_dt(0) = dq_b_dt(1)
            dq_b_dt(-1) = dq_b_dt(0)
            dq_b_dt(n_steps + 1) = dq_b_dt(n_steps)
            dq_b_dt(n_steps + 2) = dq_b_dt(n_steps + 1)

            !---------------------------------------------
            !The third step of the Runge-Kutta method
            do n = 1, n_steps
                
                !Calculates the updated parameters at the new position
                a_1_1 = a_1(n) + dt * da1_dt_2(n)/2.0d0
                a_1_east = a_1(n + 1) + dt * da1_dt_2(n + 1)/2.0d0
                a_1_west = a_1(n - 1) + dt * da1_dt_2(n - 1)/2.0d0
                
                a_2_1 = a_2(n) + dt * da2_dt_2(n)/2.0d0
                a_2_east = a_2(n + 1) + dt * da2_dt_2(n + 1)/2.0d0
                a_2_west = a_2(n - 1) + dt * da2_dt_2(n - 1)/2.0d0
                
                a_3_1 = a_3(n) + dt * da3_dt_2(n)/2.0d0
                a_3_east = a_3(n + 1) + dt * da3_dt_2(n + 1)/2.0d0
                a_3_west = a_3(n - 1) + dt * da3_dt_2(n - 1)/2.0d0
                
                z_b_east = z_b(n + 1) + dt * dzb_dt_2(n + 1)/2.0d0
                z_b_west = z_b(n - 1) + dt * dzb_dt_2(n - 1)/2.0d0
                
                h_1 = (a_1_1 - (1.0d0 - rho_w(n)/rho_s) * a_3_1)/rho_w(n)
                h_east = (a_1_east - (1.0d0 - rho_w(n)/rho_s) * a_3_east)/rho_w(n)
                h_west = (a_1_west - (1.0d0 - rho_w(n)/rho_s) * a_3_west)/rho_w(n)
                
                !Calculates the sediment concentration from the model parameters
                conc(n) = a_3_1/(rho_s * h_1)
                
                !Calculates the sediment-water mixture concentration
                rho(n) = rho_w(n) * (1.0d0 - conc(n)) + rho_s * conc(n)
                
                !Calculates the friction coefficient
                call friction_coefficient(bed_friction_type,n,n_steps,c_f,c,c_b,n_mannings,h,g,c_b_input,&
                                        c_input,n_mannings_input,h_input,h_lim,d_median, von_karman_constant, z_0s)
                
                !Recalculates velocity
                u_1(n) = a_2_1/(a * a_1_1)
                
                !Calculates the bed shear stress from the friction coefficient
                tau_bx_2(n) = rho(n) * abs(u_1(n))*u_1(n)*c_f(n)
                
                !Calls the particle transport equations subroutine
                call particle_transport_calcs(n, d_median_equation, use_deposition_equation, is_bed_fixed,&
                                            s,rho_s, g, d_median_input, rho_w_clear, &
                                            x, tau_bx_2, d_median, bed_slope_angle, &
                                            d_dimensionless, R_n, u_1, h,D,En, conc,bed_porosity, &
                                            Re_p, m_coefficient, dimensionless_critical_shields_parameter,rho_w, &
                                            Fr, shields_parameter,q, nu, rho, w_s_0, w_s_h, modified_critical_shields_parameter, &
                                            alpha_en, alpha_d, D_input, En_input, &
                                            u_shear, von_karman_constant, m_d, nu_input, angle_of_repose,&
                                            hyperconcentration_type, hyperconcentration_conc)
                
                !Calls the subroutine which calculates the change in bed elevation
                call bed_elevation_change(upwinding_type, n, n_steps, is_bed_fixed,&
                                        q_b,dzb_dt_3,D,En,bed_porosity, bedrock_elevation, &
                                        minmod_a, minmod_b, minmod_ab, dx) 
                
                !Ensures no entrainment occurs if a fixed bed is in place
                if (bedrock_elevation_equation .ne. 0 .and. dzb_dt_3(n) == 0.0d0) then
                
                    En(n) = 0.0d0
                    
                end if
                
                !Saves the values of the fluxes to a variable so they can be averaged later on
                En_3(n) = En(n)
                D_3(n) = D(n)
                
                !Updates the end values of the dimensionless critical shields parameter
                dimensionless_critical_shields_parameter(0) = dimensionless_critical_shields_parameter(1)
                dimensionless_critical_shields_parameter(-1) = dimensionless_critical_shields_parameter(1)
                dimensionless_critical_shields_parameter(n_steps + 1) = dimensionless_critical_shields_parameter(n_steps)
                dimensionless_critical_shields_parameter(n_steps + 2) = dimensionless_critical_shields_parameter(n_steps + 1)
                
                !Calculates the derivative of a_1
                da1_dt_3(n) = -(a_2_east - a_2_west)/(2.0d0 * dx) - rho_0(n) * dzb_dt_3(n) - rho_s*(q_b(n + 1) - q_b(n - 1))/(2.0d0 * dx)
                
                !Calculates the derivative of a_2
                da2_dt_3(n) =   -((b * a_2_east**2.0d0)/(a**2.0d0 * a_1_east)+(g*a_1_east * h_east)/(2.0d0) &
                                -(b * a_2_west**2.0d0)/(a**2.0d0 * a_1_west)-(g*a_1_west * h_west)/(2.0d0))/(2.0d0 * dx) &
                                - (rho(n) * g * h_1)*(z_b_east - z_b_west)/(2.0d0 * dx) - tau_bx_2(n) &
                              - rho_s * dq_b_dt(n) - rho_s * (u_b(n + 1) * q_b(n + 1) - u_b(n - 1) * q_b(n - 1))/(2.0d0 * dx)
                
                !Calculates the derivative of a_3
                da3_dt_3(n) = -((a_3_east * a_2_east)/(a_1_east)-(a_3_west * a_2_west)/(a_1_west))/(2.0d0 * dx) - rho_s * (D(n) - En(n))
                
            end do
            
            !Handles the western boundary being a uniform flow
            if (west_bc == 5) then
                
                da1_dt_3(0) = 0.0d0
                da2_dt_3(0) = da2_dt_3(1)
                da3_dt_3(0) = 0.0d0
                dzb_dt_3(0) = dzb_dt_3(1)
            
            !Handles the western boundary being a transmissive flow
            else if (west_bc == 2) then
                
                da1_dt_3(0) = da1_dt_3(1)
                da2_dt_3(0) = da2_dt_3(1)
                da3_dt_3(0) = da3_dt_3(1)
                dzb_dt_3(0) = dzb_dt_3(1)
                
                da1_dt_3(-1) = da1_dt_3(0)
                da2_dt_3(-1) = da2_dt_3(0)
                da3_dt_3(-1) = da3_dt_3(0)
                dzb_dt_3(-1) = dzb_dt_3(0)
                
            end if

            !Handles the eastern boundary being a uniform flow
            if (east_bc == 5) then
            
                da1_dt_3(n_steps + 1) = 0.0d0
                da2_dt_3(n_steps + 1) = da2_dt_3(n_steps)
                da3_dt_3(n_steps + 1) = 0.0d0
                dzb_dt_3(n_steps + 1) = dzb_dt_3(n_steps)
                  
            !Applies boundary conditions if the eastern boundary is a transmissive flow
            else if (east_bc == 2) then
                
                da1_dt_3(n_steps + 1) = da1_dt_3(n_steps)
                da2_dt_3(n_steps + 1) = da2_dt_3(n_steps)
                da3_dt_3(n_steps + 1) = da3_dt_3(n_steps)
                dzb_dt_3(n_steps + 1) = dzb_dt_3(n_steps)
                
                da1_dt_3(n_steps + 2) = da1_dt_3(n_steps)
                da2_dt_3(n_steps + 2) = da2_dt_3(n_steps)
                da3_dt_3(n_steps + 2) = da3_dt_3(n_steps)
                dzb_dt_3(n_steps + 2) = dzb_dt_3(n_steps)
                
            end if
                
            !Updates the parameters for the full range of steps including the end steps
            do n = -1, n_steps + 2
                
                a_1_1 = a_1(n) + dt * da1_dt_3(n)/2.0d0
                a_2_1 = a_2(n) + dt * da2_dt_3(n)/2.0d0
                u_1(n) = a_2_1/(a * a_1_1)

                !Calculates the bedload flux q_b and bed particle velocity u_b
                call bedload_intermediate_calc(bedload_method, q_b, n, n_steps, u_b, coeff_1, coeff_2, s, rho_s,&
                                                rho_w,rho, conc, tau_bx, u_1, c_f, shields_parameter, g, d_median,&
                                                dimensionless_critical_shields_parameter,dq_b_dt, q_b_0, dt, coeff_3, &
                                                is_bed_fixed, end_program_at_mobilisation,simulation_batch_name, &
                                                hyperconcentration_conc)
                
            end do  
            
            !Updates the boundary conditions for the fluxes
            q_b(0) = q_b(1)
            q_b(-1) = q_b(0)
            u_b(0) = u_b(1)
            u_b(-1) = u_b(0)
            q_b(n_steps + 1)= q_b(n_steps)
            q_b(n_steps + 2)= q_b(n_steps + 1)
            u_b(n_steps + 1)= u_b(n_steps)
            u_b(n_steps + 2)= u_b(n_steps + 1)
            dq_b_dt(0) = dq_b_dt(1)
            dq_b_dt(-1) = dq_b_dt(0)
            dq_b_dt(n_steps + 1) = dq_b_dt(n_steps)
            dq_b_dt(n_steps + 2) = dq_b_dt(n_steps + 1)
            
            !---------------------------------------------
            !The FOURTH step of the Runge-Kutta method
            do n = 1, n_steps
                
                a_1_1 = a_1(n) + dt * da1_dt_3(n)
                a_1_east = a_1(n + 1) + dt * da1_dt_3(n + 1)
                a_1_west = a_1(n - 1) + dt * da1_dt_3(n - 1)
                
                a_2_1 = a_2(n) + dt * da2_dt_3(n)
                a_2_east = a_2(n + 1) + dt * da2_dt_3(n + 1)
                a_2_west = a_2(n - 1) + dt * da2_dt_3(n - 1)
                
                a_3_1 = a_3(n) + dt * da3_dt_3(n)
                a_3_east = a_3(n + 1) + dt * da3_dt_3(n + 1)
                a_3_west = a_3(n - 1) + dt * da3_dt_3(n - 1)
                
                z_b_east = z_b(n + 1) + dt * dzb_dt_3(n + 1)
                z_b_west = z_b(n - 1) + dt * dzb_dt_3(n - 1)
                
                h_1 = (a_1_1 - (1.0d0 - rho_w(n)/rho_s) * a_3_1)/rho_w(n)
                h_east = (a_1_east - (1.0d0 - rho_w(n)/rho_s) * a_3_east)/rho_w(n)
                h_west = (a_1_west - (1.0d0 - rho_w(n)/rho_s) * a_3_west)/rho_w(n)
                
                !Calculates the sediment concentration from the problem coefficients
                conc(n) = a_3_1/(rho_s * h_1)
                
                !Calculates the sediment-water mixture concentration
                rho(n) = rho_w(n) * (1.0d0 - conc(n)) + rho_s * conc(n)
                
                !Calculates the friction coefficient
                call friction_coefficient(bed_friction_type,n,n_steps,c_f,c,c_b,n_mannings,h,g,c_b_input,&
                                        c_input,n_mannings_input,h_input,h_lim,d_median, von_karman_constant, z_0s)
                
                !Recalculates velocity
                u_1(n) = a_2_1/(a * a_1_1)
                
                !Calculates the friction coefficient
                tau_bx_3(n) = rho(n) * abs(u_1(n))*u_1(n) * c_f(n)
                
                !Calls the particle transport equations
                call particle_transport_calcs(n, d_median_equation, use_deposition_equation, is_bed_fixed,&
                                            s,rho_s, g, d_median_input, rho_w_clear, &
                                            x, tau_bx_3, d_median, bed_slope_angle, &
                                            d_dimensionless, R_n, u_1, h, D, En, conc,bed_porosity, &
                                            Re_p, m_coefficient, dimensionless_critical_shields_parameter,rho_w, &
                                            Fr, shields_parameter,q, nu, rho, w_s_0, w_s_h, modified_critical_shields_parameter, &
                                            alpha_en, alpha_d, D_input, En_input, &
                                            u_shear, von_karman_constant, m_d, nu_input, angle_of_repose,&
                                            hyperconcentration_type, hyperconcentration_conc)
                
                !Calculates the difference in bed elevation in each cell
                call bed_elevation_change(upwinding_type, n, n_steps, is_bed_fixed,&
                                        q_b,dzb_dt_4,D,En,bed_porosity, bedrock_elevation, &
                                        minmod_a, minmod_b, minmod_ab, dx)  
                
                !Ensures no entrainment occurs if a fixed bed is in place
                if (bedrock_elevation_equation .ne. 0 .and. dzb_dt_3(n) == 0.0d0) then
                
                    En(n) = 0.0d0
                    
                end if
   
                !Saves the values of the fluxes to a variable so they can be averaged later on
                En_3(n) = En(n)
                D_3(n) = D(n)
                En_4(n) = En(n)
                D_4(n) = D(n)
                dimensionless_critical_shields_parameter(0) = dimensionless_critical_shields_parameter(1)
                dimensionless_critical_shields_parameter(-1) = dimensionless_critical_shields_parameter(1)
                dimensionless_critical_shields_parameter(n_steps + 1) = dimensionless_critical_shields_parameter(n_steps)
                dimensionless_critical_shields_parameter(n_steps + 2) = dimensionless_critical_shields_parameter(n_steps + 1)
                
                !Calculates the derivative of a_1
                da1_dt_4(n) = -(a_2_east - a_2_west)/(2.0d0 * dx) - rho_0(n) * dzb_dt_4(n) - rho_s*(q_b(n + 1) - q_b(n - 1))/(2.0d0 * dx)
                
                !Calculates the derivative of a_2
                da2_dt_4(n) =   -((b * a_2_east**2.0d0)/(a**2.0d0 * a_1_east)+(g*a_1_east * h_east)/(2.0d0) &
                                -(b * a_2_west**2.0d0)/(a**2.0d0 * a_1_west)-(g*a_1_west * h_west)/(2.0d0))/(2.0d0 * dx) &
                                - (rho(n) * g * h_1)*(z_b_east - z_b_west)/(2.0d0 * dx) - tau_bx_3(n) &
                              - rho_s * dq_b_dt(n) - rho_s * (u_b(n + 1) * q_b(n + 1) - u_b(n - 1) * q_b(n - 1))/(2.0d0 * dx)
                
                !Calculates the derivative of a_3
                da3_dt_4(n) = -((a_3_east * a_2_east)/(a_1_east)-(a_3_west * a_2_west)/(a_1_west))/(2.0d0 * dx) - rho_s * (D(n) - En(n))
                
                !Carries out the actual Runge Kutta method by adding over multiple steps
                a_1_new(n) = a_1(n) + dt * (da1_dt_1(n) + 2.0d0 * da1_dt_2(n) +2.0d0 * da1_dt_3(n) + da1_dt_4(n))/6.0d0
                a_2_new(n) = a_2(n) + dt * (da2_dt_1(n) + 2.0d0 * da2_dt_2(n) +2.0d0 * da2_dt_3(n) + da2_dt_4(n))/6.0d0
                a_3_new(n) = a_3(n) + dt * (da3_dt_1(n) + 2.0d0 * da3_dt_2(n) +2.0d0 * da3_dt_3(n) + da3_dt_4(n))/6.0d0
                z_b_new(n) = z_b(n) + dt * (dzb_dt_1(n) + 2.0d0 * dzb_dt_2(n) +2.0d0 * dzb_dt_3(n) + dzb_dt_4(n))/6.0d0
                
                !print*,n,dzb_dt_1(n),dzb_dt_2(n),dzb_dt_4(n),dzb_dt_3(n)
                
                !Calculates the averaged entrainment and deposition over these four steps
                En(n) = (En_1(n) + 2.0d0 * En_2(n) +2.0d0 * En_3(n) + En_4(n))/6.0d0
                D(n) = (D_1(n) + 2.0d0 * D_2(n) +2.0d0 * D_3(n) + D_4(n))/6.0d0
                
                if (z_b_new(n) <= bedrock_elevation(n) .AND. bedrock_elevation(n) .ne. -100.0d0) then
                
                    z_b_new(n) = bedrock_elevation(n)
                    
                end if
                
            end do
            
            !Updates the boundary conditions
            !Handles the western boundary condition
            if (west_bc == 5) then
                
                a_1_new(0) = a_1_new(1)
                a_1_new(-1) = a_1_new(0)
                a_2_new(0) = a_2_new(1)
                a_2_new(-1) = a_2_new(0)
                a_3_new(0) = a_3_new(1)
                a_3_new(-1) = a_3_new(0)
                z_b_new(0) = z_b_new(1)
                z_b_new(-1) = z_b_new(0)


            else if (west_bc == 2) then
                
                !a_1_new(0) = a_1_new(1)
                !a_1_new(-1) = a_1_new(0)
                
                !a_2_new(0) = a_2_new(1) 
                !a_2_new(-1) = a_2_new(0) 

                !a_3_new(0) = a_3_new(1)
                !a_3_new(-1) = a_3_new(0) 
                
                
                    
                a_2_new(0) = 2.0d0 * a_2_new(1) - a_2_new(2)
                a_2_new(-1) = 2.0d0 * a_2_new(0) - a_2_new(1)
            
                
                
                
                !a_1_new(1) = 2.0d0 * a_1_new(2) - a_1_new(3)
                a_1_new(0) = 2.0d0 * a_1_new(1) - a_1_new(2)
                a_1_new(-1) = 2.0d0 * a_1_new(0) - a_1_new(1)
                
                !a_2_new(1) = 2.0d0 * a_2_new(2) - a_2_new(3)
                !a_2_new(0) = 2.0d0 * a_2_new(1) - a_2_new(2)
                !a_2_new(-1) = 2.0d0 * a_2_new(0) - a_2_new(1)

                a_3_new(1) = 2.0d0 * a_3_new(2) - a_3_new(3)
                a_3_new(0) = 2.0d0 * a_3_new(1) - a_3_new(2)
                !a_3_new(-1) = 2.0d0 * a_3_new(0) - a_3_new(1)
            
                !z_b_new(1) = 2.0d0 * z_b_new(2) - z_b_new(3)
                z_b_new(0) = 2.0d0 * z_b_new(1) - z_b_new(2)
                z_b_new(-1) = 2.0d0 * z_b_new(0) - z_b_new(1)

                
            else
                
                return
                
            end if
            
            !Handles the eastern boundary condition
            if (east_bc == 5) then
                
                a_1_new(n_steps + 1) = a_1_new(n_steps)
                a_1_new(n_steps + 2) = a_1_new(n_steps + 1)
                
                a_2_new(n_steps + 1) = a_2_new(n_steps)
                a_2_new(n_steps + 2) = a_2_new(n_steps + 1)
                
                a_3_new(n_steps + 1) = a_3_new(n_steps)
                a_3_new(n_steps + 2) = a_3_new(n_steps + 1)
                
                z_b_new(n_steps + 1) = z_b_new(n_steps)
                z_b_new(n_steps + 2) = z_b_new(n_steps + 1)
                
            else if (east_bc == 2) then
                
                !a_1_new(n_steps + 1) = a_1_new(n_steps)
                !a_1_new(n_steps + 2) = a_1_new(n_steps + 1)
                
                !a_2_new(n_steps + 1) = a_2_new(n_steps)
                !a_2_new(n_steps + 2) = a_2_new(n_steps + 1)
                
                !a_3_new(n_steps + 1) = a_3_new(n_steps)
                !a_3_new(n_steps + 2) = a_3_new(n_steps + 1)
                
                a_1_new(n_steps + 1) = 2.0d0 * a_1_new(n_steps) - a_1_new(n_steps - 1)
                a_1_new(n_steps + 2) = 2.0d0 * a_1_new(n_steps + 1) - a_1_new(n_steps)
                
                a_2_new(n_steps + 1) = 2.0d0 * a_2_new(n_steps) - a_2_new(n_steps - 1)
                a_2_new(n_steps + 2) = 2.0d0 * a_2_new(n_steps + 1) - a_2_new(n_steps)
                
                a_3_new(n_steps + 1) = 2.0d0 * a_3_new(n_steps) - a_3_new(n_steps - 1)
                a_3_new(n_steps + 2) = 2.0d0 * a_3_new(n_steps + 1) - a_3_new(n_steps)
                
                z_b_new(n_steps + 1) = 2.0d0 * z_b_new(n_steps) - z_b_new(n_steps - 1)
                z_b_new(n_steps + 2) = 2.0d0 * z_b_new(n_steps + 1) - z_b_new(n_steps)
                
                
            end if
             
            !Calculates the depth of erosion/deposition for the step for use in visualisation of results
            do n = -1, n_steps + 2
                
                if (z_b_new(n) - z_b(n) >= 0.0d0) then
                    
                    Deposition_depth(n) = z_b_new(n) - z_b(n)
                    Erosion_depth(n) = 0.0d0
                
                else
                    
                    Erosion_depth(n) = z_b(n) - z_b_new(n)
                    Deposition_depth(n) = 0.0d0
                    
                end if
                    
            end do
                
            !Updates the output parameters not used in the next stepconc(n_steps + 1) = 2.0d0 * conc(n_steps) - conc(n_steps - 1)
            conc(0) = 2.0d0 * conc(1) - conc(2)
            conc(-1) = 2.0d0 * conc(1) - conc(0)
            conc(n_steps + 1) = 2.0d0 * conc(n_steps) - conc(n_steps - 1)
            conc(n_steps + 2) = 2.0d0 * conc(n_steps + 1) - conc(n_steps)
            
            R_n(0) = 2.0d0 * R_n(1) - R_n(2)
            R_n(-1) = 2.0d0 * R_n(1) - R_n(0)
            R_n(n_steps + 1) = 2.0d0 * R_n(n_steps) - R_n(n_steps - 1)
            R_n(n_steps + 2) = 2.0d0 * R_n(n_steps + 1) - R_n(n_steps)
            
            En(0) = 2.0d0 * En(1) - En(2)
            En(-1) = 2.0d0 * En(1) - En(0)
            En(n_steps + 1) = 2.0d0 * En(n_steps) - En(n_steps - 1)
            En(n_steps + 2) = 2.0d0 * En(n_steps + 1) - En(n_steps)
            
            D(0) = 2.0d0 * D(1) - D(2)
            D(-1) = 2.0d0 * D(1) - D(0)
            D(n_steps + 1) = 2.0d0 * D(n_steps) - D(n_steps - 1)
            D(n_steps + 2) = 2.0d0 * D(n_steps + 1) - D(n_steps)
            
            !Updates the results for the next timestep
            do n = -1, n_steps + 2  
                
                u_1(n) = a_2(n)/(a*a_1(n))
                
                call bedload_intermediate_calc(bedload_method, q_b, n, n_steps, u_b, coeff_1, coeff_2, s, rho_s,&
                                                rho_w,rho, conc, tau_bx, u_1, c_f, shields_parameter, g, d_median,&
                                                dimensionless_critical_shields_parameter,dq_b_dt, q_b_0, dt, coeff_3, &
                                                is_bed_fixed, end_program_at_mobilisation,simulation_batch_name, &
                                                hyperconcentration_conc)

                q_b_0(n) = q_b(n)
                a_1(n) = a_1_new(n)
                a_2(n) = a_2_new(n)
                a_3(n) = a_3_new(n)
                z_b(n) = z_b_new(n)
                h(n) = (a_1(n) - (1.0d0 - rho_w(n)/rho_s)*a_3(n))/rho_w(n)
                u(n) = (a_2(n))/(a*a_1(n))
                q(n) = h(n) * u(n)
                conc(n) = a_3(n)/(rho_s * h(n))
                free_surface_elevation(n) = h(n) + z_b(n)
                Fr(n) = u(n)/sqrt(g * h(n))
                
                !print*,n,h(n),u(n)
                !print*,n,h(n),z_b(n)
                
            end do
            
            !Updates the boundary conditions for the bedload transport
            q_b_0(0)= q_b_0(1)
            q_b_0(-1) = q_b_0(0)
            q_b_0(n_steps + 1)= q_b_0(n_steps)
            q_b_0(n_steps + 2)= q_b_0(n_steps + 1)
            u_b(0) = u_b(1)
            u_b(-1) = u_b(0)
            u_b(n_steps + 1)= u_b(n_steps)
            u_b(n_steps + 2)= u_b(n_steps + 1)
            dq_b_dt(0) = dq_b_dt(1)
            dq_b_dt(-1) = dq_b_dt(0)
            dq_b_dt(n_steps + 1) = dq_b_dt(n_steps)
            dq_b_dt(n_steps + 2) = dq_b_dt(n_steps + 1)

            
            if (update_grain_size == 0) then
            
                continue
                
            else if (update_grain_size == 1) then
                    
                    call alter_grain_size(d_median,n_steps,n,update_grain_size,z_b_new,z_b,&
                                            z_b_initial,d_median_initial,angle_of_repose,bed_slope_angle, q_b, dx,d_median_new)
    
                
            end if

            !Sets the bed sediment and total sediment masses to an initial value of zero
            bed_sediment_area = 0.0d0
            total_sediment_mass = 0.0d0
            
            !Sums the total area of the bed sediment and the total sediment mass of the suspended sediment over all spatial steps to check mass conservation
            do n = 1, n_steps + 1
                
                bed_sediment_area = bed_sediment_area + dx * (z_b(n) + 0.5d0 * (z_b(n) - z_b(n -1)))       
                total_sediment_mass  = total_sediment_mass + (a_3(n) + 0.5d0 * (a_3(n) - a_3(n -1)))*dx
                
            end do
            
            bed_sediment_area = bed_sediment_area - 0.5 * dx * (z_b(1) + 0.5d0 * (z_b(1) - z_b(0))) -  dx * 0.5 * (z_b(n_steps + 1) + 0.5d0 * (z_b(n_steps + 1) - z_b(n_steps)))
            total_sediment_mass = total_sediment_mass  - 0.5d0 * ((a_3(1) - a_3(0)) + a_3(1))*dx - 0.50d0 * ((a_3(n_steps+1) - a_3(n_steps)) + a_3(n_steps+1))*dx + bed_sediment_area * rho_s * (1.0d0-bed_porosity(n_steps)) 
            
        end subroutine runge_kutta_coupled_calc
                                    
                                    
    
end module Runge_Kutta_coupled_calculation