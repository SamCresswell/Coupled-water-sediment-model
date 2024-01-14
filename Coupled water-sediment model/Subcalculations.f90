module subcalculations
    
    use initialisation
    use variables
    
    contains 
    
        !-------------------------------------------------------------------------
        !This subroutine calculates the spatial step length
        subroutine step_length(dx,x_tot,n_steps)
    
            implicit none
        
            double precision :: dx,x_tot
            integer :: n_steps
        
            dx = x_tot/n_steps
    
        end subroutine step_length
        
        !-------------------------------------------------------------------------
        !This subroutine calculates the time step length
        subroutine time_step(n_timesteps,dt,t_0,t_end)
    
            implicit none
        
            double precision :: dt,t_0,t_end
            integer :: n_timesteps
        
            dt = (t_end - t_0)/n_timesteps
    
        end subroutine time_step
        
        !-------------------------------------------------------------------------
        !This subroutine calculates the initial bed porosity
        subroutine bed_porosity_calculation(bed_porosity,d_median,n,n_steps,bed_porosity_input,calculate_bed_porosity)
    
            implicit none
        
            double precision :: bed_porosity_input
            integer :: n_steps, n
            logical :: calculate_bed_porosity
            double precision,dimension(-1:n_steps + 2) :: d_median,bed_porosity
        
            if (calculate_bed_porosity == .FALSE.) then
                
                bed_porosity(n) = bed_porosity_input
                
            else if (calculate_bed_porosity == .TRUE.) then
                
                bed_porosity(n) = 0.13d0 + 0.21d0/((100.0d0 * d_median(n) + 0.002d0) ** 0.21d0)
                
            end if
    
        end subroutine bed_porosity_calculation
        
        !-------------------------------------------------------------------------
        !This subroutine calculates the median grain size
        subroutine median_diameter(d_median,n,n_steps,d_median_input,d_median_equation,x)
    
            implicit none
        
            double precision :: d_median_input
            integer :: n_steps, n
            integer :: d_median_equation
            double precision,dimension(-1:n_steps + 2) :: d_median,x
        
            !Calculates the value of d_median at this point
            if (d_median_equation == 0) then
                
                d_median(n) = d_median_input
             
            !For the three-slope Karnali case
            else if (d_median_equation == 1) then  
                
                if (x(n) < 26000.0d0) then
                                    
                    d_median(n) = 0.0450d0
                                    
                else if (x(n) >= 26000.d0 .and. x(n) < 31000.d0) then
                                    
                    d_median(n) = 0.045d0 + (0.0035d0 - 0.045d0) * (x(n) - 26000.d0)/(5000.0d0)
                                    
                else if (x(n) >= 31000.d0) then
                                    
                    d_median(n) = 0.0035d0
                                
                end if
                
            !For the two-slope Karnali case
            else if (d_median_equation == 2) then  
                
                if (x(n) < 28500) then
                                    
                    d_median(n) = 0.0450d0
                                    
                else if (x(n) >= 28500) then
                                    
                    d_median(n) = 0.010d0
                                
                end if

                
            !For the three-slope Karnali case with no grain transition
            else if (d_median_equation == 3) then  
                
                if (x(n) < 26000.0d0) then
                                    
                    d_median(n) = 0.0450d0
                                    
                else if (x(n) >= 26000.d0 .and. x(n) < 31000.d0) then
                                    
                    d_median(n) = 0.0450d0
                    !d_median(n) = 0.0035 
                                    
                else if (x(n) >= 31000.d0) then
                                    
                    d_median(n) = 0.0035d0
                                
                end if
                
            !For the two-slope Karnali case with a grain size transition
            else if (d_median_equation == 4) then  
                
                if (x(n) < 28500.0d0) then
                                    
                    d_median(n) = 0.0450d0
                                    
                else if (x(n) >= 28500.d0) then
                                    
                    d_median(n) = 0.0035d0
                                
                end if
                
            else
                
                print*,'Invalid grain size equation! Stopping program.'
                stop

            end if

        end subroutine median_diameter
        
        !-------------------------------------------------------------------------
        !This subroutine calculates the saturated bed density
        subroutine saturated_bed_density(rho_0,rho_w,bed_porosity,rho_s,n,n_steps)
    
            implicit none
        
            double precision :: rho_s
            integer :: n, n_steps
            double precision,dimension(-1:n_steps + 2) :: bed_porosity, rho_0, rho_w
            
            rho_0(n) = rho_w(n) * bed_porosity(n) + rho_s * (1.0d0 - bed_porosity(n))
    
        end subroutine saturated_bed_density
        
        !-------------------------------------------------------------------------
        !This subroutine calculates the change in grain size if enabled
        subroutine alter_grain_size(d_median,n_steps,n,update_grain_size,z_b_new,z_b,&
                    z_b_initial,d_median_initial,angle_of_repose,bed_slope_angle, q_b, dx,d_median_new)
    
            implicit none
        
            !double precision :: rho_w,rho_s
            integer :: n, n_steps, update_grain_size
            double precision :: angle_of_repose, dx
            double precision,dimension(-1:n_steps + 2) :: d_median,z_b_new,z_b, z_b_initial,d_median_initial,bed_slope_angle, q_b,d_median_new
            
            !Stores the current value of all grain sizes so cells can be updated in a loop independently
            do n = -1, n_steps + 2
                
                d_median_new(n) = d_median(n)
                
            end do
            
            do n = 0, n_steps + 2
            
                !If the upstream bedload flux and deposition in a cell exceed the upstream grain diameter, updates it to that value
                if (q_b(n - 1) + D(n) > (d_median(n - 1))) then
                
                    if (d_median(n).ne.d_median(n - 1)) then
                        
                        d_median_new(n) = d_median(n - 1)
                        
                        !Includes an update message in the console for that cell
                        print*,'Cell at ',n,' was updated from median grain size of ',d_median(n),' to ',d_median_new(n)
                        
                    end if
                    
                end if
             
            end do 
            
            !Updates all cells to their new value
            do n = -1, n_steps + 2
                
                d_median(n) = d_median_new(n)
                 
            end do
    
        end subroutine alter_grain_size
        
        !-------------------------------------------------------------------------
        !This subroutine ramps up the density of water concentration based on the flowrate in each cell
        subroutine ramp_concentration(q, n, q_ramp_0, rho_w, rho_w_clear, rho_w_multiplier, pi, q_ramp_end)
    
            implicit none
        
            integer :: n 
            double precision :: q_ramp_0, rho_w_clear, rho_w_multiplier, pi, q_ramp_end
            double precision,dimension(-1:n_steps + 2) :: q, rho_w
            
            !Ramps the density of water to simulate a concentrated flow with a washload based on the flowrate
            if (q(n) < 2.0d0 * q_ramp_0) then
                    
                rho_w(n) = rho_w_clear
            else 
                    
                rho_w(n) = rho_w_clear + rho_w_clear * (rho_w_multiplier - 1.0d0) * (sin((pi * q(n))/(2 * q_ramp_end)))**2.0d0
                         
            end if
                            
        end subroutine ramp_concentration
    
        !-------------------------------------------------------------------------
        !This subroutine calculates the particle settling velocity and the relative density of water and sediment
        subroutine particle_transport_calcs(n, d_median_equation, use_deposition_equation, is_bed_fixed,&
                                            s,rho_s, g, d_median_input, rho_w_clear, &
                                            x, tau_bx, d_median, bed_slope_angle, &
                                            d_dimensionless, R_n, u, h,D,En, conc,bed_porosity, &
                                            Re_p, m_coefficient, dimensionless_critical_shields_parameter,rho_w, &
                                            Fr, shields_parameter,q, nu, rho, w_s_0, w_s_h, modified_critical_shields_parameter, &
                                            alpha_en, alpha_d, D_input, En_input, &
                                            u_shear, von_karman_constant, m_d, nu_input, angle_of_repose,&
                                            hyperconcentration_type, hyperconcentration_conc)
    
            implicit none
        
            integer :: n, d_median_equation, use_deposition_equation, is_bed_fixed, hyperconcentration_type
            double precision :: rho_s, g, d_median_input, rho_w_clear
            double precision,dimension(-1:n_steps + 2) :: x, tau_bx, d_median, bed_slope_angle
            double precision,dimension(-1:n_steps + 2) :: d_dimensionless, R_n, u, h,D,En, conc,bed_porosity, s
            double precision,dimension(-1:n_steps + 2) :: Re_p, m_coefficient, dimensionless_critical_shields_parameter,rho_w
            double precision,dimension(-1:n_steps + 2) :: Fr, shields_parameter,q, nu, rho, w_s_0, w_s_h, modified_critical_shields_parameter
            double precision :: alpha_en, alpha_d, D_input, En_input
            double precision :: u_shear, von_karman_constant, m_d, nu_input, angle_of_repose, hyperconcentration_conc
            
            !Calculates the relative density of the sediment to the sediment-water mixture
            s(n) = rho_s/rho_w(n) - 1.0d0   
              
            !Writes the dynamic viscosity of the fluid to an array
            nu(n) = nu_input
            
            !First, the particle settling velocity is computed
            w_s_0(n) = (s(n) * g * d_median(n) ** 2.0d0)/(18 * nu(n) + sqrt(0.75d0 * s(n) * g * d_median(n) ** 3.0d0))
        
            !Calculates the particle Reynolds number 
            Re_p(n) = w_s_0(n) * d_median(n) / nu(n)
            
            !Calculates the m exponent in the modified Richardson and Zaki equation depending on the flow regime, from the method in Baas et al. (2022)
            if (Re_p(n) <= 0.1d0) then
                
                m_coefficient(n) = 5.037d0
                
            else if (Re_p(n) > 0.1d0 .and. Re_p(n) <= 10.0d0) then
                
                m_coefficient(n) = 4.565d0 - 0.472d0 * log10(Re_p(n))
                
            else if (Re_p(n) > 10.0d0 .and. Re_p(n) <= 500.0d0) then
                
                m_coefficient(n) = 5.074d0 - 0.981d0 * log10(Re_p(n))
                
            else if (Re_p(n) >= 500.0d0) then
                
                m_coefficient(n) = 2.436
                
            end if
            
            !Calculates the hindered particle settling velocity
            w_s_h(n) = w_s_0(n) * (1.0d0 - conc(n)) ** m_coefficient(n)
            
            if (hyperconcentration_type == 1) then
                
                w_s_h(n) = w_s_0(n) * (1.0d0 - (conc(n) + hyperconcentration_conc*0.01)) ** m_coefficient(n)
                
            else if (hyperconcentration_type == 2) then
                
                w_s_h(n) = w_s_0(n) * (1.0d0 - (conc(n) + (rho(n) - rho_w(n))/(rho_s - rho_w(n)))) ** m_coefficient(n)
                
            end if
            
            !Calculates the dimensionless grain diameter
            d_dimensionless(n) = d_median(n) * (g * s(n)/(nu(n)**2.0d0))**(1.0d0/3.0d0)
            
            !Dimensionless critical particle settling bed stress
            dimensionless_critical_shields_parameter(n) = 0.3d0/(1.0d0 + 1.2d0 * d_dimensionless(n)) + 0.055d0 * (1 - exp(-0.02d0 * d_dimensionless(n)))
            
            !Calculates the modified critical Shield's parameter
            modified_critical_shields_parameter(n) = dimensionless_critical_shields_parameter(n) * (sin(bed_slope_angle(n) + angle_of_repose))/(sin(angle_of_repose))
            dimensionless_critical_shields_parameter(n) = modified_critical_shields_parameter(n)
            
            !Calculates the shear velocity
            u_shear = sqrt(tau_bx(n)/rho(n))
            
            !Calculates the Shields parameter
            shields_parameter(n) = tau_bx(n)/(g * d_median(n) * (rho_s - rho_w(n)))
            
            !Calculates the Rouse number
            R_n(n) = w_s_h(n)/(von_karman_constant * u_shear)
            
            !Sets the deposition value, using a constant or equation
            if (use_deposition_equation == 0) then
                
                D(n) = D_input

            else if (use_deposition_equation == 1) then
                
                D(n) = conc(n) * w_s_h(n)
                     
            else if (use_deposition_equation == 2) then
                
                !Calculates the depth-averaged concentration of suspended sediment to the bear-bed concentration
                alpha_d = min(2.0d0*d_median(n),(1 - bed_porosity(n))/conc(n))
                
                !Calculates the deposition from the unhindered settling velocity
                D(n) = w_s_0(n) * alpha_d * conc(n) * (1.0d0 - alpha_d * conc(n)) ** m_d
            
            !Calculates the deposition as a function of the unhindered settling velocity and concentration
            else if (use_deposition_equation == 3) then
                
                D(n) = w_s_0(n) * conc(n) 
           
            !Calculates deposition using a simplified version of the Richardson and Zaki equation
            else if (use_deposition_equation == 4) then
                
                m_coefficient(n) = (4.45d0 + 18 * d_median(n)/h(n)) * Re_p(n) ** (-0.1d0)
                
                D(n) = conc(n) * w_s_h(n) * (1 - conc(n)) ** m_coefficient(n) 
                
            end if
            
            !Sets the entrainment value, using a constant or equation
            !For simply using a fixed constant value
            if (use_entrainment_equation == 0) then
                
                En(n) = En_input
            
            !For the equation used in the horizontal evenly-mixed tank case
            else if (use_entrainment_equation == 1) then
                
                En(n) = En_input * ((1.0d0 - 0.2d0)/(0.2d0))
           
            !Calculates entrainment using the full method in Cao et al. (2006)
            else if (use_entrainment_equation == 2) then
                    
                if (shields_parameter(n) >= dimensionless_critical_shields_parameter(n)) then
                           
                    En(n) = 160.0d0/((sqrt(s(n) * g) /nu(n))**0.8d0) * (1.0d0 - bed_porosity(n))/dimensionless_critical_shields_parameter(n) * (shields_parameter(n) - dimensionless_critical_shields_parameter(n))/h(n) * d_median(n)** (-0.2d0) * (7.0d0 * u(n))/6.0d0 
                    
                else 
                        
                    En(n) = 0.0d0
                        
                end if
                
            !Calculates entrainment using the simplified method in Cao et al. (2006)
            else if (use_entrainment_equation == 3) then
                    
                if (shields_parameter(n) >= dimensionless_critical_shields_parameter(n)) then
                           
                    En(n) = alpha_en * (shields_parameter(n) - dimensionless_critical_shields_parameter(n)) * u(n) * h(n)**(-1.0d0) * d_median(n)**(-0.2d0)
                    
                else 
                        
                    En(n) = 0.0d0
                        
                end if
                
            !Calculates the entrainment as a function of the deposition coefficient alpha_d
            else if (use_entrainment_equation == 4) then
                
                if (shields_parameter(n) - dimensionless_critical_shields_parameter(n) <= 0.0d0) then
                
                    En(n) = 0.0d0
                
                else
                 
                    !Calculates the depth-averaged concentration of suspended sediment to the bear-bed concentration
                    alpha_d = h(n) / min(2.0d0 ,9.0d0 * d_median(n) * shields_parameter(n))
                        
                    En(n) = alpha_d * w_s_h(n) * conc(n)
                
                end if
           
            !Calculates the entrainment as a function of the Froude number
            else if (use_entrainment_equation == 6) then
                    
                !Calculates the depth-averaged concentration of suspended sediment to the bear-bed concentration
                
                if (shields_parameter(n) - dimensionless_critical_shields_parameter(n) <= 0) then
                
                    En(n) = 0.0d0
                
                else
                
                    Fr(n) = u(n)/sqrt(g * h(n))
                    En(n) = 0.000474d0 * (u_shear/w_s_0(n)) ** (1.77d0) * Fr(n) ** (1.18d0)
    
                end if
                    
            else
                        
                return        
                
            end if
            
            !If the bed is fixed, ensures that both bedload fluxes are set to zero
            if (is_bed_fixed == 1) then
                
                En(n) = 0.0d0
                D(n) = 0.0d0
            
            end if
                
        end subroutine particle_transport_calcs
        
        !-------------------------------------------------------------------------
        !This subroutine calculates the bed slope for a simple single-slope model
        subroutine bed_slope(s_0,z_b_end,z_b_0,x_tot)
    
            implicit none
        
            double precision :: s_0,z_b_end,z_b_0,x_tot 
        
            s_0 = (z_b_0 - z_b_end)/x_tot
    
        end subroutine bed_slope
        
        !-------------------------------------------------------------------------
        !This subroutine calculates the change in bed elevation at the start of each step in
        !the sediment-water Runge-Kutta model
        subroutine bed_elevation_change(upwinding_type, n, n_steps, is_bed_fixed,&
                                        q_b,dzb_dt,D,En,bed_porosity, bedrock_elevation, &
                                        minmod_a, minmod_b, minmod_ab, dx)
            implicit none
            
            integer :: upwinding_type, n, n_steps, is_bed_fixed
            double precision, dimension(-1:n_steps + 2) :: q_b,dzb_dt,D,En,bed_porosity, bedrock_elevation
            double precision, dimension(-1:n_steps + 2) :: minmod_a, minmod_b, minmod_ab
            double precision :: dx
     
            !Calculates the bed elevation using a simple central difference formula
            if (upwinding_type == 1) then
                    
                    dzb_dt(n) = (D(n) - En(n) - (q_b(n + 1) - q_b(n - 1))/(2.0d0 * dx)) * 1.0d0/(1.0d0 - bed_porosity(n))
                
            !Calculates the bed elevation using upwinding only (unstable)
            else if (upwinding_type == 2) then
                    
                dzb_dt(n) = (D(n) - En(n) - (3.0d0 * q_b(n) - 4.0d0 * q_b(n - 1) + q_b(n - 2))/(2.0d0 * dx)) * 1.0d0/(1.0d0 - bed_porosity(n))
              
            !Calculates the bed elevation using a combination of upwinding and the central difference only (unstable for a non-horizontal bed slope)
            else if (upwinding_type == 3) then
                    

                dzb_dt(n) = (D(n) - En(n) - 0.5d0 * (( q_b(n + 1) - q_b(n - 1))/(2.0d0 * dx) + (3.0d0 * q_b(n) - 4.0d0 * q_b(n - 1) + q_b(n - 2))/(2.0d0 * dx)))* 1.0d0/(1.0d0 - bed_porosity(n))
             
            !Calculates the bed elevation by applying the minmod function and a central difference
            else if (upwinding_type == 4) then
                    
                    !Calculates the minmod function slopes 'a' and 'b' around the point of in
                    minmod_a(n) = (q_b(n) - q_b(n - 1))/(dx)
                    minmod_b(n) = (q_b(n + 1) - q_b(n - 1))/(2.0d0 * dx)
                    
                    !Applies the minmod function, by checking if the two slopes are in the same direction
                    if (minmod_a(n) > 0.0d0 .and. minmod_b(n) > 0.0d0) then
                        
                        minmod_ab(n) = min(minmod_a(n),minmod_b(n))
                        
                    else if (minmod_a(n) < 0.0d0 .and. minmod_b(n) < 0.0d0) then
                        
                        minmod_ab(n) = max(minmod_a(n),minmod_b(n))
                        
                    else
                        minmod_ab(n) = 0.0d0

                    end if
                     
                    !Calculates the bed elevation change as a linear combination of the central difference and the limited slope. Adjust the coefficients as necessary.
                    dzb_dt(n) = (D(n) - En(n) - (0.3d0 * minmod_ab(n) + 0.6d0 * ( q_b(n + 1) - q_b(n - 1))/(2.0d0 *dx))) * 1.0d0/(1.0d0 - bed_porosity(n))
                    
                end if 
            
                !If bedrock excavation is enabled, ensures that there is no further entrainment
                !If the depth of excavation is limited, sets the new bed elevation to the elevation of the top of the bedrock if this value is exceeded
                if (z_b(n) + dzb_dt(n) <= bedrock_elevation(n) .AND. bedrock_elevation(n) .ne. -100.0d0) then
                
                    dzb_dt(n) = 0.0d0
                    
                else 
                
                    continue
                
                end if
                
                !Ensures that no bed movement occurs if a fixed bed has been selected
                if (is_bed_fixed == 1) then
                    
                    dzb_dt(n) = 0.0d0

                end if
                
        end subroutine bed_elevation_change
        !-------------------------------------------------------------------------
        !This subroutine calculates the geometry for the channel
        subroutine channel_geometry(s_0,z_b_0,z_b_end,x_tot,bedrock_elevation_west,bedrock_elevation_east,&
                                    x,z_b,bedrock_elevation, &
                                    n,n_steps,west_bc,east_bc,calculation_type, &
                                    bed_equation, bedrock_elevation_equation)
    
            implicit none
        
            double precision :: s_0,z_b_0,z_b_end,x_tot,bedrock_elevation_west,bedrock_elevation_east
            double precision, dimension(-1:n_steps + 2) :: x,z_b,bedrock_elevation
            integer :: n,n_steps,west_bc,east_bc,calculation_type
            integer :: bed_equation, bedrock_elevation_equation
            
            !Calculates the bed slope
            if (bed_equation == 0) then
                
                call bed_slope(s_0,z_b_end,z_b_0,x_tot)
                
            end if
            
            !Calculates the x-position for each step first
            !For the Adams-Bashford method
            if (calculation_type == 1) then
                
                do n = 1, n_steps+1
                
                    x(n) = dble(n - 1) * x_tot/(n_steps)
     
                end do
                
            !If the Runge-Kutta method is used for water only    
            else if (calculation_type == 2 .or. calculation_type == 3) then
                
                do n = -1,n_steps + 2
                
                    x(n) = dx * ((n) - 0.50d0)
            
                end do
                
            end if
                
            !Calculates the bed elevation from the step x-positions
            do n = -1,n_steps + 2
                
                !Handles the different initial bed elevation cases
                if (bed_equation == 0) then
                
                    z_b(n) = z_b_0 - (z_b_0 - z_b_end) * x(n)/x_tot
                        
                else if (bed_equation == 1) then
            
                    
                            
                    if (x(n) < 300.0d0 .OR. x(n) > 500.0d0) then
                                
                        z_b(n) = 0.0d0
                                
                    else 
                                    
                        z_b(n) = 1.0d0 * (sin((pi * (x(n) - 300.0d0))/(500.0d0-300.0d0)))**2.0d0
                                    
                    end if
       
                !For the deposition and entrainment test cases
                else if (bed_equation == 2) then
                            
                    z_b(n) = 1.0d0
                            
                !For the Karnali case
                else if (bed_equation == 3) then
                            
                    if (x(n) < 26000) then
                                    
                        z_b(n) = 78.5d0  - x(n)/26000.0d0 * (78.5d0 - 26.5d0)
                                    
                    else if (x(n) >= 26000 .and. x(n) < 31000) then
                                    
                        z_b(n) = 26.5d0  - (x(n) - 26000.0d0)/5000.0d0 * (26.5d0 - 19.0d0)
                                    
                    else if (x(n) >= 31000) then
                                    
                        z_b(n) = 19.0d0  - (x(n) - 31000.0d0)/19000.0d0 * 19.0d0
                                
                    end if
                                
                !For the two-slope case
                else if (bed_equation == 4) then
                            
                    if (x(n) < 28500.0d0) then
                                    
                        z_b(n) = 78.5d0  - x(n)/28500.0d0 * (57.0d0)
                                    
                    else if (x(n) >= 28500.0d0) then
                                    
                        z_b(n) = 21.5d0  - (x(n) - 28500.0d0)/21500.0d0 * (21.5d0)
                                
                    end if
                                
                else if (bed_equation == 5) then
                            
                    if (x(n) < 20000.0d0) then
                                
                        z_b(n) = 40.0d0 - x(n)/20000.0d0 * 30.0d0
                                
                    else 
                                    
                        z_b(n) = 10.0d0 - (x(n) - 20000.0d0) * 10.0d0/10000.0d0
                                    
                    end if
                    
                else
                        
                    print*,'Invalid bed equation! Stopping program.'
                            
                    stop
                            
                end if
                 
            end do

            !Calculates the bedrock equation
            do n = -1,n_steps + 2
                
                if (bedrock_elevation_equation == 0) then
                    
                    bedrock_elevation(n) = -100.0d0
                    
                    continue
                    
                !For bedrock of constant height
                else if (bedrock_elevation_equation == 1) then
                    
                    bedrock_elevation(n) = bedrock_elevation_west

                !For a linear bedrock slope
                else if (bedrock_elevation_equation == 2) then
                    
                    bedrock_elevation(n) = bedrock_elevation_west - x(n)/x_tot * (bedrock_elevation_west - bedrock_elevation_east)
                    
                end if
                
            end do
            
            !Applies the upstream boundary conditions to the bed elevation
            if (west_bc == 5) then
            
                    z_b(0) = 2.0d0 * z_b(1) - z_b(2)
                    z_b(-1) = 2.0d0 * z_b(0) - z_b(1)
                    
            else if (west_bc == 2) then
            
                    z_b(0) = 2.0d0 * z_b(1) - z_b(2)
                    z_b(-1) = 2.0d0 * z_b(0) - z_b(1)
                    
            end if
            
            !Applies the downstream boundary conditions to the bed elevation
            if (east_bc == 5) then
                
                    z_b(n_steps + 1) = 2.0d0 * z_b(n_steps) - z_b(n_steps - 1)
                    z_b(n_steps + 2) = 2.0d0 * z_b(n_steps + 1) - z_b(n_steps)
                    
            else if (east_bc == 2) then
                
                    z_b(n_steps + 1) = 2.0d0 * z_b(n_steps) - z_b(n_steps - 1)
                    z_b(n_steps + 2) = 2.0d0 * z_b(n_steps + 1) - z_b(n_steps)
                    
            end if
    
        end subroutine channel_geometry
        
        !-------------------------------------------------------------------------
        !This subroutine handles initial parameters
        subroutine initial_parameters(n,n_steps,west_bc,east_bc,calculation_type, bed_equation, &
                                        h,q,q_b_0,conc,free_surface_elevation,rho,rho_w, &
                                        dh_dt_1,dh_dt_2,dh_dt_3,dh_dt_4,dq_dt_1,dq_dt_2,dq_dt_3,dq_dt_4,dh_dt_a,dq_dt_a, &
                                        a_1,a_2,a_3,z_b, &
                                        h_input,q_input,dq_dt_b,dh_dt_b,q_b_0_input,rho_s,conc_0)
        
            implicit none
            
            integer :: n,n_steps,west_bc,east_bc,calculation_type, bed_equation
            double precision,dimension(-1:n_steps+2) :: h,q,q_b_0,conc,free_surface_elevation,rho,rho_w
            double precision,dimension(-1:n_steps+2) :: dh_dt_1,dh_dt_2,dh_dt_3,dh_dt_4,dq_dt_1,dq_dt_2,dq_dt_3,dq_dt_4,dh_dt_a,dq_dt_a
            double precision,dimension(-1:n_steps+2) :: a_1,a_2,a_3,z_b
            double precision :: h_input,q_input,dq_dt_b,dh_dt_b,q_b_0_input,rho_s,conc_0

            if (calculation_type == 1) then
                
                do n = 1, n_steps
                    
                    !Sets the initial specific flowrate
                    q(n) = q_input
                    
                    !Sets the input flow depth
                    if (bed_equation == 1) then
                            
                            h(n) = 10.0d0 - z_b(n)
                        
                    else 
                    
                        h(n) = h_input
                        
                    end if
                    
                    !Calculates the flow velocity for all cells
                    u(n) = q(n)/h(n)
                
                end do
            
                !Adds the values outside the domain for the case with uniform flow at the upstream boundary
                !Handles the uniform flow case
                if (west_bc == 5) then
                
                    h(0) = 2.0d0 * h(1) - h(2)
                    h(-1) = 2.0d0 * h(0) - h(1)
                
                    q(0) = q(1)
                    q(-1) = q(0)
                
                    z_b(0) = 2.0d0 * z_b(1) - z_b(2)
                    z_b(-1) = 2.0d0 * z_b(0) - z_b(1)   
                
                !Handles the permissive flow case
                else if (west_bc == 2) then
                    
                    h(1) = 2.0d0 * h(2) - h(3)
                    q(1) = 2.0d0 * q(2) - q(3)
        
                    else
                
                    end if
        
                !Adds the values outside the domain for the downstream boundary
                !case with uniform flow at the upstream boundary
                if (east_bc == 5) then
                
                    h(n_steps + 1) = 2.0d0 * h(n_steps) - h(n_steps - 1)
                    h(n_steps + 2) = 2.0d0 * h(n_steps + 1) - h(n_steps)
                
                    q(n_steps + 1) = q(n_steps)
                    q(n_steps + 2) = q(n_steps + 1)
                
                    z_b(n_steps + 1) = 2.0d0 * z_b(n_steps) - z_b(n_steps - 1)
                    z_b(n_steps + 2) = 2.0d0 * z_b(n_steps + 1) - z_b(n_steps)
                
                !Handles the permissive flow case
                else if (east_bc == 2) then
                
                    h(n_steps) = 2.0d0 * h(n_steps - 1) - h(n_steps + 2) 
                    q(n_steps) = 2.0d0 * q(n_steps - 1) - q(n_steps + 2) 
                
                end if
            
                !Sets the initial values of the derivatives
                do n = 0,n_steps + 1
                
                        dh_dt_a(n) = 0.0d0
                        dh_dt_b = 0.0d0
                        dq_dt_a(n) = 0.0d0
                        dq_dt_b = 0.0d0
                    
                end do
            
            end if
            
            !For the Runge-Kutta method without sediment transport
            if (calculation_type == 2) then
                
                do n = -1, n_steps + 2
                q(n) = q_input

                
                !Sets the input flow depth
                if (bed_equation == 1) then
                            
                    h(n) = 10.0d0 - z_b(n)
                        
                else 
                    
                    h(n) = h_input
                        
                end if
                  
                !Sets the value of the velocity for all cells
                u(n) = q(n)/h(n)
                
                end do
                
                !Adds the values outside the domain for the case with uniform flow at the upstream boundary
                !Handles the uniform flow case
                if (west_bc == 5) then
                
                    h(0) = 2.0d0 * h(1) - h(2)
                    h(-1) = 2.0d0 * h(0) - h(1)
                
                    q(0) = q(1)
                    q(-1) = q(0)
                
                    z_b(0) = 2.0d0 * z_b(1) - z_b(2)
                    z_b(-1) = 2.0d0 * z_b(0) - z_b(1)   
                
                !Handles the permissive flow case
                else if (west_bc == 2) then
                    
                    h(1) = 2.0d0 * h(2) - h(3)
                    h(n_steps) = 2.0d0 * h(n_steps - 1) - h(n_steps + 2) 
                    
                    q(0) = q(1)
                    q(-1) = q(0)
                    
                    z_b(n_steps + 1) = 2.0d0 * z_b(n_steps) - z_b(n_steps - 1)
                    z_b(n_steps + 2) = 2.0d0 * z_b(n_steps + 1) - z_b(n_steps)
                    
                else
                
                end if
        
                !Adds the values outside the domain for the downstream boundary
                !case with uniform flow at the upstream boundary
                if (east_bc == 5) then
                
                    h(n_steps + 1) = 2.0d0 * h(n_steps) - h(n_steps - 1)
                    h(n_steps + 2) = 2.0d0 * h(n_steps + 1) - h(n_steps)
                
                    q(n_steps + 1) = q(n_steps)
                    q(n_steps + 2) = q(n_steps + 1)
                
                    z_b(n_steps + 1) = 2.0d0 * z_b(n_steps) - z_b(n_steps - 1)
                    z_b(n_steps + 2) = 2.0d0 * z_b(n_steps + 1) - z_b(n_steps)
                
                !Handles the permissive flow case
                else if (east_bc == 2) then
                
                    h(n_steps + 1) = 2.0d0 * h(n_steps) - h(n_steps - 1)
                    h(n_steps + 2) = 2.0d0 * h(n_steps + 1) - h(n_steps)
                
                    q(n_steps + 1) = q(n_steps)
                    q(n_steps + 2) = q(n_steps + 1)
                
                    z_b(n_steps + 1) = 2.0d0 * z_b(n_steps) - z_b(n_steps - 1)
                    z_b(n_steps + 2) = 2.0d0 * z_b(n_steps + 1) - z_b(n_steps)
                
                end if
                
                !Sets the initial fluxes to zero
                do n = -1,n_steps + 2
                
                    dh_dt_1(n) = 0.0d0
                    dh_dt_2(n) = 0.0d0
                    dh_dt_3(n) = 0.0d0
                    dh_dt_4(n) = 0.0d0
                    dq_dt_1(n) = 0.0d0
                    dq_dt_2(n) = 0.0d0
                    dq_dt_3(n) = 0.0d0
                    dq_dt_4(n) = 0.0d0
                
                end do
                
            end if
            
            !Handles the a-parameters if the sediment Runge-Kutta model is deployed
            if (calculation_type == 3) then
                
                !Sets the initial density of the water-sediment mixture
                do n = -1, n_steps + 2
                    
                    rho(n) = rho_w(n) + conc_0 * (rho_s - rho_w(n))
                    
                end do
               
                !Iterates over all spatial steps to set the parameters to their initial values
                do n = -1, n_steps + 2
                    
                    !Sets the initial specific flowrate
                    q(n) = q_input
                    
                    !Sets the bedload transport rate over all cells
                    q_b_0(n) = q_b_0_input
                    
                    !Sets the input flow depth
                    !For the sandbar under steady flow case:
                    if (bed_equation == 1) then
                            
                        h(n) = 10.0d0 - z_b(n)
                            
                    else 
                        
                        h(n) = h_input
  
                    end if
                    
                    !Calculates the flow velocity
                    u(n) = q(n)/h(n)

                    !Sets the value of coefficient 'a_1'
                    a_1(n) = rho(n) * h(n)
                    !Sets the value of the coefficient 'a_2' to an initial zero value
                    a_2 = a * rho(n) * h(n) * u(n)
                    !Sets the value of the coefficient 'a_3' to an initial zero value
                    a_3(n) = rho_s * h(n) * conc_0
                    !Sets the concentration to an initial constant
                    conc(n) = conc_0

                end do
                
                !Now sets the boundary conditions
                !Handles uniform flow
                if (west_bc == 5) then
                    
                    a_1(0) = 2.0d0 * a_1(1) - a_1(2)
                    a_1(-1) = 2.0d0 * a_1(0) - a_1(1)
                    a_2(0) = a_2(1)
                    a_2(-1) = a_2(0)
                    a_3(0) = 2.0d0 * a_3(1) - a_3(2)
                    a_3(-1) = 2.0d0 * a_3(0) - a_3(1)
                    h(0) = 2.0d0 * h(1) - h(2)
                    h(-1) = 2.0d0 * h(0) - h(1)
                    q_b_0(0) = q_b_0(1)
                    q_b_0(-1) = q_b_0(0)
                    q(0) = q(1)
                    q(-1) = q(0)
                
                !Handles the transmissive condition
                else if (west_bc == 2) then
                    
                    a_1(0) = 2.0d0 * a_1(1) - a_1(2)
                    a_1(-1) = 2.0d0 * a_1(0) - a_1(1)
                    a_2(0) = 2.0d0 * a_2(1) - a_2(2)
                    a_2(-1) = 2.0d0 * a_2(0) - a_2(1)
                    a_3(0) = 2.0d0 * a_3(1) - a_3(2)
                    a_3(-1) = 2.0d0 * a_3(0) - a_3(1)
                    
                    h(0) = 2.0d0 * h(1) - h(2)
                    h(-1) = 2.0d0 * h(0) - h(1)
                    
                    q_b_0(0) = q_b_0(1)
                    q_b_0(-1) = q_b_0(0)
                    
                    q(0) = q(1)
                    q(-1) = q(0)
                    
                else
                    
                    return
                    
                end if
                
                if (east_bc == 5) then
                    
                    a_1(n_steps + 1) = 2.0d0 * a_1(1) - a_1(2)
                    a_1(n_steps + 2) = 2.0d0 * a_1(0) - a_1(1)
                    a_2(n_steps + 1) = a_2(n_steps)
                    a_2(n_steps + 2) = a_2(n_steps + 1)
                    a_3(n_steps + 1) = 2.0d0 * a_3(n_steps) - a_3(n_steps - 1)
                    a_3(n_steps + 2) = 2.0d0 * a_3(n_steps + 1) - a_3(n_steps)
                    h(n_steps + 1) = 2.0d0 * h(n_steps) - h(n_steps - 1)
                    h(n_steps + 2) = 2.0d0 * h(n_steps + 1) - h(n_steps)
                    q_b_0(n_steps + 1) = q_b_0(n_steps)
                    q_b_0(n_steps + 2) = q_b_0(n_steps + 1)
                    q(n_steps + 1) = q(n_steps)
                    q(n_steps + 2) = q(n_steps + 1)
                    
                else if (east_bc == 2) then
                    
                    a_1(n_steps + 1) = 2.0d0 * a_1(n_steps) - a_1(n_steps - 1)
                    a_1(n_steps + 2) = 2.0d0 * a_1(n_steps + 1) - a_1(n_steps)
                    a_2(n_steps + 1) = 2.0d0 * a_2(n_steps) - a_2(n_steps - 1)
                    a_2(n_steps + 2) = 2.0d0 * a_2(n_steps + 1) - a_2(n_steps)
                    a_3(n_steps + 1) = 2.0d0 * a_3(n_steps) - a_3(n_steps - 1)
                    a_3(n_steps + 2) = 2.0d0 * a_3(n_steps + 1) - a_3(n_steps)
                    
                    h(n_steps + 1) = 2.0d0 * h(n_steps) - h(n_steps - 1)
                    h(n_steps + 2) = 2.0d0 * h(n_steps + 1) - h(n_steps)
                    
                    q_b_0(n_steps + 1) = 2.0d0 * q_b_0(n_steps) - q_b_0(n_steps - 1)
                    q_b_0(n_steps + 2) = 2.0d0 * q_b_0(n_steps + 1) - q_b_0(n_steps)
                    
                    q(n_steps + 1) = q(n_steps)
                    q(n_steps + 2) = q(n_steps + 1)
                    
                end if
                
                !Updates the output parameters not used in the next step
                conc(0) = 2.0d0 * conc(1) - conc(2)
                conc(-1) = 2.0d0 * conc(1) - conc(0)
                conc(n_steps + 1) = 2.0d0 * conc(n_steps) - conc(n_steps - 1)
                conc(n_steps + 2) = 2.0d0 * conc(n_steps + 1) - conc(n_steps)
            
                R_n(0) = 2.0d0 * R_n(1) - R_n(2)
                R_n(-1) = 2.0d0 * R_n(1) - R_n(0)
                R_n(n_steps + 1) = 2.0d0 * R_n(n_steps) - R_n(n_steps - 1)
                R_n(n_steps + 2) = 2.0d0 * R_n(n_steps + 1) - R_n(n_steps)
                
            end if
            
            !Calculates the free surface elevation
            do n = -1, n_steps+2
                
                free_surface_elevation = h(n) + z_b(n)
                
            end do
            
        end subroutine initial_parameters
        
        !---------------------------------------------------------------------------
        !This subroutine calculates the friction coefficient using one of four methods
        subroutine friction_coefficient(bed_friction_type,n,n_steps,c_f,c,c_b,n_mannings,h,g,c_b_input,&
                                        c_input,n_mannings_input,h_input,h_lim,d_median, von_karman_constant, z_0s)

            implicit none
    
            integer :: bed_friction_type,n,n_steps 
            double precision,dimension(-1:n_steps + 2) :: c_f,c,c_b,n_mannings,h, z_0s, d_median  
            double precision :: g,c_b_input,c_input,n_mannings_input,h_input,h_lim, von_karman_constant
    
                if (bed_friction_type == 1) then
        
                    c_b(n) = c_b_input
        
                    c_f(n) = c_b(n)
        
                else if (bed_friction_type == 2) then
        
                    c(n) = c_input
            
                    c_f(n) = g/c(n)**2.0d0
        
                !Uses Manning's equation
                else if (bed_friction_type == 3.and.h(n) >= h_lim) then
                    
                    n_mannings(n) = n_mannings_input
        
                    c_f(n) = g * n_mannings(n)**2.0d0/(h(n)**(1.0d0/3.0d0))
                    
                else if (bed_friction_type == 3.and.h(n) < h_lim) then
        
                    c_f(n) = 0.0d0
                    
                !Calculates the friction coefficient using roughness length from particle median diameter
                else if (bed_friction_type == 4) then
                    
                    if (d_median(n) < 0.004d0) then
                    
                        z_0s(n) = (2.5d0 * d_median(n))/(30.0d0)
                        
                    else if (d_median(n) >= 0.004d0) then
                        
                        z_0s(n) = (6.8d0 * d_median(n))/(30.0d0)
                    
                    end if
                    
                    c_f(n) = (von_karman_constant/(1 + log(z_0s(n)/h(n))))**2.0d0
                    
                else 
                    
                end if
    
        end subroutine friction_coefficient
        
        !-------------------------------------------------------------------------
        !This subroutine is a modified version of the previous one which does not operate via an array
        subroutine friction_coefficient_1D(bed_friction_type,n,n_steps,c_f_1,c,c_b,n_mannings,h,g,c_b_input,c_input,n_mannings_input,h_lim,h_f)

            implicit none
    
            integer :: bed_friction_type,n,n_steps
            
            double precision,dimension(-1:n_steps + 2) :: c,c_b,n_mannings,h
            
            double precision :: g,c_b_input,c_input,n_mannings_input,h_lim,c_f_1,h_f
            
                h_f = h(n)

                if (bed_friction_type == 1) then
        
                    c_b(n) = c_b_input
                    c_f_1 = c_b(n)
        
                else if (bed_friction_type == 2) then
        
                    c(n)= c_input
                    c_f_1 = g/c(n)**2.0d0
        
                else if (bed_friction_type == 3.and.h_f >= h_lim) then
            
                    n_mannings(n) = n_mannings_input
                    c_f = g * n_mannings_input**2.0d0/(h**(4.0d0/3.0d0))
                    
                else if (bed_friction_type == 3.and.h_f < h_lim) then
        
                    n_mannings(n) = n_mannings_input
                    c_f_1 = 0.0d0
                    
                else 
                    
                    print*,'Invalid bed friction condition! h_f is', h_f
            
                    call exit()
    
                end if
    
        end subroutine friction_coefficient_1D
        
        !-------------------------------------------------------------------------
        !This subroutine allows the upstream input discharge to be 'ramped up' to simulate a flood 
        subroutine ramp_up(n, n_steps,ramp_up_type,hyperconcentration_type, &
                            t,t_ramp,q_ramp_end,q_ramp_0,q_input, rho_w_clear, q_fixed, rho_w_multiplier, &
                            q, h, rho_w)

            implicit none
    
            integer :: n, n_steps,ramp_up_type,hyperconcentration_type
            double precision :: t,t_ramp,q_ramp_end,q_ramp_0,q_input, rho_w_clear, q_fixed, rho_w_multiplier
            double precision,dimension(-1:n_steps + 2) :: q, h, rho_w
    
            !Behaviour changes depending on whether the ramp-up time considers wind only, inflows only or neither
            q_fixed = q_input
    
            !Considers the case where there is ramping for the inflow
            if (ramp_up_type == 1) then
            
                if (t < t_ramp_0) then
                    
                    q_fixed = q_ramp_0
                    
                else if (t >= t_ramp_0 + t_ramp) then
                    
                    q_fixed = q_ramp_end
                    
                else 
                   
                    q_fixed = q_ramp_0 + (q_ramp_end - q_ramp_0) * (t - t_ramp_0)/t_ramp
                    
                end if
                
                
            !Deploys a custom ramping equation
            else if (ramp_up_type == 2) then
                
                     if (t < t_ramp_0) then
                         
                        q_fixed = q_ramp_0
                            
                     else if (t >= t_ramp_0 .and. t < t_ramp_0 + t_ramp) then
                    
                        q_fixed = q_ramp_0 +  (q_ramp_end - q_ramp_0) * (sin((pi * (t - t_ramp_0))/((t_ramp_0 + t_ramp) - t_ramp_0)))**2.0d0

                    else if (t >= t_ramp_0 + t_ramp) then
                        
                        q_fixed = q_input
   
                    end if   

            end if
    
        end subroutine ramp_up
                            
        !-------------------------------------------------------------------------
        !This subroutine is a modified version of the previous one which does not operate via an array
        subroutine bedload_intermediate_calc(bedload_method, q_b, n, n_steps, u_b, coeff_1, coeff_2, s, rho_s,&
                                                rho_w,rho, conc, tau_bx, u_1, c_f, shields_parameter, g, d_median,&
                                                dimensionless_critical_shields_parameter,dq_b_dt, q_b_0, dt, coeff_3, &
                                                is_bed_fixed, end_program_at_mobilisation,simulation_batch_name, &
                                                hyperconcentration_conc)
        
            integer :: bedload_method, n , n_steps, is_bed_fixed, end_program_at_mobilisation
            double precision,dimension(-1:n_steps + 2) :: q_b, u_b, rho, rho_w, s, u_1, tau_bx, conc, shields_parameter
            double precision,dimension(-1:n_steps + 2) :: dimensionless_critical_shields_parameter, q_b_0, dq_b_dt, c_f, d_median
            double precision :: coeff_1, coeff_2, coeff_3, rho_s, g, dt, hyperconcentration_conc
            character(len = 5) :: simulation_batch_name
            
                !Uses the grass formula 
                if (bedload_method == 1) then
                    
                    q_b(n) = coeff_1 * u_1(n) **3.0d0
                    u_b(n) = coeff_2 * u_1(n)
                    
                !Uses the Meyer-Peter-Mueller formula
                else if (bedload_method == 2) then
                    
                !Calculates the 's' parameter for simplification of writing equations
                s(n) = rho_s/rho_w(n) - 1.0d0  
              
                !Calculates the friction coefficient
                rho(n) = rho_w(n) * (1.0d0 - conc(n)) + rho_s * conc(n)
                tau_bx(n) = rho(n) * abs(u_1(n)) * u_1(n) * c_f(n)  
                shields_parameter(n) = tau_bx(n)/(g * d_median(n) * (rho_s - rho_w(n)))
                    
                    !Ensures no bed movement occurs if the Shields parameter does not exceed the critical value
                    if (shields_parameter(n) - dimensionless_critical_shields_parameter(n) <= 0 .or. is_bed_fixed == 1) then
                        
                        q_b(n) = 0.0d0
                        u_b(n) = 0.0d0
                        
                    else
                        
                        !Ends the program at the first bed mobilisation if applicable and plays a sound
                        if (end_program_at_mobilisation == 1) then
                            
                            print*,'Threshold reached!!'
                            print*,'n,q(n), u(n),u_shear(n):,hyp_conc(n\)'
                            print*,n,q(n),u(n),sqrt(tau_bx(n)/rho(n)),hyperconcentration_conc
                            write(23, '(F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10,",",A)') q(n),u(n),sqrt(tau_bx(n)/rho(n)),tau_bx(n),shields_parameter(n),hyperconcentration_conc,simulation_batch_name
                            call system('powershell (New-Object Media.SoundPlayer "..\Sounds\Splash.wav").PlaySync()')
                            stop
                            
                        end if
                    
                        !Calculates the bedload flux and bed particle velocity
                        q_b(n) = coeff_3 * 8.0d0 * sqrt((s(n)) * g * d_median(n)**3.0d0) * (shields_parameter(n) - dimensionless_critical_shields_parameter(n)) ** 1.5d0 
                        u_b(n) = q_b(n)/(10.0d0 * shields_parameter(n) * d_median(n))
                        
                    end if
            
                end if
                
                !Calculates the change in the bedload flux
                dq_b_dt(n) = (q_b(n) - q_b_0(n))/(dt * 2.0d0)
                
        end subroutine bedload_intermediate_calc
            
        !-------------------------------------------------------------------------
        !This subroutine calculates the Courant number and tests convergence
        
        subroutine courant_number_check(courant_number,u,dt,dx,t,mean_courant_number,has_exceeded_value)
        
        implicit none
        
        double precision, dimension(-1:n_steps + 2) :: courant_number, u 
        double precision :: dt,dx,t
        double precision, dimension(0:n_timesteps) :: mean_courant_number
        integer :: has_exceeded_value
        
        !Calculates the Courant number for all cells and prints a console warnign if it exceeds 1.0
        do n = -1, n_steps + 2
    
            courant_number(n) = u(n) * dt/dx
            
            if (courant_number(n) >= 1.0d0 .and. has_exceeded_value .ne. 1) then
                
                print '(A)'
                print*, 'WARNING: COURANT NUMBER HAS EXCEEDED 1!'
                print '(A)'
                
            end if
            
        end do
        
        end subroutine courant_number_check
        
        !-------------------------------------------------------------------------
        !This subroutine calculates the mean value of each parameter over each timestep
        
        subroutine mean_values(timestep, n_steps, n_timesteps, &
                                q, h, u, conc, q_b, u_b, En, D, R_n, courant_number, &
                                mean_q, mean_h, mean_u, mean_conc, mean_q_b, mean_u_b, mean_En, mean_D, mean_R_n, mean_courant_number)
        
        implicit none
        
        integer :: n_steps, timestep, n_timesteps
        double precision, dimension(-1:n_steps + 2) :: q, h, u, conc, q_b, u_b, En, D, R_n, courant_number
        double precision, dimension(0:n_timesteps) :: mean_q, mean_h, mean_u, mean_conc, mean_q_b, mean_u_b, mean_En, mean_D, mean_R_n, mean_courant_number
        
        !Calculates the mean value of each parameter for all cells in the domain
        mean_q(timestep) = sum(q) / dble(n_steps + 4)
        mean_h(timestep) = sum(h) / dble(n_steps + 4)
        mean_u(timestep) = sum(u) / dble(n_steps + 4)
        mean_conc(timestep) = sum(conc) / dble(n_steps + 4)
        mean_q_b(timestep) = sum(q_b) / dble(n_steps + 4)
        mean_u_b(timestep) = sum(u_b) / dble(n_steps + 4)
        mean_En(timestep) = sum(En) / dble(n_steps + 4)
        mean_D(timestep) = sum(D) / dble(n_steps + 4)
        mean_R_n(timestep) = sum(R_n) / dble(n_steps + 4)
        mean_courant_number(timestep) = sum(courant_number) / dble(n_steps + 4)
        
        end subroutine mean_values
        
    end module subcalculations