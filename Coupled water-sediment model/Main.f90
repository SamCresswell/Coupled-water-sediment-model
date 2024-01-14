program main
    
    !This program uses four modules for ease of reading
    use initialisation
    use variables
    use subcalculations
    use Runge_Kutta_calculation
    use Runge_Kutta_coupled_calculation
    use adams_bashford
    
    implicit none
    
    !Calls the initial CPU time to allow later assessment of program efficiency
    call cpu_time(t_cpu_0)
    
    !Reads the input data csv file, in which ALL parameters except the number of steps are stored
    call read_input_data
    
    !Calculates dx and dt
    call step_length(dx,x_tot,n_steps)
    call time_step(n_timesteps,dt,t_0,t_end)

    
    !Converts the angle of repose to a value in radians
    angle_of_repose = pi/180.0d0 * angle_of_repose
    
    !Calls the subroutine which deals with the channel geometry, generating the x-positions, bed
    !slope and elevations
    call channel_geometry(s_0,z_b_0,z_b_end,x_tot,bedrock_elevation_west,bedrock_elevation_east,&
                                    x,z_b,bedrock_elevation, &
                                    n,n_steps,west_bc,east_bc,calculation_type, &
                                    bed_equation, bedrock_elevation_equation)
    
    !Calculates the bed median diameters, the bed porosity and the saturated bed density
    do n = -1, n_steps + 2
        
        !Calculates the density of water or its multiplier if a hyperconcentrated flow is being used
        if (hyperconcentration_type == 0) then
            
            rho_w(n) = rho_w_clear
            
        else if (hyperconcentration_type == 1) then
            
            rho_w(n) = rho_w_clear * (1.0d0 - hyperconcentration_conc * 0.01d0) + rho_s * hyperconcentration_conc * 0.01d0
            
        else if (hyperconcentration_type == 2) then
            
            rho_w_multiplier = (rho_w_clear * (1.0d0 - hyperconcentration_conc * 0.01d0) + rho_s * hyperconcentration_conc * 0.01d0)/(rho_w_clear)
            rho_w(n) = rho_w_clear
            
        end if
           
        !Calls subroutines to calculate the median diameter in each cell, saturated bed density and the bed porosity
        call median_diameter(d_median,n,n_steps,d_median_input,d_median_equation,x)
        call bed_porosity_calculation(bed_porosity,d_median,n,n_steps,bed_porosity_input,calculate_bed_porosity)
        call saturated_bed_density(rho_0,rho_w,bed_porosity,rho_s,n,n_steps)
        
    end do
    
    !Saves the initial grain size value
    d_median_initial = d_median
    
    !Sets the initial flow and depth for all cells to the global input value, and sets the initial surface stress 
    call initial_parameters(n,n_steps,west_bc,east_bc,calculation_type, bed_equation, &
                                        h,q,q_b_0,conc,free_surface_elevation,rho,rho_w, &
                                        dh_dt_1,dh_dt_2,dh_dt_3,dh_dt_4,dq_dt_1,dq_dt_2,dq_dt_3,dq_dt_4,dh_dt_a,dq_dt_a, &
                                        a_1,a_2,a_3,z_b, &
                                        h_input,q_input,dq_dt_b,dh_dt_b,q_b_0_input,rho_s,conc_0)
    
    !Calculates the mean value of each parameter over all spatial steps for that timestep
    call mean_values(timestep, n_steps, n_timesteps, &
                            q, h, u, conc, q_b, u_b, En, D, R_n, courant_number, &
                            mean_q, mean_h, mean_u, mean_conc, mean_q_b, mean_u_b, mean_En, mean_D, mean_R_n, mean_courant_number)
    
    !Opens an output file for the flowrate
    open (1, file='output_flowrate.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the flow depth
    open (2, file='output_depth.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the fluid velocity
    open (3, file='output_velocity.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output fil for exporting data about the simulation which may be read by the python visualisation
    open (4, file='output_data.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the Courant number
    open (5, file='output_courant_number.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the bed elevation
    open (6, file='output_bed_elevation.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the volumetric concentration
    open (7, file='output_concentration.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the deposition
    open (8, file='output_deposition.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the entrainment flux
    open (9, file='output_entrainment.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the Rouse number
    open (10, file='output_Rouse_number.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the free surface elevation
    open (11, file='output_free_surface_elevation.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the bedload flux
    open (12, file='output_bed_flowrate.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the bed particle velocity
    open (13, file='output_bed_particle_velocity.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the initial and final bed placements
    open (14, file='output_bed_elevation_comparison.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the change in depth
    open (15, file='output_bed_elevation_change.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the median diameter
    open (16, file='output_median_diameter.csv', status ='replace',action = 'write',access = 'append')
    !Opens a tenth output file which contains the erosion
    open (17, file='output_erosion_depth.csv', status ='replace',action = 'write',access = 'append')
    !Opens a tenth output file which contains the deposition depth
    open (18, file='output_deposition_depth.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the 'water' density
    open (19, file='output_water_density.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the mean values of each parameter
    open (20, file='output_mean_values.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the instantaneous value of each parameter
    open (21, file='output_instantaneous_values.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for the simulation averaged results
    open (22, file ='output_simulation_results.csv', action ='write', status='old', access ='append')
    !Opens an output file for the simulation averaged results
    open (23, file ='output_bed_mobilisation_results.csv', action ='write', status='old', access ='append')
    !Opens an output file for Froude number
    open (24, file='output_Froude_number.csv', status ='replace',action = 'write',access = 'append')
    !Opens an output file for free surface comparison
    open (25, file='output_free_surface_comparison.csv', status ='replace',action = 'write',access = 'append')
    
    !Creates an initial line in the file for the time-series for each output file
    write(1, '(1000(F25.10,","))') t,x
    write(2, '(1000(F25.10,","))') t,x
    write(3, '(1000(F25.10,","))') t,x
    write(5, '(1000(F25.10,","))') t,x
    write(6, '(1000(F25.10,","))') t,x
    write(7, '(1000(F25.10,","))') t,x
    write(8, '(1000(F25.10,","))') t,x
    write(9, '(1000(F25.10,","))') t,x
    write(10, '(1000(F25.10,","))') t,x
    write(11, '(1000(F25.10,","))') t,x
    write(12, '(1000(F25.10,","))') t,x
    write(13, '(1000(F25.10,","))') t,x
    write(14, '(1000(F25.10,","))') t,x
    write(15, '(1000(F25.10,","))') t,x
    write(16, '(1000(F25.10,","))') t,x
    write(24, '(1000(F25.10,","))') t,x
    write(25, '(1000(F25.10,","))') t,x
    
    !Prints console output for the program parameters
    print*,'The parameter case being run is: ',parameter_case
    print*,'The calculation type is',calculation_type
    
    if (calculation_type == 1) then
        print*,'(Adams-Bashford method)'
    else if (calculation_type == 2) then
        print*, '(Runge-Kutta method)'
    else if (calculation_type == 3) then
        print*,'(Runge-Kutta method, with sediment incorporated)'
    else
        print*,'Invalid calculation type, stopping program.'
        stop
    end if
    
    print*,'The bedload method is',bedload_method
    
    if (bed_equation .ne. 0) then   
        print*,'Custom bed equation used: ',bed_equation
    else 
        print*,'Default bed elevation profile used'
    end if
    
    if (hyperconcentration_type == 0) then   
        print*,'Hyperconcentration disabled.'
    else if (hyperconcentration_type == 1) then
        print*,'Constant hyperconcentration used: ',hyperconcentration_conc, '%'
    else if (hyperconcentration_type == 2) then
        print*,'Ramping hyperconcentration used, multiplier: ',rho_w_multiplier
    end if
    
    print*,'The number of spatial steps is',n_steps
    print*,'The time step is ',dt
    print*,'The spatial step is ',dx
    print*,'The bed-slope upwinding type is',upwinding_type
    print*,'The bed slope is ',s_0
    print*,'The upstream boundary condition is ',west_bc
    print*,'The downstream boundary condition is ',east_bc
    print*,'The input specific flowrate is ',q_input
    print*,'The input flow depth is ',h_input
    print*,'The friction type is ',bed_friction_type
    print*,'The ramping type is',ramp_up_type
    print*,'Initial q_ramp is ',q_ramp_0
    print*,'Final q_ramp is ',q_ramp_end
    print*,'End time is ',t_end
    print*,'The number of time steps is ',n_timesteps
    print*,'z_b_0 is',z_b_0
    print*,'Profile factor a is ',a
    print*,'Profile factor b is ',b
    print*,'The initial bed porosity is ',bed_porosity(1)*100,'%'
    print*,'Sediment transport coefficient 1 is ',coeff_1_input
    print*,'Sediment transport coefficient 2 is ',coeff_2_input
    print*,'Sediment transport coefficient 3 is ',coeff_3_input
    print*,'The input deposition rate D is ',D_input
    print*,'The entrainment rate E is ',En_input 
    print*,'The initial saturated bed density is ',rho_0(1)
    print*,'The input bedload transport rate is ',q_b_0_input
    print*,'The input concentration is ',conc_0
    print*,'The upwinding type is ',upwinding_type
    
    if (bedrock_elevation_equation == 0) then
        
        print*,'Excavation depth is not limited.'
        
    else 
        
        print*,'Excavation depth limited, equation used: ',bedrock_elevation_equation

    end if
    
    if (d_median_equation == 0) then
        
        print*,'Particle median diameter taken as constant: .',d_median_input
        
    else 
        
        print*,'Particle median diameter equation used: ',d_median_equation

    end if
    
    if (calculation_type == 3 .AND. fixed_bed == 0) then
        
        print*,'Bed is movable.'
        
    else if (calculation_type == 3 .AND. fixed_bed == 1) then
        
        print*,'Bed is fixed.'
        
    else if (calculation_type == 3 .AND. fixed_bed == 2) then
        
        print*,'Bed is initially fixed.'

    end if
    
    print*,'Angle of repose is', angle_of_repose * 180.0d0/pi
    print*,'The critical Shield''s parameter is',critical_shields_parameter
    print '(A)'
    
    !Stores the value of the initial bed elevation for later comparison
    z_b_initial = z_b
    
    !Sets the number of iterations to zero
    iters = 1
    
    !Writes the three calculated arrays to a separate output file each, with the value for time of the current timestep in column 1
    write(1, '(F25.10, ",", 1000(F25.10,","))') t, q
    write(2, '(F25.10, ",", 1000(F25.10,","))') t,h
    write(3, '(F25.10, ",", 1000(F25.10,","))') t,u
    write(6, '(F25.10, ",", 1000(F25.10,","))') t,z_b
    write(7, '(F25.10, ",", 1000(F25.10,","))') t,conc
    write(11, '(F25.10, ",", 1000(F25.10,","))') t,free_surface_elevation
    write(12, '(F25.10, ",", 1000(F25.10,","))') t,q_b_0
    write(14, '(F25.10, ",", 1000(F25.10,","))') t,z_b
    write(16, '(F25.10, ",", 1000(F25.10,","))') t, d_median
    
    !Sets an initial value of the new bed elevation to ensure that there is initially no change in the bed elevation
    z_b_new = 0.0d0
    
    write(15, '(F25.10, ",", 1000(F25.10,","))') t,z_b_new
    write(17, '(F25.10, ",", 1000(F25.10,","))') t, z_b_new
    write(18, '(F25.10, ",", 1000(F25.10,","))') t, z_b_new
    write(19, '(F25.10, ",", 1000(F25.10,","))') t, rho_w
    
    !Sets the initial time to zero, and the initial iteration number to 1
    t = 0.0d0
    iters = 1

    print*,'Performing calculations:'
    
    !Runs the model until uniform flow is established
    if (fixed_bed == 2) then
        
        print '(A)'
        print*,'ESTABLISHING UNIFORM FLOW: ',iters
        print '(A)'
          
        !Sets sediment transport processes to zero
        coeff_1 = 0.0d0
        coeff_2 = 0.0d0
        is_bed_fixed = 1
        
        !Calculates the number of uniform flow timesteps
        n_uniform_flow_timesteps = uniform_flow_time/dt
        
        do timestep = 1,n_uniform_flow_timesteps
        
            !Calculates the current time within the simulation to achieve steady-state flow
            t = dt * (timestep)
            
            !For running a very large number of iterations, provides some visual output
            if (mod(iters,1000) == 0) then
            
                print*,'Iteration: ',timestep,'Time: ',t+dt,'The bed sediment area is : ',bed_sediment_area
                print*,'The total sediment mass is : ',total_sediment_mass
                print*,'h(1) is ',h(1),'h(n_steps) is ',h(n_steps)
                print*,'q(1) is ',q(1),'q(n_steps) is ',q(n_steps)
             
            end if
        
            !Calls the couple-Runge-Kutta calculation to run to uniform flow if applicable
            call runge_kutta_coupled_calc(n,n_steps, west_bc,east_bc,&
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
            
        !Adds one to the number of iterations
        iters = iters + 1
        
        end do
            
    end if
    
    !Resets the initial number of iterations to 1 and the initial time to zero
    iters = 1
    t = 0.0d0
    
    !Calculates the initial bed sediment area and total sediment mass
    bed_sediment_area = 0.0d0
    total_sediment_mass = 0.0d0    
            
    !Calculates the bed sedument area from all cells
    do n = 1, n_steps + 1
                
        bed_sediment_area = bed_sediment_area + dx * (z_b(n) + 0.5d0 * (z_b(n) - z_b(n -1))) 
        total_sediment_mass  = total_sediment_mass + (a_3(n) + 0.5d0 * (a_3(n) - a_3(n -1)))*dx

    end do
        
    !Gets the correct bed sediment area and bed sediment mass by removing values outside the domain boundary
    bed_sediment_area = bed_sediment_area - 0.5d0 * dx * (z_b(1) + 0.5d0 * (z_b(1) - z_b(0))) -  dx * 0.5d0 * (z_b(n_steps + 1) + 0.5d0 * (z_b(n_steps + 1) - z_b(n_steps)))
    total_sediment_mass = total_sediment_mass  - 0.5d0 * ((a_3(1) - a_3(0)) + a_3(1))*dx - 0.50d0 * ((a_3(n_steps+1) - a_3(n_steps)) + a_3(n_steps+1))*dx !+ bed_sediment_area * rho_s * (1.0d0-bed_porosity) 
    
    !Sets the bed movement coefficients depending on the fixed-bed status
    if (fixed_bed == 0) then
        coeff_1 = coeff_1_input
        coeff_2 = coeff_2_input
        coeff_3 = coeff_3_input
        
        is_bed_fixed = 0
            
    else if (fixed_bed == 1) then
        
        coeff_1 = 0.0d0
        coeff_2 = 0.0d0
        coeff_3 = 0.0d0
        is_bed_fixed = 1
        
    else if (fixed_bed == 2) then
        
        coeff_1 = coeff_1_input
        coeff_2 = coeff_2_input
        coeff_3 = coeff_3_input
        is_bed_fixed = 0
            
    end if
        
    print '(A)'
    print*,'MAIN CALCULATION: ',iters
    print '(A)'

    !The main loop; re-runs the calculations for every time-step
    do timestep = 1,n_timesteps 
            
        !Prints the iteration number for debugging
        if (end_prematurely == 'yes') then
            
            print '(A)'
            print*,'ITERATION: ',iters
            print*,'Time is ',t, 's'
            print '(A)'
            
        end if
        
        !For running a very large number of iterations, provides some visual output
        if (n_timesteps.ge.10000.and.end_prematurely=='no'.and.mod(iters,1000) == 0) then
            
            print*,'Iteration: ',timestep,'Time: ',t+dt,'The bed sediment area is : ',bed_sediment_area
            print*,'The total sediment mass is : ',total_sediment_mass
            print*,'h(1) is ',h(1),'h(n_steps) is ',h(n_steps)
            print*,'q(1) is ',q(1),'q(n_steps) is ',q(n_steps)
            print*,'conc(1) is',conc(1),'rho_w(1) is',rho_w(1)
               
        end if
     
        !Updates the time by adding the timestep
        t = dt * (timestep)
        
         !Calls the Adam-Bashford method
        if (calculation_type == 1) then
            
            call adams_bashford_method(q,h,h_new,q_new,c,c_b,n_mannings,dh_dt_a,dq_dt_a,dh_dt_b,&
                                        dq_dt_b,c_f_1,h_1,g,c_b_input,c_input,n_mannings_input,h_lim,&
                                        n, n_steps,z_b,dt, d_median, von_karman_constant, z_0s, h_input,&
                                        q_fixed,drive_flow)
        
        !Calls the main Runge-Kutta calculation
        else if (calculation_type == 2) then
            
            call main_calculation(n,n_steps,dx,dt,q,dh_dt_1,ramp_up_type,t,t_ramp,q_ramp_end,q_ramp_0,&
                                    bed_friction_type,c_f,c,c_b,n_mannings,h,g,c_b_input,c_input,n_mannings_input,h_input,h_lim,&
                                    u,rho_w,dq_dt_1,tau_bx,west_bc,east_bc,h_1,h_2,h_3,q_1,q_2,q_3,h_east,h_west,q_east,q_west, &
                                    tau_bx_1,tau_bx_2,tau_bx_3,u_1,u_2,u_3,z_b,i, drive_flow)
            
        !Calls the coupled Runge-Kutta method    
        else if (calculation_type == 3) then
            
            call runge_kutta_coupled_calc(n,n_steps, west_bc,east_bc,&
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
        
        end if
            
        !Calls the subroutine measuring the Courant number
        call courant_number_check(courant_number,u,dt,dx,t,mean_courant_number,has_exceeded_value)

        !Writes the three calculated arrays to a separate output file each, with the value for time of the current timestep in column 1
        write(1, '(F25.10, ",", 1000(F25.10,","))') t, q
        write(2, '(F25.10, ",", 1000(F25.10,","))') t, h
        write(3, '(F25.10, ",", 1000(F25.10,","))') t, u
        
        !If the sediment-water model is being used, writes additional parameters to the output files
        write(5, '(F25.10, ",", 1000(F25.10,","))') t, courant_number
        write(6, '(F25.10, ",", 1000(F25.10,","))') t, z_b
        write(7, '(F25.10, ",", 1000(F25.10,","))') t, conc
        write(8, '(F25.10, ",", 1000(F25.10,","))') t, D
        write(9, '(F25.10, ",", 1000(F25.10,","))') t, En
        write(10, '(F25.10, ",", 1000(F25.10,","))') t, R_n
        write(11, '(F25.10, ",", 1000(F25.10,","))') t, free_surface_elevation
        write(12, '(F25.10, ",", 1000(F25.10,","))') t, q_b_0
        write(13, '(F25.10, ",", 1000(F25.10,","))') t, u_b
        write(15, '(F25.10, ",", 1000(F25.10,","))') t,z_b - z_b_initial
        write(16, '(F25.10, ",", 1000(F25.10,","))') t, d_median
        write(17, '(F25.10, ",", 1000(F25.10,","))') t, erosion_depth
        write(18, '(F25.10, ",", 1000(F25.10,","))') t, deposition_depth
        write(19, '(F25.10, ",", 1000(F25.10,","))') t, rho_w
        write(24, '(F25.10, ",", 1000(F25.10,","))') t, Fr
        
        !Calculates the mean value of each parameter over all spatial steps for that timestep
        call mean_values(timestep, n_steps, n_timesteps, &
                                q, h, u, conc, q_b, u_b, En, D, R_n, courant_number, &
                                mean_q, mean_h, mean_u, mean_conc, mean_q_b, mean_u_b, mean_En, mean_D, mean_R_n, mean_courant_number)
        
        !Writes the mean values to a file        
        write(20, '(F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10)')t,mean_q(timestep), mean_u(timestep), mean_h(timestep), mean_conc(timestep), mean_q_b(timestep), mean_u_b(timestep), mean_En(timestep), mean_D(timestep), mean_R_n(timestep), mean_Courant_number(timestep)
        
        !Writes the values at position x = 50 to a file
        write(21, '(F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10)')t,q(int(n_steps/2)), u(int(n_steps/2)), h(int(n_steps/2)), conc(int(n_steps/2)), q_b(int(n_steps/2)), u_b(int(n_steps/2)), En(int(n_steps/2)), D(int(n_steps/2)), R_n(int(n_steps/2)), Courant_number(int(n_steps/2))
        
        !This if block ends the program early after running only a specified number of iterations
        !Use it for testing purposes
        if (end_prematurely == 'yes'.and.iters >= n_iters_to_run) then
            
            print*,'Program ended after ',n_iters_to_run,' iterations.'
            stop
            
        end if
        
        !Adds the free surface elevation at different points
        if (iters == 1 .or. iters == 148750 .or. iters == n_timesteps) then
            
            write(25, '(F25.10, ",", 1000(F25.10,","))') t, free_surface_elevation
            
        end if
       
        !Adds one to the number of iterations
        iters = iters + 1
            
    end do
    
    !Calculates the total CPU time taken
    call cpu_time(t_cpu_end)
    total_cpu_time = t_cpu_end - t_cpu_0
    print '(A)'
    print*,'The total CPU time was ',total_cpu_time,' s'
    
    !Prints the final calculated values at the end of the output
    print '(A)'
    print*,'Final values (n, h(n), u(n)):'
    do n = 1,n_steps
        
        print*,n,h(n),u(n)
        
    end do
    
    print '(A)'
    print*,'Final values (n, z_b(n)),eta(n):'
    do n = 1,n_steps
        
        print*,n,z_b(n),free_surface_elevation(n)
        
    end do
    
    print '(A)'
    print*,'Final values (n, q_b(n),u_b(n):'
    do n = 1,n_steps
        
        print*,n,q_b(n),u_b(n)
        
    end do
    
    !Displays the mean values after the program has finished running if a simple case is being run
    if (display_mean_values == 1) then
        
        print '(A)'
        print*,'Final mean values:'
        print*,'Mean depth:',sum(h) / dble(n_steps + 4)
        print*,'Mean velocity:',sum(u) / dble(n_steps + 4)
        print*,'Mean flowrate:',sum(q) / dble(n_steps + 4)
        print*,'Mean concentration:',sum(conc) / dble(n_steps + 4)
        print*,'Mean bed flowrate:',sum(q_b) / dble(n_steps + 4)
        print*,'Mean bed velocity:',sum(u_b) / dble(n_steps + 4)
        print*,'Mean entrainment:',sum(En) / dble(n_steps + 4)
        print*,'Mean deposition:',sum(D) / dble(n_steps + 4)
        print*,'Mean Rouse number:',sum(R_n) / dble(n_steps + 4)
        print*,'Mean Courant number:',sum(Courant_number) / dble(n_steps + 4)
        
    end if
    
    !Writes the averaged final values to a file with the relevant parameters
    write(22, '(F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",F25.10, ",",A, ",",F25.10)')sum(q) / dble(n_steps + 4),sum(h) / dble(n_steps + 4),sum(u) / dble(n_steps + 4), sum(conc) / dble(n_steps + 4), sum(q_b) / dble(n_steps + 4), sum(u_b) / dble(n_steps + 4), sum(En) / dble(n_steps + 4), sum(D) / dble(n_steps + 4), sum(R_n) / dble(n_steps + 4),sum(Courant_number) / dble(n_steps + 4),hyperconcentration_conc,coeff_3,total_cpu_time, simulation_batch_name,d_median_input
    
    !Saves the program performance output to a file
    write(4, '(A,",",A)') 'Parameter:','Value:'
    write(4, '(A,",",F25.10)') 'Total CPU time',total_cpu_time
    write(4, '(A,",",I6)') 'Number of spatial steps', n_steps
    write(4, '(A,",",I6)') 'Number of time steps', n_timesteps 
    write(4, '(A,",",F25.10)') 'Start time', t_0
    write(4, '(A,",",F25.10)') 'Final time', t_end
    write(4, '(A,",",F25.10)') 'Timestep', dt
    write(4, '(A,",",F25.10)') 'Spatial step', dx
    write(4, '(A,",",I10)') 'Parameter case', parameter_case
    write(4, '(A,",",F25.10)') 'Total length', x_tot
    
    !Writes the final bed elevation values to a file
    write(14, '(F25.10, ",", 1000(F8.5,","))') t,z_b
    
    !Closes all the output files
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(6)
    close(7)
    close(8)
    close(9)
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
    close(23)
    close(24)
    close(25)

    !Plays a sound to alert the user that the simulation has finished running
    call system('powershell (New-Object Media.SoundPlayer "..\Sounds\Music.wav").PlaySync()')
    
end program main