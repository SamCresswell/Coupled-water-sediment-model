module variables
    
    !Uses the step number module to be able to ensure arrays are the correct size
    use initialisation
    
    implicit none
    
    !Links to the filepath where the re-usable parameters file is saved
    character(*), parameter :: file_path = '..\Parameters\'
    
    !-------------------------------------------------------------------------
    !All variables used by the program are declared here
    
    !Constants are declared here
    !Accleration due to gravity (m s^-2)
    double precision :: g
    !The mathematical constant pi
    double precision, parameter :: pi = 4.d0 * datan(1.d0)
    !The von Karman constant
    double precision :: von_karman_constant

    !Material properties are declared here
    !The density of clear water
    double precision :: rho_w_clear 
    !Sediment density (kg m^-3)
    double precision :: rho_s
    !Water density (kg m^-3)
    double precision,dimension(-1:n_steps + 2) :: rho_w
    !The density of the saturated bed (kgm^-3)
    double precision,dimension(-1:n_steps + 2) :: rho_0
    !The density of the combined sediment-water mixture (kg m^-3)
    double precision,dimension(-1:n_steps + 2) :: rho
    !The input constant bed porosity (unitless)
    double precision :: bed_porosity_input
    !The porosity of the bed (unitless)
    double precision,dimension(-1:n_steps + 2) :: bed_porosity
    !Whether to calculate the bed porosity
    logical :: calculate_bed_porosity
    
    !Program parameters are defined here
    !Current step
    integer :: n
    !Current number of iterations
    integer :: iters
    !Extra integration variable
    integer :: i
    !Calculation type
    !Set to 1 for Adams-Bashford, 2 for fourth-order Runge-Kutta method and 3 for the Runge-Kutta sediment model
    integer :: calculation_type
    !The upwinding type
    !Set to 1 for no upwinding, 2 for upwinding, 3 for a linear combination of a central difference and upwinding
    !and 4 to use a central difference with a slope limiter
    integer :: upwinding_type
    !This is the dummy variable, which is used to read data to be ignored
    character(len=0) :: e
    !This is a debugging variable
    integer :: e2
    !This is the parameters filepath variable
    character(len = 100) :: parameter_file
    
    !Geometry variables are declared here
    !Channel spatial grid of x-positions
    double precision,dimension(-1:n_steps+2) :: x
    !Total length of domain (m)
    double precision :: x_tot
    !Which bed elevation equation to use
    !Set to 0 for no custom bed equation and a constant slope, 1 for a the humped bed case, 
    !2 for the 1m elevation constant slope, 3 for the three-slope Karnali case,
    !4 for the 2-slope Karnali case, 5 for an experimental Karnali case
    integer :: bed_equation
    !Bed elevation (m):
    double precision,dimension(-1:n_steps + 2) :: z_b
    !Initial bed elevation (m)
    double precision :: z_b_0
    !Final bed elevation (m)
    double precision :: z_b_end
    !The initial bed profile
    double precision,dimension(-1:n_steps + 2) :: z_b_initial
    !Length of each step (m)
    double precision :: dx
    !Bed slope
    double precision :: s_0

    !Flow variables are defined here
    !Depth of flow (m)
    double precision,dimension(-1:n_steps+2) :: h
    !Free surface elevation (m)
    double precision,dimension(-1:n_steps+2) :: free_surface_elevation
    !Specific flowrate (m^2 s^-1)
    double precision,dimension(-1:n_steps+2) :: q
    !Flow velocity (m s^-1)
    double precision,dimension(-1:n_steps+2) :: u
    !Current flow velocity (m s^-1)
    double precision,dimension(-1:n_steps+2) :: u_current
    !New iteration of depth of flow (m)
    double precision,dimension(-1:n_steps+2) :: h_new
    !New iteration of specific flowrate (m^2 s^-1)
    double precision,dimension(-1:n_steps+2) :: q_new
    !New iteration of flow velocity (m s^-1)
    double precision,dimension(-1:n_steps+2) :: u_new
    !The first calculated height (m)
    double precision :: h_1
    !The second calculated height (m)
    double precision :: h_2
    !The third calculated height (m)
    double precision :: h_3
    !The calculated upstream height (m)
    double precision :: h_west
    !The calculated downstream height (m)
    double precision :: h_east
    !The first calculated specific flowrate (m^2s^-1)
    double precision :: q_1
    !The second calculated specific flowrate (m^2s^-1)
    double precision :: q_2
    !The third calculated specific flowrate (m^2s^-1)
    double precision :: q_3
    !The calculated upstream specific flowrate (m^2s^-1)
    double precision :: q_west
    !The calculated downstream specific flowrate (m^2s^-1)
    double precision :: q_east
    !The first calculated flow velocity (m^s^-1)
    double precision,dimension(-1:n_steps+2) :: u_1
    !The second calculated flow velocity (m^s^-1)
    double precision :: u_2
    !The third calculated flow velocity (m^s^-1)
    double precision :: u_3
     !Extra iteration variable for the friction coefficient for the Adams-Bashford method
    integer :: n_f
    !Extra height variable for the friction calculation for the Adams-Bashford method
    double precision :: h_f
    
    !Time variables are declared here
    !Current time (s)
    double precision :: t
    !Initial time (s)
    double precision :: t_0
    !Final time (s)
    double precision :: t_end
    !Timestep (s)
    double precision :: dt
    !Current timestep
    integer :: timestep
    !The initial running time to establish uniform flow
    double precision :: uniform_flow_time
    
    !Boundary conditions are declared here
    !Western boundary condition type
    !Set to 5 for uniform flow and 2 for transmissive flow
    integer :: west_bc
    !Eastern boundary condition type
    !Set to 5 for uniform flow and 2 for transmissive flow
    integer :: east_bc
    
    !Friction conditions are declared here
    !The type of bed friction
    !Set to 1 for cf = cb, 2 for cf = g/cb^2, 3 for cf = g*cb^2/h^(1/3)
    !and set to 4 to use the bed roughness length
    integer :: bed_friction_type
    !The calculated friction coefficient
    double precision,dimension(-1:n_steps+2) :: c_f
    !The simple bed friction coefficient
    double precision,dimension(-1:n_steps+2) :: c_b
    !Chezy coefficient
    double precision,dimension(-1:n_steps+2) :: c
    !Manning's roughness
    double precision,dimension(-1:n_steps+2) :: n_mannings
    !Input simple bed friction coefficient
    double precision :: c_b_input
    !character(len=10) :: c_b_input
    !Input Chezy coefficient
    double precision :: c_input
    !Input Manning's roughness
    double precision:: n_mannings_input
    !The bed friction coefficient
    double precision,dimension(-1:n_steps+2) :: tau_bx 
    !The first bed friction coefficient
    double precision,dimension(-1:n_steps+2) :: tau_bx_1 
    !The second bed friction coefficient
    double precision,dimension(-1:n_steps+2) :: tau_bx_2 
    !The third bed friction coefficient
    double precision,dimension(-1:n_steps+2) :: tau_bx_3 
    !The one-dimensional friction coefficient
    double precision :: c_f_1
    !The bed roughness length
    double precision,dimension(-1:n_steps+2) :: z_0s
    
    !Sediment transport parameters are defined here
    !Input load transport rate (m^2s^-1)
    double precision :: q_b_0_input
    !The bed load transport rate (m^2s^-1)
    double precision,dimension(-1:n_steps+2) :: q_b_0
    !The sediment transport rate (m^2s^-1)
    double precision,dimension(-2:n_steps+2) :: q_b
    !The bed velocity (ms^-1)
    double precision,dimension(-1:n_steps+2) :: u_b
    !Whether a custom deposition equation will be used, and which one if so
    !Set to zero for a constant value
    integer :: use_deposition_equation
    !Whether a custom entrainment equation will be used, and which one if so
    !Set to zero for a constant value
    integer :: use_entrainment_equation
    !Input deposition rate (m^2s^-1)
    double precision :: D_input
    !double precision :: D
    double precision,dimension(-1:n_steps + 2) :: D
    !The input entrainment rate (m^2s^-1)
    double precision :: En_input
    !Entrainment rate (m^2s^-1)
    double precision,dimension(-1:n_steps + 2) :: En
    !Entrainment rate (m^2s^-1) value 1
    double precision,dimension(-1:n_steps + 2) :: En_1
    !Entrainment rate (m^2s^-1) value 2
    double precision,dimension(-1:n_steps + 2) :: En_2
    !Entrainment rate (m^2s^-1) value 3
    double precision,dimension(-1:n_steps + 2) :: En_3
    !Entrainment rate (m^2s^-1) value 1
    double precision,dimension(-1:n_steps + 2) :: En_4
    !Deposition rate (m^2s^-1) value 1
    double precision,dimension(-1:n_steps + 2) :: D_1
    !Deposition rate (m^2s^-1) value 2
    double precision,dimension(-1:n_steps + 2) :: D_2
    !Deposition rate (m^2s^-1) value 3
    double precision,dimension(-1:n_steps + 2) :: D_3
    !Deposition rate (m^2s^-1) value 1
    double precision,dimension(-1:n_steps + 2) :: D_4
    !Coefficient 1 linking bed transport rate and flow velocity
    double precision :: coeff_1
    !Coefficient 2 linking bed velocity and flow velocity
    double precision :: coeff_2
    !Coefficient 3 for the Mueller-Peter-Meyer equation
    double precision :: coeff_3
    !Input coefficient 3 for the Mueller-Peter-Meyer equation
    double precision :: coeff_3_input
    !Input for coefficient 1 linking bed transport rate and flow velocity
    double precision :: coeff_1_input
    !Input for coefficient 2 linking bed velocity and flow velocity
    double precision :: coeff_2_input
    !Which method to use to calculated the bedload transport rate
    integer :: bedload_method
    !The sediment concentration
    double precision,dimension(-1:n_steps + 2) :: conc
    !The initial suspended sediment concentration
    double precision :: conc_0
    !Which equation to use for d_median
    !Use 0 for a constant value
    integer :: d_median_equation
    !The input median particle diameter (m)
    double precision :: d_median_input
    !The median particle diameter as a function of space (m)
    double precision,dimension(-1:n_steps + 2) :: d_median
    !The initial median particle diameter as a function of space (m)
    double precision,dimension(-1:n_steps + 2) :: d_median_initial
    !The new median particle diameter as a function of space (m)
    double precision,dimension(-1:n_steps + 2) :: d_median_new
    !The relative density of sediment to water (rho_s/rho_w)
    double precision,dimension(-1:n_steps + 2) :: s
    !The input dynamic viscosity
    double precision :: nu_input
    !The kinematic viscosity of the water (m^2s^-1)
    double precision,dimension(-1:n_steps + 2) :: nu
    !The particle settling velocity (ms^-1)
    double precision,dimension(-1:n_steps + 2) :: w_s_0
    !The hindered particle settling velocity (ms^-1)
    double precision,dimension(-1:n_steps + 2) :: w_s_h
    !The dimensionless particle diameter parameter
    double precision,dimension(-1:n_steps + 2) :: d_dimensionless 
    !The particle Reynolds number'
    double precision,dimension(-1:n_steps + 2) :: Re_p
    !The m coefficient
    double precision,dimension(-1:n_steps + 2) :: m_coefficient
    !Whether the bed is fixed
    integer :: fixed_bed
    !Whether the bed is currently fixed
    integer :: is_bed_fixed
    !Uniform flow timestep
    integer :: n_uniform_flow_timesteps
    !The deposition exponent
    double precision :: m_d
    !The entrainment coefficient
    double precision :: alpha_en = 0.015
    !The erosion (m)
    double precision,dimension(-1:n_steps + 2) :: erosion_depth
    !The deposition depth (m)
    double precision,dimension(-1:n_steps + 2) :: deposition_depth
    !The deposition coefficient
    double precision :: alpha_d
    !The Rouse number
    double precision,dimension(-1:n_steps + 2) :: R_n
    !The shear velocity (ms^-1)
    double precision :: u_shear
    !The Froude number for the flow
    double precision,dimension(-1:n_steps + 2) :: Fr
    !Whether to update the grain diameter
    integer :: update_grain_size
    
    !Input parameters are declared here
    !Initial depth (m)
    double precision :: h_input
    !Initial flow (m^2s^-1)
    double precision :: q_input
    !Upstream flowrate (m^2s^-1)
    double precision :: q_fixed
    
    !'A' parameters for the sediment-water model are declared here
    !Parameter 'a_1', equal to rho_w * h_input
    double precision,dimension(-1:n_steps+2) :: a_1
    !Parameter 'a_2', equal to a * rho * h * u
    double precision,dimension(-1:n_steps+2) :: a_2 
    !Parameter 'a_3', equal to rho_s * c * h
    double precision,dimension(-1:n_steps+2) :: a_3
    
    !Ramping conditionbs are declared here
    !The ramping up condition
    !Set to 1 for no ramp up, set to 2 for ramped up wind and 3 for ramped up inflows and 5 for custom ramping equation
    integer :: ramp_up_type
    !Ramped up inlet flow discharge over total time
    double precision :: q_ramp_0 
    !Ramped up inlet flow discharge at time t
    double precision,dimension(n_timesteps) :: q_ramp
    !Ramped up inlet flow discharge at final time
    double precision :: q_ramp_end
    !Ramping start time
    double precision :: t_ramp_0
    !Total ramping up time
    double precision :: t_ramp
    !Whether to drive the flowrate
    integer :: drive_flow
    !The hyperconcentration type
    !Set to 0 for no hyperconcentration, 1 for a constant value and 3 for matching the flow ramping
    integer :: hyperconcentration_type
    !The volumetric concentration of a simulated hyperconcentrated flow (%)
    double precision :: hyperconcentration_conc
    !The multiplier for the clear water density and the maximum hyperconcentrated water density
    double precision :: rho_w_multiplier
    !The profile factor for velocity
    double precision :: a
    !The profile factor for the non-uniform suspended sediment concentration and moment distribution
    double precision :: b
    
    !Limiting conditions are declared here
    !The limiting depth (m)
    double precision :: h_lim
    
    !Differential terms are declared here
    !Derivative of head with respect to time at first position
    double precision,dimension(-1: n_steps +2) :: dh_dt_1
    !Derivative of flowrate with respect to time at first position
    double precision,dimension(-1: n_steps +2) :: dq_dt_1
    !Derivative of head with respect to time at second position
    double precision,dimension(-1: n_steps +2) :: dh_dt_2
    !Derivative of flowrate with respect to time at second position
    double precision,dimension(-1: n_steps +2) :: dq_dt_2
    !Derivative of head with respect to time at third position
    double precision,dimension(-1: n_steps +2) :: dh_dt_3
    !Derivative of flowrate with respect to time at third position
    double precision,dimension(-1: n_steps +2) :: dq_dt_3
    !Derivative of head with respect to time at fourth position
    double precision,dimension(-1: n_steps +2) :: dh_dt_4
    !Derivative of flowrate with respect to time at fourth position
    double precision,dimension(-1: n_steps +2) :: dq_dt_4
    
    !Derivatives for Adam-Bashford method
    !Derivative of head with respect to time at first position
    double precision,dimension(n_steps) :: dh_dt_a
    !Derivative of flowrate with respect to time at first position
    double precision,dimension(n_steps) :: dq_dt_a
    !Derivative of head with respect to time at second position
    double precision :: dh_dt_b
    !Derivative of flowrate with respect to time at second position
    double precision :: dq_dt_b
    
    !Derivatives for the sediment-water Runge Kutta method are declared here
    !The change in the bedload transport rate with respect to time
    double precision,dimension(-1: n_steps +2) :: dq_b_dt
    !The first calculated derivative of the change in bed gradient with respect to time
    double precision,dimension(-1: n_steps +2) :: dzb_dt_1
    !The second calculated derivative of the change in bed gradient with respect to time
    double precision,dimension(-1: n_steps +2) :: dzb_dt_2
    !The third calculated derivative of the change in bed gradient with respect to time
    double precision,dimension(-1: n_steps +2) :: dzb_dt_3
    !The fourth calculated derivative of the change in bed gradient with respect to time
    double precision,dimension(-1: n_steps +2) :: dzb_dt_4
    !The first calculated derivative of the change in 'a_1' with respect to time
    double precision,dimension(-1: n_steps +2) :: da1_dt_1
    !The second calculated derivative of the change in 'a_1' with respect to time
    double precision,dimension(-1: n_steps +2) :: da1_dt_2
    !The third calculated derivative of the change in 'a_1' with respect to time
    double precision,dimension(-1: n_steps +2) :: da1_dt_3
    !The fourth calculated derivative of the change in 'a_1' with respect to time
    double precision,dimension(-1: n_steps +2) :: da1_dt_4
    !The first calculated derivative of the change in 'a_2' with respect to time
    double precision,dimension(-1: n_steps +2) :: da2_dt_1
    !The second calculated derivative of the change in 'a_2' with respect to time
    double precision,dimension(-1: n_steps +2) :: da2_dt_2
    !The third calculated derivative of the change in 'a_2' with respect to time
    double precision,dimension(-1: n_steps +2) :: da2_dt_3
    !The fourth calculated derivative of the change in 'a_2' with respect to time
    double precision,dimension(-1: n_steps +2) :: da2_dt_4
    !The first calculated derivative of the change in 'a_3' with respect to time
    double precision,dimension(-1: n_steps +2) :: da3_dt_1
    !The second calculated derivative of the change in 'a_3' with respect to time
    double precision,dimension(-1: n_steps +2) :: da3_dt_2
    !The third calculated derivative of the change in 'a_3' with respect to time
    double precision,dimension(-1: n_steps +2) :: da3_dt_3
    !The fourth calculated derivative of the change in 'a_3' with respect to time
    double precision,dimension(-1: n_steps +2) :: da3_dt_4
    
    !Temporary parameters for the sediment-water Runge-Kutta method are stored here
    !Middle value of the first parameter 'a_1'
    double precision :: a_1_1
    !East value of the first parameter 'a_1'
    double precision :: a_1_east
    !West value of the first parameter 'a_1'
    double precision :: a_1_west
    !Middle value of the second parameter 'a_2'
    double precision :: a_2_1
    !East value of the second parameter 'a_2'
    double precision :: a_2_east
    !West value of the second parameter 'a_2'
    double precision :: a_2_west
    !Middle value of the third parameter 'a_3'
    double precision :: a_3_1
    !East value of the third parameter 'a_3'
    double precision :: a_3_east
    !West value of the third parameter 'a_3'
    double precision :: a_3_west
    !East value of the bed elevation 
    double precision :: z_b_east
    !West value of the bed elevation 
    double precision :: z_b_west
    
    !Parameters to be updated in the sediment-water Runge-Kutta method are declared here
    !The new value of parameter 'a_1'
    double precision,dimension(-1: n_steps +2) :: a_1_new
    !The new value of parameter 'a_2'
    double precision,dimension(-1: n_steps +2) :: a_2_new
    !The new value of parameter 'a_3'
    double precision,dimension(-1: n_steps +2) :: a_3_new
    !The new value of dzb/dt
    double precision,dimension(-1: n_steps +2) :: z_b_new
    
    !Model tracking terms are declared here
    !Initial CPU time
    double precision :: t_cpu_0 
    !Final CPU time
    double precision :: t_cpu_end     
    !Total CPU time
    double precision :: total_cpu_time      
    !The Courant number
    double precision,dimension(-1: n_steps +2) :: courant_number
    !Whether the Courant number has been exceeded
    integer :: has_exceeded_value
    
    !Variables for the bedrock equation are stored here
    !Whether to use a limiting bedrock equation
    !Set to 0 for no limiting bedrock layer, 1 for a linear equation and 2 for a custom equation
    integer :: bedrock_elevation_equation
    !The bedrock elevation at the western end
    double precision :: bedrock_elevation_west
    !The bedrock elevation at the eastern end
    double precision :: bedrock_elevation_east
    !The bedrock elevation at each point
    double precision,dimension(-1: n_steps +2) :: bedrock_elevation
    
    !Slope limiting variables go here
    !The current bed slope
    double precision,dimension(-1: n_steps +2) :: S_b
    !The bed slope angle
    double precision,dimension(-1: n_steps +2) :: bed_slope_angle
    !The angle of repose of the slope (degrees)'
    double precision :: angle_of_repose
    !The shields parameter (in N m^-2)
    double precision,dimension(-1: n_steps +2) :: shields_parameter
    !The critical shields parameter
    double precision :: critical_shields_parameter
    !The dimensionless critical shields parameter
    double precision,dimension(-1: n_steps +2) :: dimensionless_critical_shields_parameter
    !The modified critical shields parameter
    double precision,dimension(-1: n_steps +2) :: modified_critical_shields_parameter
    !The first minmod slope function
    double precision,dimension(-1: n_steps +2) :: minmod_a
    !The second minmod slope function
    double precision,dimension(-1: n_steps +2) :: minmod_b
    !The combined minmod slope function
    double precision,dimension(-1: n_steps +2) :: minmod_ab
    
    !Mean values go in this block
    !The mean specific flowrate (m^2s^-1)
    double precision,dimension(0: n_timesteps) :: mean_q
    !The mean water depth (m)
    double precision,dimension(0: n_timesteps) :: mean_h
    !The mean flow velocity (ms^-1)
    double precision,dimension(0: n_timesteps) :: mean_u
    !The mean flow velocity (ms^-1)
    double precision,dimension(0: n_timesteps) :: mean_conc    
    !The mean bed flowrate (m^2s^-1)
    double precision,dimension(0: n_timesteps) :: mean_q_b
    !The mean bed velocity (ms^-1)
    double precision,dimension(0: n_timesteps) :: mean_u_b
    !The mean entrainment (ms^-1)
    double precision,dimension(0: n_timesteps) :: mean_En
    !The mean deposition (ms^-1)
    double precision,dimension(0: n_timesteps) :: mean_D
    !The mean Rouse number
    double precision,dimension(0: n_timesteps) :: mean_R_n
    !The mean Courant number
    double precision,dimension(0: n_timesteps) :: mean_Courant_number
    
    !Verification parameters are declared here
    !The area under the graph
    double precision :: bed_sediment_area 
    !The total mass of sediment
    double precision :: total_sediment_mass
    
    contains
    
    !-------------------------------------------------------------------------
    !This subroutine reads the input data from the csv file
        subroutine read_input_data
        
            !Reads the relevant parameters file
            !For the generic test case
            if (parameter_case == 0) then
                
                parameter_file = 'Parameters.csv'
            
            !For the sediment mound test case
            else if (parameter_case == 1) then
                
                parameter_file = 'Parameters_test_case_1.csv'
                
            !For the deposition test case
            else if (parameter_case == 2) then
                
                parameter_file = 'Parameters_test_case_2.csv'
            
            !For the entrainment test case
            else if (parameter_case == 3) then
                
                parameter_file = 'Parameters_test_case_3.csv'
            
            !For the triple-slope case
            else if (parameter_case == 4) then
                
                parameter_file = 'Parameters_Karnali_triple_slope.csv'
                
            else if (parameter_case == 5) then
                
                parameter_file = 'Parameters_Karnali_single_slope.csv'
                
            else if (parameter_case == 6) then
                
                parameter_file = 'Parameters_Karnali_double_slope.csv'
                
            else if (parameter_case == 7) then
                
                parameter_file = 'Parameters_Karnali_GST_evolution.csv'
                
            else if (parameter_case == 8) then
                
                parameter_file = 'Parameters_Karnali_fixed_bed.csv'

                
            else
                
                print*,'Invalid parameters set, stopping program.'
                
                stop
                
            end if
        
            !Opens the csv file, which is saved in the same directory
            open(2, file = file_path//parameter_file, status = 'old', access = 'sequential')
    
            !Skips the first line with the header
            read(2,*)
            
            !Reads the calculation parameters block
            !Reads the value for calculation type
            read(2,*)e,e,e,calculation_type
            !Reads the upwinding type
            read(2,*)e,e,e,upwinding_type
            !Reads whether to use a fixed bed
            read(2,*)e,e,e,fixed_bed
            !Whether to run the model to establish uniform flow first
            read(2,*)e,e,e,uniform_flow_time
            
            !Reads the time parameters block
            read(2,*)
            !Reads the value for initial time
            read(2,*)e,e,e,t_0
            !Reads the value for the final time
            read(2,*)e,e,e,t_end
            
            !Reads the input parameters block
            read(2,*)
            !Reads the initial input depth
            read(2,*)e,e,e,h_input
            !Reads the initial input flowrate
            read(2,*)e,e,e,q_input
            !Reads the initial input flowrate
            read(2,*)e,e,e,drive_flow
            
            !Reads the boundary conditions block
            read(2,*)
            !Reads the value for the western boundary condition
            read(2,*)e,e,e,west_bc
            !Reads the value for the eastern boundary condition
            read(2,*)e,e,e,east_bc    

            !Reads the ramping block
            read(2,*) 
            !Reads the ramping type value
            read(2,*)e,e,e,ramp_up_type
            !Reads the initial ramping time
            read(2,*)e,e,e,t_ramp_0
            !Reads the total ramping time
            read(2,*)e,e,e,t_ramp
            !Reads the initial flow ramping value
            read(2,*)e,e,e,q_ramp_0
            !Reads the final flow ramping value
            read(2,*)e,e,e,q_ramp_end
            !Reads the hyperconcentration type
            read(2,*)e,e,e,hyperconcentration_type
            !The volumetric concentration of the simulated hyperconcentration
            read(2,*)e,e,e,hyperconcentration_conc
            
            !Reads the problem geometry block
            read(2,*)
            !Reads the value for the total length 
            read(2,*)e,e,e,x_tot
            !Reads the value for whether to use a custom bed elevation
            read(2,*)e,e,e,bed_equation
            !Reads the value for the initial elevation
            read(2,*)e,e,e,z_b_0
            !Reads the value for the final elevation
            read(2,*)e,e,e,z_b_end
            !Reads whether to use an equation for the median particle diameter
            read(2,*)e,e,e,d_median_equation
            !Reads the input median particle diameter
            read(2,*)e,e,e,d_median_input
            
            !Reads the friction conditions block
            read(2,*)
            !Reads the value for the type of friction calculation used
            read(2,*)e,e,e,bed_friction_type
            !Reads the value for the input simple bed friction type
            read(2,*)e,e,e,c_b_input
            !Reads the value for the input Chezy coefficient
            read(2,*)e,e,e,C_input
            !Reads the value for the input Manning's coefficient
            read(2,*)e,e,e,n_mannings_input
            !Reads the limiting depth value
            read(2,*)e,e,e,h_lim
            
            !Reads the sediment transport and velocity profile block
            read(2,*)
            !Reads the initial suspended sediment concentration
            read(2,*)e,e,e,conc_0
            !Reads the value for whether a custom deposition equation will be used
            read(2,*)e,e,e,use_deposition_equation
            !Reads the value for whether a custom entrainment equation will be used
            read(2,*)e,e,e,use_entrainment_equation
            !Reads the input deposition rate value
            read(2,*)e,e,e,D_input
            !Reads the entrainment rate
            read(2,*)e,e,e,En_input
            !Reads the deposition exponent
            read(2,*)e,e,e,m_d
            !Reads the entrainment coefficient
            read(2,*)e,e,e,alpha_en
            !Reads the bed load transport rate
            read(2,*)e,e,e,q_b_0_input
            !Reads the first sediment transport coefficient
            read(2,*)e,e,e,bedload_method
            !Reads the first sediment transport coefficient
            read(2,*)e,e,e,coeff_1_input
            !Reads the second sediment transport coefficient
            read(2,*)e,e,e,coeff_2_input
            !Coefficient 3 for the Mueller-Peter-Meyer equation
            read(2,*)e,e,e,coeff_3_input
            !Reads the limiting depth value
            read(2,*)e,e,e,a
            !Reads the convergence determination value
            read(2,*)e,e,e,b
            
            !Reads the constants block
            read(2,*)
            !Reads the value for acceleration due to gravity
            read(2,*)e,e,e,g
            !Reads the value for acceleration due to gravity
            read(2,*)e,e,e,von_Karman_constant
            
            !Reads the material properties block
            read(2,*)
            !Reads the value for quartzite sediment density
            read(2,*)e,e,e,rho_s
            !Reads the value for water density
            read(2,*)e,e,e,rho_w_clear
            !Whether to calculate bed porosity
            read(2,*)e,e,e,calculate_bed_porosity
            !Reads the value for bed porosity
            read(2,*)e,e,e,bed_porosity_input
            !Reads the kinematic viscosity of the water
            read(2,*)e,e,e,nu_input
            !Reads the angle of repose
            read(2,*)e,e,e,angle_of_repose
            
            !Reads the bedrock parameters block
            read(2,*)
            !Reads the bedrock elevation equation
            read(2,*)e,e,e,bedrock_elevation_equation
            !Reads the bedrock elevation at the western end
            read(2,*)e,e,e,bedrock_elevation_west
            !Reads the bedrock elevation at the eastern end
            read(2,*)e,e,e,bedrock_elevation_east
            !Reads whether to update the 
            read(2,*)e,e,e,update_grain_size
        
            !Closes the input parameters file
            close ( 2, status = 'keep') 
            
        end subroutine read_input_data
        
    end module variables
