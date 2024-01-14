!This module exists solely to define the number of steps. Unfortunately due to the limitations of Fortran,
!the numbers of timesteps and spatial steps have to be defined within the program and cannot be placed in a csv file
module initialisation
    
    implicit none
    
    !Set the number of spatial and timesteps here
    integer,parameter :: n_steps = 100
    integer,parameter :: n_timesteps = 7500
    
    !Use this block for testing the loop without running all iterations
    !A parameter controlling whether the program ends prematurely for testing purposes
    !Set to 'yes' for breaking the loop after x iterations
    character(len=3) :: end_prematurely = 'no'
    !The number of iterations to run
    integer :: n_iters_to_run = 1
    
    !Use this parameter to set the model parameter set. 
    !Set to 0 for the test case, 1 for the moving sandbar test case, 2 for the deposition test case
    !3 for the entrainment test case, 4 for the three-slope Karnali case, 5 for the one-slope Karnali case
    !6 for the double-slope gravel Karnali case, 7 for the experimental GST evolution case
    !and 8 for the experimental fixed-sand bed case
    integer,parameter :: parameter_case = 0
    
    !Set this to 1 to display mean values when the program has run
    integer,parameter :: display_mean_values = 1
    
    !Set this to 1 to end the program as soon as bed mobilisation occurs and output this value
    integer,parameter :: end_program_at_mobilisation = 0
    
    !Use this parameter as the simulation batch name for exported results
    character(len = 5) :: simulation_batch_name = 'K0_40'
    
end module initialisation