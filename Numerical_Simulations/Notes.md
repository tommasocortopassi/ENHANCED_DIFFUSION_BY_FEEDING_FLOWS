## Notes
We provide the functions, scripts and data to perform numerical simulations
of tracers advected by either the synthetic population A or B as defined in
the article. For more details, check the documentation of each function.

# Functions
 - test_displacement_synth_pop      :   Simulates the trajectory of tracers advected
                                        by synthetic populations of cells 
                                        (i.e. Regularised Stokeslets). Returns 
                                        relevant info, such as MSD_xy, MSD_z    
                                        and segmented tracers trajectories, as
                                        well as info on encounter prameters
                                        and a sample of full tracers trajectories

- segmented_trajectories_msd         :  Computes the MSD (in both xy and z
                                        variables) assuming only relevant displacements 
                                        occur, i.e. displacements satisfying a 
                                        given threshold condition on their norm.
                                        Also returns the directional correlation    
                                        between successive relevant displacements.

- select_synth_population           :   Either randomises the forcing angles after
                                        each inactive period (population A) 
                                        or it doesn't (population B).

- v_NN                              :   Computes the velocity field induced on 
                                        a tracer by its nearest neighbour only.

- total_v                           :   Computes the velocity field induced on 
                                        a tracer by the whole population of cells.

- ScriptSynthPop                    :   Script run to compute the MSD, after a given 
                                        time T, of a tracer initially
                                        placed at different heights.

# Data

- options_synth_population          :   Example of data to use as input for 
                                        test_displacement_synth_pop.
                                        