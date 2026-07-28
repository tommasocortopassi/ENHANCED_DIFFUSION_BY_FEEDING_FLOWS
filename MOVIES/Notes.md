## Notes

In this folder we have the functions used for creating the movies, the 
necessary data and the movies themselves. Launch the scripts to obtain the 
corresponding movies.

# Functions & Scripts

 - movie_traj            :  Given data on tracers and Stylonychia trajectories
                            and forcing directions, creates a video.

 - fit_trajectories     :   Linearly interpolates missing NaN values in the 
                            data on Stylonychia trajectories and force orientations.

 - simulate_tracer      :   Computes the trajectory of a tracer particle subject
                            to the action of an ensemble of cells.
   Script_Movie1
 - Script_Movie2        :   Launch to recreate the corresponding movies.
   Script_Movie3

 


# IMPORTANT REMARK
In the Dataset1.csv provided (see RAW_DATA folder) the angles are measured with 
respect to the x axis, but positive angles are assumed to be clockwise. Moreover 
the y axis is assumed to point downwards. This is an important remark, 
since a typical behaviour of Stylonychia is the performance of clockwise SSRs. 
Without taking into account the convention used in the dataset the movies would
incorrectly show counterclockwise SSRs. In other functions this difference is 
immaterial since it does not affect the induced mean squared displacement, however
in movis it wuld incorrectly show counter clockwise side stepping reactions.


