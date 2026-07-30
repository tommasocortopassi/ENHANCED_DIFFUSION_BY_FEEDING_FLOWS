# ENHANCED DIFFUSION BY COLLECTIVE FEEDING FLOWS OF UNICELLULAR CILIATES

We provide the codes for:

-   Numerical Monte Carlo simulations based on synthetic populations of 
    cells defined in the paper
  
-   Theoretical model based on binary collisions between tracers and
    their nearest neighbour

More details can be found in the article.

## Repository structure

The project is organized as follows:

# FUNCTIONS
- Reg_Stokeslet_velocity_near_wall.m :  Function computing the velocity field
                                        induced on a fluid by a regularised Stokeslet 
                                        above a no slip floor;

- setup_paths                        :  Script to launch at the beginning. 
                                        See ##Setup.

# SUBFOLDERS
- Numerical_Simulations              :  Contains functions, script and data 
                                        for Monte Carlo simulations of tracer 
                                        particles advected by a synthetic population 
                                        of cells. More details on a separate Notes.md
                                        file inside this folder.

- Theoretical_Model                   : Contains the function and data for 
                                        the theoretical binary collisions model.
                                        More details on a separate Notes.md
                                        file inside this folder.
 
- RAW_DATA                            : Contains the raw Dataset of cells 
                                        trajectories and force orientations,
                                        and the functions needed to extract     
                                        and manipulate such data. More details on 
                                        a separate Notes.md file inside this folder. Dataset
                                        from "Cooperative mixing through hydrodynamic 
                                        interactions in Stylonychia lemnae" by R.Turuban, G.
                                        Noselli, A. Beran, A. De Simone.

- RESULTS_OF_EXPERIMENTS              : Stores the results of both Monte Carlo 
                                        simulations with synthetic populations 
                                        and of the binary collision model,
                                        saved in different subfolders. In the case
                                        of Monte Carlo simulations we organise
                                        subfolders in terms of the radius R
                                        and of the synthetic population used.
                                        We also have a subfolder dedicated to
                                        the case of varying initial tracer height, 
                                        used in the paper to validate
                                        the theoretical binary collision model.
                

- FIGURES                             : Figures used in the paper and in the 
                                        Supporting Material, divided by
                                        the section they were used at.

- MOVIES                              : Movies showing the typical motion of  
                                        Stylonychia and of a tracer advected 
                                        by an ensemble of Stylonychia . 
                                        We report also the data and the 
                                        functions used for making these videos.

## Requirements

The code was developed and tested in MATLAB R2024b.

## Setup

Before running the code, execute: setup_paths

## Notes

AI tools have been used to improve the readability of the documentation and 
of the comments in the codes. See the dedicated section in the article. 
