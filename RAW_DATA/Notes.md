## Notes
We provide the raw data on Stylonychia trajectories as "Dataset1.csv". 
The data is taken from the Supporting Information of 
"Cooperative mixing through hydrodynamic interactions in Stylonychia lemnae"
by R.Turuban, G.Noselli, A.Beran, A.De Simone.

## Functions

More details are contained in the documentation of each function.

 - Data_collect             :   Extracts data from the dataset so that it 
                                is usable in numerical experiments.
 
 - extract_SSR_angles       :   Extracts a probability distribution
                                for angles of SSRs of Stylonychia.
                                REMARK: When used with the dataset provided, 
                                change the sign of cells angles. In the dataset,
                                angles rotating clockwise are considered positive, 
                                while in the code we use the convention of 
                                considering counterclockwise angles as positive.  

- compute_pmf_t             :   Computes the probability distributions of active
                                and inactive times for Stylonychia from 
                                the available data.

- apply_SSR_rotation_shift  :   Applies a fixed angle rotation to the forcing 
                                directions of the cells after every SSR event.