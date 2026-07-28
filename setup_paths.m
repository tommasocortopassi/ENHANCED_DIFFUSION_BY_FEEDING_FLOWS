function setup_paths()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SETUP_PATHS 
    %   Adds the project directory and all sub-folders to the MATLAB path.
    %
    %   SETUP_PATHS() determines the directory containing this function and
    %   recursively adds it, together with all its sub-folders, to the MATLAB
    %   search path using GENPATH and ADDPATH.
    %
    %   This function should be called once at the beginning of a MATLAB
    %   session to ensure that all project scripts, functions, and utilities
    %   are accessible regardless of the current working directory.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        root = fileparts(mfilename('fullpath'));
        addpath(genpath(root));
end