%Script running a numerical simulation of tracers displaced by 
%actuators whose trajectories are sampled from a given dataset. 
% The parameters are stored in options_synth_population



if ~exist("options_synth_population","var")
    load("options_synth_population.mat");
end


if isfield(options_synth_population,"starting_heights") && ~isempty(options_synth_population.starting_heights)
    starting_heights=options_synth_population.starting_heights;
else
    starting_heights=options_synth_population.z0;
end

for i=1:1:length(start_heights)
options_synth_population.starting_height=start_heights(i);
[times, ~, ~, ~, ~, ~, ~, MSDxy_aux, MSDz_aux]=...
    test_displacement_synth_pop(options_synth_population);
MSD_xy(i)=MSDxy_aux(end);
MSD_z(i)= MSDz_aux(end);
i
end

if isscalar(start_heights)

    s=num2str(start_heights);
    s=strrep(s, '.', '_');
    
    filename = "Test_synth_pop_Np="+ num2str(options_synth_population.Np) ...
        + "_R=" + num2str(options_synth_population.R) + ...
        "_T=" + num2str(options_synth_population.T)+ ...
        "_only_NN=" + num2str(options_synth_population.only_NN)+ ...
        "_f=" + num2str(options_synth_population.f)+ ...
        "_randomise=" + num2str(options_synth_population.randomise_angle) + s+".mat";
    clear s;
else
     filename = "Test_synth_pop_Np="+ num2str(options_synth_population.Np) ...
        + "_R=" + num2str(options_synth_population.R) + ...
        "_T=" + num2str(options_synth_population.T)+ ...
        "_only_NN=" + num2str(options_synth_population.only_NN)+ ...
        "_f=" + num2str(options_synth_population.f)+ ...
        "_randomise=" + num2str(options_synth_population.randomise_angle) + ...
        "variable_height"+".mat";
end
project_root = fileparts(which('setup_paths'));
out_dir = fullfile(project_root, 'RESULTS_OF_EXPERIMENTS', 'NUMERICAL_SIMULATIONS_RESULTS');

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

% Initial name
full_filename = fullfile(out_dir, filename);

% If it alredy exists add _1, _2, ...
counter = 1;
while exist(full_filename, 'file')
    full_filename = fullfile(out_dir, ...
        erase(filename, ".mat") + "_" + num2str(counter) + ".mat");
    counter = counter + 1;
end

save(full_filename);

clearvars -except options_synth_population;