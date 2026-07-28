function [times, particle] = Data_collect(filename, noise_threshold)
    %% DATA_COLLECT  Read microorganism trajectories from a CSV file.
    %
    % The input file is assumed to contain one time row followed by data rows
    % arranged in blocks of three for each microorganism:
    %   1) x-position,
    %   2) y-position,
    %   3) force orientation.
    %
    % Microorganisms containing too many missing values are discarded. The
    % admissible trajectories are returned as fields p1, p2, ... of the output
    % structure.
    %
    % INPUT
    %   filename            :   Name of the CSV file, with or without the '.csv'
    %                           extension.
    %   noise_threshold     :   Scalar between 0 and 1. Maximum allowed fraction
    %                           of NaN entries in the 3 x Nt data block
    %                           associated with one microorganism. For example,
    %                           if noise_threshold = 0.5, then a microorganism
    %                           is discarded when more than 50% of its stored
    %                           entries are NaN.
    %
    % OUTPUT
    %   times               :   1 x Nt vector of time samples, given by the first row of the
    %                           CSV file.
    %   particle            :   Structure whose fields p1, p2, ... contain the admissible
    %                           microorganism trajectories. Each field is a 3 x Nt matrix
    %                           whose first two rows are the x- and y-positions, and whose
    %                           third row is the force orientation in radians.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Validate the noise tolerance
    if ~isscalar(noise_threshold) || noise_threshold < 0 || noise_threshold > 1
        error('noise_threshold must be a scalar between 0 and 1.');
    end
    
    % Accept filenames provided either with or without the .csv extension
    filename = string(filename);
    if ~endsWith(filename, ".csv")
        filename = filename + ".csv";
    end
    
    % Read the data matrix from file
    data = readmatrix(filename);
    
    % Basic format checks
    if isempty(data) || size(data, 1) < 4
        error('The input file must contain at least one time row and one 3-row trajectory block.');
    end
    
    times = data(1, :);
    coordinates = data(2:end, :);
    
    if mod(size(coordinates, 1), 3) ~= 0
        error('The number of data rows after the time row must be a multiple of 3.');
    end
    
    N = size(coordinates, 1) / 3;
    
    particle = struct();
    j = 1;
    
    for i = 1:N
        rows = 3*(i-1) + (1:3);
        block = coordinates(rows, :);
    
        % Compute the fraction of missing entries over the full 3 x Nt block
        nan_fraction = sum(isnan(block), 'all') / numel(block);
    
        if nan_fraction <= noise_threshold
            name = "p" + num2str(j);
            particle.(name) = block;
            j = j + 1;
        end
    end
end