function p = select_synth_population(particles, randomise_angle)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SELECT_SYNTH_POPULATION
    %
    % The function takes a struct of trajectories in input and returns a
    % modified struct. 
    %
    % INPUT
    %   - particles                     : Structure containing the trajectories
    %   - randomise_angle               : Parameter whose value is either 0 or 1.
    %
    %               * randomise_angle=0 : Does not modify force orientations.
    %                                     Corresponds to Population B in the article.
    %               * randomise_angle=1 : Randomises forcing directions. 
    %                                     Corresponds to Population A in the article.
    % OUTPUT
    %   - p                             : modified structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 || isempty(randomise_angle)
        randomise_angle = 0;
    end
    
    
    p = particles;
    
    prefix = "p";
    all_names = fieldnames(particles);
    mask = startsWith(string(all_names), prefix);
    names = all_names(mask);
    
    for i = 1:numel(names)
        name = names{i};
        A = particles.(name);
    
        % Active frame = column without NaN
        is_active_frame = ~any(isnan(A), 1);
    
        % Find blocks of active frames
        diff_mask = diff([0, is_active_frame, 0]);
        starts = find(diff_mask == 1);
        ends = find(diff_mask == -1) - 1;
    
        for b = 1:numel(starts)
            idx_start = starts(b);
            idx_end = ends(b);
            block_length = idx_end - idx_start + 1;
    
    
            % Eventually fix orientation
            if randomise_angle == 0 
                angle = A(3, idx_start:idx_end);
            elseif randomise_angle == 1
                angle = A(3, idx_start:idx_end) + 2 * pi *rand ;
                angle = mod(angle + pi, 2*pi) - pi;
            else
                error("Invalid values for randomise_angle: select either 0 or 1")
            end
            p.(name)(3, idx_start:idx_end) = angle;
        end
    end
end