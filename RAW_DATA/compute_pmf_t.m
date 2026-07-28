function [t, pmf_t_active, pmf_t_inactive] = compute_pmf_t(particle, dt)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% COMPUTE_PMF_T
    %
    % The function takes as inputs trajectories of actuators and the time 
    % step and outputs the PMFs (probability mass functions) of the activation and waiting times.
    %
    % INPUT
    %   particle        :   Struct with fields .p1, .p2, ... containing 3xN matrices. 
    %                       The first 2 rows represent [x,y] trajectories of actuators.
    %                       The third row represents the direction the force is directed to.
    %   dt              :   Duration of a time frame
    %
    % OUTPUT
    %   t               :   Vector [dt, 2*dt, ..., K*dt] of times
    %   pmf_t_active    :   PMF of activation times
    %   pmf_t_inactive  :   PMF of waiting times
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if nargin < 2 || isempty(dt) || ~isscalar(dt) || dt <= 0
        error('dt must be a positive real number.');
    end
    
    %% Extract field names and initialize global variables
    fields = fieldnames(particle);
    all_active_durations = [];
    all_inactive_durations = [];
    tot_nans = 0;
    tot_frames = 0;
    
    %% Loop over actuators
    for i = 1:numel(fields)
        
        A = particle.(fields{i});
        if isempty(A)
            continue;
        end
        if size(A,1) < 1
            continue;
        end
        
        %Since I only care about entries of columns of A 
        %being NaN or not, I extract the first row only
        p = A(1,:);
    
        if isempty(p)
            continue;
        end
        
        tot_nans   = tot_nans + sum(isnan(p));
        tot_frames = tot_frames + numel(p);
    
        %% Active segments (not NaN) 
        active = ~isnan(p);
        d_active = diff([false, active, false]);
        
        %Find start and end positions of active segments
        start_active = find(d_active == 1);
        end_active   = find(d_active == -1) - 1;
        
        %Find lengths of active segments
        active_lengths = end_active - start_active + 1;
    
        if ~isempty(active_lengths)
            all_active_durations = [all_active_durations, active_lengths]; 
        end
    
        %% Inactive segments (NaN)
    
        inactive = isnan(p);
        d_inactive = diff([false, inactive, false]);
    
        %Find start and end positions of inactive segments
        start_inactive = find(d_inactive == 1);
        end_inactive   = find(d_inactive == -1) - 1;
        
        %Find lengths of inactive segments
        inactive_lengths = end_inactive - start_inactive + 1;
    
        if ~isempty(inactive_lengths)
            all_inactive_durations = [all_inactive_durations, inactive_lengths]; 
        end
    
    end
    
    %Maximum duration of active/inactive times
    max_time=max([all_active_durations,all_inactive_durations]);
    
    %% PMF active times
    if isempty(all_active_durations)
        pmf_t_active = [];
    else
        counts_active = histcounts(all_active_durations,...
            0.5:1:(max_time + 0.5));
        pmf_t_active = counts_active / sum(counts_active);
    end
    
    %% PMF inactive times
    
    if isempty(all_inactive_durations)
        pmf_t_inactive = [];
    else
        counts_inactive = histcounts(all_inactive_durations, ...
            0.5:1:(max_time + 0.5));
        pmf_t_inactive = counts_inactive / sum(counts_inactive);
    end
    
    %Compute possible active/inactive times
    t = (1:max_time) * dt;
end