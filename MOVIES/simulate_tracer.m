function tracer_pos = simulate_tracer(times_data, particle, x0, f, h, delta)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SIMULATE_TRACER  
    % 
    % Integrates a tracer particle driven by all Stokeslet actuators.
    %
    % The coordinates in `particle` are in micrometers and are divided
    % by 200 internally to match the rescaled length unit used in the article.
    % The tracer dynamics are solved with ode45 using the total_v velocity field
    % (sum of contributions from every actuator at each instant).
    %
    % -------------------------------------------------------------------------
    % INPUTS
    %   times_data  - (1 x Nt) or (Nt x 1) time vector [s] corresponding to the
    %                 columns of each actuator trajectory.
    %
    %   particle    - Struct of actuator data. Fields must be named 'p1', 'p2', ...
    %                 Each field is a (3 x Nt) matrix:
    %                   rows 1-2 : x, y positions [raw units, divided by 200 here]
    %                   row  3   : orientation angle theta [rad]
    %
    %   x0          - (3 x 1) initial tracer position [rescaled units]: [z; x; y].
    %
    %   f           - Scalar. Actuator force/strength parameter passed to total_v.
    %
    %   h           - Scalar. Fixed actuator height (z-coordinate) [rescaled units].
    %
    %   delta       - Scalar. Radius of actuators (i.e. of regularised Stokeslets).
    %
    % -------------------------------------------------------------------------
    % OUTPUTS
    %
    %   tracer_pos  - (3 x Nt) integrated tracer trajectory [rescaled units]:
    %                   row 1 : z coordinate
    %                   rows 2-3 : x, y coordinates
    %
    % -------------------------------------------------------------------------
    % EXAMPLE
    %   [times, pos] = simulate_tracer(t, particle, [0.5; 0; 0], 1.0, 0.75, 0.1);
    %
    % -------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % =========================================================================
    % 1.  BUILD ACTUATOR POSITION AND ORIENTATION CELL ARRAYS
    % =========================================================================
    
    names    = fieldnames(particle);
    act_names = names(startsWith(string(names), 'p'));   % fields p1, p2, ...
    Np       = numel(act_names);
    Nt       = numel(times_data);
    dt       = times_data(2) - times_data(1);            % assumes uniform spacing
    
    pos_cell = cell(1, Np);   % pos_cell{j} : 3 x Nt  [z; x; y] in rescaled units
    ori_cell = cell(1, Np);   % ori_cell{j} : 1 x Nt  orientation angle [rad]
    
    for j = 1:Np
        A          = particle.(act_names{j});
        pos_cell{j} = [h * ones(1, Nt); ...   % z is fixed at actuator height
                       A(1:2, :) / 200];      % xy divided by 200 → rescaled units
        ori_cell{j} = A(3, :);                % orientation row (already in rad)
    end
    
    % =========================================================================
    % 2.  INTEGRATE TRACER WITH ODE45
    % =========================================================================
    
    options = odeset('AbsTol', 1e-6, 'RelTol', 1e-4);
    
    [~, X_ode] = ode45( ...
        @(t, x) total_v(t, x, f, h, delta, ori_cell, pos_cell, times_data(1), dt), ...
        times_data, x0, options);
    
    % ode45 returns (Nt x 3); transpose to (3 x Nt) to match column-per-timestep convention
    tracer_pos = X_ode.';

end  
